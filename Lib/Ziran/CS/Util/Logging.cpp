#include <algorithm>
#include <iostream>
#include <signal.h>
#include <string>
#include <vector>

#include <Ziran/CS/Util/Active.h>
#include <Ziran/CS/Util/ErrorContext.h>
#include <Ziran/CS/Util/PrettyPrinting.h>
#include <Ziran/CS/Util/StackTrace.h>
#include <Ziran/CS/Util/Logging.h>

#ifndef _WIN32
#include <unistd.h>
#else
#include <WinSock2.h>
#include <Windows.h>
#include <io.h>
#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define isatty _isatty
#endif

namespace ZIRAN {

std::ostream& TimingResult::print(std::ostream& os, bool print_self, int time_width) const
{
    os << std::fixed
       << std::setw(time_width) << total_seconds.count() << "s";
    if (print_self)
        os << " (" << std::setw(time_width) << self_seconds.count() << "s)";
    os << " elapsed for " << name;
    return os;
}

std::string LogMessage::toString()
{
    std::stringstream is(entry);
    std::stringstream os;
    std::string line;
    std::getline(is, line);
    Color c = (level < Warning) ? Color::GREY : (level < Error) ? Color::YELLOW : Color::RED;
    std::string indent;
    if (level >= Fatal) {
        os << "\n"
           << c << level.text << Color::RESET
           << "  Received fatal signal: "
           << exitReasonName(signal)
           << "(" << signal.number << ")\n";
    }
    else {
        os << c << level.text << Color::RESET << "  " << line << '\n';
        indent = "         ";
    }
    while (std::getline(is, line)) {
        os << indent << line << '\n';
    }
    return os.str();
}

void LogWorker::bgSave(LogMessage&& msg)
{
    static bool overloaded = false;
    static unsigned count = 0;
    count++;
    if (msg.level >= Warning || !overloaded)
        write(std::move(msg), file);

    if (!overloaded && count == 1024) {
        count = 0;
        unsigned backlog = bg->size();
        if (backlog > 1024 * 1024) {
            overloaded = true;
            std::cerr << Color::YELLOW << "WARNING" << Color::RESET << "  Logging overloaded" << std::endl;
            file << "WARNING  Logging overloaded\n";
        }
    }
    if (overloaded) {
        unsigned backlog = bg->size();
        if (backlog < 1024) {
            overloaded = false;
            std::cerr << Color::YELLOW << "WARNING" << Color::RESET << "  Resuming logging" << std::endl;
            file << "WARNING  Resuming logging\n";
            count = 0;
        }
    }
}

void LogWorker::bgRecordTiming(TimingResult&& tr, std::thread::id id)
{
    auto search = timings.emplace(id, StdVector<TimingResult>());
    if (search.second) {
        thread_ids.emplace_back(id);
    }
    auto& thread_timings = search.first->second;
    bool found = false;
    for (auto& recorded : thread_timings)
        if (recorded.name == tr.name) {
            found = true;
            recorded.total_seconds += tr.total_seconds;
            recorded.self_seconds += tr.self_seconds;
        }
    if (!found) {
        thread_timings.emplace_back(std::move(tr));
    }
}

void LogWorker::bgPrintTimings()
{
    std::stringstream ss;
    for (size_t i = 0; i < thread_ids.size(); i++) {
        ss << "Total timing for thread: " << i << '\n';
        for (const auto& recorded : timings[thread_ids[i]]) {
            recorded.print(ss, true, 12);
            ss << '\n';
        }
    }
    if (logLevelEnabled(Timing))
        write(LogMessage{ Timing, { SIGABRT }, ss.str() }, file);
}

void LogWorker::bgOpenLogFile(const std::string& filename, bool append)
{
    std::ios_base::iostate exceptionMask = file.exceptions() | std::ios::failbit;
    file.exceptions(exceptionMask);
    auto mode = std::fstream::out;
    if (append)
        mode |= std::fstream::app;
    else
        mode |= std::fstream::trunc;
    file.open(filename, mode);
}

std::weak_ptr<LogWorker>& LogWorker::getGlobal()
{
    static std::weak_ptr<LogWorker> worker;
    return worker;
}

void LogWorker::write(LogMessage&& msg, std::ofstream& file)
{
    std::string s = msg.toString();
    std::string stripped;
    if (msg.level >= Level::getConsoleLevel()) {
        std::ostream& out = (msg.level >= Error) ? std::cerr : std::cout;
        if (isatty(STDOUT_FILENO)) {
            out << s;
        }
        else {
            stripped = stripColors(s);
            out << stripped;
        }
        out.flush();
    }

    if (file.is_open()) {
        if (stripped.empty())
            stripped = stripColors(s);
        file << stripped;
    }

    if (msg.level >= Fatal) {
        finishWriting(file);
        exitWithDefaultSignalHandler(msg.signal);
        // should never reach this point
        perror("Exited after receiving FATAL trigger. Flush message status: ");
        exit(msg.signal.number);
    }
}

void LogWorker::finishWriting(std::ofstream& file)
{
    std::cout.flush();
    std::cerr.flush();
    file.flush();
    file.close();
}

LogWorker::LogWorker(bool should_exit)
    : bg(Active::create([should_exit]() {
        bool should_stop = interuptedTwice();
        if (should_stop && should_exit) {
            Signal s{ getExitSignalNumber() };
            exitWithDefaultSignalHandler(s);
        }
        return should_stop;
    }))
{
}

void LogWorker::save(LogMessage&& msg)
{
    bg->send([this, m(std::move(msg))]() mutable { bgSave(std::move(m)); });
}

void LogWorker::save(TimingResult&& msg, std::thread::id id)
{
    bg->send([this, m(std::move(msg)), id]() mutable { bgRecordTiming(std::move(m), id); });
}

void LogWorker::openLogFile(const std::string& filename, bool append)
{
    bg->send([this, filename, append]() mutable { bgOpenLogFile(filename, append); });
}

void LogWorker::printTimings()
{
    bg->send([this]() { bgPrintTimings(); });
}

/** Write to the global background logger if it exists otherwise write synchronously
      */
void LogWorker::log(LogMessage&& msg)
{
    if (auto global_logger = getGlobal().lock()) {
        global_logger->save(std::move(msg));
        //Stop any futher messages after this
        if (msg.level >= Fatal) {
            getGlobal().reset();
            while (!isExitingThread()) {
                //Wait for the logging thread to exit the program
                std::this_thread::sleep_for(std::chrono::seconds(1));
            }
        }
    }
    else {
        std::ofstream fakeFile;
        write(std::move(msg), fakeFile);
    }
}

std::shared_ptr<LogWorker> LogWorker::initializeLogging()
{
    initializeColors();
    // Block SIGINT from the logging thread
    SignalBlocker(SIGINT);
    if (auto global_logger = getGlobal().lock())
        return global_logger;
    else {
        auto new_logger = std::make_shared<LogWorker>(true);
        getGlobal() = new_logger;
        new_logger->bg->send([]() { exitingThreadId(); });
        return new_logger;
    }
}

void LogWorker::recordTiming(const TimingResult& tr, const std::thread::id& id)
{
    if (auto global_logger = getGlobal().lock()) {
        global_logger->save(TimingResult{ tr }, id);
    }
}

LogCapture::LogCapture(const Level& level, bool in_signal_handler)
    : stream(getStream())
    , level(level)
    , signal{ SIGABRT }
    , in_signal_handler(in_signal_handler)
    , send(true)
{
    // Check if a signal has been caught
    // If so wait for program to end and
    // don't add to the log
    if (interupted() && !in_signal_handler && !isExitingThread())
        while (true) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
}

LogCapture::~LogCapture()
{
    int count = std::min(ZIRAN_MAX_LOG_ENTRY, (int)stream.tellp());
    std::string result(count, ' ');
    stream.read(&result[0], count);
    if ((level >= Fatal) && (signal.number != SIGINT)) {
        result += "\n\n----------------------------------------STACKTRACE-----------------------------------------\n\n";
        result += stackTrace(2);
        result += "\n----------------------------------------CONTEXT--------------------------------------------\n\n";
        ErrorContext::addContext(result);
        result += "\n-------------------------------------------------------------------------------------------\n";
    }
    if (send)
        LogWorker::log(LogMessage{ level, signal, result });
}

void LogCapture::capture()
{
}

void LogCapture::capture(TimingResult& tr, bool print_now)
{
    LogWorker::recordTiming(tr, std::this_thread::get_id());
    if (print_now)
        tr.print(stream);
    else
        send = false;
}

void LogCapture::capture(const TimingResult& tr, bool print_now)
{
    LogWorker::recordTiming(tr, std::this_thread::get_id());
    if (print_now)
        tr.print(stream);
    else
        send = false;
}

// Stringstreams need a global mutex to be created
// so only do create them one per thread
std::stringstream& LogCapture::getStream()
{
    const static thread_local std::stringstream default_format;
    static thread_local std::stringstream thread_stream;

    thread_stream.str(std::string());
    thread_stream.clear();
    thread_stream.copyfmt(default_format);
    return thread_stream;
}
} // namespace ZIRAN
