#include <Ziran/CS/Util/Timer.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/PrettyPrinting.h>
#include <ctime>
#include <iomanip>

namespace ZIRAN {
Timer::Timer(bool start_now)
{
    if (start_now)
        start();
}
void Timer::start()
{
    start_time = std::chrono::steady_clock::now();
}

std::chrono::duration<double> Timer::click(bool reset)
{
    to_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = to_time - start_time;
    if (reset)
        start_time = to_time;
    return elapsed_seconds;
}

std::chrono::duration<double> Timer::dump(const std::string& s, bool reset)
{
    std::chrono::duration<double> time = click(reset);
    ZIRAN_TIME(std::fixed, time.count(), "s elapsed for ", s);
    return time;
}

ScopedTimer::ScopedTimer(const std::string& function_name, const std::string& file, int line_number, bool print_now)
    : print_now(print_now)
{
    auto& call_stack = callStack();
    if (!call_stack.empty()) {
        call_stack.top().get().enteringSubCall();
    }
    std::string ss(function_name);
    size_t found = file.find_last_of("/\\");
    ss += " " + file.substr(found + 1) + ":" + std::to_string(line_number);
    result = TimingResult{ ss, (int)call_stack.size(), std::chrono::seconds(0), std::chrono::seconds(0) };
    call_stack.emplace(std::ref(*this));
    ZIRAN_ASSERT(this == &call_stack.top().get());
}

ScopedTimer::~ScopedTimer()
{
    auto& call_stack = callStack();
    auto elapsed_seconds = timer.click();
    result.total_seconds += elapsed_seconds;
    result.self_seconds += elapsed_seconds;
    ZIRAN_ASSERT(!call_stack.empty());
    ZIRAN_ASSERT(this == &call_stack.top().get());
    call_stack.pop();
    if (!call_stack.empty()) {
        call_stack.top().get().leavingSubCall();
    }
    ZIRAN_TIME(result, print_now);
}

void ScopedTimer::enteringSubCall()
{
    auto elapsed_seconds = timer.click();
    result.total_seconds += elapsed_seconds;
    result.self_seconds += elapsed_seconds;
}

void ScopedTimer::leavingSubCall()
{
    auto elapsed_seconds = timer.click();
    result.total_seconds += elapsed_seconds;
}

void ScopedTimer::printETA(double n)
{
    using namespace std::chrono;

    auto elapsed_seconds = timer.click();
    result.total_seconds += elapsed_seconds;
    result.self_seconds += elapsed_seconds;

    auto now = system_clock::now();
    auto estimate = duration_cast<system_clock::duration>(n * result.total_seconds);
    std::time_t finish = system_clock::to_time_t(now + estimate);
    ZIRAN_TIME("ETA in ", PrettyPrintDuration(estimate), '\n', "    on ", std::put_time(std::localtime(&finish), "%a %b %d \n    at %r"));
}

std::stack<std::reference_wrapper<ScopedTimer>>& ScopedTimer::callStack()
{
    static thread_local std::stack<std::reference_wrapper<ScopedTimer>> stack;
    return stack;
}
} // namespace ZIRAN
