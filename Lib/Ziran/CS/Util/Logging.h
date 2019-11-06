#ifndef LOGGING_H
#define LOGGING_H
#include <chrono>
#include <fstream>
#include <unordered_map>
#include <iomanip>

#include <Ziran/CS/Util/Signals.h>

namespace ZIRAN {

class Active;

/** The Ziran logging design is based on that of g3log.  It has a background thread which handles all the io.  */

#define ZIRAN_LOG(level, ...)                       \
    do {                                            \
        if (logLevelEnabled(level))                 \
            LogCapture(level).capture(__VA_ARGS__); \
    } while (false)

#define ZIRAN_LOG_IF(level, boolean_expression, ...)    \
    do {                                                \
        if (logLevelEnabled(level))                     \
            if ((boolean_expression))                   \
                LogCapture(level).capture(__VA_ARGS__); \
    } while (false)

#define ZIRAN_MAX_LOG_ENTRY 16384

struct Level {

    const int value;
    const char text[8];

    static int& getLogLevel()
    {
        static int logLevel = 0;
        return logLevel;
    }

    static int& getConsoleLevel()
    {
        static int consoleLevel = 0;
        return consoleLevel;
    }

    bool operator<(const Level& other) const
    {
        return value < other.value;
    }
    bool operator<=(const Level& other) const
    {
        return value <= other.value;
    }
    bool operator>(const Level& other) const
    {
        return value > other.value;
    }
    bool operator>=(const Level& other) const
    {
        return value >= other.value;
    }
    bool operator==(const Level& other) const
    {
        return value == other.value;
    }
    bool operator!=(const Level& other) const
    {
        return value != other.value;
    }
    bool operator<(int other) const
    {
        return value < other;
    }
    bool operator<=(int other) const
    {
        return value <= other;
    }
    bool operator>(int other) const
    {
        return value > other;
    }
    bool operator>=(int other) const
    {
        return value >= other;
    }
    bool operator==(int other) const
    {
        return value == other;
    }
    bool operator!=(int other) const
    {
        return value != other;
    }
};

const Level Debug{ 0, "  DEBUG" },
    Verbose{ 50, "VERBOSE" },
    Info{ 100, "   INFO" },
    Timing{ 200, " TIMING" },
    Event{ 400, "  EVENT" },
    Warning{ 500, "WARNING" },
    Error{ 1000, "  ERROR" },
    Fatal{ 9000, "  FATAL" };

struct TimingResult {
    std::string name;
    int level;
    std::chrono::duration<double> total_seconds;
    std::chrono::duration<double> self_seconds;

    std::ostream& print(std::ostream& os, bool print_self = false, int time_width = 8) const;
};

#define ZIRAN_DBUG(...) ZIRAN_LOG(Debug, ##__VA_ARGS__)
#define ZIRAN_DEBUG(X) ZIRAN_DBUG(#X, ": \n", X)
#define ZIRAN_VERB(...) ZIRAN_LOG(Verbose, ##__VA_ARGS__)
#define ZIRAN_INFO(...) ZIRAN_LOG(Info, ##__VA_ARGS__)
#define ZIRAN_TIME(...) ZIRAN_LOG(Timing, ##__VA_ARGS__)
#define ZIRAN_EVENT(...) ZIRAN_LOG(Event, ##__VA_ARGS__)
#define ZIRAN_WARN(...) ZIRAN_LOG(Warning, ##__VA_ARGS__)
#define ZIRAN_ERR(...) ZIRAN_LOG(Error, ##__VA_ARGS__)
#define ZIRAN_FATAL(...) ZIRAN_LOG(Fatal, ##__VA_ARGS__)

#define ZIRAN_DBUG_IF(...) ZIRAN_LOG_IF(Debug, ##__VA_ARGS__)
#define ZIRAN_DEBUG_IF(condition, X) ZIRAN_DBUG_IF(condition, #X, ": \n", X)
#define ZIRAN_VERB_IF(...) ZIRAN_LOG_IF(Verbose, ##__VA_ARGS__)
#define ZIRAN_INFO_IF(...) ZIRAN_LOG_IF(Info, ##__VA_ARGS__)
#define ZIRAN_TIME_IF(...) ZIRAN_LOG_IF(Timing, ##__VA_ARGS__)
#define ZIRAN_EVENT_IF(...) ZIRAN_LOG_IF(Timing, ##__VA_ARGS__)
#define ZIRAN_WARN_IF(...) ZIRAN_LOG_IF(Warning, ##__VA_ARGS__)
#define ZIRAN_ERR_IF(...) ZIRAN_LOG_IF(Error, ##__VA_ARGS__)
#define ZIRAN_FATAL_IF(...) ZIRAN_LOG_IF(Fatal, ##__VA_ARGS__)

inline bool logLevelEnabled(const Level& level)
{
    return level >= Level::getLogLevel();
}

struct LogMessage {

    std::string toString();

    Level level;
    Signal signal;
    std::string entry;
};

class LogWorker final {
    std::ofstream file;
    std::unordered_map<std::thread::id, StdVector<TimingResult>> timings;
    StdVector<std::thread::id> thread_ids;
    std::unique_ptr<Active> bg;

    void bgSave(LogMessage&& msg);

    void bgRecordTiming(TimingResult&& tr, std::thread::id id);

    void bgPrintTimings();

    void bgOpenLogFile(const std::string& filename, bool append);

    static std::weak_ptr<LogWorker>& getGlobal();

    static void write(LogMessage&& msg, std::ofstream& file);

    static void finishWriting(std::ofstream& file);

public:
    LogWorker(bool should_exit = false);

    LogWorker(const LogWorker&) = delete;
    LogWorker& operator=(const LogWorker&) = delete;

    void save(LogMessage&& msg);

    void save(TimingResult&& msg, std::thread::id id);

    void openLogFile(const std::string& filename, bool append);

    void printTimings();

    /** Write to the global background logger if it exists otherwise write synchronously
      */
    static void log(LogMessage&& msg);

    static std::shared_ptr<LogWorker> initializeLogging();

    static void recordTiming(const TimingResult& tr, const std::thread::id& id);
};

class LogCapture {
public:
    LogCapture(const Level& level, bool in_signal_handler = false);

    ~LogCapture();

    void capture();

    void capture(TimingResult& tr, bool print_now);

    void capture(const TimingResult& tr, bool print_now);

    template <class... Ts>
    void capture(Signal s, Ts&&... rest)
    {
        signal = s;
        capture(std::forward<Ts>(rest)...);
    }

    template <class T, class... Ts>
    void capture(T&& message, Ts&&... rest)
    {
        if (stream.tellp() < ZIRAN_MAX_LOG_ENTRY) {
            stream << std::forward<T>(message);
            capture(std::forward<Ts>(rest)...);
        }
    }

    // Stringstreams need a global mutex to be created
    // so only do create them one per thread
    std::stringstream& getStream();

    std::stringstream& stream;
    const Level& level;
    Signal signal;
    bool in_signal_handler;
    bool send;
};
} // namespace ZIRAN
#endif
