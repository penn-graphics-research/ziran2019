#ifndef TIMER_H
#define TIMER_H
#include <Ziran/CS/Util/Logging.h>
#include <chrono>
#include <stack>
#include <string>

namespace ZIRAN {

/**
    Timer. We can use either system timer or stready timer
*/
class Timer {
public:
    Timer(bool start_now = true);

    /**
 	  \brief Start timing
    */
    void start();

    /**
  	 \return time elapsed since last click in seconds
    */
    std::chrono::duration<double> click(bool reset = true);

    std::chrono::duration<double> dump(const std::string& s, bool reset = true);

private:
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> to_time;
};

class ScopedTimer {
public:
    ScopedTimer(const std::string& function_name, const std::string& file, int line_number, bool print_now = true);

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

    ~ScopedTimer();

    void enteringSubCall();

    void leavingSubCall();

    void printETA(double n);

private:
    TimingResult result;
    Timer timer;
    bool print_now;

    static std::stack<std::reference_wrapper<ScopedTimer>>& callStack();
};

#define ZIRAN_TIMER() ScopedTimer ziran_timer(__func__, __FILE__, __LINE__, true)
#define ZIRAN_QUIET_TIMER() ScopedTimer ziran_timer(__func__, __FILE__, __LINE__, false)

#define ZIRAN_STRING_HELPER(s) #s
#define ZIRAN_STRING(s) ZIRAN_STRING_HELPER(s)

#define ZIRAN_ETA(n) ziran_timer.printETA(n)
} // namespace ZIRAN

#endif
