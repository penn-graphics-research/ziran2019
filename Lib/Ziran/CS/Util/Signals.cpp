#include <Ziran/CS/Util/Signals.h>

namespace ZIRAN {

std::atomic<uint64_t>& firstExit()
{
    static std::atomic<uint64_t> first_exit{ 0 };
    return first_exit;
}

std::atomic<int>& exitSignalNumber()
{
    static std::atomic<int> exit_signal{ SIGINT };
    return exit_signal;
}

std::thread::id exitingThreadId()
{
    const static std::thread::id exiting_thread_id{ std::this_thread::get_id() };
    return exiting_thread_id;
}

bool isExitingThread()
{
    return exitingThreadId() == std::this_thread::get_id();
}

void setExitSignalNumber(int number)
{
    exitSignalNumber().store(number, std::memory_order_release);
}

int getExitSignalNumber()
{
    return exitSignalNumber().load(std::memory_order_acquire);
}

bool interupted()
{
    auto const count = firstExit().load(std::memory_order_relaxed);
    return (count > 0);
}

bool interuptedTwice()
{
    auto const count = firstExit().load(std::memory_order_relaxed);
    return (count > 1);
}

bool shouldDoExit()
{
    auto const count = firstExit().fetch_add(1, std::memory_order_relaxed);
    return (0 == count);
}

void restoreSignalHandler(int signal_number)
{
    signal(signal_number, SIG_DFL);
}

void exitWithDefaultSignalHandler(Signal signal)
{
    const int signal_number = signal.number;
    restoreSignalHandler(signal_number);
    raise(signal_number);
}

/// string representation of signal ID
std::string exitReasonName(Signal signal)
{
    switch (signal.number) {
    case SIGINT:
        return "SIGINT";
        break;
    case SIGABRT:
        return "SIGABRT";
        break;
    case SIGFPE:
        return "SIGFPE";
        break;
    case SIGSEGV:
        return "SIGSEGV";
        break;
    case SIGILL:
        return "SIGILL";
        break;
    case SIGTERM:
        return "SIGTERM";
        break;
    default:
        std::string s;
        s.reserve(18);
        s = "UNKNOWN SIGNAL(" + std::to_string(signal.number) + ")";
        return s;
    }
}
} // namespace ZIRAN
