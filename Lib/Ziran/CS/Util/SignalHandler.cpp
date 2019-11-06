#include <Ziran/CS/Util/SignalHandler.h>
#include <Ziran/CS/Util/Signals.h>
#include <Ziran/CS/Util/Logging.h>
#include <Ziran/CS/Util/StackTrace.h>

#include <mutex>

namespace ZIRAN {

const std::map<int, std::string> kSignals = {
    { SIGINT, "SIGINT" },
    { SIGABRT, "SIGABRT" },
    { SIGFPE, "SIGFPE" },
    { SIGILL, "SIGILL" },
    { SIGSEGV, "SIGSEGV" },
    { SIGTERM, "SIGTERM" },
};

std::map<int, std::string>& gSignals()
{
    static std::map<int, std::string> gSignals = kSignals;
    return gSignals;
}

// Dump of stack then exit through background worker
// ALL thanks to this thread at StackOverflow. Pretty much borrowed from:
// Ref: http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes
void signalHandler(int signal_number)
{
    // Only one signal will be allowed past this point
    if (false == shouldDoExit()) {
        if (isExitingThread())
            return;
        while (true) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    }

    setExitSignalNumber(signal_number);

    {
        Signal signal{ signal_number };
        LogCapture(Fatal, true).capture(signal);
    }
}

//
// Installs FATAL signal handler that is enough to handle most fatal events
//  on *NIX systems
void installSignalHandler()
{
    stackTraceString();
    // do it verbose style - install all signal actions
    const auto& signals = gSignals();
    for (const auto& sig_pair : signals) {
        signal(sig_pair.first, &signalHandler);
    }
}

// This will override the default signal handler setup and instead
// install a custom set of signals to handle
void overrideSetupSignals(const std::map<int, std::string> overrideSignals)
{
    stackTraceString();
    static std::mutex signalLock;
    std::lock_guard<std::mutex> guard(signalLock);
    auto& signals = gSignals();
    for (const auto& sig : signals) {
        restoreSignalHandler(sig.first);
    }

    signals = overrideSignals;
    installSignalHandler(); // installs all the signal handling for gSignals
}

// restores the signal handler back to default
void restoreSignalHandlerToDefault()
{
    overrideSetupSignals(kSignals);
}
} // namespace ZIRAN
