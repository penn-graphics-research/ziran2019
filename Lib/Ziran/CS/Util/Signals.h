#ifndef SIGNALS_H
#define SIGNALS_H
#include <atomic>
#include <csignal>
#include <thread>

#include <Ziran/CS/Util/Meta.h>

namespace ZIRAN {
struct Signal {
    int number;
};

std::atomic<uint64_t>& firstExit();

std::atomic<int>& exitSignalNumber();

std::thread::id exitingThreadId();

bool isExitingThread();

void setExitSignalNumber(int number);

int getExitSignalNumber();

bool interupted();

bool interuptedTwice();

bool shouldDoExit();

void restoreSignalHandler(int signal_number);

void exitWithDefaultSignalHandler(Signal signal);

/// string representation of signal ID
std::string exitReasonName(Signal signal);

struct SignalBlocker {
    sigset_t set;

    template <class... SignalNumber>
    SignalBlocker(SignalNumber... s)
    {
        sigemptyset(&set);
        Call{ sigaddset(&set, s)... };
        pthread_sigmask(SIG_BLOCK, &set, NULL);
    }

    ~SignalBlocker()
    {
        pthread_sigmask(SIG_UNBLOCK, &set, NULL);
    }
};
} // namespace ZIRAN
#endif
