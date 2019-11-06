#ifndef SIGNAL_HANDLER_H
#define SIGNAL_HANDLER_H
#include <map>

namespace ZIRAN {

extern const std::map<int, std::string> kSignals;

std::map<int, std::string>& gSignals();

// Dump of stack then exit through background worker
// ALL thanks to this thread at StackOverflow. Pretty much borrowed from:
// Ref: http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes
void signalHandler(int signal_number);

//
// Installs FATAL signal handler that is enough to handle most fatal events
//  on *NIX systems
void installSignalHandler();

// This will override the default signal handler setup and instead
// install a custom set of signals to handle
void overrideSetupSignals(const std::map<int, std::string> overrideSignals);

// restores the signal handler back to default
void restoreSignalHandlerToDefault();
} // namespace ZIRAN
#endif
