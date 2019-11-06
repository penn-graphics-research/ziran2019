#ifndef STACK_TRACE_H
#define STACK_TRACE_H

#include <string>

namespace ZIRAN {

/*
   To avoid allocating in a signal handler the string
   can be allocated before hand
*/
std::string& stackTraceString();
std::string& stackTrace(int frames_to_skip = 1);
} // namespace ZIRAN
#endif
