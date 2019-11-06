#ifndef DEBUG_UTIL_H
#define DEBUG_UTIL_H

#include <iostream>
#include <sstream>
#include <string>
#include <Ziran/CS/Util/Logging.h>

#define ZIRAN_DEBUG_FUNCTION_NAME ((const char*)__FUNCTION__) // cast to const char* to work around error in noreturn
/**
   Throws a std::runtime_error when the condition is false.
   First parameter is an expression that evaluates to a boolean.
   Second parameter is optional. It can be a string description of the error.
   Use it as follows:
      int t=3;
      ZIRAN_ASSERT(t>4,"t should be greater than 4");
*/

#define ZIRAN_ASSERT(condition, ...) \
    if (condition) {                 \
    }                                \
    else                             \
        ZIRAN::DEBUG_UTILITIES::assertionFailed(ZIRAN_DEBUG_FUNCTION_NAME, __FILE__, __LINE__, #condition, ##__VA_ARGS__)

namespace ZIRAN {

namespace DEBUG_UTILITIES {

template <class... Ts>
[[noreturn]] inline void assertionFailed(const char* function, const char* file, unsigned int line, const char* condition, Ts&&... messages)
{
    ZIRAN_ERR("In file: ", file, ":", line, ":\n",
        "In function ", function, ":\n",
        std::forward<Ts>(messages)...,
        "\ncondition = ", condition);
    throw std::runtime_error("Assertion failed");
}

[[noreturn]] inline void assertionFailed(const char* function, const char* file, unsigned int line, const char* condition)
{
    assertionFailed(function, file, line, condition, "Assertion failed");
}
}
} // namespace ZIRAN::DEBUG_UTILITIES
#endif
