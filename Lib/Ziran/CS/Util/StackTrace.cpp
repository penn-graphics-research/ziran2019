#include <Ziran/CS/Util/StackTrace.h>

namespace ZIRAN {

/*
   To avoid allocating in a signal handler the string
   can be allocated before hand
*/
std::string& stackTraceString()
{
    static std::string s;
    s.reserve(600);
    s.clear();
    return s;
}
} // namespace ZIRAN

#ifndef _WIN32
#include <cxxabi.h>
#include <execinfo.h>
#include <string>
#include <cstring>
#include <cstdlib>

namespace ZIRAN {
std::string& stackTrace(int frames_to_skip)
{
    const size_t max_dump_size = 50;
    void* dump[max_dump_size];
    size_t size = backtrace(dump, max_dump_size);
    char** messages = backtrace_symbols(dump, size); // overwrite sigaction with caller's address

    // dump stack: skip first frame, since that is here
    std::string& oss = stackTraceString();
    for (size_t idx = frames_to_skip; idx < size && messages != nullptr; ++idx) {
        char *mangled_name = 0, *offset_begin = 0;
#ifdef __APPLE__
        // OSX style stack trace
        for (char* p = messages[idx]; *p; ++p) {
            if ((*p == '_') && (*(p - 1) == ' '))
                mangled_name = p - 1;
            else if (*p == '+')
                offset_begin = p - 1;
        }

        if (mangled_name && offset_begin && (mangled_name < offset_begin)) {
            *mangled_name++ = '\0';
            *offset_begin++ = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply
            // __cxa_demangle():
            int status;
            char* real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);
            if (status == 0) {
                oss += messages[idx];
                oss += ' ';
                oss += real_name;
                oss += offset_begin;
                oss += '\n';
            }
            else {
                oss += messages[idx];
                oss += mangled_name;
                oss += offset_begin;
                oss += '\n';
            }
            free(real_name); // mallocated by abi::__cxa_demangle(...)
        }
        else {
            // no demangling done -- just dump the whole line
            oss += messages[idx];
            oss += '\n';
        }
#else
        char* offset_end = 0;
        // find parantheses and +address offset surrounding mangled name
        for (char* p = messages[idx]; *p; ++p) {
            if (*p == '(') {
                mangled_name = p;
            }
            else if (*p == '+') {
                offset_begin = p;
            }
            else if (*p == ')') {
                offset_end = p;
                break;
            }
        }
        // if the line could be processed, attempt to demangle the symbol
        if (mangled_name && offset_begin && offset_end && mangled_name < offset_begin) {
            *mangled_name++ = '\0';
            *offset_begin++ = '\0';
            *offset_end++ = '\0';

            int status;
            char* real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);
            // if demangling is successful, output the demangled function name
            int frame = idx - frames_to_skip + 1;
            int first_digit = frame / 10;
            oss += (first_digit == 0) ? ' ' : '0' + first_digit;
            oss += '0' + frame % 10;
            oss += "  ";
            oss += messages[idx];
            int extra = 34 - strlen(messages[idx]);
            for (int i = 0; i < extra; i++)
                oss += ' ';
            if (status == 0) {
                oss += "  ";
                oss += real_name;
                oss += '+';
            } // otherwise, output the mangled function name
            else {
                oss += "  ";
                oss += mangled_name;
                oss += '+';
            }
            oss += offset_begin;
            oss += offset_end;
            oss += '\n';
            free(real_name); // mallocated by abi::__cxa_demangle(...)
        }
        else {
            // no demangling done -- just dump the whole line
            int first_digit = idx / 10;
            oss += (first_digit == 0) ? ' ' : '0' + first_digit;
            oss += '0' + idx % 10;
            oss += "  ";
            oss += messages[idx];
            oss += '\n';
        }
#endif
    } // END: for(size_t idx = 1; idx < size && messages != nullptr; ++idx)
    free(messages);
    return oss;
}
} // namespace ZIRAN
#else //Windows
namespace ZIRAN {
std::string stackTrace(int frames_to_skip)
{
    //TODO Figure out how do get a stack trace on windows
    return "";
}
} // namespace ZIRAN
#endif
