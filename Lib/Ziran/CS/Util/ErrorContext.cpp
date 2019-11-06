#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/ErrorContext.h>
namespace ZIRAN {

ErrorContext::ErrorContext(const char* function_name, const char* file, int line_number, const char* var_name)
    : stream(getStream())
{
    std::string ss(file);
    size_t found = ss.find_last_of("/\\");
    stream << ss.substr(found + 1) << ':' << line_number << ' ' << function_name << ' ' << var_name << ":\t";
}

ErrorContext::~ErrorContext()
{
    auto& call_stack = callStack();
    ZIRAN_ASSERT(!call_stack.empty());
    call_stack.pop_back();
}

std::vector<std::string>& ErrorContext::callStack()
{
    static thread_local std::vector<std::string> call_stack;
    return call_stack;
}

void ErrorContext::addContext(std::string& result)
{
    auto& call_stack = callStack();
    for (const auto& s : call_stack)
        result += s;
}

// Stringstreams need a global mutex to be created
// so only do create them one per thread
std::stringstream& ErrorContext::getStream()
{
    const static thread_local std::stringstream default_format;
    static thread_local std::stringstream thread_stream;

    thread_stream.str(std::string());
    thread_stream.clear();
    thread_stream.copyfmt(default_format);
    return thread_stream;
}
} // namespace ZIRAN
