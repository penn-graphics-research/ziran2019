#ifndef ERROR_CONTEXT_H
#define ERROR_CONTEXT_H
#include <Ziran/CS/Util/Forward.h>
#include <functional>
#include <vector>
#include <string>
#include <sstream>

namespace ZIRAN {

class ErrorContext {
public:
    ErrorContext(const char* function_name, const char* file, int line_number, const char* var_name);

    ErrorContext(const ErrorContext&) = delete;
    ErrorContext& operator=(const ErrorContext&) = delete;

    ~ErrorContext();

    template <typename Derived>
    void capture(const Eigen::DenseBase<Derived>& var)
    {
        Eigen::IOFormat fmt(4, Eigen::DontAlignCols, ", ", ";", "", "", "[", "]");
        stream << var.format(fmt) << '\n';
        auto& call_stack = callStack();
        call_stack.emplace_back(stream.str());
    }

    template <typename T>
    void capture(const T& var)
    {
        stream << var << '\n';
        auto& call_stack = callStack();
        call_stack.emplace_back(stream.str());
    }

    static void addContext(std::string& result);

private:
    // Stringstreams need a global mutex to be created
    // so only create them one per thread
    std::stringstream& getStream();

    std::stringstream& stream;

    static std::vector<std::string>& callStack();
};

#define ZIRAN_COMBINE1(X, Y) X##Y
#define ZIRAN_COMBINE(X, Y) ZIRAN_COMBINE1(X, Y)

#define ZIRAN_GET_MACRO(_1, _2, NAME, ...) NAME

#define ZIRAN_CONTEXT2(string, var)                                                        \
    ErrorContext ZIRAN_COMBINE(context, __LINE__){ __func__, __FILE__, __LINE__, string }; \
    ZIRAN_COMBINE(context, __LINE__)                                                       \
        .capture(var)

#define ZIRAN_CONTEXT1(var) ZIRAN_CONTEXT2(#var, var)

#define ZIRAN_CONTEXT(...)                                       \
    ZIRAN_GET_MACRO(__VA_ARGS__, ZIRAN_CONTEXT2, ZIRAN_CONTEXT1) \
    (__VA_ARGS__)

#ifdef NDEBUG
#define ZIRAN_DEBUG_CONTEXT(var)
#else
#define ZIRAN_DEBUG_CONTEXT(...) ZIRAN_CONTEXT(__VA_ARGS__)
#endif
} // namespace ZIRAN

#endif
