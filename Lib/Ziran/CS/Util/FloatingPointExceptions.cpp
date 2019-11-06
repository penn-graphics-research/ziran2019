#include <Ziran/CS/Util/FloatingPointExceptions.h>

namespace ZIRAN {

namespace FPE {
std::ostream& operator<<(std::ostream& os, State s)
{
    std::string sep = "";
    auto print_cause = [&](State check, const std::string& name) mutable {
        if (State::None != (s & check)) {
            os << sep << name;
            sep = ", ";
        }
    };
    print_cause(State::Invalid, "Invalid Operation");
    print_cause(State::DivZero, "Divide by Zero");
    print_cause(State::Denorm, "Denormal Number");
    print_cause(State::Overflow, "Overflow");
    print_cause(State::Underflow, "Underflow");
    print_cause(State::Inexact, "Inexact");
    return os;
}
WatchedScope::WatchedScope(Mask m)
    : initial(getMask())
{
    setMask(Mask::DisableAll & ~m);
}

WatchedScope::~WatchedScope()
{
    setMask(initial);
}
}
} // namespace ZIRAN::FPE
