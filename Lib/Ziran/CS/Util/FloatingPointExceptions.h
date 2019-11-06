#ifndef FLOATING_POINT_EXCEPTIONS_H
#define FLOATING_POINT_EXCEPTIONS_H
#include <Ziran/CS/Util/PlatformSpecific.h>
#include <type_traits>
#include <xmmintrin.h>
#include <iostream>

namespace ZIRAN {

namespace FPE {

enum class State : unsigned int {
    None = 0,
    Invalid = _MM_EXCEPT_INVALID,
    DivZero = _MM_EXCEPT_DIV_ZERO,
    Denorm = _MM_EXCEPT_DENORM,
    Overflow = _MM_EXCEPT_OVERFLOW,
    Underflow = _MM_EXCEPT_UNDERFLOW,
    Inexact = _MM_EXCEPT_INEXACT
};

inline State operator|(State a, State b)
{
    using T = std::underlying_type_t<State>;
    return static_cast<State>(static_cast<T>(a) | static_cast<T>(b));
}

inline State operator^(State a, State b)
{
    using T = std::underlying_type_t<State>;
    return static_cast<State>(static_cast<T>(a) ^ static_cast<T>(b));
}

inline State operator&(State a, State b)
{
    using T = std::underlying_type_t<State>;
    return static_cast<State>(static_cast<T>(a) & static_cast<T>(b));
}

inline State operator~(State a)
{
    using T = std::underlying_type_t<State>;
    return static_cast<State>(~static_cast<T>(a));
}

inline State getState()
{
    return static_cast<State>(_MM_GET_EXCEPTION_STATE());
}

inline void setState(State a)
{
    using T = std::underlying_type_t<State>;
    _MM_SET_EXCEPTION_STATE(static_cast<T>(a));
}

std::ostream& operator<<(std::ostream& os, State s);

/**
  A mask is used to turn off floating point exceptions
  triggering SIGFPE.  They can be combined with bitwise |
  */
enum class Mask : unsigned int {
    Invalid = _MM_MASK_INVALID,
    DivZero = _MM_MASK_DIV_ZERO,
    Denorm = _MM_MASK_DENORM,
    Overflow = _MM_MASK_OVERFLOW,
    Underflow = _MM_MASK_UNDERFLOW,
    Inexact = _MM_MASK_INEXACT,
    DisableAll = Invalid | DivZero | Denorm | Overflow | Underflow | Inexact
};

ZIRAN_FORCE_INLINE inline Mask operator|(Mask a, Mask b)
{
    using T = std::underlying_type_t<Mask>;
    return static_cast<Mask>(static_cast<T>(a) | static_cast<T>(b));
}

ZIRAN_FORCE_INLINE inline Mask operator^(Mask a, Mask b)
{
    using T = std::underlying_type_t<Mask>;
    return static_cast<Mask>(static_cast<T>(a) ^ static_cast<T>(b));
}

ZIRAN_FORCE_INLINE inline Mask operator&(Mask a, Mask b)
{
    using T = std::underlying_type_t<Mask>;
    return static_cast<Mask>(static_cast<T>(a) & static_cast<T>(b));
}

ZIRAN_FORCE_INLINE inline Mask operator~(Mask a)
{
    using T = std::underlying_type_t<Mask>;
    return static_cast<Mask>(~static_cast<T>(a));
}

ZIRAN_FORCE_INLINE inline Mask getMask()
{
    return static_cast<Mask>(_MM_GET_EXCEPTION_MASK());
}

/* Turn off exceptions for the masked bits */
ZIRAN_FORCE_INLINE inline void setMask(Mask a)
{
    using T = std::underlying_type_t<Mask>;
    _MM_SET_EXCEPTION_MASK(static_cast<T>(a));
}

/**
  During the lifetime of an instance a WatchedScope
  any floating point exceptions matching the 
  the bitflags passed in will result in a SIGFPE.
  (Unless of course another function changes the global
  signal mask within the scope)
  For example to enable SIGFPE for overflow, divide by zero
  and invalid operations you would use

    FPE::WatchedScope w(FPE::Mask::Overflow | FPE::Mask::Invalid | FPE::Mask::DivZero);

  */
class WatchedScope {
    //Store the
    Mask initial;

public:
    WatchedScope(Mask m);

    ~WatchedScope();
};

enum class Rounding {
    Nearest = _MM_ROUND_NEAREST,
    Down = _MM_ROUND_DOWN,
    Up = _MM_ROUND_UP,
    Zero = _MM_ROUND_TOWARD_ZERO
};

ZIRAN_FORCE_INLINE inline void setRounding(Rounding x)
{
    using T = std::underlying_type_t<Rounding>;
    _MM_SET_ROUNDING_MODE(static_cast<T>(x));
}

ZIRAN_FORCE_INLINE inline Rounding getRounding()
{
    return static_cast<Rounding>(_MM_GET_ROUNDING_MODE());
}

enum class DenormalToZero {
    Enable = _MM_FLUSH_ZERO_ON,
    Disable = _MM_FLUSH_ZERO_OFF
};

ZIRAN_FORCE_INLINE inline void setFlush(DenormalToZero x)
{
    using T = std::underlying_type_t<DenormalToZero>;
    _MM_SET_FLUSH_ZERO_MODE(static_cast<T>(x));
}

ZIRAN_FORCE_INLINE inline DenormalToZero getFlush()
{
    return static_cast<DenormalToZero>(_MM_GET_FLUSH_ZERO_MODE());
}
}
} // namespace ZIRAN::FPE
#endif
