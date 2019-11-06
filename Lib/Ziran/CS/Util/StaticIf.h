#pragma once

#include <utility>
#include <type_traits>

namespace ZIRAN {
namespace STATIC_IF {
struct identity {
    template <typename T>
    T operator()(T&& x) const
    {
        return std::forward<T>(x);
    }
};

template <bool Cond>
struct statement {
    template <typename F>
    void then(const F& f)
    {
        f(identity());
    }

    template <typename F>
    void else_(const F&)
    {
    }
};

template <>
struct statement<false> {
    template <typename F>
    void then(const F&)
    {
    }

    template <typename F>
    void else_(const F& f)
    {
        f(identity());
    }
};

template <bool Cond, typename F>
inline statement<Cond> static_if(F const& f)
{
    statement<Cond> if_;
    if_.then(f);
    return if_;
}
} // namespace STATIC_IF

using STATIC_IF::static_if;
#define ZIRAN_STATIC_IF(x) static_if<(x)>([&,this](const auto& _____) -> void {
#define ZIRAN_STATIC_ELSE \
    }).else_([&,this](const auto &_____) -> void {
#define ZIRAN_STATIC_END_IF \
    });

// After we switch to C++17, we should use
// (Note the the behaviour of 'return' is still different.)

/*
#define ZIRAN_STATIC_IF(x) if constexpr(x) {
#define ZIRAN_STATIC_ELSE \
    } else {
#define ZIRAN_STATIC_END_IF \
    }
*/
} // namespace ZIRAN
