#ifndef PLATFORM_SPECIFIC_H
#define PLATFORM_SPECIFIC_H

// Microsoft Visual Studio

#if defined(_MSC_VER)

#define ZIRAN_FORCE_INLINE __forceinline

#define ZIRAN_BIG_CONSTANT(x) (x)

#define ZIRAN_EXPORT __declspec(dllexport)
#define ZIRAN_LOCAL

// Other compilers

#else // defined(_MSC_VER)

#define ZIRAN_FORCE_INLINE __attribute__((always_inline))

#define ZIRAN_BIG_CONSTANT(x) (x##LLU)

#if __GNUC__ >= 4

#define ZIRAN_EXPORT __attribute__((visibility("default")))

#define ZIRAN_LOCAL __attribute__((visibility("hidden")))

#endif

#endif // !defined(_MSC_VER)

#ifdef __has_cpp_attribute
#if __has_cpp_attribute(fallthrough)
#define ZIRAN_FALLTHROUGH [[fallthrough]]
#endif
#endif
#ifndef ZIRAN_FALLTHROUGH
#define ZIRAN_FALLTHROUGH
#endif

#endif // ifndef PLATFORM_SPECIFIC_H
