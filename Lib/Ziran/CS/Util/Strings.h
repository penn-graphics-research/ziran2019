#ifndef STRINGS_H
#define STRINGS_H
#include <Ziran/CS/Util/Meta.h>
#include <memory>
#include <string>
#include <type_traits>

namespace ZIRAN {

template <typename T>
inline std::enable_if_t<std::is_same<T, std::string>::value, std::string>
fromStr(const std::string& argv)
{
    return argv;
}

template <typename T>
inline std::enable_if_t<std::is_same<T, int>::value, int>
fromStr(const std::string& argv)
{
    return std::stoi(argv);
}

template <typename T>
inline std::enable_if_t<std::is_same<T, long>::value, long>
fromStr(const std::string& argv)
{
    return std::stol(argv);
}

template <typename T>
inline std::enable_if_t<std::is_same<T, long long>::value, long long>
fromStr(const std::string& argv)
{
    return std::stoll(argv);
}

template <typename T>
inline std::enable_if_t<std::is_same<T, unsigned long>::value, unsigned long>
fromStr(const std::string& argv)
{
    return std::stoul(argv);
}

template <typename T>
inline std::enable_if_t<std::is_same<T, unsigned long long>::value, unsigned long long>
fromStr(const std::string& argv)
{
    return std::stoull(argv);
}

template <typename T>
inline std::enable_if_t<std::is_same<T, float>::value, float>
fromStr(const std::string& argv)
{
    return std::stof(argv);
}

template <typename T>
inline std::enable_if_t<std::is_same<T, double>::value, double>
fromStr(const std::string& argv)
{
    return std::stod(argv);
}

template <typename TP>
inline std::enable_if_t<IsUniquePtr<TP>::value, TP>
fromStr(const std::string& arg)
{
    using T = typename TP::element_type;
    return std::make_unique<T>(fromStr<T>(arg));
}

template <typename TP>
inline std::enable_if_t<IsSharedPtr<TP>::value, TP>
fromStr(const std::string& arg)
{
    using T = typename TP::element_type;
    return std::make_shared<T>(fromStr<T>(arg));
}
} // namespace ZIRAN
#endif
