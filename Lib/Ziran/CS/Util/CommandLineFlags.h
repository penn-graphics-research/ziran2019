#ifndef COMMAND_LINE_FLAGS_H
#define COMMAND_LINE_FLAGS_H
#include <cassert>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include <Ziran/CS/Util/Strings.h>

namespace ZIRAN {
namespace FLAGS {

template <typename T>
inline int parse(int argc, char** argv, T& value)
{
    if (argc < 2)
        throw std::runtime_error(std::string("Not enough arguments to ") + argv[0]);
    value = fromStr<T>(argv[1]);
    return 2;
}

inline int parse(int argc, char** argv, bool& value)
{
    value = true;
    return 1;
}

/**
  Class for command line flags
*/
class Flag {
public:
    using Callback = std::function<int(int, char**)>;
    /// Name of the flag
    const char* name;
    /// Help message for flag
    const char* help;
    /** callback to parse the arguments of the flag
      passed in the number of arguments left
      and a pointer to argument list
      starting with the flag itself
      return number of arguments parsed
    */
    Callback callback;

    /// The callback constructor is the most general
    Flag(const char* name, const char* help, Callback&& callback)
        : name(name)
        , help(help)
        , callback(std::move(callback))
    {
    }

    Flag(const char* name, const char* help, const Callback& callback)
        : name(name)
        , help(help)
        , callback(callback)
    {
    }

    /// The std::vector constructor allows for repeated arguments to be stored
    template <class T, class Alloc>
    Flag(const char* name, const char* help, std::vector<T, Alloc>& vec)
        : name(name)
        , help(help)
    {

        callback = [&vec](int argc, char** argv) -> int {
            T value;
            int r = parse(argc, argv, value);
            vec.push_back(value);
            return r;
        };
    }

    // This constructor stores the parsed argument in the passed value
    template <class T>
    Flag(const char* name, const char* help, T& value)
        : name(name)
        , help(help)
    {
        callback = [&value](int argc, char** argv) -> int { return parse(argc, argv, value); };
    }

    /// This constructor is used for flags that don't take arguments
    /// it simply calls the passed callback when the named argument is passed
    Flag(const char* name, const char* help, std::function<void(void)>&& f)
        : name(name)
        , help(help)
    {
        callback = [stored_f = std::move(f)] // C++14 generalized lamba capture
            (int argc, char** argv)
            -> int {
            stored_f();
            return 1;
        };
    }

    /// This constructor parses a single argument of type T
    /// and then calls the callback on that value
    template <class T>
    Flag(const char* name, const char* help, std::function<void(T)>&& f)
        : name(name)
        , help(help)
    {
        callback =
            [stored_f = std::move(f)] // C++14 generalized lamba capture
            (int argc, char** argv)
            -> int {
            T value;
            int r = parse(argc, argv, value);
            stored_f(value);
            return r;
        };
    }
};

/**
  Functor to compare char* as strings
*/
struct StringCompare {
    bool operator()(const char* s, const char* t) const
    {
        return std::strcmp(s, t) < 0;
    }
};

/**
  Registry for commandline flags
*/
class Registry {
    /// Map from flag string to Flag
    using Callback = std::function<int(int, char**)>;
    std::map<const char*, std::shared_ptr<Flag>, StringCompare> flag_map;
    StdVector<Callback> non_flag_callbacks;

public:
    /// Returns reference to the global registry
    static Registry& global()
    {
        static Registry global_registry;
        return global_registry;
    }

    /// Add a flag to the registry
    void add(const Flag& flag)
    {
        bool inserted = flag_map.emplace(flag.name, std::make_shared<Flag>(flag)).second;
        assert(inserted);
        (void)inserted;
    }

    void add(Flag&& flag)
    {
        bool inserted = flag_map.emplace(flag.name, std::make_shared<Flag>(std::move(flag))).second;
        assert(inserted);
        (void)inserted;
    }

    /// Perfect forwarding create flag
    template <typename... Args>
    void add(Args&&... args)
    {
        add(Flag(std::forward<Args>(args)...));
    }

    template <class Func>
    void addNonFlagCallback(Func&& c)
    {
        non_flag_callbacks.emplace_back(std::forward<Func>(c));
    }

    /// Find a flag in the registry
    std::shared_ptr<Flag> find(const char* flag)
    {
        auto i = flag_map.find(flag);
        if (i == flag_map.end())
            return nullptr;
        else
            return i->second;
    }

    void parseFlags(int argc, char** argv, bool skip_unknown = false)
    {
        for (int a = 1; a < argc;) {
            auto i = find(argv[a]);
            if (i == nullptr) {
                // flag not found in registry
                int p = 0;
                for (size_t j = 0; j < non_flag_callbacks.size(); j++) {
                    p = non_flag_callbacks[j](argc - a, argv + a);
                    if (p > 0) {
                        a += p;
                        break;
                    }
                }
                if (p > 0 || skip_unknown)
                    continue;
                else
                    throw std::runtime_error(std::string("Unknown flag ") + argv[a]);
            }
            int p = i->callback(argc - a, argv + a);
            assert(p > 0);
            a += p;
        }
    }

    void printUsage(std::ostream& os)
    {
        for (const auto& kv : flag_map) {
            os << "  " << kv.first << "\n"
               << "      " << kv.second->help << std::endl;
        }
    }
};

inline void PrintUsage(std::ostream& os)
{
    Registry::global().printUsage(os);
}

inline void ParseFlags(int argc, char** argv, bool skip_unknown = false)
{
    Registry::global().parseFlags(argc, argv);
}

template <class Func>
inline void addNonFlagCallback(Func&& c)
{
    Registry::global().addNonFlagCallback(std::forward<Func>(c));
}

template <class T>
inline void addMandatoryArgument(T& value)
{
    bool set = false;
    addNonFlagCallback(
        [&value, set](int argc, char** argv) mutable -> int {
            if (set)
                return 0;
            set = true;
            value = fromStr<T>(argv[0]);
            return 1;
        });
}

/**
  This class is used to register flags in the global registry
  The reason it's a class and not a function is so that you can register flags
  at global initialization time (i.e. before main)
  */
class Register {
public:
    /// Perfect forwarding
    template <typename... Args>
    Register(Args&&... args) noexcept
    {
        Registry::global().add(std::forward<Args>(args)...);
    }
};
}
} // namespace ZIRAN::FLAGS
#endif
