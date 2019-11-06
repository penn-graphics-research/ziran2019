#ifndef FACTORY_H
#define FACTORY_H value
#include <memory>
#include <Ziran/CS/Util/PlatformSpecific.h>
namespace ZIRAN {
class ZIRAN_EXPORT FactoryBase {
public:
    virtual ~FactoryBase() {}
    virtual const std::type_info& createdTypeInfo() const = 0;
    virtual bool supported(const char* sim_name, bool use_double, int dimension) = 0;
};

template <class Interface>
class ZIRAN_EXPORT IFactory : public FactoryBase {
public:
    virtual ~IFactory() {}
};

template <class Interface, class... Args>
class ZIRAN_EXPORT AFactory : public virtual IFactory<Interface> {
public:
    virtual ~AFactory() {}

    virtual std::unique_ptr<Interface> create(Args... args) const = 0;
};

template <class Derived, class Interface, class... Args>
class ZIRAN_EXPORT Factory : public AFactory<Interface, Args...> {
public:
    virtual ~Factory() {}

    virtual std::unique_ptr<Interface> create(Args... args) const override
    {
        return std::make_unique<Derived>(args...);
    }

    const std::type_info& createdTypeInfo() const override
    {
        return typeid(Derived);
    }
};
} // namespace ZIRAN
#endif
