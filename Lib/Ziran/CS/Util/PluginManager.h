#ifndef PLUGIN_MANAGER_H
#define PLUGIN_MANAGER_H
#include <unordered_map>
#include <typeindex>
#include <vector>
#include <cassert>
#include <Ziran/CS/Util/Factory.h>
#include <Ziran/CS/Util/Plugin.h>
#include <Ziran/CS/Util/SharedLibrary.h>
#include <Ziran/CS/Util/Debug.h>
namespace ZIRAN {

class ZIRAN_EXPORT PluginManager {
    using TypeMap = std::unordered_multimap<std::type_index, std::unique_ptr<FactoryBase>>;
    using TypeMapIter = typename TypeMap::const_iterator;
    std::vector<SharedLibrary> shared_libraries;
    std::vector<PluginBase*> plugins;
    std::vector<PluginDetails*> plugin_details;
    TypeMap data;

public:
    template <class Interface>
    void registerFactory(std::unique_ptr<IFactory<Interface>>&& factory)
    {
        data.emplace(typeid(Interface), std::move(factory));
    }

    int numPlugins() const
    {
        return plugins.size();
    }

    const PluginDetails& getPluginDetails(int i) const
    {
        assert(i < int(plugin_details.size()));
        return *plugin_details[i];
    }

    PluginBase* getPluginBase(int i) const
    {
        assert(i < int(plugins.size()));
        return plugins[i];
    }

    template <class Interface>
    class FactoryIter {
        TypeMapIter current;

    public:
        FactoryIter(const TypeMapIter& t)
            : current(t)
        {
        }

        operator IFactory<Interface>*()
        {
            auto p = dynamic_cast<IFactory<Interface>*>(current->second.get());
            ZIRAN_ASSERT(p != nullptr);
            return p;
        }

        template <class... Args>
        operator AFactory<Interface, Args...>*()
        {
            return dynamic_cast<AFactory<Interface, Args...>*>(current->second.get());
        }

        IFactory<Interface>& operator*()
        {
            auto p = dynamic_cast<IFactory<Interface>*>(current->second.get());
            ZIRAN_ASSERT(p != nullptr);
            return *p;
        }

        FactoryIter<Interface>& operator++()
        {
            ++current;
            return *this;
        }

        bool operator!=(const FactoryIter<Interface>& other) const
        {
            return current != other.current;
        }

        bool operator==(const FactoryIter<Interface>& other) const
        {
            return current == other.current;
        }
    };

    template <class Interface>
    class FactoryRange {
        std::pair<TypeMapIter, TypeMapIter> p;

    public:
        FactoryRange(const std::pair<TypeMapIter, TypeMapIter>& p)
            : p(p)
        {
        }

        FactoryIter<Interface> begin() const
        {
            return p.first;
        }

        FactoryIter<Interface> end() const
        {
            return p.second;
        }
    };

    void loadAllPlugins(const char* path = ZIRAN_PLUGIN_DIR);

    void loadPlugin(const char* path);

    template <class Interface>
    FactoryRange<Interface> getAll() const
    {
        return FactoryRange<Interface>(data.equal_range(typeid(Interface)));
    }
};
} // namespace ZIRAN
#endif /* ifndef PLUGIN_MANAGER_H */
