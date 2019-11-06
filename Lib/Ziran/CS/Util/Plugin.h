#ifndef PLUGIN_H
#define PLUGIN_H
#include <Ziran/CS/Util/PlatformSpecific.h>
namespace ZIRAN {

class PluginManager;

// Define the API version.
// This value is incremented whenever there are ABI breaking changes.
// i.e. the plugin base class changes
#define ZIRAN_PLUGIN_API_VERSION 2

class ZIRAN_EXPORT PluginBase {
public:
    PluginBase() {}
    virtual ~PluginBase() {}
    virtual void registerFactories(PluginManager& manager) = 0;
};

// Plugin details structure that's exposed to the application.
struct PluginDetails {
    int apiVersion;
    const char* fileName;
    const char* className;
    const char* pluginVersion;
    PluginBase* (*initializeFunc)();
};

#define ZIRAN_PLUGIN(classType, pluginVersion)       \
    extern "C" {                                     \
    ZIRAN_EXPORT ZIRAN::PluginBase* get##classType() \
    {                                                \
        static classType singleton;                  \
        return &singleton;                           \
    }                                                \
    ZIRAN_EXPORT ZIRAN::PluginDetails exports = {    \
        ZIRAN_PLUGIN_API_VERSION,                    \
        __FILE__,                                    \
        #classType,                                  \
        pluginVersion,                               \
        get##classType                               \
    };                                               \
    }
} // namespace ZIRAN
#endif
