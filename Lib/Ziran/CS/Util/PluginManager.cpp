#include <Ziran/CS/Util/Logging.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/PluginManager.h>
#include <dirent.h>
namespace ZIRAN {

void PluginManager::loadPlugin(const char* path)
{
    ZIRAN_INFO("Loading: ", path);
    shared_libraries.emplace_back(path);
    PluginDetails* info;
    shared_libraries.back().sym("exports", reinterpret_cast<void**>(&info));
    plugin_details.emplace_back(info);
    ZIRAN_INFO("Plugin Loaded: ");
    ZIRAN_INFO("\tAPI Version: ", info->apiVersion);
    ZIRAN_INFO("\tFile Name: ", info->fileName);
    ZIRAN_INFO("\tClass Name: ", info->className);
    ZIRAN_INFO("\tPlugin Version: ", info->pluginVersion);

    // API Version checking
    ZIRAN_ASSERT(info->apiVersion == ZIRAN_PLUGIN_API_VERSION,
        "Plugin ABI version mismatch. Expected ", ZIRAN_PLUGIN_API_VERSION, " got ", info->apiVersion);

    // Instantiate the plugin
    plugins.emplace_back(info->initializeFunc());
    plugins.back()->registerFactories(*this);
}

void PluginManager::loadAllPlugins(const char* path)
{
    ZIRAN_INFO("Loading all plugins in ", path);
    const char* extension = ".so";
    int len = strlen(extension);
    DIR* dirp = opendir(path);
    ZIRAN_ASSERT(dirp != nullptr, "Could not open directory ", path);
    while (dirent* dp = readdir(dirp)) {
        int file_len = strlen(dp->d_name);
        if (len > file_len)
            continue;
        if (strncmp(&dp->d_name[file_len - len], extension, len) == 0) {
            std::string s(path);
            s += '/';
            s += dp->d_name;
            loadPlugin(s.c_str());
        }
    }
    closedir(dirp);
}
} // namespace ZIRAN
