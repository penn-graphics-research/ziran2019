#include <Ziran/CS/Util/SharedLibrary.h>
#include <dlfcn.h>
#include <Ziran/CS/Util/Debug.h>
namespace ZIRAN {
SharedLibrary::SharedLibrary(const char* path)
{
    m_shared_library = dlopen(path, RTLD_GLOBAL | RTLD_NOW);
    ZIRAN_ASSERT(m_shared_library != nullptr, "Could not load ", path, "\n error ", std::string(dlerror()));
}

SharedLibrary::~SharedLibrary()
{
    if (m_shared_library != nullptr)
        dlclose(m_shared_library);
}

void SharedLibrary::sym(const char* name, void** symbol)
{
    *symbol = dlsym(m_shared_library, name);
    char* error = dlerror();
    ZIRAN_ASSERT(error == nullptr, "Could not load symbol\n error ", std::string(error));
}
} // namespace ZIRAN
