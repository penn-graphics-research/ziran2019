#ifndef SHARED_LIBRARY_H
#define SHARED_LIBRARY_H
namespace ZIRAN {
class SharedLibrary {
    void* m_shared_library;

public:
    SharedLibrary(const char* path);

    ~SharedLibrary();

    void sym(const char* name, void** symbol);
};
} // namespace ZIRAN
#endif
