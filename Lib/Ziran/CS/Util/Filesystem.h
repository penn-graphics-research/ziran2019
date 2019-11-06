#ifndef FILESYSTEM_H
#define FILESYSTEM_H
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <system_error>
#ifdef _WIN32
#include <direct.h>
#endif

namespace ZIRAN {

namespace FILESYSTEM {
/**
  Creates a directory if it does not exist
*/
inline void createDir(const std::string& dir)
{
    if (dir.empty()) {
        throw(std::runtime_error("Can't create directory with no name"));
    }
    struct stat info;
    if (stat(dir.c_str(), &info) != 0) {
        int status;
#ifdef _WIN32
        status = _mkdir(dir.c_str());
#else
        status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
        if (status == 0)
            return;
        if (errno != EEXIST)
            throw(std::system_error(errno, std::system_category()));
    }
}
/**
  Creates a whole path of directories if any don't exist
  like mkdir -p
*/
inline void createPath(const std::string& dir)
{
    size_t pos = dir.find('/', 1);
    for (; pos != std::string::npos; pos = dir.find('/', pos + 1)) {
        createDir(dir.substr(0, pos));
    }
    createDir(dir);
}

/**
  Read file into an ostream
  reads whole file into memory so should
  only be used for small files
*/
inline void readFile(std::ostream& contents, const std::string& filename)
{
    std::ifstream in;
    std::ios_base::iostate exceptionMask = in.exceptions() | std::ios::failbit;
    in.exceptions(exceptionMask);
    in.open(filename, std::ios::in | std::ios::binary);
    contents << in.rdbuf();
    in.close();
}

inline std::string fileToString(std::istream& contents)
{
    std::string str;
    contents.seekg(0, std::ios::end);
    str.reserve(contents.tellg());
    contents.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(contents)),
        std::istreambuf_iterator<char>());
    return str;
}
}
} // namespace ZIRAN::FILESYSTEM
#endif
