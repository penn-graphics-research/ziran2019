#ifndef USC_HAIR_DATA_IO_H
#define USC_HAIR_DATA_IO_H
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {

template <class T, int dim>
void readUscHairData(std::istream& in, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments, const T percentage_to_keep);

template <class T, int dim>
void readUscHairData(const std::string& filename, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments, const T percentage_to_keep = 1);
} // namespace ZIRAN
#endif
