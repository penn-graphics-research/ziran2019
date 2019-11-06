#ifndef BINARY_CURVE_COLLECTION_IO_H

#define BINARY_CURVE_COLLECTION_IO_H

#include <string>
#include <Ziran/CS/Util/Forward.h>
namespace ZIRAN {

template <class T, int dim>
void readBinaryCurveCollectionFile(const std::string& filename, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments);
}

#endif
