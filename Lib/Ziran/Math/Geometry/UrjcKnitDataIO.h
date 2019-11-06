#ifndef URJC_KNIT_DATA_IO_H
#define URJC_KNIT_DATA_IO_H
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {

template <class T, int dim>
void readUrjcKnitData(const std::string& filename, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments);
}
#endif
