#ifndef TETGEN_IO
#define TETGEN_IO
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {
template <class T, int dim>
void writeNode(std::ostream& out, const StdVector<Vector<T, dim>>& X);
}
#endif
