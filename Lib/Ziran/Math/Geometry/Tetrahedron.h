#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {
template <class T>
void barycentricWeights(
    const Vector<T, 3>& p,
    const Vector<T, 3>& a,
    const Vector<T, 3>& b,
    const Vector<T, 3>& c,
    const Vector<T, 3>& d,
    Vector<T, 4>& weights);
}
#endif
