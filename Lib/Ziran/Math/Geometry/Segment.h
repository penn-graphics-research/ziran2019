#ifndef SEGMENT_H
#define SEGMENT_H
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {
template <class T, int dim>
void barycentricWeights(
    const Vector<T, dim>& p,
    const Vector<T, dim>& a,
    const Vector<T, dim>& b,
    Vector<T, 2>& weights);

/**
  Computes the barycentric weights of the closest point on the segment ab to p

*/
template <class T, int dim>
void interiorBarycentricWeights(
    const Vector<T, dim>& p,
    const Vector<T, dim>& a,
    const Vector<T, dim>& b,
    Vector<T, 2>& weights);
} // namespace ZIRAN
#endif
