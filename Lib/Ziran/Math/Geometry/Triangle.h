#ifndef TRIANGLE_H
#define TRIANGLE_H
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {
template <class T, int dim>
void barycentricWeights(
    const Vector<T, dim>& p,
    const Vector<T, dim>& a,
    const Vector<T, dim>& b,
    const Vector<T, dim>& c,
    Vector<T, 3>& weights);

/**
  Computes the barycentric weights of the closest point on the triangle abc to p

*/
template <class T, int dim>
void interiorBarycentricWeights(
    const Vector<T, dim>& p,
    const Vector<T, dim>& a,
    const Vector<T, dim>& b,
    const Vector<T, dim>& c,
    Vector<T, 3>& weights);
} // namespace ZIRAN

#endif
