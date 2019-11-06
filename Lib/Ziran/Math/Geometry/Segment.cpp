#include "Segment.h"

namespace ZIRAN {

template <class T, int dim>
void barycentricWeights(
    const Vector<T, dim>& p,
    const Vector<T, dim>& a,
    const Vector<T, dim>& b,
    Vector<T, 2>& weights)
{
    Vector<T, dim> v0, v1;
    v0 = b - a;
    v1 = p - a;
    T d00 = v0.dot(v0);
    T d01 = v0.dot(v1);
    T inv_d00 = d00 ? (1 / d00) : T(0.5);
    T w = d01 * inv_d00;
    weights << 1 - w, w;
}

/**
  Computes the barycentric weights of the closest point on the segment ab to p

*/

template <class T, int dim>
void interiorBarycentricWeights(
    const Vector<T, dim>& p,
    const Vector<T, dim>& a,
    const Vector<T, dim>& b,
    Vector<T, 2>& weights)
{
    assert(p == p);
    assert(a == a);
    assert(b == b);

    Vector<T, 2> interior;
    barycentricWeights<T, dim>(p, a, b, interior);
    if (interior[0] < 0)
        weights << 0, 1;
    else if (interior[1] < 0)
        weights << 1, 0;
    else
        weights = interior;
}

template void barycentricWeights<double, 3>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template void barycentricWeights<float, 3>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1>&);

template void barycentricWeights<double, 2>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template void barycentricWeights<float, 2>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1>&);
template void interiorBarycentricWeights<double, 3>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template void interiorBarycentricWeights<float, 3>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1>&);

template void interiorBarycentricWeights<double, 2>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template void interiorBarycentricWeights<float, 2>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1>&);
} // namespace ZIRAN
