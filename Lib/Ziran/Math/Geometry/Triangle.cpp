#include "Triangle.h"

namespace ZIRAN {

template <class T, int dim>
void barycentricWeights(
    const Vector<T, dim>& p,
    const Vector<T, dim>& a,
    const Vector<T, dim>& b,
    const Vector<T, dim>& c,
    Vector<T, 3>& weights)
{
    Vector<T, dim> v0, v1, v2;
    v0 = b - a;
    v1 = c - a;
    v2 = p - a;
    T d00 = v0.dot(v0);
    T d01 = v0.dot(v1);
    T d11 = v1.dot(v1);
    T d20 = v2.dot(v0);
    T d21 = v2.dot(v1);
    T inv_denom = (T)1 / (d00 * d11 - d01 * d01);
    weights[1] = (d11 * d20 - d01 * d21) * inv_denom;
    weights[2] = (d00 * d21 - d01 * d20) * inv_denom;
    weights[0] = (T)1 - weights[1] - weights[2];
}

/**
  Computes the barycentric weights of the closest point on the triangle abc to p

*/
template <class T, int dim>
void interiorBarycentricWeights(
    const Vector<T, dim>& p,
    const Vector<T, dim>& a,
    const Vector<T, dim>& b,
    const Vector<T, dim>& c,
    Vector<T, 3>& weights)
{
    assert(p == p);
    assert(a == a);
    assert(b == b);
    assert(c == c);

    using TV = Vector<T, dim>;
    TV v0, v1, v2;
    v0 = b - a;
    v1 = c - a;
    v2 = p - a;

    T d20 = v2.dot(v0);
    T d21 = v2.dot(v1);
    if (d20 <= 0 && d21 <= 0) {
        weights << 1, 0, 0;
        return;
    }

    TV v3 = p - b;
    T d30 = v3.dot(v0);
    T d31 = v3.dot(v1);
    if (d30 >= 0 && d31 <= d30) {
        weights << 0, 1, 0;
        return;
    }

    T vC = d20 * d31 - d30 * d21;
    if (vC < 0 && d20 >= 0 && d30 <= 0) {
        T diff = (d20 - d30);
        T v = diff ? (d20 / diff) : 0;
        weights << (T)1 - v, v, 0;
        return;
    }

    TV v4 = p - c;
    T d40 = v4.dot(v0);
    T d41 = v4.dot(v1);
    if (d41 >= 0 && d40 <= d41) {
        weights << 0, 0, 1;
        return;
    }

    T vB = d40 * d21 - d20 * d41;
    if (vB < 0 && d21 >= 0 && d41 <= 0) {
        T diff = d21 - d41;
        T w = diff ? (d21 / diff) : 0;
        weights << 1 - w, 0, w;
        return;
    }

    T vA = d30 * d41 - d40 * d31;
    T d310 = d31 - d30;
    T d410 = d41 - d40;
    if (vA < 0 && d310 >= 0 && d410 <= 0) {
        T diff = d310 - d410;
        T w = diff ? (d310 / diff) : 0;
        weights << 0, 1 - w, w;
        return;
    }

    T total_volume = (vA + vB + vC);
    if (total_volume > 0) {
        T inv = 1 / total_volume;
        weights << vA * inv, vB * inv, vC * inv;
    }
    else {
        T diff0 = d20 - d30;
        T diff1 = d21 - d41;
        if (diff0 > diff1) {
            T v = diff0 ? (d20 / diff0) : 0;
            weights << (T)1 - v, v, 0;
        }
        else {
            T w = diff1 ? (d21 / diff1) : 0;
            weights << 1 - w, 0, w;
        }
    }
}
template void interiorBarycentricWeights<double, 3>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template void interiorBarycentricWeights<float, 3>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
template void barycentricWeights<double, 3>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template void barycentricWeights<float, 3>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1>&);

template void interiorBarycentricWeights<double, 2>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template void interiorBarycentricWeights<float, 2>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
template void barycentricWeights<double, 2>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template void barycentricWeights<float, 2>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
} // namespace ZIRAN
