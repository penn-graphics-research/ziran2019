#include "Tetrahedron.h"
#include <Ziran/CS/Util/Debug.h>
#include <Eigen/Dense>

namespace ZIRAN {

template <class T>
void barycentricWeights(
    const Vector<T, 3>& p,
    const Vector<T, 3>& a,
    const Vector<T, 3>& b,
    const Vector<T, 3>& c,
    const Vector<T, 3>& d,
    Vector<T, 4>& weights)
{
    Matrix<T, 3, 3> Dm;
    Dm.col(0) = a - d;
    Dm.col(1) = b - d;
    Dm.col(2) = c - d;

    Vector<T, 3> weights_for_abc = Dm.partialPivLu().solve(p - d);
    weights.template head<3>() = weights_for_abc;
    weights(3) = (T)1 - weights_for_abc.array().sum();
    ZIRAN_ASSERT((weights.array() >= (T)0).all(), weights.transpose());
    ZIRAN_ASSERT((weights.array() <= (T)1).all(), weights.transpose());
    ZIRAN_ASSERT(std::abs(weights.array().sum() - (T)1) < 1e-7, weights.transpose());

    Vector<T, 3> test_p = Vector<T, 3>::Constant((T)0);
    test_p += weights(0) * a;
    test_p += weights(1) * b;
    test_p += weights(2) * c;
    test_p += weights(3) * d;
    ZIRAN_ASSERT((p - test_p).squaredNorm() < 1e-7);
}

template void barycentricWeights<double>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 4, 1, 0, 4, 1>&);
template void barycentricWeights<float>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1>&);
} // namespace ZIRAN
