/**
SVD based on Eigen's built-in Jacobi SVD
*/
#ifndef EIGEN_BUILTIN_JACOBI_SVD_H
#define EIGEN_BUILTIN_JACOBI_SVD_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>

namespace ZIRAN {

template <class T>
void singularValueDecompositionEigenBuiltinJacobi(const Matrix<T, 3, 3>& A,
    Matrix<T, 3, 3>& U,
    Vector<T, 3>& sigma,
    Matrix<T, 3, 3>& V)
{
    Eigen::JacobiSVD<Matrix<T, 3, 3>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    sigma = svd.singularValues();
    U = svd.matrixU();
    V = svd.matrixV();
    if (U.determinant() < 0) {
        U.col(2) *= -1;
        sigma(2) *= -1;
    }
    if (V.determinant() < 0) {
        V.col(2) *= -1;
        sigma(2) *= -1;
    }
}
} // namespace ZIRAN
#endif
