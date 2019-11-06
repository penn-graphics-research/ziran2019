#ifndef GIVENS_QR_H
#define GIVENS_QR_H

#include <Ziran/Math/Linear/Givens.h>
#include <iostream>

namespace ZIRAN {

template <class T, int m, int n>
inline void inplaceGivensR(Matrix<T, m, n>& A)
{
    for (int j = 0; j < A.cols(); j++) {
        for (int i = A.rows() - 1; i > j; i--) {
            GivensRotation<T> r(A(i - 1, j), A(i, j), i - 1, i);
            r.rowRotation(A);
        }
    }
}

template <class T, int m, int n, int l>
inline void simultaneousGivensQR(Matrix<T, m, n>& A, Matrix<T, m, l>& M)
{
    for (int j = 0; j < A.cols(); j++) {
        for (int i = A.rows() - 1; i > j; i--) {
            GivensRotation<T> r(A(i - 1, j), A(i, j), i - 1, i);
            r.rowRotation(A);
            r.rowRotation(M);
        }
    }
    // A <- Q.transpose() * A = R
    // M <- Q.transpose() * M
    // thin or thick A, M not yet tested
}

template <class T, int m, int n>
inline void inplaceGivensQR(Matrix<T, m, n>& A, Matrix<T, m, m>& Q)
{
    Q.setIdentity();
    simultaneousGivensQR(A, Q);
    Q.transposeInPlace();
}

template <class T, int m, int n>
inline void GivensQR(const Matrix<T, m, n>& A, Matrix<T, m, m>& Q, Matrix<T, m, n>& R)
{
    R = A;
    inplaceGivensQR(R, Q);
}

template <class T, int m, int n>
inline void thinGivensQR(const Eigen::Matrix<T, m, n, 0, m, n>& A, Eigen::Matrix<T, m, n, 0, m, n>& Q0, Eigen::Matrix<T, n, n, 0, n, n>& R0)
{
    Matrix<T, m, m> Q(A.rows(), A.rows());
    Matrix<T, m, n> R = A;
    inplaceGivensQR(R, Q);
    Q0 = Q.topLeftCorner(A.rows(), A.cols());
    R0 = R.topLeftCorner(A.cols(), A.cols());
}
} // namespace ZIRAN
#endif
