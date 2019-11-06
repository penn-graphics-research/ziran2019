/**
   A collection of rountines for small matricies.
*/
#ifndef DENSE_EXT_H
#define DENSE_EXT_H
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/Meta.h>
#include <Ziran/CS/Util/Strings.h>
#include <Ziran/Math/MathTools.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <tick/requires.h>

namespace ZIRAN {

template <class TM, class Hessian, class TV1, class TV2>
static void dPdFContractTwoVectors(TM& result, const Hessian& A, const TV1& u, const TV2& v)
{
    result = TM::Zero();
    constexpr int rows = TM::RowsAtCompileTime;
    constexpr int cols = TM::ColsAtCompileTime;
    for (int delta = 0; delta < cols; delta++) {
        for (int gamma = 0; gamma < rows; gamma++) {
            result += A.template block<rows, cols>(rows * gamma, cols * delta) * u(gamma) * v(delta);
        }
    }
}

namespace EIGEN_EXT {

template <class T, int dim>
Vector<int, dim> vectorFloor(const Vector<T, dim>& x)
{
    Vector<int, dim> v;
    for (size_t d = 0; d < dim; d++)
        v(d) = std::floor(x(d));
    return v;
}

/**
   \brief frobenius norm
   \param[in] x matrix.

   return the Frobenius norm of x
*/
template <class Derived>
auto norm(const Eigen::MatrixBase<Derived>& x)
{
    return x.norm();
}

/**
   \brief frobenius norm
   \param[in] x scalar

   return the absolute value of x
*/
template <class T, TICK_REQUIRES(std::is_floating_point<T>::value)>
T norm(const T& x)
{
    using std::fabs;
    return fabs(x);
}

/**
   \brief squared frobenius norm
   \param[in] x matrix.

   return the square of the Frobenius norm of x
*/
template <class Derived>
auto squaredNorm(const Eigen::MatrixBase<Derived>& x)
{
    return x.squaredNorm();
}

/**
   \brief squared frobenius norm
   \param[in] x scalar

   return x*x
*/
template <class T, TICK_REQUIRES(std::is_floating_point<T>::value)>
T squaredNorm(const T& x)
{
    return x * x;
}

/**
   \brief approxEqual
   \param[in] x matrix.
   \param[in] y matrix.
   \param[in] eps tolerance

   return if x and y are relatively close in the Frobenius norm
*/
template <class T1, class T2, class T3>
bool approxEqual(const T1& x, const T2& y, T3 eps)
{
    return norm(x - y) < std::max((T3)1, squaredNorm(y)) * eps;
}

constexpr int vecRowsAtCompileTime(int m, int n)
{
    return (m == Eigen::Dynamic) ? -1 : ((n == Eigen::Dynamic) ? -1 : m * n);
}

/**
   \brief Returns a map into a matrix as a vector
   \param[in] x matrix
   \returns map into x as a vector by stacking its columns.

   The map shares the same memory as the input matrix.
   Only supports column major matrices.
*/
template <class T, int m, int n, int flags, TICK_REQUIRES((flags & Eigen::ColMajor) == Eigen::ColMajor)>
Eigen::Map<const Eigen::Matrix<T, vecRowsAtCompileTime(m, n), 1>>
vec(const Eigen::Matrix<T, m, n, flags, m, n>& x)
{
    using VectorXT = Eigen::Matrix<T, vecRowsAtCompileTime(m, n), 1>;
    Eigen::Map<const VectorXT> x_vec(x.data(), x.rows() * x.cols(), 1);
    return x_vec;
}

template <class T, int m, int n, int flags, TICK_REQUIRES((flags & Eigen::ColMajor) == Eigen::ColMajor)>
Eigen::Map<Eigen::Matrix<T, vecRowsAtCompileTime(m, n), 1>>
vec(Eigen::Matrix<T, m, n, flags, m, n>& x)
{
    using VectorXT = Eigen::Matrix<T, vecRowsAtCompileTime(m, n), 1>;
    Eigen::Map<VectorXT> x_vec(x.data(), x.rows() * x.cols(), 1);
    return x_vec;
}

/**
   \brief Contract derivative with differential
   \param[out] df matrix differential
   \param[in] dfdx Jacobian of f
   \param[in] dx matrix differential

   Computes  \f$\delta f_{ij} = \frac{\partial f_{ij}}{\partial x_{kl}} \delta x_{kl}\f$
   Assumes the storage order of dfdx is such that this contraction is dfdx*vec(dx) where
   where vec(dx) is dx as a vector in its storage order.
*/
template <class T, int k, int l, int m, int n, int kl, int mn, int flags, TICK_REQUIRES(kl == k * l && mn == m * n)>
void contract(Eigen::Matrix<T, k, l, flags, k, l>& df, const Eigen::Matrix<T, kl, mn, flags, kl, mn>& dfdx, const Eigen::Matrix<T, m, n, flags, m, n>& dx)
{
    vec(df).noalias() = dfdx * vec(dx);
}

template <class T, int k, int l, int m, int n, int options, int flags>
void contract(Eigen::Matrix<T, k, l, flags, k, l>& df, const Eigen::SparseMatrix<T, options>& dfdx, const Eigen::Matrix<T, m, n, flags, m, n>& dx)
{
    vec(df).noalias() = dfdx * vec(dx);
}

/**
   \brief Contract derivative with differential
   \param[out] df scalar differential
   \param[in] dfdx Gradient of f
   \param[in] dx matrix differential

   Computes  \f$\delta f = \frac{\partial f}{\partial x_{kl}} \delta x_{kl}\f$
*/
template <class OutT, class Derived, TICK_REQUIRES(std::is_arithmetic<OutT>::value)>
void contract(OutT& df,
    const Eigen::MatrixBase<Derived>& dfdx,
    const Eigen::MatrixBase<Derived>& dx)
{
    df = dfdx.cwiseProduct(dx).sum();
}

/**
   \brief Contract derivative with differential
   \param[out] df differential
   \param[in] dfdx Jacobian of f
   \param[in] dx scalar differential

   Computes  \f$\delta f_{ij} = \frac{\partial f_{ij}}{\partial x} \delta x\f$
*/
template <class InT, class OutT, class DerT, TICK_REQUIRES(std::is_arithmetic<InT>::value)>
void contract(OutT& df, const DerT& dfdx, const InT& dx)
{
    df = dfdx * dx;
}

/**
   \brief Contract derivative with differential
   \param[out] df scalar differential
   \param[in] dfdx Gradient of f
   \param[in] dx vector differential

   Computes  \f$\delta f = \frac{\partial f}{\partial x_k} \delta x_k\f$
*/
template <class InT, class OutT, class DerT, TICK_REQUIRES(std::is_arithmetic<OutT>::value&& InT::ColsAtCompileTime == 1)>
void contract(OutT& df, const DerT& dfdx, const InT& dx)
{
    df = dfdx.dot(dx);
}

/**
   \brief  firstInvariant of F'*F
   \param[in] F matrix.
   \return Frobenius norm squared of F
*/
template <class Mat>
ScalarType<Mat> firstInvariant(const Mat& F)
{
    return F.squaredNorm();
}

/**
   \brief 1X1 JF^(-T)
   \param[in] F matrix
   \param[out] A the cofactor matrix of F
*/
template <class TF, class TA, TICK_REQUIRES(isSize<TF>(1, 1) && isSize<TA>(1, 1))>
inline void cofactorMatrix(const TF& F, TA& A)
{
    A(0, 0) = F(0, 0);
}
/**
   \brief 2X2 JF^(-T)
   \param[in] F matrix
   \param[out] A the cofactor matrix of F
*/
template <class TF, class TA, TICK_REQUIRES(isSize<TF>(2, 2) && isSize<TA>(2, 2))>
inline void cofactorMatrix(const TF& F, TA& A)
{
    A(0, 0) = F(1, 1);
    A(1, 0) = -F(0, 1);
    A(0, 1) = -F(1, 0);
    A(1, 1) = F(0, 0);
}
/**
   \brief 3X3 JF^(-T)
   \param[in] F input matrix
   \param[out] A the cofactor matrix of F
*/
template <class TF, class TA, TICK_REQUIRES(isSize<TF>(3, 3) && isSize<TA>(3, 3))>
inline void cofactorMatrix(const TF& F, TA& A)
{
    A(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
    A(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
    A(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
    A(1, 0) = F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2);
    A(1, 1) = F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0);
    A(1, 2) = F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1);
    A(2, 0) = F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1);
    A(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
    A(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
}
/**
   \brief take tensor product of two matrices A and B, store the result in M
   M_{ijkl} += scale*A_{ik}B_{jl}
*/
template <class T, int dim>
inline void addScaledTensorProduct_ik_jl(const Matrix<T, dim, dim>& A, const Matrix<T, dim, dim>& B, const T& scale, Eigen::Matrix<T, dim * dim, dim * dim>& M)
{
    for (int j2 = 0; j2 < dim; j2++)
        for (int i2 = 0; i2 < dim; i2++)
            for (int j1 = 0; j1 < dim; j1++)
                for (int i1 = 0; i1 < dim; i1++) {
                    int index_1 = j1 * dim + i1;
                    int index_2 = j2 * dim + i2;
                    M(index_1, index_2) += scale * A(i1, i2) * B(j1, j2);
                }
}
template <class T, int dim>
inline void tensorProduct_ik_jl(const Matrix<T, dim, dim>& A, const Matrix<T, dim, dim>& B, Eigen::Matrix<T, dim * dim, dim * dim>& M)
{
    for (int j2 = 0; j2 < dim; j2++)
        for (int i2 = 0; i2 < dim; i2++)
            for (int j1 = 0; j1 < dim; j1++)
                for (int i1 = 0; i1 < dim; i1++) {
                    int index_1 = j1 * dim + i1;
                    int index_2 = j2 * dim + i2;
                    M(index_1, index_2) = A(i1, i2) * B(j1, j2);
                }
}

/**
   \brief take tensor product of two matrices A and B, store the result in M
   M_{ijkl} += scale*A_{il}B_{kj}
*/
template <class T, int dim>
inline void addScaledTensorProduct_il_kj(const Matrix<T, dim, dim>& A, const Matrix<T, dim, dim>& B, const T& scale, Eigen::Matrix<T, dim * dim, dim * dim>& M)
{
    for (int j2 = 0; j2 < dim; j2++)
        for (int i2 = 0; i2 < dim; i2++)
            for (int j1 = 0; j1 < dim; j1++)
                for (int i1 = 0; i1 < dim; i1++) {
                    int index_1 = j1 * dim + i1;
                    int index_2 = j2 * dim + i2;
                    M(index_1, index_2) += scale * A(i1, j2) * B(i2, j1);
                }
}
template <class T, int dim>
inline void tensorProduct_il_kj(const Matrix<T, dim, dim>& A, const Matrix<T, dim, dim>& B, Eigen::Matrix<T, dim * dim, dim * dim>& M)
{
    for (int j2 = 0; j2 < dim; j2++)
        for (int i2 = 0; i2 < dim; i2++)
            for (int j1 = 0; j1 < dim; j1++)
                for (int i1 = 0; i1 < dim; i1++) {
                    int index_1 = j1 * dim + i1;
                    int index_2 = j2 * dim + i2;
                    M(index_1, index_2) = A(i1, j2) * B(i2, j1);
                }
}

template <class T>
void rotationalDifferential(const Vector<T, 1>& R, const Vector<T, 1>& S, const Vector<T, 1>& dF, Vector<T, 1>& dR)
{
    dR(0, 0) = (T)0;
}
template <class T>
void addScaledRotationalDifferential(const Vector<T, 1>& R, const Vector<T, 1>& S, const Vector<T, 1>& dF, T scale, Vector<T, 1>& M)
{
}
template <class T>
void rotationalDifferential(const Matrix<T, 2, 2>& R, const Matrix<T, 2, 2>& S, const Matrix<T, 2, 2>& dF, Matrix<T, 2, 2>& dR)
{
    //We dont need the result of R^T dF after the omega stage, so we use dR to store the temporary matrix
    dR.noalias() = R.transpose() * dF;
    T trace_s = S.trace();
    ZIRAN_ASSERT(trace_s != 0, "dR computation encountered division by zero");
    T omega = (dR(0, 1) - dR(1, 0)) / trace_s;
    dR(0, 0) = -omega * R(0, 1);
    dR(1, 0) = -omega * R(1, 1);
    dR(0, 1) = omega * R(0, 0);
    dR(1, 1) = omega * R(1, 0);
}
template <class T>
void addScaledRotationalDifferential(const Matrix<T, 2, 2>& R, const Matrix<T, 2, 2>& S, const Matrix<T, 2, 2>& dF, T scale, Matrix<T, 2, 2>& M)
{
    T trace_s = S.trace();
    ZIRAN_ASSERT(trace_s != 0, "dR computation encountered division by zero");
    //Majorly exploited but basically derived from "rotationalDifferential'
    T omega = (R(0, 0) * dF(0, 1) + R(1, 0) * dF(1, 1) - R(0, 1) * dF(0, 0) - R(1, 1) * dF(1, 0)) / trace_s * scale;
    M(0, 0) -= omega * R(0, 1);
    M(1, 0) -= omega * R(1, 1);
    M(0, 1) += omega * R(0, 0);
    M(1, 1) += omega * R(1, 0);
}

/**
   \brief add scaled 1X1 rotational Derivative (dRdF)
   \param[in] R Rotation matrix of F
   \param[in] S Symmetric matrix of F
   \param[in] scale The scale factor
   \param[out] M The matrix to add scale * dRdF to
*/
template <class T>
void addScaledRotationalDerivative(const Vector<T, 1>& R, const Vector<T, 1>& S, T scale, Vector<T, 1>& M)
{
}
/**
   \brief 2X2 rotational Derivative (dRdF)
   \param[in] R Rotation matrix of F
   \param[in] S Symmetric matrix of F
   \param[out] dRdF the derivative of the rotation matrix.
*/
template <class T>
void rotationalDerivative(const Matrix<T, 2, 2>& R, const Matrix<T, 2, 2>& S, Matrix<T, 4, 4>& dRdF)
{
    T trace_s = S.trace();
    ZIRAN_ASSERT(trace_s != 0, "dRdF computation encountered division by zero");
    T one_over_trace = (T)1 / trace_s;
    Matrix<T, 2, 2> RE;
    RE << -R(0, 1), R(0, 0), -R(1, 1), R(1, 0);
    Vector<T, 4> vec_RE = vec(RE);
    dRdF.noalias() = one_over_trace * (vec_RE * (vec_RE.transpose()));
}

/**
   \brief add scaled 2X2 rotational Derivative (dRdF)
   \param[in] R Rotation matrix of F
   \param[in] S Symmetric matrix of F
   \param[in] scale The scale factor
   \param[out] M The matrix to add scale * dRdF to
*/
template <class T>
void addScaledRotationalDerivative(const Matrix<T, 2, 2>& R, const Matrix<T, 2, 2>& S, T scale, Matrix<T, 4, 4>& M)
{
    T trace_s = S.trace();
    ZIRAN_ASSERT(trace_s != 0, "dRdF computation encountered division by zero");
    T scale_over_trace = scale / trace_s;
    Matrix<T, 2, 2> RE;
    RE << -R(0, 1), R(0, 0), -R(1, 1), R(1, 0);
    Vector<T, 4> vec_RE = vec(RE);
    M.noalias() += scale_over_trace * (vec_RE * (vec_RE.transpose()));
}

/**
   \brief 3x3 rotational Differential (dR)
   \param[in] R Rotation matrix of F
   \param[in] S Symmetric matrix of F
   \param[in] dF amount of change in F
   \param[out] dR the differential that is returned by the function

   See the [notes](../Notes/derivative_of_R.pdf) for details.
*/
/**
   \brief 3x3 rotational Derivative (dRdF)
   \param[in] R Rotation matrix of F
   \param[in] S Symmetric matrix of F
   \param[out] dRdF the derivative of the rotation matrix.
*/
template <class T>
void rotationalDerivative(const Matrix<T, 3, 3>& R, const Matrix<T, 3, 3>& S, Matrix<T, 9, 9>& dRdF)
{
    using TM = Matrix<T, 3, 3>;
    TM S_hat = -S;
    S_hat.diagonal().array() += S.trace();
    T b = S_hat.determinant();
    ZIRAN_ASSERT(b != 0, "dRdF computation encountered division by zero");
    T a = 1 / b;
    TM RS_hat = R * S_hat;
    TM aRS_hat = a * RS_hat;
    TM aRS_hatRT = aRS_hat * R.transpose();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int row_index = 3 * j + i;
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    int column_index = 3 * l + k;
                    dRdF(row_index, column_index) = aRS_hatRT(i, k) * S_hat(j, l) - aRS_hat(i, l) * RS_hat(k, j);
                }
            }
        }
    }
}
template <class T>
void rotationalDifferential(const Matrix<T, 3, 3>& R, const Matrix<T, 3, 3>& S, const Matrix<T, 3, 3>& dF, Matrix<T, 3, 3>& dR)
{
    using TM = Matrix<T, 3, 3>;

    TM S_hat = -S;
    S_hat.diagonal().array() += S.trace();
    T b = S_hat.determinant();
    ZIRAN_ASSERT(b != 0, "dRdF computation encountered division by zero");

    TM A = R.transpose() * dF;
    dR.noalias() = (1 / b) * R * S_hat * (A - A.transpose()) * S_hat;
}
template <class T>
void addScaledRotationalDifferential(const Matrix<T, 3, 3>& R, const Matrix<T, 3, 3>& S, const Matrix<T, 3, 3>& dF, T scale, Matrix<T, 3, 3>& M)
{
    using TM = Matrix<T, 3, 3>;

    TM S_hat = -S;
    S_hat.diagonal().array() += S.trace();
    T b = S_hat.determinant();
    ZIRAN_ASSERT(b != 0, "dRdF computation encountered division by zero");
    TM A = R.transpose() * dF;
    M.noalias() += (scale / b) * R * S_hat * (A - A.transpose()) * S_hat;
}
/**
   \brief add scaled 3X3 rotational Derivative (dRdF)
   \param[in] R Rotation matrix of F
   \param[in] S Symmetric matrix of F
   \param[in] scale The scale factor
   \param[in,out] M The matrix to add scale * dRdF to
*/
template <class T>
void addScaledRotationalDerivative(const Matrix<T, 3, 3>& R, const Matrix<T, 3, 3>& S, T scale, Matrix<T, 9, 9>& M)
{
    using TM = Matrix<T, 3, 3>;
    TM S_hat = -S;
    S_hat.diagonal().array() += S.trace();
    T b = S_hat.determinant();
    ZIRAN_ASSERT(b != 0, "dRdF computation encountered division by zero");
    T a = scale / b;
    TM RS_hat = R * S_hat;
    TM aRS_hat = a * RS_hat;
    TM aRS_hatRT = aRS_hat * R.transpose();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int row_index = 3 * j + i;
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    int column_index = 3 * l + k;
                    M(row_index, column_index) += aRS_hatRT(i, k) * S_hat(j, l) - aRS_hat(i, l) * RS_hat(k, j);
                }
            }
        }
    }
}

template <class T>
void QRDifferential(const Vector<T, 1>& Q, const Vector<T, 1>& R, const Vector<T, 1>& dF, Vector<T, 1>& dQ, Vector<T, 1>& dR)
{
    dQ(0) = (T)0;
    dR(0) = dF(0);
}

template <class T>
void QRDifferential(const Matrix<T, 2, 2>& Q, const Matrix<T, 2, 2>& R, const Matrix<T, 2, 2>& dF, Matrix<T, 2, 2>& dQ, Matrix<T, 2, 2>& dR)
{
    using TM = Matrix<T, 2, 2>;
    ZIRAN_ASSERT(R(0, 0) != 0, "dQ computation encountered division by zero");
    TM QtdF = Q.transpose() * dF;
    T a = -QtdF(1, 0) / R(0, 0);
    TM QtdQ;
    QtdQ << 0, a, -a, 0;
    dQ = Q * QtdQ;
    dR = Q.transpose() * dF - QtdQ * R;
}

template <class T>
void QRDifferential(const Matrix<T, 3, 3>& Q, const Matrix<T, 3, 3>& R, const Matrix<T, 3, 3>& dF, Matrix<T, 3, 3>& dQ, Matrix<T, 3, 3>& dR)
{
    using TM = Matrix<T, 3, 3>;
    ZIRAN_ASSERT(R(0, 0) != 0 && R(1, 1) != 0, "dQ computation encountered division by zero");
    TM QtdF = Q.transpose() * dF;
    T w3 = QtdF(1, 0) / R(0, 0);
    T w2 = -QtdF(2, 0) / R(0, 0);
    T w1 = (QtdF(2, 1) + w2 * R(0, 1)) / R(1, 1);
    TM QtdQ;
    QtdQ << 0, -w3, w2, w3, 0, -w1, -w2, w1, 0;
    dQ = Q * QtdQ;
    dR = Q.transpose() * dF - QtdQ * R;
}

template <class T>
bool nearKink(const Matrix<T, 1, 1>& S, T tol)
{
    return false;
}

template <class T>
bool nearKink(const Matrix<T, 2, 2>& S, T tol)
{
    T trace_s = S.trace();
    return std::abs(trace_s) < tol;
}

template <class T>
bool nearKink(const Matrix<T, 3, 3>& S, T tol)
{
    using TM = Matrix<T, 3, 3>;
    TM S_hat = -S;
    S_hat.diagonal().array() += S.trace();
    T b = S_hat.determinant();
    return std::abs(b) < tol;
}

/**
   \brief 1X1 d( JF^(-T) )
   \param[in] F input matrix
   \param[in] dF d(F)
   \param[out] result d( JF^(-T) )
*/
template <class T>
void cofactorMatrixDifferential(const Vector<T, 1>& F, const Vector<T, 1>& dF, Vector<T, 1>& result)
{
    result = dF;
}
template <class T>
void addScaledCofactorMatrixDifferential(const Vector<T, 1>& F, const Vector<T, 1>& dF, T scale, Vector<T, 1>& M)
{
    M += scale * dF;
}

/**
   \brief 2X2 d( JF^(-T) )
   \param[in] F input matrix
   \param[in] dF d(F)
   \param[out] result d( JF^(-T) )
*/
template <class T>
void cofactorMatrixDifferential(const Matrix<T, 2, 2>& F, const Matrix<T, 2, 2>& dF, Matrix<T, 2, 2>& result)
{
    cofactorMatrix(dF, result);
}

template <class T>
void addScaledCofactorMatrixDifferential(const Matrix<T, 2, 2>& F, const Matrix<T, 2, 2>& dF, T scale, Matrix<T, 2, 2>& M)
{
    M(0, 0) += scale * dF(1, 1);
    M(1, 0) -= scale * dF(0, 1);
    M(0, 1) -= scale * dF(1, 0);
    M(1, 1) += scale * dF(0, 0);
}

/**
   \brief 1X1 add scaled d( JF^(-T) )/dF
   \param[in] F input matrix
   \param[in] T scale
   \param[in,out] result d( JF^(-T) )/dF
*/
template <class T>
void addScaledCofactorMatrixDerivative(const Vector<T, 1>& F, T scale, Vector<T, 1>& result)
{
    result(0, 0) += scale;
}

/**
   \brief nxn dense matrix contraction
   \param[in] F1 input matrix
   \param[in] F2 input matrix
   \param[in,out] result F1:F2
*/
template <class T, int n>
T contractMatrices(const Matrix<T, n, n>& F1, const Matrix<T, n, n>& F2)
{
    T result = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            result += F1(i, j) * F2(i, j);
    return result;
}

/**
   \brief 1X1 d( JF^(-T) )/dF
   \param[in] F input matrix
   \param[out] result d( JF^(-T) )/dF << dP11/dF11, dP11/dF21, dP11/dF12, dP11/dF22, dP21/dF11, ...
*/
template <class T>
void cofactorMatrixDerivative(const Vector<T, 1>& F, Vector<T, 1>& result)
{
    result(0, 0) = 1;
}

/**
   \brief 2X2 add scaled d( JF^(-T) )/dF
   \param[in] F input matrix
   \param[in] T scale
   \param[in,out] result d( JF^(-T) )/dF
*/
template <class T>
void addScaledCofactorMatrixDerivative(const Matrix<T, 2, 2>& F, T scale, Matrix<T, 4, 4>& result)
{
    result(3, 0) += scale;
    result(2, 1) -= scale;
    result(1, 2) -= scale;
    result(0, 3) += scale;
}

/**
   \brief 2X2 d( JF^(-T) )/dF
   \param[in] F input matrix
   \param[out] result d( JF^(-T) )/dF << dP11/dF11, dP11/dF21, dP11/dF12, dP11/dF22, dP21/dF11, ...
*/
template <class T>
void cofactorMatrixDerivative(const Matrix<T, 2, 2>& F, Matrix<T, 4, 4>& result)
{
    result(0, 0) = 0;
    result(1, 0) = 0;
    result(2, 0) = 0;
    result(3, 0) = 1;
    result(0, 1) = 0;
    result(1, 1) = 0;
    result(2, 1) = -1;
    result(3, 1) = 0;
    result(0, 2) = 0;
    result(1, 2) = -1;
    result(2, 2) = 0;
    result(3, 2) = 0;
    result(0, 3) = 1;
    result(1, 3) = 0;
    result(2, 3) = 0;
    result(3, 3) = 0;
}

/**
   \brief 3X3 d( JF^(-T) )
   \param[in] F input matrix
   \param[in] dF d(F)
   \param[out] result d( JF^(-T) )
*/
template <class T>
void cofactorMatrixDifferential(const Matrix<T, 3, 3>& F, const Matrix<T, 3, 3>& dF, Matrix<T, 3, 3>& result)
{
    result(0, 0) = dF(1, 1) * F(2, 2) + F(1, 1) * dF(2, 2) - dF(2, 1) * F(1, 2) - F(2, 1) * dF(1, 2);
    result(1, 0) = dF(2, 1) * F(0, 2) + F(2, 1) * dF(0, 2) - dF(0, 1) * F(2, 2) - F(0, 1) * dF(2, 2);
    result(2, 0) = dF(0, 1) * F(1, 2) + F(0, 1) * dF(1, 2) - dF(1, 1) * F(0, 2) - F(1, 1) * dF(0, 2);
    result(0, 1) = dF(2, 0) * F(1, 2) + F(2, 0) * dF(1, 2) - dF(1, 0) * F(2, 2) - F(1, 0) * dF(2, 2);
    result(1, 1) = dF(0, 0) * F(2, 2) + F(0, 0) * dF(2, 2) - dF(2, 0) * F(0, 2) - F(2, 0) * dF(0, 2);
    result(2, 1) = dF(1, 0) * F(0, 2) + F(1, 0) * dF(0, 2) - dF(0, 0) * F(1, 2) - F(0, 0) * dF(1, 2);
    result(0, 2) = dF(1, 0) * F(2, 1) + F(1, 0) * dF(2, 1) - dF(2, 0) * F(1, 1) - F(2, 0) * dF(1, 1);
    result(1, 2) = dF(2, 0) * F(0, 1) + F(2, 0) * dF(0, 1) - dF(0, 0) * F(2, 1) - F(0, 0) * dF(2, 1);
    result(2, 2) = dF(0, 0) * F(1, 1) + F(0, 0) * dF(1, 1) - dF(1, 0) * F(0, 1) - F(1, 0) * dF(0, 1);
}
template <class T>
void addScaledCofactorMatrixDifferential(const Matrix<T, 3, 3>& F, const Matrix<T, 3, 3>& dF, T scale, Matrix<T, 3, 3>& M)
{
    M(0, 0) += scale * (dF(1, 1) * F(2, 2) + F(1, 1) * dF(2, 2) - dF(2, 1) * F(1, 2) - F(2, 1) * dF(1, 2));
    M(1, 0) += scale * (dF(2, 1) * F(0, 2) + F(2, 1) * dF(0, 2) - dF(0, 1) * F(2, 2) - F(0, 1) * dF(2, 2));
    M(2, 0) += scale * (dF(0, 1) * F(1, 2) + F(0, 1) * dF(1, 2) - dF(1, 1) * F(0, 2) - F(1, 1) * dF(0, 2));
    M(0, 1) += scale * (dF(2, 0) * F(1, 2) + F(2, 0) * dF(1, 2) - dF(1, 0) * F(2, 2) - F(1, 0) * dF(2, 2));
    M(1, 1) += scale * (dF(0, 0) * F(2, 2) + F(0, 0) * dF(2, 2) - dF(2, 0) * F(0, 2) - F(2, 0) * dF(0, 2));
    M(2, 1) += scale * (dF(1, 0) * F(0, 2) + F(1, 0) * dF(0, 2) - dF(0, 0) * F(1, 2) - F(0, 0) * dF(1, 2));
    M(0, 2) += scale * (dF(1, 0) * F(2, 1) + F(1, 0) * dF(2, 1) - dF(2, 0) * F(1, 1) - F(2, 0) * dF(1, 1));
    M(1, 2) += scale * (dF(2, 0) * F(0, 1) + F(2, 0) * dF(0, 1) - dF(0, 0) * F(2, 1) - F(0, 0) * dF(2, 1));
    M(2, 2) += scale * (dF(0, 0) * F(1, 1) + F(0, 0) * dF(1, 1) - dF(1, 0) * F(0, 1) - F(1, 0) * dF(0, 1));
}

/**
   \brief 3X3 add scaled d( JF^(-T) )/dF
   \param[in] F input matrix
   \param[in, out] result  += scale * d( JF^(-T) )/dF
*/
template <class T>
void addScaledCofactorMatrixDerivative(const Matrix<T, 3, 3>& F, T scale, Matrix<T, 9, 9>& result)
{
    Matrix<T, 3, 3> A = scale * F;
    result(4, 0) += A(2, 2);
    result(5, 0) += -A(1, 2);
    result(7, 0) += -A(2, 1);
    result(8, 0) += A(1, 1);
    result(3, 1) += -A(2, 2);
    result(5, 1) += A(0, 2);
    result(6, 1) += A(2, 1);
    result(8, 1) += -A(0, 1);
    result(3, 2) += A(1, 2);
    result(4, 2) += -A(0, 2);
    result(6, 2) += -A(1, 1);
    result(7, 2) += A(0, 1);
    result(1, 3) += -A(2, 2);
    result(2, 3) += A(1, 2);
    result(7, 3) += A(2, 0);
    result(8, 3) += -A(1, 0);
    result(0, 4) += A(2, 2);
    result(2, 4) += -A(0, 2);
    result(6, 4) += -A(2, 0);
    result(8, 4) += A(0, 0);
    result(0, 5) += -A(1, 2);
    result(1, 5) += A(0, 2);
    result(6, 5) += A(1, 0);
    result(7, 5) += -A(0, 0);
    result(1, 6) += A(2, 1);
    result(2, 6) += -A(1, 1);
    result(4, 6) += -A(2, 0);
    result(5, 6) += A(1, 0);
    result(0, 7) += -A(2, 1);
    result(2, 7) += A(0, 1);
    result(3, 7) += A(2, 0);
    result(5, 7) += -A(0, 0);
    result(0, 8) += A(1, 1);
    result(1, 8) += -A(0, 1);
    result(3, 8) += -A(1, 0);
    result(4, 8) += A(0, 0);
}

/**
   \brief 3X3 d( JF^(-T) )/dF
   \param[in] F input matrix
   \param[out] result d( JF^(-T) )/dF << dP11/dF11, dP11/dF21, ..., dP11/dF33, dP21/dF11, ...
*/
template <class T>
void cofactorMatrixDerivative(const Matrix<T, 3, 3>& F, Matrix<T, 9, 9>& result)
{
    result.setZero();
    addScaledCofactorMatrixDerivative(F, (T)1, result);
}

// for fracture project
template <class T, int dim>
Vector<T, dim> deviatoric(const Vector<T, dim>& input)
{
    T sum = 0;
    for (int i = 0; i < dim; ++i)
        sum += input(i);
    sum /= (T)dim;

    Vector<T, dim> output;
    for (int i = 0; i < dim; ++i)
        output(i) = input(i) - sum;

    return output;
}

// end for fracture project

} // namespace EIGEN_EXT
namespace FLAGS {
template <class T, int m, int n>
inline int parse(int argc, char** argv, Matrix<T, m, n>& value)
{
    if (argc < m * n + 1)
        throw std::runtime_error(std::string("Not enough arguments to ") + argv[0]);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            value(i, j) = fromStr<T>(argv[1 + i * m + j]);
    return m * n + 1;
}
} // namespace FLAGS
} // namespace ZIRAN
#endif
