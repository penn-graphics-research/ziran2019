/**
Givens rotation
*/
#ifndef GIVENS_H
#define GIVENS_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Ziran/Math/MathTools.h>
#include <Ziran/CS/Util/PlatformSpecific.h>

namespace ZIRAN {

/**
    Class for givens rotation.
    Row rotation G*A corresponds to something like
    c -s  0
    ( s  c  0 ) A
    0  0  1
    Column rotation A G' corresponds to something like
    c -s  0
    A ( s  c  0 )
    0  0  1

    c and s are always computed so that
    ( c -s ) ( a )  =  ( * )
    s  c     b       ( 0 )

    Assume rowi<rowk.
    */
template <class T>
class GivensRotation {
public:
    int rowi;
    int rowk;
    T c;
    T s;

    ZIRAN_FORCE_INLINE GivensRotation(int rowi_in, int rowk_in)
        : rowi(rowi_in)
        , rowk(rowk_in)
        , c(1)
        , s(0)
    {
    }

    ZIRAN_FORCE_INLINE GivensRotation(T a, T b, int rowi_in, int rowk_in)
        : rowi(rowi_in)
        , rowk(rowk_in)
    {
        compute(a, b);
    }

    ~GivensRotation() {}

    ZIRAN_FORCE_INLINE void setIdentity()
    {
        c = 1;
        s = 0;
    }

    ZIRAN_FORCE_INLINE void transposeInPlace()
    {
        s = -s;
    }

    /**
        Compute c and s from a and b so that
        ( c -s ) ( a )  =  ( * )
        s  c     b       ( 0 )
        */
    ZIRAN_FORCE_INLINE void compute(const T a, const T b)
    {
        using std::sqrt;

        T d = a * a + b * b;
        c = 1;
        s = 0;
        T sqrtd = sqrt(d);
        //T t = MATH_TOOLS::rsqrt(d);
        if (sqrtd) {
            T t = 1 / sqrtd;
            c = a * t;
            s = -b * t;
        }
    }

    /**
        This function computes c and s so that
        ( c -s ) ( a )  =  ( 0 )
        s  c     b       ( * )
        */
    ZIRAN_FORCE_INLINE void computeUnconventional(const T a, const T b)
    {
        using std::sqrt;

        T d = a * a + b * b;
        c = 0;
        s = 1;
        T sqrtd = sqrt(d);
        //T t = MATH_TOOLS::rsqrt(d);
        if (sqrtd) {
            T t = 1 / sqrtd;
            s = a * t;
            c = b * t;
        }
    }

    /**
      Fill the R with the entries of this rotation
        */
    template <class MatrixType>
    ZIRAN_FORCE_INLINE void fill(const MatrixType& R) const
    {
        MatrixType& A = const_cast<MatrixType&>(R);
        A = MatrixType::Identity();
        A(rowi, rowi) = c;
        A(rowk, rowi) = -s;
        A(rowi, rowk) = s;
        A(rowk, rowk) = c;
    }

    /**
        This function does something like Q^T A -> A 
        [ c -s  0 ]
        [ s  c  0 ] A -> A
        [ 0  0  1 ]
        It only affects row i and row k of A.
        */
    template <class MatrixType>
    ZIRAN_FORCE_INLINE void rowRotation(MatrixType& A) const
    {
        for (int j = 0; j < A.cols(); j++) {
            T tau1 = A(rowi, j);
            T tau2 = A(rowk, j);
            A(rowi, j) = c * tau1 - s * tau2;
            A(rowk, j) = s * tau1 + c * tau2;
        }
        //not type safe :/
    }

    /**
        This function does something like A Q -> A
           [ c  s  0 ]
        A  [-s  c  0 ]  -> A
           [ 0  0  1 ]
        It only affects column i and column k of A.
        */
    template <class MatrixType>
    ZIRAN_FORCE_INLINE void columnRotation(MatrixType& A) const
    {
        for (int j = 0; j < A.rows(); j++) {
            T tau1 = A(j, rowi);
            T tau2 = A(j, rowk);
            A(j, rowi) = c * tau1 - s * tau2;
            A(j, rowk) = s * tau1 + c * tau2;
        }
        //not type safe :/
    }

    /**
      Multiply givens must be for same row and column
      **/
    ZIRAN_FORCE_INLINE void operator*=(const GivensRotation<T>& A)
    {
        T new_c = c * A.c - s * A.s;
        T new_s = s * A.c + c * A.s;
        c = new_c;
        s = new_s;
    }

    /**
      Multiply givens must be for same row and column
      **/
    ZIRAN_FORCE_INLINE GivensRotation<T> operator*(const GivensRotation<T>& A) const
    {
        GivensRotation<T> r(*this);
        r *= A;
        return r;
    }
};

/**
    \brief zero chasing the 3X3 matrix to bidiagonal form
    original form of H:
    x x 0
    x x x
    0 0 x
    after zero chase:
    x x 0
    0 x x
    0 0 x
    */
template <class T>
inline ZIRAN_FORCE_INLINE void zeroChase(Matrix<T, 3, 3>& H, Matrix<T, 3, 3>& U, Matrix<T, 3, 3>& V)
{

    /**
        Reduce H to of form
        x x +
        0 x x
        0 0 x
        */
    GivensRotation<T> r1(H(0, 0), H(1, 0), 0, 1);
    /**
        Reduce H to of form
        x x 0
        0 x x
        0 + x
        Can calculate r2 without multiplying by r1 since both entries are in first two
        rows thus no need to divide by sqrt(a^2+b^2)
        */
    GivensRotation<T> r2(1, 2);
    if (H(1, 0) != 0)
        r2.compute(H(0, 0) * H(0, 1) + H(1, 0) * H(1, 1), H(0, 0) * H(0, 2) + H(1, 0) * H(1, 2));
    else
        r2.compute(H(0, 1), H(0, 2));

    r1.rowRotation(H);

    /* GivensRotation<T> r2(H(0, 1), H(0, 2), 1, 2); */
    r2.columnRotation(H);
    r2.columnRotation(V);

    /**
        Reduce H to of form
        x x 0
        0 x x
        0 0 x
        */
    GivensRotation<T> r3(H(1, 1), H(2, 1), 1, 2);
    r3.rowRotation(H);

    // Save this till end for better cache coherency
    // r1.rowRotation(u_transpose);
    // r3.rowRotation(u_transpose);
    r1.columnRotation(U);
    r3.columnRotation(U);
}

/**
     \brief make a 3X3 matrix to upper bidiagonal form
     original form of H:   x x x
                           x x x
                           x x x
     after zero chase:
                           x x 0
                           0 x x
                           0 0 x
  */
template <class T>
inline ZIRAN_FORCE_INLINE void makeUpperBidiag(Matrix<T, 3, 3>& H, Matrix<T, 3, 3>& U, Matrix<T, 3, 3>& V)
{
    U = Matrix<T, 3, 3>::Identity();
    V = Matrix<T, 3, 3>::Identity();

    /**
      Reduce H to of form
                          x x x
                          x x x
                          0 x x
    */

    GivensRotation<T> r(H(1, 0), H(2, 0), 1, 2);
    r.rowRotation(H);
    // r.rowRotation(u_transpose);
    r.columnRotation(U);
    // zeroChase(H, u_transpose, V);
    zeroChase(H, U, V);
}

/**
     \brief make a 3X3 matrix to lambda shape
     original form of H:   x x x
     *                     x x x
     *                     x x x
     after :
     *                     x 0 0
     *                     x x 0
     *                     x 0 x
  */
template <class T>
inline ZIRAN_FORCE_INLINE void makeLambdaShape(Matrix<T, 3, 3>& H, Matrix<T, 3, 3>& U, Matrix<T, 3, 3>& V)
{
    U = Matrix<T, 3, 3>::Identity();
    V = Matrix<T, 3, 3>::Identity();

    /**
      Reduce H to of form
      *                    x x 0
      *                    x x x
      *                    x x x
      */

    GivensRotation<T> r1(H(0, 1), H(0, 2), 1, 2);
    r1.columnRotation(H);
    r1.columnRotation(V);

    /**
      Reduce H to of form
      *                    x x 0
      *                    x x 0
      *                    x x x
      */

    r1.computeUnconventional(H(1, 2), H(2, 2));
    r1.rowRotation(H);
    r1.columnRotation(U);

    /**
      Reduce H to of form
      *                    x x 0
      *                    x x 0
      *                    x 0 x
      */

    GivensRotation<T> r2(H(2, 0), H(2, 1), 0, 1);
    r2.columnRotation(H);
    r2.columnRotation(V);

    /**
      Reduce H to of form
      *                    x 0 0
      *                    x x 0
      *                    x 0 x
      */
    r2.computeUnconventional(H(0, 1), H(1, 1));
    r2.rowRotation(H);
    r2.columnRotation(U);
}
} // namespace ZIRAN

#endif
