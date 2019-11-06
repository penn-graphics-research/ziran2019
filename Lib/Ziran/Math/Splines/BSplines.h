#ifndef B_SPLINES_H
#define B_SPLINES_H
#include <Eigen/Core>
#include <Ziran/Math/MathTools.h>
#include <tick/requires.h>
#include <iostream>

namespace ZIRAN {

template <int interpolation_degree, class T, TICK_REQUIRES(interpolation_degree == 1)>
inline int baseNode(const T& x)
{
    return MATH_TOOLS::int_floor(x);
}

template <int interpolation_degree, class T, TICK_REQUIRES(interpolation_degree != 1)>
inline int baseNode(const T& x)
{
    return MATH_TOOLS::int_floor(x - (T)0.5 * (interpolation_degree - 1));
}

template <int interpolation_degree, class T, int dim>
inline Vector<int, dim> baseNode(const Vector<T, dim>& x)
{
    Vector<int, dim> base;
    for (int d = 0; d < dim; d++)
        base(d) = baseNode<interpolation_degree, T>(x(d));
    return base;
}

/**
\brief Compute linear interpolation weights
\param x the point in a h=1 grid
\param w the weights.
\param dw the weight gradients.
*/
template <class T>
inline void computeBSplineWeights(const T x, int& base_node, Vector<T, 2>& w, Vector<T, 2>* dw = 0, Vector<T, 2>* ddw = 0)
{
    base_node = baseNode<1>(x);
    T dx = x - base_node;
    w << 1 - dx, dx;
    if (dw)
        (*dw) << -1, 1;
    if (ddw)
        (*ddw) << 0, 0;
}

/**;p-
\brief Compute quadratic B-spline interpolation weights
\param x the point in a h=1 grid
\param w the weights.
\param dw the weight gradients.
*/
template <class T>
inline void computeBSplineWeights(const T x, int& base_node, Vector<T, 3>& w, Vector<T, 3>* dw = 0, Vector<T, 3>* ddw = 0)
{
    base_node = baseNode<2>(x);
    T d0 = x - base_node;
    T z = ((T)1.5 - d0);
    T z2 = z * z;
    w(0) = (T)0.5 * z2;
    T d1 = d0 - 1;
    w(1) = (T)0.75 - d1 * d1;
    T d2 = 1 - d1;
    T zz = (T)1.5 - d2;
    T zz2 = zz * zz;
    w(2) = (T)0.5 * zz2;

    if (dw) {
        (*dw)(0) = -z;
        (*dw)(1) = -(T)2 * d1;
        (*dw)(2) = zz;
    }

    if (ddw) {
        (*ddw)(0) = 1;
        (*ddw)(1) = -2;
        (*ddw)(2) = 1;
    }
}

/**
\brief Compute cubic B-spline interpolation weights
\param x the point in a h=1 grid
\param w the weights.
\param dw the weight gradients.
*/
template <class T>
inline void computeBSplineWeights(const T x, int& base_node, Vector<T, 4>& w, Vector<T, 4>* dw = 0, Vector<T, 4>* ddw = 0)
{
    base_node = baseNode<3>(x);
    T d0 = x - base_node;
    T z = 2 - d0;
    T z3 = z * z * z;
    w(0) = ((T)1 / (T)6) * z3;
    T d1 = d0 - 1;
    T zz2 = d1 * d1;
    w(1) = ((T)0.5 * d1 - 1) * zz2 + (T)2 / (T)3;
    T d2 = 1 - d1;
    T zzz2 = d2 * d2;
    w(2) = ((T)0.5 * d2 - 1) * zzz2 + (T)2 / (T)3;
    T d3 = 1 + d2;
    T zzzz = 2 - d3;
    T zzzz3 = zzzz * zzzz * zzzz;
    w(3) = ((T)1 / (T)6) * zzzz3;

    if (dw) {
        (*dw)(0) = -(T)0.5 * z * z;
        (*dw)(1) = ((T)1.5 * d1 - (T)2) * d1;
        (*dw)(2) = (-(T)1.5 * d2 + (T)2) * d2;
        (*dw)(3) = (T)0.5 * zzzz * zzzz;
    }

    if (ddw) {
        (*ddw)(0) = (T)2 + base_node - x;
        (*ddw)(1) = -(T)2 + (T)3 * (-(T)1 - base_node + x);
        (*ddw)(2) = -(T)2 + (T)3 * ((T)2 + base_node - x);
        (*ddw)(3) = -(T)1 - base_node + x;
    }
}
template <class T>
inline void computeBSplineWeightsBoundary(const T x, const int size, int& base_node, Vector<T, 2>& w, Vector<T, 2>* dw = 0, Vector<T, 2>* ddw = 0)
{
    base_node = baseNode<1>(x);
    // std::cout<<" X is "<<x<<" size is "<<size<<" base node is "<<base_node<<std::endl;
    if (base_node >= size - 3) {
        // std::cout << " X is " << x << " size is " << size << " base node is " << base_node << std::endl;
        w << 1, 0;
        if (dw)
            (*dw) << 1, 1;
        if (ddw)
            (*ddw) << 0, 0;
        return;
    }
    T dx = x - base_node;
    w << 1 - dx, dx;
    if (dw)
        (*dw) << -1, 1;
    if (ddw)
        (*ddw) << 0, 0;
}
template <class T>
inline void computeBSplineWeightsBoundary(const T x, const int size, int& base_node, Vector<T, 3>& w, Vector<T, 3>* dw = 0, Vector<T, 3>* ddw = 0)
{
    base_node = baseNode<2>(x);
    T d0 = x - base_node;
    if (base_node == -1) {
        T z = -(T)1.5 + d0;
        w(0) = 4 * z * z;
        w(1) = -(T)28 / 3 + (T)44 / 3 * d0 - (T)16 / 3 * d0 * d0;
        T zz = 1 - d0;
        w(2) = (T)4 / 3 * zz * zz;
        if (dw) {
            (*dw)(0) = (T)8 * z;
            (*dw)(1) = (T)44 / 3 - (T)32 / 3 * d0;
            (*dw)(2) = -(T)8 / 3 * zz;
        }
        if (ddw) {
            (*ddw)(0) = 8;
            (*ddw)(1) = -(T)32 / 3;
            (*ddw)(2) = (T)8 / 3;
        }
    }
    else if (base_node == 0) {
        T d02 = d0 * d0;
        w(0) = 1.5 - 2 * d0 + (T)2 / 3 * d02;
        w(1) = -(T).625 + 2.5 * d0 - (T)7 / 6 * d02;
        T z = (T).5 - d0;
        w(2) = (T).5 * z * z;
        ;
        if (dw) {
            (*dw)(0) = -(T)2 + (T)4 / 3 * d0;
            (*dw)(1) = (T)2.5 - (T)7 / 3 * d0;
            (*dw)(2) = -z;
        }
        if (ddw) {
            (*ddw)(0) = (T)4 / 3;
            (*ddw)(1) = -(T)7 / 3;
            (*ddw)(2) = 1;
        }
    }
    else if (base_node >= 1 && base_node <= size - 6) {
        T z = ((T)1.5 - d0);
        T z2 = z * z;
        w(0) = (T)0.5 * z2;
        T d1 = d0 - 1;
        w(1) = (T)0.75 - d1 * d1;
        T d2 = 1 - d1;
        T zz = (T)1.5 - d2;
        T zz2 = zz * zz;
        w(2) = (T)0.5 * zz2;

        if (dw) {
            (*dw)(0) = -z;
            (*dw)(1) = -(T)2 * d1;
            (*dw)(2) = zz;
        }

        if (ddw) {
            (*ddw)(0) = 1;
            (*ddw)(1) = -2;
            (*ddw)(2) = 1;
        }
    }
    else if (base_node == size - 4) {
        T z = -(T)1 + d0;
        T d02 = d0 * d0;
        w(0) = (T)4 / 3 * z * z;
        w(1) = -(T)4 / 3 + (T)20 / 3 * d0 - (T)16 / 3 * d02;
        T zz = -(T).5 + d0;
        w(2) = (T)4 * zz * zz;
        ;
        if (dw) {
            (*dw)(0) = (T)8 / 3 * z;
            (*dw)(1) = (T)20 / 3 - (T)32 / 3 * d0;
            (*dw)(2) = 8 * zz;
        }
        if (ddw) {
            (*ddw)(0) = (T)8 / 3;
            (*ddw)(1) = -(T)32 / 3;
            (*ddw)(2) = 8;
        }
    }
    else if (base_node == size - 5) {
        T z = (T)1.5 - d0;
        w(0) = (T).5 * z * z;
        w(1) = -(T)7 / 24 + (T)13 / 6 * d0 - (T)7 / 6 * d0 * d0;
        T zz = (T).5 - d0;
        w(2) = (T)2 / 3 * zz * zz;
        ;
        if (dw) {
            (*dw)(0) = -z;
            (*dw)(1) = (T)13 / 6 - (T)7 / 3 * d0;
            (*dw)(2) = -(T)4 / 3 * zz;
        }
        if (ddw) {
            (*ddw)(0) = 1;
            (*ddw)(1) = -(T)7 / 3;
            (*ddw)(2) = (T)4 / 3;
        }
    }
}
template <class T>
inline void computeBSplineWeightsBoundary(const T x, const int size, int& base_node, Vector<T, 4>& w, Vector<T, 4>* dw = 0, Vector<T, 4>* ddw = 0)
{

    T base_node_x = floor(x - (T)0.5 * (3 - 1));
    base_node = (int)base_node_x;
    T x2 = x * x;
    T x3 = x * x2;
    //  std::cout<<"base node is "<<base_node<<" size is "<<size<<std::endl<<std::flush;
    if (base_node >= 1 && base_node < size - 6) {
        // T d0 = x - base_node_x;
        // w(0) = -((T)1 / (T)6) * (d0-2)*(d0-2)*(d0-2);
        // w(1) = ((T)1/(T)6)*(-5+21*d0-15*d0*d0+3*d0*d0*d0);
        // w(2) = ((T)1/(T)6)*(4-12*d0+12*d0*d0-3*d0*d0*d0);
        // w(3) = ((T)1 / (T)6) * (-1+d0)* (-1+d0)* (-1+d0);
        //
        // if (dw) {
        //     (*dw)(0) = -(T)0.5 * (d0-2) * (d0-2);
        //     (*dw)(1) = ((T)1.5 *(1-d0) - (T)2) * (1-d0);
        //     (*dw)(2) = -(T)1.5 * d0 +4*d0 -2;
        //     (*dw)(3) = (T)0.5 * (d0-1)*(d0-1);
        // }
        //
        // if (ddw) {
        //     (*ddw)(0) = (T)2 -d0;
        //     (*ddw)(1) = -(T)2 + (T)3 * (-(T)1 +d0);
        //     (*ddw)(2) = -(T)2 + (T)3 * ((T)2 -d0);
        //     (*ddw)(3) = -(T)1 +d0;
        // }
        // }
        T d0 = x - base_node_x;
        T z = 2 - d0;
        T z3 = z * z * z;
        w(0) = ((T)1 / (T)6) * z3;
        T d1 = d0 - 1;
        T zz2 = d1 * d1;
        w(1) = ((T)0.5 * d1 - 1) * zz2 + (T)2 / (T)3;
        T d2 = 1 - d1;
        T zzz2 = d2 * d2;
        w(2) = ((T)0.5 * d2 - 1) * zzz2 + (T)2 / (T)3;
        T d3 = 1 + d2;
        T zzzz = 2 - d3;
        T zzzz3 = zzzz * zzzz * zzzz;
        w(3) = ((T)1 / (T)6) * zzzz3;

        if (dw) {
            (*dw)(0) = -(T)0.5 * z * z;
            (*dw)(1) = ((T)1.5 * d1 - (T)2) * d1;
            (*dw)(2) = (-(T)1.5 * d2 + (T)2) * d2;
            (*dw)(3) = (T)0.5 * zzzz * zzzz;
        }

        if (ddw) {

            (*ddw)(0) = (T)2 + base_node_x - x;
            (*ddw)(1) = -(T)2 + (T)3 * (-(T)1 - base_node_x + x);
            (*ddw)(2) = -(T)2 + (T)3 * ((T)2 + base_node_x - x);
            (*ddw)(3) = -(T)1 - base_node_x + x;
        }
    }
    else if (base_node == 0) {
        T z = 2 - x;
        T z3 = z * z * z;
        w(0) = (T)0.25 * z3;
        w(1) = -(T)1.5 + ((T)1 / (T)12) * x * ((T)54 + x * (-(T)36 + (T)7 * x));
        w(2) = (T)2 / (T)3 - (T)0.5 * z * z * x;
        T zzzz = (-1 + x);
        T zzzz3 = zzzz * zzzz * zzzz;
        w(3) = (T)1 / (T)6 * zzzz3;
        if (dw) {
            (*dw)(0) = -(T)0.75 * z * z;
            (*dw)(1) = (T)4.5 - (T)6 * x + (T)21 / (T)12 * x2;
            (*dw)(2) = (T)-2 + (T)4 * x - (T)1.5 * x2;
            (*dw)(3) = (T)0.5 * zzzz * zzzz;
        }

        if (ddw) {

            (*ddw)(0) = (T)1.5 * z;
            (*ddw)(1) = -(T)6 + (T)42 / (T)12 * x;
            (*ddw)(2) = (T)4 - (T)3 * x;
            (*ddw)(3) = -(T)1 + x;
        }
    }
    else if (base_node == -1) {
        T z = x - 1;
        T z3 = z * z * z;
        w(0) = -z3;
        w(1) = (T)0.25 * x * ((T)12 - (T)18 * x + (T)7 * x2);
        w(2) = (T)1 / (T)12 * ((T)18 - (T)11 * x) * x2;
        w(3) = (T)1 / (T)6 * x3;
        if (dw) {
            (*dw)(0) = -(T)3 * z * z;
            (*dw)(1) = (T)3 - 9 * x + (T)21 / (T)4 * x2;
            (*dw)(2) = (T)3 * x - (T)33 / (T)12 * x2;
            (*dw)(3) = (T)0.5 * x2;
        }

        if (ddw) {

            (*ddw)(0) = -(T)6 * z;

            (*ddw)(1) = -(T)9 + (T)21 / (T)2 * x;

            (*ddw)(2) = (T)3 - (T)33 / (T)6 * x;

            (*ddw)(3) = (T)x;
        }
    }
    else if (base_node == size - 6) {
        //  std::cout<<"size-3 is "<< x<< "base node is "<< base_node<<std::endl<<std::flush;
        T z = x - base_node_x - (T)1;
        T z2 = z * z;
        T z3 = z * z * z;
        T y = -1 + z;
        w(0) = -(T)1 / (T)6 * y * y * y;
        w(1) = (T)1 / (T)6 * ((T)4 - (T)6 * z2 + (T)3 * z3);
        w(2) = (T)1 / (T)12 * ((T)2 + (T)6 * z + (T)6 * z2 - (T)7 * z3);
        w(3) = (T)0.25 * z3;

        if (dw) {
            (*dw)(0) = -(T)0.5 * y * y;
            (*dw)(1) = -(T)2 * z + (T)1.5 * z2;
            (*dw)(2) = (T).5 + z - (T)21 / (T)12 * z2;
            (*dw)(3) = (T)0.75 * z2;
        }
        if (ddw) {
            (*ddw)(0) = -y;
            (*ddw)(1) = -(T)2 + (T)3 * z;
            (*ddw)(2) = (T)1 - (T)7 / (T)2 * z;
            (*ddw)(3) = (T)1.5 * z;
        }
    }
    else if (base_node == size - 5) {

        //   std::cout<<"size-2 is "<< x<< "base node is "<< base_node<<std::endl<<std::flush;
        T z = x - base_node_x - (T)1;
        T z2 = z * z;
        T z3 = z * z * z;
        T y = -1 + z;
        w(0) = -(T)1 / (T)6 * y * y * y;
        w(1) = (T)1 / (T)12 * y * y * ((T)7 + (T)11 * z);
        w(2) = (T)0.25 * ((T)1 + (T)3 * z + (T)3 * z2 - (T)7 * z3);
        w(3) = z3;

        if (dw) {
            (*dw)(0) = -(T)0.5 * y * y;
            (*dw)(1) = -(T)0.25 - (T)2.5 * z + (T)11 / (T)4 * z2;
            (*dw)(2) = (T).75 + (T)1.5 * z - (T)21 / (T)4 * z2;
            (*dw)(3) = (T)3 * z2;
        }
        if (ddw) {
            (*ddw)(0) = -y;
            (*ddw)(1) = -(T)2.5 + (T)11 / (T)2 * z;
            (*ddw)(2) = (T)1.5 - (T)21 / (T)2 * z;
            (*ddw)(3) = (T)6 * z;
        }
    }
    else if (x == size - 3) {
        base_node -= 1;
        w(0) = (T)0;
        w(1) = (T)0;
        w(2) = (T)0;
        w(3) = (T)1;
        if (dw) {
            (*dw)(0) = (T)0;
            (*dw)(1) = (T)0;
            (*dw)(2) = (T)-3;
            (*dw)(3) = (T)3;
        }
        if (ddw) {
            (*ddw)(0) = (T)0;
            (*ddw)(1) = (T)3;
            (*ddw)(2) = (T)-9;
            (*ddw)(3) = (T)6;
        }
    }
    else {
        std::cout << "Base node is" << base_node << " size is " << size << " X is " << x << std::endl
                  << std::flush;
    }
}

} // namespace ZIRAN

#endif
