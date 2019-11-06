#ifndef SPH_KERNELS_H
#define SPH_KERNELS_H
#include <Ziran/Math/MathTools.h>

namespace ZIRAN {
namespace SPH_KERNELS {

//############################################################################
// Wendland kernel
//############################################################################

template <class T>
T wendlandKernelQ(const T q)
{
    using MATH_TOOLS::to_the_fourth;
    T M = (q < 2) ? ((1 + 2 * q) * to_the_fourth(1 - q / 2)) : 0;
    return M;
}

template <class T>
T wendlandKernelGradientQ(const T q)
{
    using MATH_TOOLS::cube;
    using MATH_TOOLS::to_the_fourth;
    return (q < 2) ? (to_the_fourth(q - (T)2) / (T)8 + cube(q - 2) * (1 + 2 * q) / (T)4) : 0;
}

// This computes W_b(x_a)
template <class TV>
typename TV::Scalar wendlandKernel(const TV& xa, const TV& xb, const typename TV::Scalar h)
{
    using MATH_TOOLS::cube;
    using T = typename TV::Scalar;
    static const int dim = TV::RowsAtCompileTime;
    static const T Cd = (dim == 2) ? ((T)7 / ((T)4 * M_PI)) : ((T)21 / ((T)16 * M_PI));
    T r = (xa - xb).norm();
    T q = r / h;
    T M = wendlandKernelQ(q);
    return Cd * M / std::pow(h, (T)dim);
}

// This computes grad W_b(x_a)
template <class TV>
TV wendlandKernelGradient(const TV& xa, const TV& xb, const typename TV::Scalar h)
{
    static const int dim = TV::RowsAtCompileTime;
    using T = typename TV::Scalar;
    static const T Cd = (dim == 2) ? ((T)7 / ((T)4 * M_PI)) : ((T)21 / ((T)16 * M_PI));
    T r = (xa - xb).norm();
    T q = r / h;
    T dMdq = wendlandKernelGradientQ(q);
    return (Cd * dMdq / (std::pow(h, (T)(dim + 1)) * r)) * (xa - xb);
}

//############################################################################
// Spiky kernel
//############################################################################

template <class T>
T spikyKernelGradientQ(const T q)
{
    using MATH_TOOLS::sqr;

    T dMdq = 0;
    if (q <= 2)
        dMdq = -3 * sqr(2 - q);
    else
        dMdq = 0;
    return dMdq;
}

// This computes grad W_b(x_a)
template <class TV>
TV spikyKernelGradient(const TV& xa, const TV& xb, const typename TV::Scalar h)
{
    using T = typename TV::Scalar;
    static const int dim = TV::RowsAtCompileTime;
    static const T Cd = (dim == 2) ? ((T)5 / ((T)16 * M_PI)) : ((T)15 / ((T)64 * M_PI));
    T r = (xa - xb).norm();
    T q = r / h;
    T dMdq = spikyKernelGradientQ(q);
    return (Cd * dMdq / (std::pow(h, (T)(dim + 1)) * r)) * (xa - xb);
}

//############################################################################
// Cubic B-spline kernel
//############################################################################

template <class T>
T cubicBSplineKernelQ(const T q)
{
    using MATH_TOOLS::cube;
    T M = 0;
    if (q < 1)
        M = (T)1 / (T)6 * (-4 * cube(1 - q) + cube(2 - q));
    else if (q < 2)
        M = (T)1 / (T)6 * cube(2 - q);
    else
        M = 0;
    return M;
}

template <class T>
T cubicBSplineKernelGradientQ(const T q)
{
    using MATH_TOOLS::sqr;

    T dMdq = 0;
    if (q <= 1)
        dMdq = (T)1 / (T)6 * (-(T)12 * q + (T)9 * q * q);
    else if (q <= 2)
        dMdq = -(T)0.5 * sqr(q - (T)2);
    else
        dMdq = 0;
    return dMdq;
}

template <class T>
T cubicBSplineKernelLaplacianQ(const T q)
{
    using MATH_TOOLS::sqr;

    T d2Mdq2 = 0;
    if (q <= 1)
        d2Mdq2 = (T)1 / (T)2 * ((T)6 * q - (T)4);
    else if (q <= 2)
        d2Mdq2 = (T)2 - q;
    else
        d2Mdq2 = 0;
    return d2Mdq2;
}
// This computes W_b(x_a)
template <class TV>
typename TV::Scalar cubicBSplineKernel(const TV& xa, const TV& xb, const typename TV::Scalar h)
{
    using MATH_TOOLS::cube;
    using T = typename TV::Scalar;
    static const int dim = TV::RowsAtCompileTime;
    static const T Cd = (dim == 2) ? (15.0 / (7.0 * M_PI)) : (3.0 / (2.0 * M_PI));
    T r = (xa - xb).norm();
    T q = r / h;
    T M = cubicBSplineKernelQ(q);
    return Cd * M / std::pow(h, (T)dim);
}

// This computes grad W_b(x_a)
template <class TV>
TV cubicBSplineKernelGradient(const TV& xa, const TV& xb, const typename TV::Scalar h)
{
    using T = typename TV::Scalar;
    static const int dim = TV::RowsAtCompileTime;
    static const T Cd = (dim == 2) ? (15.0 / (7.0 * M_PI)) : (3.0 / (2.0 * M_PI));
    T r = (xa - xb).norm();
    T q = r / h;
    T dMdq = cubicBSplineKernelGradientQ(q);
    return (Cd * dMdq / (std::pow(h, (T)(dim + 1)) * r)) * (xa - xb);
}

// This compute Grad W_i(x_j)
template <class TV>
typename TV::Scalar cubicBSplineKernelLaplacian(const TV& xj, const TV& xi, const typename TV::Scalar h)
{
    static const int dim = TV::RowsAtCompileTime;
    using T = typename TV::Scalar;
    static const T Cd = (dim == 2) ? (15.0 / (7.0 * M_PI)) : (3.0 / (2.0 * M_PI));
    T r = (xi - xj).norm();
    T q = r / h;
    T dMdq = cubicBSplineKernelGradientQ(q);
    T d2Mdq2 = cubicBSplineKernelLaplacianQ(q);
    T result = 0;
    result += Cd * dMdq / (std::pow(h, (T)(dim + 1)) * r);
    result += Cd * d2Mdq2 / (std::pow(h, (T)(dim + 2)));
    return result;
}
}
} // namespace ZIRAN::SPH_KERNELS
#endif
