#include <Ziran/Physics/ConstitutiveModel/StvkWithHencky.h>
#include <Ziran/Physics/ConstitutiveModel/HyperelasticConstitutiveModel.h>
#include <Ziran/Math/Linear/DenseExt.h>
#include <Ziran/Math/Linear/ImplicitQRSVD.h>
#include <Ziran/CS/Util/BinaryIO.h>

namespace ZIRAN {

template <class T, int _dim>
StvkWithHencky<T, _dim>::StvkWithHencky(const T E, const T nu)
{
    setLameParameters(E, nu);
}

template <class T, int _dim>
void StvkWithHencky<T, _dim>::setLameParameters(const T E, const T nu)
{
    lambda = E * nu / (((T)1 + nu) * ((T)1 - (T)2 * nu));
    mu = E / ((T)2 * ((T)1 + nu));
}

template <class T, int _dim>
void StvkWithHencky<T, _dim>::updateScratch(const TM& new_F, Scratch& scratch) const
{
    using namespace EIGEN_EXT;
    scratch.F = new_F;
    singularValueDecomposition(scratch.F, scratch.U, scratch.sigma, scratch.V);
    scratch.log_sigma = scratch.sigma.array().abs().log();
}

template <class T, int _dim>
Matrix<T, 2, 2> StvkWithHencky<T, _dim>::Bij(const TV& sigma, int i, int j, T clamp_value) const
{
    auto mingClampMagnitude = [&](const T input) {
        T magnitude = input > 0 ? input : -input;
        T sign = input > 0 ? 1.f : -1.f;
        T output = magnitude > clamp_value ? magnitude : clamp_value;
        return output * sign;
    };
    T eps = (T)1e-8;
    if constexpr (dim == 2) {
        TV dE = firstPiolaDiagonal(sigma);
        T B_Pij = 0.5 * (dE[i] + dE[j]) / mingClampMagnitude(sigma[i] + sigma[j]);

        TV logS = sigma.array().abs().log();
        T prod = sigma(0) * sigma(1);
        T q = std::max(sigma(0) / sigma(1) - 1, -1 + eps);
        T h = (std::fabs(q) < eps) ? 1 : (std::log1p(q) / q);
        T t = h / sigma(1);
        T z = logS(1) - t * sigma(1);
        T B_Mij = -0.5 * (lambda * (logS(0) + logS(1)) + 2 * mu * z) / prod;

        Matrix<T, 2, 2> B_P_Const;
        B_P_Const << 1, 1, 1, 1;
        Matrix<T, 2, 2> B_M_Const;
        B_M_Const << 1, -1, -1, 1;
        return B_M_Const * B_Pij + B_P_Const * B_Mij;
    }
    else {
        TV dE = firstPiolaDiagonal(sigma);
        T B_Pij = 0.5 * (dE[i] + dE[j]) / mingClampMagnitude(sigma[i] + sigma[j]);

        TV logS = sigma.array().abs().log();
        T sum_log = logS(0) + logS(1) + logS(2);
        T prodij = sigma(i) * sigma(j);
        T B_Mij = -0.5 * (lambda * sum_log + 2 * mu * MATH_TOOLS::diff_interlock_log_over_diff(sigma(i), sigma(j), logS(j), eps)) / prodij;

        Matrix<T, 2, 2> B_P_Const;
        B_P_Const << 1, 1, 1, 1;
        Matrix<T, 2, 2> B_M_Const;
        B_M_Const << 1, -1, -1, 1;
        return B_M_Const * B_Pij + B_P_Const * B_Mij;
    }
}

template <class T, int _dim>
Matrix<T, 2, 2> StvkWithHencky<T, _dim>::BijFull(const TV& zi, const TV& sigma, const T& second_term, int i, int j, T clamp_value) const
{
    if constexpr (dim == 2) {
        auto mingClampMagnitude = [&](const T input) {
            T magnitude = input > 0 ? input : -input;
            T sign = input > 0 ? 1.f : -1.f;
            T output = magnitude > clamp_value ? magnitude : clamp_value;
            return output * sign;
        };
        TV dE = firstPiolaDiagonal(zi);
        T B_Pij = 0.5 * (dE[i] + dE[j]) / mingClampMagnitude(sigma[i] + sigma[j]);

        T eps = (T)1e-8;
        TV logS = zi.array().abs().log();
        T prod = zi(0) * zi(1);
        T q = std::max(zi(0) / zi(1) - 1, -1 + eps);
        T h = (std::fabs(q) < eps) ? 1 : (std::log1p(q) / q);
        T t = h / zi(1);
        T z = logS(1) - t * zi(1);
        T B_Mij = -0.5 * (lambda * (logS(0) + logS(1)) + 2 * mu * z) / prod;
        B_Mij *= second_term;

        Matrix<T, 2, 2> B_P_Const;
        B_P_Const << 1, 1, 1, 1;
        Matrix<T, 2, 2> B_M_Const;
        B_M_Const << 1, -1, -1, 1;
        return B_M_Const * B_Pij + B_P_Const * B_Mij;
    }
    else {
        ZIRAN_ASSERT(false, "not implemented");
    }
}

/**
       psi = mu tr((log S)^2) + 1/2 lambda (tr(log S))^2
     */
template <class T, int _dim>
T StvkWithHencky<T, _dim>::psi(const Scratch& s) const
{
    TV log_sigma_squared = s.log_sigma.array().square();
    T trace_log_sigma = s.log_sigma.array().sum();
    return mu * log_sigma_squared.array().sum() + (T).5 * lambda * trace_log_sigma * trace_log_sigma;
}

/**
       P = U (2 mu S^{-1} (log S) + lambda tr(log S) S^{-1}) V^T
     */
template <class T, int _dim>
void StvkWithHencky<T, _dim>::firstPiola(const Scratch& s, TM& P) const
{
    TV sigma_inverse = s.sigma.array().inverse();
    TM mu_term = sigma_inverse.asDiagonal();
    mu_term *= s.log_sigma.asDiagonal();
    TM lambda_term = s.log_sigma.array().sum() * sigma_inverse.asDiagonal();
    TM Phat = (T)2 * mu * mu_term + lambda * lambda_term;
    P = s.U * Phat * (s.V.transpose());
}

template <class T, int _dim>
void StvkWithHencky<T, _dim>::firstPiolaDerivative(const Scratch& s, Hessian& dPdF) const
{
    // TODO: implement, and modify hessianImplemented() function.
    ZIRAN_ASSERT(false, "not implemented");
}

template <class T, int _dim>
T StvkWithHencky<T, _dim>::psiDiagonal(const TV& sigma) const
{
    TV logS = sigma.array().abs().log();
    TV logS_squared = logS.array().square();
    T trace_logS = logS.array().sum();
    return mu * logS_squared.array().sum() + (T).5 * lambda * trace_logS * trace_logS;
}

template <class T, int _dim>
typename StvkWithHencky<T, _dim>::TV StvkWithHencky<T, _dim>::firstPiolaDiagonal(const TV& sigma) const
{
    TV dE;
    if constexpr (dim == 1) {
        // TODO: implement, and modify hessianImplemented() function.
        ZIRAN_ASSERT(false, "not implemented");
    }
    else if constexpr (dim == 2) {
        TV logS = sigma.array().abs().log();
        T g = 2 * mu + lambda;
        dE[0] = (g * logS(0) + lambda * logS(1)) / sigma(0);
        dE[1] = (g * logS(1) + lambda * logS(0)) / sigma(1);
    }
    else {
        TV logS = sigma.array().abs().log();
        T sum_log = logS(0) + logS(1) + logS(2);
        dE[0] = (2 * mu * logS(0) + lambda * sum_log) / sigma(0);
        dE[1] = (2 * mu * logS(1) + lambda * sum_log) / sigma(1);
        dE[2] = (2 * mu * logS(2) + lambda * sum_log) / sigma(2);
    }
    return dE;
}

template <class T, int _dim>
typename StvkWithHencky<T, _dim>::TM StvkWithHencky<T, _dim>::firstPiolaDerivativeDiagonal(const TV& sigma) const
{
    TM ddE;
    if constexpr (dim == 1) {
        // TODO: implement, and modify hessianImplemented() function.
        ZIRAN_ASSERT(false, "not implemented");
    }
    else if constexpr (dim == 2) {
        TV logS = sigma.array().abs().log();
        T g = 2 * mu + lambda;
        T prod = sigma(0) * sigma(1);
        ddE(0, 0) = (g * (1 - logS(0)) - lambda * logS(1)) / (sigma(0) * sigma(0));
        ddE(1, 1) = (g * (1 - logS(1)) - lambda * logS(0)) / (sigma(1) * sigma(1));
        ddE(0, 1) = ddE(1, 0) = lambda / prod;
    }
    else {
        TV logS = sigma.array().abs().log();
        T g = 2 * mu + lambda;
        ddE(0, 0) = (g * (1 - logS(0)) - lambda * (logS(1) + logS(2))) / (sigma(0) * sigma(0));
        ddE(1, 1) = (g * (1 - logS(1)) - lambda * (logS(0) + logS(2))) / (sigma(1) * sigma(1));
        ddE(2, 2) = (g * (1 - logS(2)) - lambda * (logS(0) + logS(1))) / (sigma(2) * sigma(2));
        ddE(0, 1) = ddE(1, 0) = lambda / (sigma(0) * sigma(1));
        ddE(0, 2) = ddE(2, 0) = lambda / (sigma(0) * sigma(2));
        ddE(1, 2) = ddE(2, 1) = lambda / (sigma(1) * sigma(2));
    }
    return ddE;
}

template <class T, int _dim>
void StvkWithHencky<T, _dim>::firstPiolaDifferential(const Scratch& s, const TM& dF, TM& dP) const
{
    // TODO: implement, and modify hessianImplemented() function.
    ZIRAN_ASSERT(false, "not implemented");
}

template <class T, int _dim>
bool StvkWithHencky<T, _dim>::isC2(const Scratch& s, T tolerance) const
{
    return s.sigma.prod() > tolerance; // due to the log sigma term
}

/**
       Returns whether dP (or dPdF) is implemented
    */
template <class T, int _dim>
bool StvkWithHencky<T, _dim>::hessianImplemented() const
{
    return false;
}

template <class T, int _dim>
void StvkWithHencky<T, _dim>::write(std::ostream& out) const
{
    writeEntry(out, mu);
    writeEntry(out, lambda);
}

template <class T, int _dim>
StvkWithHencky<T, _dim> StvkWithHencky<T, _dim>::read(std::istream& in)
{
    StvkWithHencky<T, _dim> model;
    model.mu = readEntry<T>(in);
    model.lambda = readEntry<T>(in);
    return model;
}

template class StvkWithHencky<double, 1>;
template class StvkWithHencky<double, 2>;
template class StvkWithHencky<double, 3>;
template class StvkWithHencky<float, 1>;
template class StvkWithHencky<float, 2>;
template class StvkWithHencky<float, 3>;
} // namespace ZIRAN
