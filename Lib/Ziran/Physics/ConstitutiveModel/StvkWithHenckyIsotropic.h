#ifndef STVK_WITH_HENCKY_ISOTROPIC_H
#define STVK_WITH_HENCKY_ISOTROPIC_H
#include <Ziran/Physics/ConstitutiveModel/SvdBasedIsotropicHelper.h>
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/Util/BinaryIO.h>
#include <Ziran/Math/MathTools.h>
#include <tick/requires.h>

namespace ZIRAN {

template <class Derived>
class HyperelasticConstitutiveModel;

template <typename Derived>
struct ScratchTrait;

template <class T, int _dim>
class StvkWithHenckyIsotropic;

// scratch (non-state) variables for the consitutive model
template <class T, int dim>
struct StvkWithHenckyIsotropicScratch {
    using TM = Matrix<T, dim, dim>;
    using TV = Vector<T, dim>;
    TM F, U, V;
    TV sigma;
    TV logS;
    SvdBasedIsotropicHelper<T, dim> isotropic;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    StvkWithHenckyIsotropicScratch()
        : isotropic(0)
    {
    }

    static const char* name()
    {
        return "StvkWithHenckyIsotropicScratch";
    }
};

template <class T, int _dim>
class StvkWithHenckyIsotropic : public HyperelasticConstitutiveModel<StvkWithHenckyIsotropic<T, _dim>> {
public:
    static const int dim = _dim;
    static constexpr T eps = (T)1e-6;
    using Base = HyperelasticConstitutiveModel<StvkWithHenckyIsotropic<T, dim>>;
    using TM = typename Base::TM;
    using TV = typename Base::TV;
    using Strain = TM;
    using Hessian = typename Base::Hessian;
    using Scalar = typename Base::Scalar;
    using Scratch = typename HyperelasticTraits<StvkWithHenckyIsotropic<T, dim>>::ScratchType;
    using Base::firstPiolaDerivative; // TODO: a more efficient version
    using Vec = Vector<T, Eigen::Dynamic>;
    using VecBlock = Eigen::VectorBlock<Vec>;

    T mu, lambda;

    StvkWithHenckyIsotropic(const T E = (T)1, const T nu = (T)0.3)
    {
        setLameParameters(E, nu);
    }

    void setLameParameters(const T E, const T nu)
    {
        lambda = E * nu / (((T)1 + nu) * ((T)1 - (T)2 * nu));
        mu = E / ((T)2 * ((T)1 + nu));
    }

    void updateScratchSVD(const TM& new_F, Scratch& s) const // private
    {
        s.F = new_F;
        singularValueDecomposition(s.F, s.U, s.sigma, s.V);
        s.logS = s.sigma.array().abs().log();
    }

    TICK_MEMBER_REQUIRES(dim == 1)
    void updateScratch(const TM& new_F, Scratch& s) const
    {
        using namespace MATH_TOOLS;
        updateScratchSVD(new_F, s);
        T g = 2 * mu + lambda;
        T one_over_F = 1 / new_F(0, 0);
        s.isotropic.psi0 = g * one_over_F;
        s.isotropic.psi00 = g * sqr(one_over_F);
    }

    TICK_MEMBER_REQUIRES(dim == 2)
    void updateScratch(const TM& new_F, Scratch& s) const
    {
        using namespace MATH_TOOLS;
        updateScratchSVD(new_F, s);
        T g = 2 * mu + lambda;
        T prod = s.sigma(0) * s.sigma(1);
        s.isotropic.psi0 = (g * s.logS(0) + lambda * s.logS(1)) / s.sigma(0);
        s.isotropic.psi1 = (g * s.logS(1) + lambda * s.logS(0)) / s.sigma(1);
        s.isotropic.psi00 = (g * (1 - s.logS(0)) - lambda * s.logS(1)) / sqr(s.sigma(0));
        s.isotropic.psi11 = (g * (1 - s.logS(1)) - lambda * s.logS(0)) / sqr(s.sigma(1));
        s.isotropic.psi01 = lambda / prod;

        // (psi0-psi1)/(sigma0-sigma1)
        T q = std::max(s.sigma(0) / s.sigma(1) - 1, -1 + eps);
        T h = (std::fabs(q) < eps) ? 1 : (std::log1p(q) / q);
        T t = h / s.sigma(1);
        T z = s.logS(1) - t * s.sigma(1);
        s.isotropic.m01 = -(lambda * (s.logS(0) + s.logS(1)) + 2 * mu * z) / prod;

        // (psi0+psi1)/(sigma0+sigma1)
        s.isotropic.p01 = (s.isotropic.psi0 + s.isotropic.psi1) / clamp_small_magnitude(s.sigma(0) + s.sigma(1), eps);
    }

    TICK_MEMBER_REQUIRES(dim == 3)
    void updateScratch(const TM& new_F, Scratch& s) const
    {
        using namespace MATH_TOOLS;
        updateScratchSVD(new_F, s);
        T g = 2 * mu + lambda;
        T sum_log = s.logS(0) + s.logS(1) + s.logS(2);
        T prod01 = s.sigma(0) * s.sigma(1);
        T prod02 = s.sigma(0) * s.sigma(2);
        T prod12 = s.sigma(1) * s.sigma(2);
        s.isotropic.psi0 = (2 * mu * s.logS(0) + lambda * sum_log) / s.sigma(0);
        s.isotropic.psi1 = (2 * mu * s.logS(1) + lambda * sum_log) / s.sigma(1);
        s.isotropic.psi2 = (2 * mu * s.logS(2) + lambda * sum_log) / s.sigma(2);
        s.isotropic.psi00 = (g * (1 - s.logS(0)) - lambda * (s.logS(1) + s.logS(2))) / sqr(s.sigma(0));
        s.isotropic.psi11 = (g * (1 - s.logS(1)) - lambda * (s.logS(0) + s.logS(2))) / sqr(s.sigma(1));
        s.isotropic.psi22 = (g * (1 - s.logS(2)) - lambda * (s.logS(0) + s.logS(1))) / sqr(s.sigma(2));
        s.isotropic.psi01 = lambda / (s.sigma(0) * s.sigma(1));
        s.isotropic.psi02 = lambda / (s.sigma(0) * s.sigma(2));
        s.isotropic.psi12 = lambda / (s.sigma(1) * s.sigma(2));

        // (psiA-psiB)/(sigmaA-sigmaB)
        s.isotropic.m01 = -(lambda * sum_log + 2 * mu * diff_interlock_log_over_diff(s.sigma(0), s.sigma(1), s.logS(1), eps)) / prod01;
        s.isotropic.m02 = -(lambda * sum_log + 2 * mu * diff_interlock_log_over_diff(s.sigma(0), s.sigma(2), s.logS(2), eps)) / prod02;
        s.isotropic.m12 = -(lambda * sum_log + 2 * mu * diff_interlock_log_over_diff(s.sigma(1), s.sigma(2), s.logS(2), eps)) / prod12;

        // (psiA+psiB)/(sigmaA+sigmaB)
        s.isotropic.p01 = (s.isotropic.psi0 + s.isotropic.psi1) / clamp_small_magnitude(s.sigma(0) + s.sigma(1), eps);
        s.isotropic.p02 = (s.isotropic.psi0 + s.isotropic.psi2) / clamp_small_magnitude(s.sigma(0) + s.sigma(2), eps);
        s.isotropic.p12 = (s.isotropic.psi1 + s.isotropic.psi2) / clamp_small_magnitude(s.sigma(1) + s.sigma(2), eps);
    }

    static constexpr bool diagonalDifferentiable()
    {
        return false; //TODO implement diagonal functions
    }

    /**
       psi = mu tr((log S)^2) + 1/2 lambda (tr(log S))^2
     */
    T psi(const Scratch& s) const
    {
        TV logS_squared = s.logS.array().square();
        T trace_logS = s.logS.array().sum();
        return mu * logS_squared.array().sum() + (T).5 * lambda * trace_logS * trace_logS;
    }

    void firstPiola(const Scratch& s, TM& P) const
    {
        TV P_hat;
        s.isotropic.computePHat(P_hat);
        P = s.U * P_hat.asDiagonal() * s.V.transpose();
    }

    void firstPiolaDifferential(const Scratch& s, const TM& dF, TM& dP) const
    {
        TM D = s.U.transpose() * dF * s.V;
        TM K;
        s.isotropic.dPdFOfSigmaContract(D, K);
        dP = s.U * K * s.V.transpose();
    }

    Matrix<T, 2, 2> Bij(const TV& sigma, int i, int j, T clamp_value) const
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

    TV firstPiolaDiagonal(const TV& sigma) const
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

    TM firstPiolaDerivativeDiagonal(const TV& sigma) const
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

    bool isC2(const Scratch& s, T tolerance) const
    {
        return s.sigma.prod() > tolerance; // due to the log sigma term
    }

    /**
       Returns whether dP (or dPdF) is implemented
    */
    bool hessianImplemented() const
    {
        return true;
    }

    void write(std::ostream& out) const
    {
        writeEntry(out, mu);
        writeEntry(out, lambda);
    }

    static StvkWithHenckyIsotropic<T, dim> read(std::istream& in)
    {
        StvkWithHenckyIsotropic<T, dim> model;
        model.mu = readEntry<T>(in);
        model.lambda = readEntry<T>(in);
        return model;
    }

    static const char* name()
    {
        return "StvkWithHenckyIsotropic";
    }

    inline static AttributeName<StvkWithHenckyIsotropic<T, dim>> attributeName()
    {
        return AttributeName<StvkWithHenckyIsotropic<T, dim>>("StvkWithHenckyIsotropic");
    }

    static const char* scratch_name()
    {
        return Scratch::name();
    }
};

template <class T, int dim>
struct HyperelasticTraits<StvkWithHenckyIsotropic<T, dim>> {
    using ScratchType = StvkWithHenckyIsotropicScratch<T, dim>;
};

template <class T, int dim>
struct RW<StvkWithHenckyIsotropicScratch<T, dim>> {
    using Tag = NoWriteTag<StvkWithHenckyIsotropicScratch<T, dim>>;
};
} // namespace ZIRAN

#endif
