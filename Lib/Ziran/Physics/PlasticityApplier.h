#ifndef PLASTICITY_APPLIER_H
#define PLASTICITY_APPLIER_H

#include <Ziran/CS/DataStructure/DataManager.h>
#include <Ziran/CS/DataStructure/DisjointRanges.h>
#include <Ziran/Math/Geometry/Rotation.h>
#include <Ziran/Math/Linear/Decomposition.h>
#include <Ziran/Math/MathTools.h>
#include <Ziran/Physics/ConstitutiveModel/ConstitutiveModel.h>

#include <algorithm>
#include <cmath>

namespace ZIRAN {

template <class T, int dim>
class CotangentOrthotropic;

class PlasticityApplierBase {
public:
    virtual ~PlasticityApplierBase()
    {
    }

    virtual void applyPlasticity(const DisjointRanges& subrange, DataManager& data_manager) = 0; // parallel
};

template <class TConst, class TPlastic, class TStrain>
class PlasticityApplier : public PlasticityApplierBase {
public:
    AttributeName<TStrain> strain_name;

    PlasticityApplier(const std::string& strain_name_in)
        : strain_name(strain_name_in)
    {
    }

    // This is assumed to be called in a parallel loop that split all particles/elements to subranges.
    void applyPlasticity(const DisjointRanges& subrange, DataManager& data_manager) override
    {
        auto constitutive_model_name = AttributeName<TConst>(TConst::name());
        auto plastic_name = AttributeName<TPlastic>(TPlastic::name());
        if (!data_manager.exist(constitutive_model_name) || !data_manager.exist(plastic_name))
            return;

        DisjointRanges subset(subrange,
            data_manager.commonRanges(constitutive_model_name,
                strain_name,
                plastic_name));
        for (auto it = data_manager.subsetIter(subset, constitutive_model_name, strain_name, plastic_name); it; ++it) {
            auto& c = it.template get<0>();
            auto& s = it.template get<1>();
            auto& p = it.template get<2>();
            p.projectStrain(c, s);
        }
    }
};

template <class T>
class DummyPlasticity {
public:
    template <class TConst>
    bool projectStrain(TConst& c, Matrix<T, TConst::dim, TConst::dim>& strain)
    {
        return false;
    }

    template <class TConst>
    Vector<T, TConst::dim> projectSigma(TConst& c, const Vector<T, TConst::dim>& sigma)
    {
        return sigma;
    }

    template <class TConst>
    Matrix<T, TConst::dim, TConst::dim> projectSigmaDerivative(TConst& c, const Vector<T, TConst::dim>& sigma)
    {
        return Matrix<T, TConst::dim, TConst::dim>::Identity();
    }

    template <class TConst>
    void projectSigmaAndDerivative(TConst& c, const Vector<T, TConst::dim>& sigma, Vector<T, TConst::dim>& projectedSigma, Matrix<T, TConst::dim, TConst::dim>& projectedSigmaDerivative)
    {
        projectedSigma = sigma;
        projectedSigmaDerivative = Matrix<T, TConst::dim, TConst::dim>::Identity();
    }

    template <int dim>
    T secondTerm(const StvkWithHencky<T, dim>& c, const Vector<T, dim>& sigma, int i, int j)
    {
        return (T)1;
    }

    static const char* name()
    {
        return "DummyPlasticity";
    }
};

template <class T>
class NonAssociativeCamClay {
public:
    T logJp,
        M, beta, xi;
    bool hardeningOn;

    NonAssociativeCamClay(T logJp = 0, T friction_angle = 0, T beta = 0, T xi = 0, int dim = 2, bool hardeningOn = true);

    template <class TConst>
    bool projectStrain(TConst& c, Matrix<T, TConst::dim, TConst::dim>& strain);

    void fillAttributesToVec3(Vector<T, 3>& data);
    static const char* name();
};

/**
   This is the Drucker Prager plasticity model from 

   Drucker-Prager Elastoplasticity for Sand Animation, 
   G. Klar, T. Gast, A. Pradhana, C. Fu, C. Schroeder, C. Jiang, J. Teran,
   ACM Transactions on Graphics (SIGGRAPH 2016).

   It assumes the StvkHencky elasticity consitutive model.
 */
template <class T>
class DruckerPragerStvkHencky {
public:
    T alpha; // friction_coeff
    T beta; // hardening coeff
    T logJp;
    T cohesion;
    bool volume_correction;

    DruckerPragerStvkHencky(const T friction_angle = 30, const T beta = 1, const T cohesion = 0, const bool volume_correction = true);

    void setParameters(const T friction_angle_in, const T beta_in, const T cohesion_in);

    // strain s is deformation F
    template <class TConst>
    bool projectStrain(TConst& c, Matrix<T, TConst::dim, TConst::dim>& strain);

    template <class TConst>
    Vector<T, TConst::dim> projectSigma(TConst& c, const Vector<T, TConst::dim>& sigma);

    template <class TConst>
    Matrix<T, TConst::dim, TConst::dim> projectSigmaDerivative(TConst& c, const Vector<T, TConst::dim>& sigma);

    template <class TConst>
    void projectSigmaAndDerivative(TConst& c, const Vector<T, TConst::dim>& sigma, Vector<T, TConst::dim>& projectedSigma, Matrix<T, TConst::dim, TConst::dim>& projectedSigmaDerivative);

    template <int dim>
    T secondTerm(const StvkWithHencky<T, dim>& c, const Vector<T, dim>& sigma, int i, int j)
    {
        typedef Vector<T, dim> TV;
        const T eps = (T)1e-6;
        TV epsilon = sigma.array().abs().max(1e-4).log() - cohesion;
        T trace_epsilon = epsilon.sum();
        TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
        T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
        if (trace_epsilon >= 0) // case II: project to tip
        {
            return 0;
        }
        T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
        T delta_gamma = epsilon_hat_norm + (dim * c.lambda + 2 * c.mu) / (2 * c.mu) * trace_epsilon * alpha;
        TV H;
        if (delta_gamma <= 0) // case I: inside yield surface
        {
            return 1;
        }
        else {
            H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat + TV::Constant(cohesion); // case III: projection
        }
        return ((T)1 - delta_gamma / epsilon_hat_norm) * MATH_TOOLS::diff_exp_over_diff(H(i), H(j), eps) * MATH_TOOLS::diff_log_over_diff(sigma(i), sigma(j), eps);
    }

    template <class TConst>
    void computeSigmaPInverse(TConst& c, const Vector<T, TConst::dim>& sigma_e, Vector<T, TConst::dim>& sigma_p_inv);

    inline static AttributeName<DruckerPragerStvkHencky<T>> attributeName()
    {
        return AttributeName<DruckerPragerStvkHencky<T>>("DruckerPragerStvkHencky");
    }

    static const char* name();
};

template <class T, int dim>
class VonMisesStvkHencky {
public:
    using TV = Vector<T, dim>;

    T yield_stress;
    T xi; //hardening scale
    T fail_stress;
    TV crack_normal;
    bool broken = false;
    bool three_way_cut = false;
    T three_way_cut_sensitivity = 1;

    VonMisesStvkHencky()
        : yield_stress(0)
        , xi(0)
        , fail_stress(0)
    {
    }

    VonMisesStvkHencky(const T yield_stress, const T fail_stress, const T xi);

    void setParameters(const T yield_stress_in, const T xi_in, const T fail_stress_in = -1);

    // strain s is deformation F
    template <class TConst>
    bool projectStrain(TConst& c, Matrix<T, TConst::dim, TConst::dim>& strain);

    template <class TConst>
    void projectStrainDiagonal(TConst& c, Vector<T, TConst::dim>& sigma);

    template <class TConst>
    Vector<T, TConst::dim> projectSigma(TConst& c, const Vector<T, TConst::dim>& sigma);

    template <class TConst>
    Matrix<T, TConst::dim, TConst::dim> projectSigmaDerivative(TConst& c, const Vector<T, TConst::dim>& sigma);

    template <class TConst>
    void projectSigmaAndDerivative(TConst& c, const Vector<T, TConst::dim>& sigma, Vector<T, TConst::dim>& projectedSigma, Matrix<T, TConst::dim, TConst::dim>& projectedSigmaDerivative);

    T secondTerm(const StvkWithHencky<T, dim>& c, const Vector<T, dim>& sigma, int i, int j)
    {
        typedef Vector<T, dim> TV;
        TV H;
        //TV epsilon = sigma.array().log();
        TV epsilon = sigma.array().max(1e-4).log();
        T trace_epsilon = epsilon.sum();
        TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
        T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
        T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
        T delta_gamma = epsilon_hat_norm - yield_stress / (2 * c.mu);
        if (delta_gamma <= 0) // case I
        {
            H = epsilon;
            return (T)1;
        }
        else {
            H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat; // case II
        }
        const T eps = (T)1e-6;
        return ((T)1 - delta_gamma / epsilon_hat_norm) * MATH_TOOLS::diff_exp_over_diff(H(i), H(j), eps) * MATH_TOOLS::diff_log_over_diff(sigma(i), sigma(j), eps);
    }

    template <class TConst>
    void computeSigmaPInverse(TConst& c, const Vector<T, TConst::dim>& sigma_e, Vector<T, TConst::dim>& sigma_p_inv);

    void write(std::ostream& out) const
    {
        writeEntry(out, yield_stress);
        writeEntry(out, xi);
        writeEntry(out, fail_stress);
        writeEntry(out, crack_normal);
        writeEntry(out, broken);
        writeEntry(out, three_way_cut);
    }

    static VonMisesStvkHencky<T, dim> read(std::istream& in)
    {
        VonMisesStvkHencky<T, dim> model;
        model.yield_stress = readEntry<T>(in);
        model.xi = readEntry<T>(in);
        model.fail_stress = readEntry<T>(in);
        model.crack_normal = readEntry<TV>(in);
        model.broken = readEntry<bool>(in);
        model.three_way_cut = readEntry<bool>(in);
        return model;
    }

    inline static AttributeName<VonMisesStvkHencky<T, dim>> attributeName()
    {
        return AttributeName<VonMisesStvkHencky<T, dim>>("VonMisesStvkHencky");
    }

    static const char* name();
};

template <class T, int dim>
struct RW<VonMisesStvkHencky<T, dim>> {
    using Tag = CustomTypeTag<VonMisesStvkHencky<T, dim>>;
};
} // namespace ZIRAN
#endif
