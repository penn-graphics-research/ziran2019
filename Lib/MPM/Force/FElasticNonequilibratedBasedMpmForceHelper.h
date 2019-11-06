#ifndef F_ELASTIC_NONEQUILIBRATED_BASED_MPM_FORCE_HELPER_H
#define F_ELASTIC_NONEQUILIBRATED_BASED_MPM_FORCE_HELPER_H

#include "FBasedMpmForceHelper.h"

#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/DataStructure/DataArray.h>
#include <Ziran/CS/DataStructure/DataManager.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Math/Linear/ImplicitQRSVD.h>

namespace ZIRAN {

template <class TCONST>
class FElasticNonequilibratedBasedMpmForceHelper : public FBasedMpmForceHelper<TCONST> {
public:
    using Base = FBasedMpmForceHelper<TCONST>;
    using T = typename TCONST::Scalar;
    using Base::dim;
    using Base::particles;
    using typename Base::Scratch;
    using typename Base::TM;
    using typename Base::TV;
    using typename Base::TVStack;

    T viscosity_d;
    T viscosity_v;

    explicit FElasticNonequilibratedBasedMpmForceHelper(Particles<T, dim>& particles)
        : Base(particles)
    {
    }

    void setParameters(T viscosity_d_input, T viscosity_v_input)
    {
        viscosity_d = viscosity_d_input;
        viscosity_v = viscosity_v_input;
    }

    void updateState(const DisjointRanges& subrange, StdVector<TM>& vtau, TVStack& fp)
    {
        DisjointRanges subset(subrange,
            particles.commonRanges(constitutive_model_name(),
                constitutive_model_scratch_name(),
                Base::element_measure_name(),
                F_name()));
        for (auto iter = particles.subsetIter(subset, constitutive_model_name(), constitutive_model_scratch_name(), Base::element_measure_name(), F_name()); iter; ++iter) {
            auto& constitutive_model = iter.template get<0>();
            auto& scratch = iter.template get<1>();
            const auto& element_measure = iter.template get<2>();
            const auto& F = iter.template get<3>();
            int p = iter.entryId();
            constitutive_model.updateScratch(F, scratch);
            TM vtau_local;
            constitutive_model.kirchhoff(scratch, vtau_local);
            vtau_local *= element_measure;
            assert(vtau_local == vtau_local);
            vtau[p] += vtau_local;
        }
    }

    void evolveStrain(const DisjointRanges& subrange, T dt, const StdVector<TM>& gradV)
    {
        DisjointRanges subset(subrange,
            particles.commonRanges(constitutive_model_name(),
                constitutive_model_scratch_name(),
                Base::element_measure_name(),
                F_name()));
        for (auto iter = particles.subsetIter(subset, constitutive_model_name(), F_name()); iter; ++iter) {
            auto& c = iter.template get<0>();
            auto& F = iter.template get<1>();
            int p = iter.entryId();
            F = (TM::Identity() + ((T)dt) * gradV[p]) * F;
            assert(gradV[p] == gradV[p]);

            T alpha = (T)2.0 * c.mu / viscosity_d;
            T beta = (T)2.0 * ((T)2.0 * c.mu + c.lambda * dim) / ((T)9.0 * viscosity_v) - (T)2.0 * c.mu / (viscosity_d * dim);

            TM U, V;
            TV sigma;
            singularValueDecomposition(F, U, sigma, V);
            TV epsilon_trial = sigma.array().abs().log();
            TV epsilon = (T)1.0 / ((T)1.0 + dt * alpha) * (epsilon_trial - dt * beta / ((T)1.0 + dt * (alpha + dim * beta)) * epsilon_trial.trace() * TV::Ones());
            TV exp_epsilon = epsilon.array().exp();
            F = U * exp_epsilon.asDiagonal() * V.transpose();
        }
    }

    inline AttributeName<TCONST> constitutive_model_name() override
    {
        return AttributeName<TCONST>(TCONST::name());
    }
    inline AttributeName<Scratch> constitutive_model_scratch_name() override
    {
        return AttributeName<Scratch>(TCONST::scratch_name());
    }

    inline AttributeName<TM> F_name() override
    {
        return AttributeName<TM>("Fe_Nonequilibrated");
    }
};
} // namespace ZIRAN
#endif /* ifndef F_ELASTIC_NONEQUILIBRATED_BASED_MPM_FORCE_HELPER_H */
