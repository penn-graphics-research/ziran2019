#ifndef MPM_FORCE_BASE_H
#define MPM_FORCE_BASE_H

#include "MpmForceHelperBase.h"

#include <MPM/MpmGrid.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Sim/Scene.h>

namespace ZIRAN {

class PlasticityApplierBase;

template <class T, int dim>
class MpmForceBase : public LagrangianForce<T, dim> {
public:
    static constexpr int interpolation_degree = ZIRAN_MPM_DEGREE;

    typedef Vector<T, 4> TV4;
    typedef Matrix<T, 4, 4> TM4;
    typedef Vector<T, dim> TV;
    typedef Vector<int, dim> IV;
    typedef Matrix<T, dim, dim> TM;

    using Base = LagrangianForce<T, dim>;
    using Base::addScaledForceDifferential;
    using Base::addScaledForces;
    using Base::reinitialize;
    using typename Base::TVStack;

    const bool& mls_mpm;
    T D_inverse;

    const T& dx;
    const T& dt;
    const TVStack& dv;
    const TVStack& vn;

    // using Base::addScaledForces;
    // using Base::addScaledForceDifferential;
    // using Base::gatherFromPadsToTVStack;

    // particle data
    Particles<T, dim>& particles;
    Scene<T, dim>& scene;

    StdVector<TV>& scratch_xp;
    TVStack& scratch_vp;
    TVStack& scratch_fp;

    StdVector<TM>& scratch_gradV;
    StdVector<TM>& scratch_stress;

    StdVector<std::unique_ptr<MpmForceHelperBase<T, dim>>> helpers;

    MpmGrid<T, dim>& grid;
    StdVector<uint64_t>& particle_base_offset;
    StdVector<int>& particle_order;
    std::vector<std::pair<int, int>>& particle_group;
    std::vector<uint64_t>& block_offset;
    int& num_nodes;
    StdVector<std::unique_ptr<PlasticityApplierBase>>& plasticity_appliers;
    bool& full_implicit;

    MpmForceBase(const bool& mls_mpm, const T& dx, const T& dt, Particles<T, dim>& particles,
        Scene<T, dim>& scene,
        StdVector<TV>& scratch_xp,
        TVStack& scratch_vp,
        TVStack& scratch_fp,
        StdVector<TM>& scratch_gradV,
        StdVector<TM>& scratch_stress,
        const TVStack& dv, const TVStack& vn,
        MpmGrid<T, dim>& grid, StdVector<uint64_t>& particle_base_offset, StdVector<int>& particle_order, std::vector<std::pair<int, int>>& particle_group, std::vector<uint64_t>& block_offset, int& num_nodes,
        StdVector<std::unique_ptr<PlasticityApplierBase>>& plasticity_appliers, bool& full_implicit);

    virtual ~MpmForceBase();

    void reinitialize() override;

    void backupStrain();

    void restoreStrain();

    void computeVAndGradV();

    void computeDvAndGradDv(const TVStack& dv) const;

    template <bool USE_MLS_MPM>
    void rasterizeForceToTVStack(const T scale, TVStack& force) const;

    void updateParticleState();

    void updateParticleImplicitState();

    void addScaledForces(const T scale, TVStack& forces) const override;

    /**
      Compute for a vector function f defined on active grid nodes, compute f and grad f on particles
      */

    template <class Func>
    void evalInterpolantAndGradient(Func&& f, TVStack& f_eval, StdVector<TM>& grad_f) const;

    template <class TCONST, class TPCONST>
    void computeStressDifferential(const StdVector<TM>& gradDv, StdVector<TM>& dstress) const;

    void addScaledForceDifferential(const T scale, const TVStack& dv, TVStack& df) const override;

    // This is only called when implicit.
    void updatePositionBasedState() override;

    void updatePositionBasedState(const StdVector<TV>& x) override;

    void evolveStrain(T dt);

    void updateStrainWithFullImplicit() override;

    T totalEnergy() const override;

    struct Proxy {
        StdVector<std::unique_ptr<MpmForceHelperBase<T, dim>>>& helpers;
        Particles<T, dim>& particles;

        Proxy(StdVector<std::unique_ptr<MpmForceHelperBase<T, dim>>>& helpers,
            Particles<T, dim>& particles)
            : helpers(helpers)
            , particles(particles)
        {
        }

        /**
           Implicit Cast to HelperT 
           Get the Force Helper with type HelperT or insert it if it's not there
           Not threadsafe
        */
        template <class HelperT>
        operator HelperT&()
        {
            for (auto& helper : helpers)
                if (HelperT* h = dynamic_cast<HelperT*>(helper.get()))
                    return *h;
            // We didn't find it
            helpers.emplace_back(std::make_unique<std::remove_const_t<HelperT>>(particles));
            return *dynamic_cast<HelperT*>(helpers.back().get());
        }
    };

    Proxy getHelper()
    {
        return Proxy(helpers, particles);
    }
};
} // namespace ZIRAN
#endif
