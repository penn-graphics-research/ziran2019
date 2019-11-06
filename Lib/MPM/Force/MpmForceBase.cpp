#include "MpmForceBase.h"

#include <Ziran/CS/Util/AttributeNamesForward.h>
#include <Ziran/CS/Util/Timer.h>
#include <Ziran/Math/Geometry/Elements.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Physics/LagrangianForce/LagrangianForce.h>
#include <Ziran/Sim/Scene.h>

namespace ZIRAN {
template <class T, int dim>
MpmForceBase<T, dim>::MpmForceBase(const bool& mls_mpm, const T& dx, const T& dt, Particles<T, dim>& particles,
    Scene<T, dim>& scene,
    StdVector<TV>& scratch_xp,
    TVStack& scratch_vp,
    TVStack& scratch_fp,
    StdVector<TM>& scratch_gradV,
    StdVector<TM>& scratch_stress,
    const TVStack& dv, const TVStack& vn,
    MpmGrid<T, dim>& grid, StdVector<uint64_t>& particle_base_offset, StdVector<int>& particle_order, std::vector<std::pair<int, int>>& particle_group, std::vector<uint64_t>& block_offset, int& num_nodes,
    StdVector<std::unique_ptr<PlasticityApplierBase>>& plasticity_appliers,
    bool& full_implicit)

    : mls_mpm(mls_mpm)
    , D_inverse(0)
    , dx(dx)
    , dt(dt)
    , dv(dv)
    , vn(vn)
    , particles(particles)
    , scene(scene)
    , scratch_xp(scratch_xp)
    , scratch_vp(scratch_vp)
    , scratch_fp(scratch_fp)
    , scratch_gradV(scratch_gradV)
    , scratch_stress(scratch_stress)
    , grid(grid)
    , particle_base_offset(particle_base_offset)
    , particle_order(particle_order)
    , particle_group(particle_group)
    , block_offset(block_offset)
    , num_nodes(num_nodes)
    , plasticity_appliers(plasticity_appliers)
    , full_implicit(full_implicit)
{
}

template <class T, int dim>
MpmForceBase<T, dim>::~MpmForceBase()
{
}

template <class T, int dim>
void MpmForceBase<T, dim>::reinitialize()
{
    if (mls_mpm && D_inverse == 0) {
        T dx2 = dx * dx;
        if (interpolation_degree == 2)
            D_inverse = 4 / dx2;
        else if (interpolation_degree == 3)
            D_inverse = 3 / dx2;
        else
            ZIRAN_ASSERT(false);
    }

    for (auto& h : helpers)
        h->reinitialize();
}

template <class T, int dim>
void MpmForceBase<T, dim>::backupStrain()
{
    ZIRAN_QUIET_TIMER();
    for (auto& h : helpers)
        h->backupStrain();
}

template <class T, int dim>
void MpmForceBase<T, dim>::restoreStrain()
{
    ZIRAN_QUIET_TIMER();
    for (auto& h : helpers)
        h->restoreStrain();
}

template <class T, int dim>
void MpmForceBase<T, dim>::computeVAndGradV()
{
    ZIRAN_QUIET_TIMER();
    evalInterpolantAndGradient([&](int node_id) -> TV { return vn.col(node_id) + dv.col(node_id); }, scratch_vp, scratch_gradV);
}

template <class T, int dim>
void MpmForceBase<T, dim>::computeDvAndGradDv(const TVStack& dv) const
{
    ZIRAN_QUIET_TIMER();
    evalInterpolantAndGradient([&](int node_id) -> TV { return dv.col(node_id); }, scratch_vp, scratch_gradV);
}

template <class T, int dim>
template <bool USE_MLS_MPM>
void MpmForceBase<T, dim>::rasterizeForceToTVStack(const T scale, TVStack& force) const
{
    ZIRAN_QUIET_TIMER();

    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = TV::Zero();
    });

    auto& Xarray = particles.X.array;
    auto grid_array = grid.grid->Get_Array();

    for (uint64_t color = 0; color < (1 << dim); ++color) {
        tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
            if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                return;
            for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                int i = particle_order[idx];
                TV& Xp = Xarray[i];

                TM4 stress_density = TM4::Zero();
                TM& stress = scratch_stress[i];
                stress_density.template topLeftCorner<dim, dim>() = stress; // stress density
                const TV& fp = scratch_fp.col(i);
                stress_density.template topRightCorner<dim, 1>() = -fp; // fp (for meshed forces).

                BSplineWeights<T, dim> spline(Xp, dx);
                grid.iterateKernel(spline, particle_base_offset[i],
                    [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                        TV4 weight = TV4::Zero();
                        weight.template topLeftCorner<dim, 1>() = dw;
                        weight(3) = w;
                        if constexpr (USE_MLS_MPM) {
                            TV4 xi_minus_xp = TV4::Zero();
                            xi_minus_xp.template topLeftCorner<dim, 1>() = node.template cast<T>() * dx - Xp; // top dim entries non-zero
                            xi_minus_xp(3) = 1;
                            TV4 delta = ((stress_density * xi_minus_xp) * (w * D_inverse));
                            g.new_v -= scale * delta.template topLeftCorner<dim, 1>(); // fi -= \sum_p (Ap (xi-xp)  - fp )w_ip Dp_inv
                        }
                        else {
                            TV4 delta = (stress_density * weight);
                            g.new_v -= scale * delta.template topLeftCorner<dim, 1>();
                        }
                    });
            }
        });
    }

    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        force.col(g.idx) += g.new_v;
    });
}

template <class T, int dim>
void MpmForceBase<T, dim>::updateParticleState()
{
    ZIRAN_QUIET_TIMER();

    // zero out scratch_fp and scratch_stress, prepare for adding lagrangian force
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particles.count),
        [&](const tbb::blocked_range<size_t>& range) {
                              for (size_t b = range.begin(), b_end = range.end(); b < b_end; ++b) {
                                  scratch_fp.col(b) = TV::Zero();
                                  scratch_stress[b] = TM::Zero(); } });
    // Lagrangian
    for (auto& lf : scene.forces)
        lf->updatePositionBasedState();

    auto& vtau = scratch_stress;
    auto& fp = scratch_fp;
    auto ranges = particles.X.ranges;
    tbb::parallel_for(ranges,
        [&](DisjointRanges& subrange) {
            for (auto& h : helpers)
                h->updateState(subrange, vtau, fp); // this adds to vtau and fp
        });

    // Lagrangian
    scene.addScaledForces(1, fp); // add forces to fp
}

template <class T, int dim>
void MpmForceBase<T, dim>::updateParticleImplicitState()
{
    ZIRAN_QUIET_TIMER();
    // zero out scratch_fp and scratch_stress, prepare for adding lagrangian force
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particles.count),
        [&](const tbb::blocked_range<size_t>& range) {
                              for (size_t b = range.begin(), b_end = range.end(); b < b_end; ++b) {
                                  scratch_fp.col(b) = TV::Zero();
                                  scratch_stress[b] = TM::Zero(); } });

    // Lagrangian
    for (auto& lf : scene.forces)
        lf->updatePositionBasedState(scratch_xp);

    auto& vPFnT = scratch_stress;
    auto& fp = scratch_fp;
    auto ranges = particles.X.ranges;
    tbb::parallel_for(ranges,
        [&](DisjointRanges& subrange) {
            for (auto& h : helpers)
                h->updateImplicitState(subrange, vPFnT, fp); // this adds to vPFnT and fp
        });

    // Lagrangian
    scene.addScaledForces(1, fp); // add forces to fp
}
/**
      Compute for a vector function f defined on active grid nodes, compute f and grad f on particles
      */
template <class T, int dim>
template <class Func>
void MpmForceBase<T, dim>::evalInterpolantAndGradient(Func&& f, TVStack& f_eval, StdVector<TM>& grad_f) const
{
    grid.iterateTouchedGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = TV::Zero();
    });
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = f((int)g.idx);
    });
    for (uint64_t color = 0; color < (1 << dim); ++color) {
        tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
            if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                return;
            for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                int i = particle_order[idx];
                TM& grad_fp = grad_f[i];
                // if (grad_fp != grad_fp)
                //     continue;

                TV& Xp = particles.X[i];
                BSplineWeights<T, dim> spline(Xp, dx);

                grad_fp = TM::Zero();
                f_eval.col(i) = TV::Zero();
                grid.iterateKernel(spline, particle_base_offset[i], [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                    grad_fp.noalias() += g.new_v * dw.transpose();
                    f_eval.col(i) += g.new_v * w;
                });
            }
        });
    }
}

template <class T, int dim>
void MpmForceBase<T, dim>::addScaledForces(const T scale, TVStack& forces) const
{
    ZIRAN_QUIET_TIMER();
    if (mls_mpm) {
        rasterizeForceToTVStack<true>(scale, forces);
    }
    else {
        rasterizeForceToTVStack<false>(scale, forces);
    }
}

template <class T, int dim>
void MpmForceBase<T, dim>::addScaledForceDifferential(const T scale, const TVStack& dv, TVStack& df) const
{
    ZIRAN_QUIET_TIMER();

    computeDvAndGradDv(dv);

    // Compute per element dF
    for (auto& lf : scene.forces)
        lf->updatePositionDifferentialBasedState(scratch_vp);

    // scratch_gradV is now temporaraly used for storing gradDV (evaluated at particles)
    // scratch_vp is now temporaraly used for storing DV (evaluated at particles)

    // zero out scratch_fp and scratch_stress, prepare for adding lagrangian force differentials
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particles.count),
        [&](const tbb::blocked_range<size_t>& range) {
                              for (size_t b = range.begin(), b_end = range.end(); b < b_end; ++b) {
                                  scratch_fp.col(b) = TV::Zero();
                                  scratch_stress[b] = TM::Zero(); } });

    auto ranges = particles.X.ranges;
    tbb::parallel_for(ranges,
        [&](DisjointRanges& subrange) {
            for (auto& h : helpers) {
                if (full_implicit) {
                    h->computeStressDifferentialWithPlasticity(subrange, scratch_gradV, scratch_stress);
                }
                else {
                    h->computeStressDifferential(subrange, scratch_gradV, scratch_stress, scratch_vp, scratch_fp);
                }
            }
            // scratch_stress is now V_p^0 dP (F_p^n)^T (dP is Ap in snow paper)
        });

    // Lagrangian
    scene.addScaledForceDifferentials(1, scratch_vp, scratch_fp);

    if (mls_mpm) {
        rasterizeForceToTVStack<true>(scale, df);
    }
    else {
        rasterizeForceToTVStack<false>(scale, df);
    }
}

// This is only called when implicit.
template <class T, int dim>
void MpmForceBase<T, dim>::updatePositionBasedState()
{
    ZIRAN_QUIET_TIMER();
    computeVAndGradV();

    restoreStrain();
    evolveStrain(dt);

    if (full_implicit) {
        using TConst = StvkWithHencky<T, dim>;
        using TPlastic = DruckerPragerStvkHencky<T>;
        auto constitutive_model_name = AttributeName<TConst>(TConst::name());
        auto plastic_name = AttributeName<TPlastic>(TPlastic::name());
        auto ranges = particles.X.ranges;
        if (!particles.exist(constitutive_model_name) || !particles.exist(plastic_name))
            return;
        tbb::parallel_for(ranges, [&](DisjointRanges& subrange) {
            DisjointRanges subset(subrange, particles.commonRanges(constitutive_model_name, plastic_name, F_name<T, dim>()));
            for (auto iter = particles.subsetIter(subset, constitutive_model_name, plastic_name, F_name<T, dim>()); iter; ++iter) {
                auto& c = iter.template get<0>();
                auto& p = iter.template get<1>();
                auto& F = iter.template get<2>();
                TM U, V;
                TV sigma;
                singularValueDecomposition(F, U, sigma, V);
                TV project_sigma = p.projectSigma(c, sigma);
                F = U * project_sigma.asDiagonal() * V.transpose();
            }
        });
    }

    // xp = xn + dt * vp
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particles.count),
        [&](const tbb::blocked_range<size_t>& range) {
                              for (size_t b = range.begin(), b_end = range.end(); b < b_end; ++b) {
                                  scratch_xp[b] = particles.X.array[b] + dt * scratch_vp.col(b); } });

    updateParticleImplicitState();
}

template <class T, int dim>
void MpmForceBase<T, dim>::updatePositionBasedState(const StdVector<TV>& x)
{
    ZIRAN_ASSERT(false, "UpdatePositionBasedState(const StdVector<TV>& x) is not supported by mpm");
}

template <class T, int dim>
void MpmForceBase<T, dim>::evolveStrain(T dt)
{
    // TODO: combinme evolve strain and update state?

    ZIRAN_QUIET_TIMER();
    auto ranges = particles.X.ranges;
    tbb::parallel_for(ranges,
        [&](DisjointRanges& subrange) {
            for (auto& h : helpers)
                h->evolveStrain(subrange, dt, scratch_gradV);
        });
}

template <class T, int dim>
void MpmForceBase<T, dim>::updateStrainWithFullImplicit()
{
    restoreStrain();
    evolveStrain(dt);
}

template <class T, int dim>
T MpmForceBase<T, dim>::totalEnergy() const
{
    ZIRAN_QUIET_TIMER();
    auto ranges = particles.X.ranges;
    T result = tbb::parallel_reduce(
        ranges, 0.0,
        [&](const DisjointRanges& subset, const double& e) -> double {
            double psi = e;
            for (auto& h : helpers)
                psi += h->totalEnergy(subset);
            return psi;
        },
        [&](const double& x, const double& y) -> double {
            return x + y;
        });

    // Lagrangian
    result += scene.totalEnergy();
    return result;
}
} // namespace ZIRAN
