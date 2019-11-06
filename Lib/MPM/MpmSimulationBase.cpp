#pragma once
#include "MpmSimulationBase.h"
#include <Ziran/Math/Geometry/SourceCollisionObject.h>

#include <math.h>

#include <math.h>

#include <tbb/concurrent_unordered_map.h>
#include <tbb/tbb.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#undef B2

#include <Ziran/CS/Util/Timer.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Math/Geometry/AnalyticLevelSet.h>
#include <Ziran/Math/Geometry/CollisionObject.h>
#include <Ziran/Math/Geometry/PartioIO.h>
#include <Ziran/Math/Splines/BSplines.h>
#include <Ziran/Physics/LagrangianForce/Inertia.h>
#include <Ziran/Physics/LagrangianForce/LagrangianForce.h>
#include <Ziran/Physics/PlasticityApplier.h>
#include <Ziran/Sim/DiffTest.h>
#include <Ziran/Sim/Scene.h>

#include <MPM/Force/MpmForceBase.h>
#include "MpmSimulationDataAnalysis.h"

namespace ZIRAN {

template <class T, int dim>
MpmSimulationBase<T, dim>::MpmSimulationBase()
    : transfer_scheme(APIC_blend_RPIC)
    , apic_rpic_ratio((T)1)
    , flip_pic_ratio((T)0)
    , flip_rpic_ratio((T)0)
    , print_stats(false)
    , write_partio(true)
    , write_meshes(true)
    , autorestart(true)
    , gravity(TV::Unit(1) * (0))
    , cfl((T).6)
    , dx(1)
    , dt((T)1.0 / 24)
    , element_partitions(0)
    , particles(scene.particles)
    , element_measure(particles.add(element_measure_name<T>()))
    , objective(*this)
    , newton(objective, (T)1e-5, 3)
    , force((forces.emplace_back(std::make_unique<MpmForceBase<T, dim>>(mls_mpm, dx, dt, particles, scene, scratch_xp, scratch_vp, scratch_fp, scratch_gradV, scratch_stress, dv, vn, grid, particle_base_offset, particle_order, particle_group, block_offset, num_nodes, plasticity_appliers, full_implicit)), forces.back()))
{
    symplectic = true;
    quasistatic = false;
    openvdb::initialize();
    inertia = std::make_unique<MassLumpedInertia<T, dim>>(mass_matrix, dv, dt);
    explicit_velocity_field = nullptr;
    ignoreCollisionObject = false;
    mls_mpm = false;
    useTrialCollision = true; //default to true
}

template <class T, int dim>
MpmSimulationBase<T, dim>::~MpmSimulationBase()
{
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::setDx(const T new_dx)
{
    dx = new_dx;
}

template <class T, int dim>
auto MpmSimulationBase<T, dim>::getMpmForce() -> MpmForceBase<T, dim>*
{
    return force.get();
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::addCollisionObject(AnalyticCollisionObject<T, dim>&& object)
{
    collision_objects.emplace_back(std::make_unique<AnalyticCollisionObject<T, dim>>(std::move(object)));
}

template <class T, int dim>
int MpmSimulationBase<T, dim>::addSourceCollisionObject(SourceCollisionObject<T, dim>&& object)
{
    collision_objects.emplace_back(std::make_unique<SourceCollisionObject<T, dim>>(std::move(object)));
    return collision_objects.size() - 1;
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::initialize()
{
    ZIRAN_TIMER();
    Base::initialize();
    ZIRAN_INFO("Initialize called");
    ZIRAN_INFO("end_frame = ", end_frame);

    if (mls_mpm) {
        ZIRAN_ASSERT(transfer_scheme == APIC_blend_RPIC && interpolation_degree != 1 && symplectic == true,
            "mls_mpm requires apic quadratic/cubic and only supports explicit for now");
        ZIRAN_INFO("Using mls_mpm.");
    }

    if ((transfer_scheme == APIC_blend_RPIC)
        && interpolation_degree != 1) {

        // use datamanager to store C
        particles.add(C_name<TM>());
        // compute D inverse
        T dx2 = dx * dx;
        if constexpr (interpolation_degree == 2)
            D_inverse = 4 / dx2;
        else if (interpolation_degree == 3)
            D_inverse = 3 / dx2;
        else
            ZIRAN_ASSERT(false);

        ZIRAN_ASSERT(apic_rpic_ratio >= 0 && apic_rpic_ratio <= 1);

        // resize C array
        const DisjointRanges& x_range = particles.X.ranges;
        particles.DataManager::get(C_name<TM>()).lazyResize(x_range, TM::Zero());
    }

    if (!symplectic) {
        objective.initialize(
            [&](TVStack& dv) {
                for (auto iter = collision_nodes.begin(); iter != collision_nodes.end(); ++iter) {
                    int node_id = (*iter).node_id;
                    (*iter).project(dv.col(node_id));
                }
            });
    }

    for (auto& em : scene.element_managers) {
        em->registerParticleReorderingCallback(particles);
    }

    writeSimulationInformation();
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::reinitialize()
{
    ZIRAN_TIMER();
    Base::reinitialize();

    if (!restarting) {
        if (element_partitions) {
            int count = 0;
            bool repartitioned = false;
            for (auto& em : scene.element_managers) {
                if (em && em->needRepartitioning()) {
                    repartitioned = true;
                    em->partitionAndReindexElements(element_partitions);

                    if (count++ == 0) {

                        StdVector<int> old_to_new;
                        StdVector<int> new_to_old;
                        if (em->renumberParticles(particles.count, old_to_new, new_to_old)) {
                            particles.reorder(new_to_old);
                        }
                    }
                }
            }

            if (repartitioned) {

                scene.segmesh_to_write.indices.clear();
                AttributeName<int> no_write("no_write");

                for (auto& em : scene.element_managers) {

                    // TODO: rebuild boundary for tet mesh

                    if (SimplexElements<T, 2, dim>* tris = dynamic_cast<SimplexElements<T, 2, dim>*>(em.get()))
                        scene.trimesh_to_write.indices = tris->indices.array;

                    if (SimplexElements<T, 1, dim>* segs = dynamic_cast<SimplexElements<T, 1, dim>*>(em.get())) {
                        if (!segs->exist(no_write)) {
                            ZIRAN_INFO("segmesh_to_write did not find no_write tag");
                            scene.segmesh_to_write.indices = segs->indices.array;
                        }
                        else {
                            ZIRAN_INFO("segmesh_to_write found no_write tag");
                            DataArray<int>& no_write_array = segs->get(no_write);
                            for (size_t z = 0; z < segs->indices.array.size(); ++z) {
                                if (no_write_array.valueId(z) == -1) // that means z is not labled no write
                                    scene.segmesh_to_write.indices.push_back(segs->indices.array[z]);
                            }
                        }
                    }
                }
            }
        }
    }

    for (size_t k = 0; k < collision_objects.size(); ++k) {
        if (step.time == 0 || restarting) {
            if (collision_objects[k]->updateState)
                collision_objects[k]->updateState(step.time, *collision_objects[k]);
        }
    }

    collision_nodes.clear();

    force->reinitialize();
    inertia->reinitialize();

    sortParticlesAndPolluteGrid();

    // TV
    scratch_xp.resize(particles.count, TV::Constant(0));
    scratch_vp.resize(dim, particles.count);
    scratch_fp.resize(dim, particles.count);
    // TM
    scratch_gradV.resize(particles.count, TM::Constant(0));
    scratch_stress.resize(particles.count, TM::Constant(0));
    // resize C array
    if ((transfer_scheme == APIC_blend_RPIC)
        && interpolation_degree != 1) {
        const DisjointRanges& x_range = particles.X.ranges;
        particles.DataManager::get(C_name<TM>()).lazyResize(x_range, TM::Zero());
    }
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::beginTimeStep(int frame, int substep, double time, double dt)
{
    this->dt = dt;
    for (auto f : begin_time_step_callbacks) {
        f(frame, substep, time, dt);
    }
}

template <class T, int dim>
T MpmSimulationBase<T, dim>::totalEnergy(std::string info)
{
    T pe = 0;
    for (auto& lf : forces) {
        pe += lf->totalEnergy();
    }

    MpmSimulationDataAnalysis<T, dim> mpm_data(*this);
    T ke = mpm_data.evalTotalKineticEnergyGrid();

    ZIRAN_INFO(info, pe, " ", ke, " ", ke + pe);
    return ke + pe;
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::advanceOneTimeStep(double dt)
{
    ZIRAN_TIMER();
    ZIRAN_INFO("Advance one time step from time ", std::setw(7), step.time, " with                     dt = ", dt);

    reinitialize();

    if (symplectic) {
        forcesUpdateParticleState();
        for (auto f : before_p2g_callbacks) {
            f();
        }
        particlesToGrid();
        for (auto f : before_euler_callbacks) {
            f();
        }
        gridVelocityExplicitUpdate(dt);
        for (size_t k = 0; k < collision_objects.size(); ++k)
            if (collision_objects[k]->updateState)
                collision_objects[k]->updateState(step.time + dt, *collision_objects[k]);
        if (useTrialCollision) {
            trialCollision();
        }
        gridToParticles(dt);
    }
    else {
        for (auto f : before_p2g_callbacks) {
            f();
        }
        particlesToGrid();
        for (auto f : before_euler_callbacks) {
            f();
        }
        backwardEulerStep();
        for (size_t k = 0; k < collision_objects.size(); ++k)
            if (collision_objects[k]->updateState)
                collision_objects[k]->updateState(step.time + dt, *collision_objects[k]);
        gridToParticles(dt);
    }
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::trialCollision()
{
    ZIRAN_TIMER();

    ZIRAN_ASSERT(symplectic);
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        if (explicit_collision_nodes[g.idx])
            return;
        TV vi = g.new_v;
        TV xi = node.template cast<T>() * dx;
        AnalyticCollisionObject<T, dim>::trialCollision(collision_objects, xi, dt, vi);
        g.new_v = vi;
    });
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::buildAndPassGradVnToMpmForceHelpers()
{
    bool need_to_compute = false;
    for (auto& h : force->helpers) {
        if (h->needGradVn()) {
            need_to_compute = true;
            break;
        }
    }
    if (need_to_compute) {
        buildGradVn(scratch_gradV);
        for (auto& h : force->helpers)
            if (h->needGradVn())
                h->getGradVn(scratch_gradV);
    }
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::buildGradVn(StdVector<TM>& grad_vn)
{
    // TODO: use real rasterized grad_vn instead of grad_v from previous time step!
    auto& Xarray = particles.X.array;
    for (int i = 0; i < particles.count; ++i) {
        TM& gradVp = grad_vn[i];
        if (gradVp != gradVp) //check for NAN/unitialized data
            continue;
        gradVp = TM::Zero();
        TV Xp = Xarray[i];
        BSplineWeights<T, dim> spline(Xp, dx);
        grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
            gradVp.noalias() += g.v * dw.transpose();
        });
    }
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::performRpicVelocityFilteringAtBeginningOfTimeStep()
{
    ZIRAN_ASSERT(!mls_mpm);
    ZIRAN_ASSERT(FLIP_blend_PIC == transfer_scheme);
    ZIRAN_ASSERT(0 == flip_pic_ratio);
    ZIRAN_ASSERT(1 != ZIRAN_MPM_DEGREE);

    reinitialize();
    forcesUpdateParticleState();

    particlesToGridWithForceHelper<false, false, true>();
    num_nodes = grid.getNumNodes();
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        if (g.m != 0) {
            g.v /= g.m;
        }
        else {
            g.v = TV::Zero();
        }
    });

    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = g.v;
    });

    auto grid_array = grid.grid->Get_Array();
    auto* C = ((transfer_scheme == APIC_blend_RPIC) && interpolation_degree != 1)
        ? &particles.DataManager::get(C_name<TM>())
        : nullptr;
    const bool USE_APIC_BLEND_RPIC = false;
    const bool USE_MPM_DEGREE_ONE = false;
    const bool USE_MLS_MPM = false;
    tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
        for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {

            int i = particle_order[idx];
            TV& Xp = particles.X[i];

            TV picV = TV::Zero();
            BSplineWeights<T, dim> spline(Xp, dx);
            if constexpr (!USE_APIC_BLEND_RPIC) {
                TV oldV = TV::Zero();
                TM& gradVp = scratch_gradV[i];
                gradVp = TM::Zero();
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    picV += w * g.new_v;
                    oldV += w * g.v;
                    gradVp.noalias() += g.new_v * dw.transpose();
                });
                particles.V[i] *= flip_pic_ratio;
                particles.V[i] += picV - flip_pic_ratio * oldV;
            }
            else if constexpr (!USE_MPM_DEGREE_ONE) {
                TM Bp = TM::Zero();
                TM& gradVp = scratch_gradV[i];
                gradVp = TM::Zero();
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    picV += w * g.new_v;
                    TV xi_minus_xp = node.template cast<T>() * dx - Xp;
                    Bp.noalias() += w * g.new_v * xi_minus_xp.transpose();
                    // !mls_mpm
                    if constexpr (!USE_MLS_MPM) {
                        gradVp.noalias() += g.new_v * dw.transpose();
                    }
                });
                particles.V[i] = picV;
                TM CC = Bp * D_inverse;
                (*C)[i] = ((apic_rpic_ratio + 1) * (T)0.5) * CC + ((apic_rpic_ratio - 1) * (T)0.5) * CC.transpose();
                // no check for NAN/unitialized data
                if constexpr (USE_MLS_MPM) {
                    gradVp = (*C)[i];
                }
            }
            else {
                TM& gradVp = scratch_gradV[i];
                gradVp = TM::Zero();
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    picV += w * g.new_v;
                    gradVp.noalias() += g.new_v * dw.transpose();
                });
                particles.V[i] = picV;
            }
        }
    });
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::particlesToGrid()
{
    if (symplectic) {
        ZIRAN_TIMER();
        if (mls_mpm) {
            if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE == 1)
                particlesToGridWithForceHelper<false, true, true>();
            else if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE != 1)
                particlesToGridWithForceHelper<false, false, true>();
            else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE == 1)
                particlesToGridWithForceHelper<true, true, true>();
            else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE != 1)
                particlesToGridWithForceHelper<true, false, true>();
        }
        else {
            if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE == 1)
                particlesToGridWithForceHelper<false, true, false>();
            else if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE != 1)
                particlesToGridWithForceHelper<false, false, false>();
            else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE == 1)
                particlesToGridWithForceHelper<true, true, false>();
            else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE != 1)
                particlesToGridWithForceHelper<true, false, false>();
        }
    }
    else {
        ZIRAN_TIMER();
        if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE == 1)
            particlesToGridHelper<false, true>();
        else if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE != 1)
            particlesToGridHelper<false, false>();
        else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE == 1)
            particlesToGridHelper<true, true>();
        else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE != 1)
            particlesToGridHelper<true, false>();
    }
    {
        ZIRAN_TIMER();
        num_nodes = grid.getNumNodes();
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            if (g.m != 0) {
                g.v /= g.m;
            }
            else {
                g.v = TV::Zero();
            }
        });
    }
}

template <class T, int dim>
template <bool USE_APIC_BLEND_RPIC, bool USE_MPM_DEGREE_ONE, bool USE_MLS_MPM>
void MpmSimulationBase<T, dim>::particlesToGridWithForceHelper()
{
    auto& Xarray = particles.X.array;
    auto& Varray = particles.V.array;
    auto& marray = particles.mass.array;
    auto grid_array = grid.grid->Get_Array();
    const StdVector<Matrix<T, dim, dim>>* Carray_pointer;
    Carray_pointer = ((transfer_scheme == APIC_blend_RPIC) && interpolation_degree != 1) ? (&(particles.DataManager::get(C_name<TM>()).array)) : NULL;

    for (uint64_t color = 0; color < (1 << dim); ++color) {
        tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
            if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                return;
            for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                int i = particle_order[idx];
                TV& Xp = Xarray[i];
                T mass = marray[i];
                TV momentum = marray[i] * Varray[i];

                TM C = TM::Zero();
                if constexpr (USE_APIC_BLEND_RPIC) {
                    if constexpr (!USE_MPM_DEGREE_ONE) {
                        C = mass * (*Carray_pointer)[i];
                    }
                    else {
                        const TM& gradVp = scratch_gradV[i];
                        C = mass * (((apic_rpic_ratio + 1) * (T)0.5) * gradVp + ((apic_rpic_ratio - 1) * (T)0.5) * gradVp.transpose());
                    }
                }
                TM4 stress_density = TM4::Zero();
                TM& stress = scratch_stress[i];
                stress_density.template topLeftCorner<dim, dim>() = stress; // stress density
                const TV& fp = scratch_fp.col(i);
                stress_density.template topRightCorner<dim, 1>() = -fp; // fp (for meshed forces).
                TM4 velocity_density = TM4::Zero();
                velocity_density.template topLeftCorner<dim, dim>() = C;
                velocity_density.template topRightCorner<dim, 1>() = momentum;
                velocity_density(3, 3) = mass;
                BSplineWeights<T, dim> spline(Xp, dx);
                if constexpr (USE_MLS_MPM) {
                    grid.iterateKernel(spline, particle_base_offset[i],
                        [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                            TV4 xi_minus_xp = TV4::Zero();
                            xi_minus_xp.template topLeftCorner<dim, 1>() = node.template cast<T>() * dx - Xp; // top dim entries non-zero
                            xi_minus_xp(3) = 1;
                            TV4 force_delta = ((stress_density * xi_minus_xp) * (w * D_inverse));
                            TV4 velocity_delta = velocity_density * xi_minus_xp * w;
                            g.m += velocity_delta(3);
                            g.v += velocity_delta.template topLeftCorner<dim, 1>();
                            g.new_v -= force_delta.template topLeftCorner<dim, 1>(); // fi -= \sum_p (Ap (xi-xp)  - fp )w_ip Dp_inv
                        });
                }
                else {
                    grid.iterateKernel(spline, particle_base_offset[i],
                        [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                            TV4 weight = TV4::Zero();
                            weight.template topLeftCorner<dim, 1>() = dw;
                            weight(3) = w;
                            TV4 xi_minus_xp = TV4::Zero();
                            xi_minus_xp.template topLeftCorner<dim, 1>() = node.template cast<T>() * dx - Xp; // top dim entries non-zero
                            xi_minus_xp(3) = 1;
                            TV4 force_delta = (stress_density * weight);
                            TV4 velocity_delta = velocity_density * xi_minus_xp * w;
                            g.m += velocity_delta(3);
                            g.v += velocity_delta.template topLeftCorner<dim, 1>();
                            g.new_v -= force_delta.template topLeftCorner<dim, 1>();
                        });
                }
            }
        });
    }
}

template <class T, int dim>
template <bool USE_APIC_BLEND_RPIC, bool USE_MPM_DEGREE_ONE>
void MpmSimulationBase<T, dim>::particlesToGridHelper()
{
    auto& Xarray = particles.X.array;
    auto& Varray = particles.V.array;
    auto& marray = particles.mass.array;
    const StdVector<Matrix<T, dim, dim>>* Carray_pointer;
    Carray_pointer = ((transfer_scheme == APIC_blend_RPIC) && interpolation_degree != 1) ? (&(particles.DataManager::get(C_name<TM>()).array)) : NULL;

    for (uint64_t color = 0; color < (1 << dim); ++color) {
        tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
            if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                return;
            for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                int i = particle_order[idx];
                TV& Xp = Xarray[i];
                T mass = marray[i];
                TV momentum = marray[i] * Varray[i];
                TM C = TM::Zero();
                if constexpr (USE_APIC_BLEND_RPIC) {
                    if constexpr (!USE_MPM_DEGREE_ONE) {
                        C = mass * (*Carray_pointer)[i];
                    }
                    else {
                        const TM& gradVp = scratch_gradV[i];
                        C = mass * (((apic_rpic_ratio + 1) * (T)0.5) * gradVp + ((apic_rpic_ratio - 1) * (T)0.5) * gradVp.transpose());
                    }
                }
                TM4 velocity_density = TM4::Zero();
                velocity_density.template topLeftCorner<dim, dim>() = C;
                velocity_density.template topRightCorner<dim, 1>() = momentum;
                velocity_density(3, 3) = mass;
                BSplineWeights<T, dim> spline(Xp, dx);
                grid.iterateKernel(spline, particle_base_offset[i], [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                    TV4 xi_minus_xp = TV4::Zero();
                    xi_minus_xp.template topLeftCorner<dim, 1>() = node.template cast<T>() * dx - Xp; // top dim entries non-zero
                    xi_minus_xp(3) = 1;
                    TV4 velocity_delta = velocity_density * xi_minus_xp * w;
                    g.m += velocity_delta(3);
                    g.v += velocity_delta.template topLeftCorner<dim, 1>();
                });
            }
        });
    }
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::startBackwardEuler()
{
    ZIRAN_TIMER();
    ZIRAN_ASSERT(objective.matrix_free, "MPM only supports matrix free");

    buildMassMatrix();
    objective.setPreconditioner(([&](const TVStack& in, TVStack& out) {
        for (int i = 0; i < num_nodes; i++) {
            for (int d = 0; d < dim; d++) {
                out(d, i) = in(d, i) / mass_matrix(i);
            }
        }
    }));

    dv.resize(dim, num_nodes);
    vn.resize(dim, num_nodes);

    buildInitialDvAndVnForNewton(); // This also builds collision_nodes
    // TODO: this should probably be done in objective reinitialize.
    // Which should be called at the beginning of newton.
    force->backupStrain();
    objective.reinitialize(); // Reinitialize matrix sparsity pattern
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::backwardEulerStep()
{
    ZIRAN_TIMER();
    startBackwardEuler();
    for (auto f : general_callbacks) {
        f(frame, BeforeNewtonSolve);
    }
    // run diff_test if we specify it
    if (diff_test)
        runDiffTest<T, dim, Objective>(num_nodes, dv, objective, diff_test_perturbation_scale);
    newton.solve(dv, verbose);
    for (auto f : general_callbacks) {
        f(frame, AfterNewtonSolve);
    }
    // run diff_test if we specify it
    if (diff_test)
        runDiffTest<T, dim, Objective>(num_nodes, dv, objective, diff_test_perturbation_scale);

    force->restoreStrain();

    constructNewVelocityFromNewtonResult();
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::addSurfaceTension()
{
    if (!st_tau)
        return;
    auto& fp = scratch_fp;
    auto& Xarray = particles.X.array;
    auto& marray = particles.mass.array;
    tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
        for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
            int i = particle_order[idx];
            TV& Xp = Xarray[i];
            T& mass = marray[i];
            BSplineWeights<T, dim> spline(Xp, dx);
            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                TV xp_minus_xi = Xp - node.template cast<T>() * dx;
                fp.col(i) += -st_tau * (g.m - mass * w) * w * xp_minus_xi;
            });
        }
    });
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::forcesUpdateParticleState()
{
    ZIRAN_TIMER();

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
            for (auto& h : force->helpers)
                h->updateState(subrange, vtau, fp); // this adds to vtau and fp
        });

    addSurfaceTension();

    // Lagrangian
    scene.addScaledForces(1, fp); // add forces to fp
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::moveNodes(const TVStack& dv_in)
{
    ZIRAN_QUIET_TIMER();
    if (dv.data() == dv_in.data())
        return;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_nodes),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(), i_end = range.end(); i < i_end; ++i) {
                dv.col(i) = dv_in.col(i);
            }
        });
}

/**
       Callback to write the simulation state

       Should be overridden to write the simulation state
    */
template <class T, int dim>
void MpmSimulationBase<T, dim>::writeState(std::ostream& out)
{
    ZIRAN_TIMER();
    ZIRAN_INFO("Write state called");

    if (write_partio) {
        std::string obj_filename = output_dir.absolutePath(outputFileName("partio", ".bgeo"));
        writePartio(obj_filename, particles);
    }

    scene.writeState(out, outputFileName("boundary", ".obj"), outputFileName("segmesh", ".poly"), output_dir, binary_ver, write_meshes);

    if (transfer_scheme == APIC_blend_RPIC && interpolation_degree == 1)
        writeSTDVector(out, scratch_gradV); // apic C matrix
}

/**
       Callback to read the simulation state

       Should be overridden to read the output of writeState
    */
template <class T, int dim>
void MpmSimulationBase<T, dim>::readState(std::istream& in)
{
    ZIRAN_QUIET_TIMER();
    ZIRAN_INFO("Read state called");
    scene.readState(in, binary_ver);

    if (transfer_scheme == APIC_blend_RPIC && interpolation_degree == 1)
        readSTDVector(in, scratch_gradV); // apic C matrix
}

/**
       Callback to calculate dt based on CFL
    */
template <class T, int dim>
double MpmSimulationBase<T, dim>::calculateDt()
{
    ZIRAN_QUIET_TIMER();
    TV p_min_corner, p_max_corner;
    Eigen::Array<T, 2, 1> max_speed = evalMaxParticleSpeed(p_min_corner, p_max_corner);
    ZIRAN_INFO("Particle Min Corner: ", p_min_corner.transpose());
    ZIRAN_INFO("Particle Max Corner: ", p_max_corner.transpose());
    p_min_corner.array() -= (interpolation_degree + 2) * dx; // expand by influenced region and buffer
    p_max_corner.array() += (interpolation_degree + 2) * dx;
    ZIRAN_INFO("Maximum particle speed: ", max_speed(0));
    ZIRAN_ASSERT(max_speed(0) == max_speed(1));

    T max_collision_object_speed = 0;
    for (size_t k = 0; k < collision_objects.size(); ++k)
        max_collision_object_speed = std::max(collision_objects[k]->evalMaxSpeed(p_min_corner, p_max_corner), max_collision_object_speed);
    ZIRAN_INFO("Maximum collision object speed: ", max_collision_object_speed);
    max_speed(1) = std::max(max_speed(1), max_collision_object_speed);

    double dt_compute = step.max_dt;
    if (max_speed(1))
        dt_compute = cfl * dx / max_speed(1);
    ZIRAN_INFO("Computed next dt (based on CFL): ", dt_compute);
    return dt_compute;
}

// Build and mass_matrix (for implicit)
template <class T, int dim>
void MpmSimulationBase<T, dim>::buildMassMatrix()
{
    ZIRAN_QUIET_TIMER();
    mass_matrix.resize(num_nodes, 1);

    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        mass_matrix(g.idx) = g.m;
    });
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::addScaledForces(const T scale, TVStack& f) // called BE Newton objective (for computing residual)
{
    for (auto& lf : forces)
        lf->addScaledForces(scale, f);
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::addScaledForceDifferentials(const T scale, const TVStack& x, TVStack& f) // called BE Newton objective (for computing residual)
{
    for (auto& lf : forces)
        lf->addScaledForceDifferential(scale, x, f);
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::gridVelocityExplicitUpdate(double dt)
{
    ZIRAN_TIMER();
    assert(symplectic);

    explicit_collision_nodes.resize((int)num_nodes, 0);

    TV dv_gravity = gravity * dt;
    if (explicit_velocity_field) {
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            T mi = g.m;
            TV xi = node.template cast<T>() * dx;
            TV vi = g.v;
            explicit_velocity_field(mi, xi, vi);
            g.new_v = vi;
        });
    }
    else {
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            T mi = g.m;
            T dt_over_mi = (T)dt / mi;
            TV xi = node.template cast<T>() * dx;
            TV old_v = g.v;
            TV force_value = g.new_v;
            TV dvijk = force_value * dt_over_mi;
            TV vi = old_v + dv_gravity + dvijk;
            TM normal_basis;

            // apply external body force field: fext(x,t)
            for (auto ff : fext) {
                TV external_force = ff(xi, (T)(step.time));
                vi += dt * external_force;
            }
            auto old_vi = vi;

            bool collided = AnalyticCollisionObject<T, dim>::
                multiObjectCollision(collision_objects, xi, vi, normal_basis);

            explicit_collision_nodes[g.idx] = (int)collided;

            if (ignoreCollisionObject)
                vi = old_vi;

            g.new_v = vi;
        });
    }
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::constructNewVelocityFromNewtonResult()
{
    ZIRAN_QUIET_TIMER();
    grid.iterateTouchedGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = TV::Zero();
    });
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = g.v + dv.col(g.idx);
    });
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::gridToParticles(double dt)
{
    if (mls_mpm) {
        ZIRAN_TIMER();
        if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE == 1)
            gridToParticlesHelper<false, true, true>(dt);
        else if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE != 1)
            gridToParticlesHelper<false, false, true>(dt);
        else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE == 1)
            gridToParticlesHelper<true, true, true>(dt);
        else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE != 1)
            gridToParticlesHelper<true, false, true>(dt);
    }
    else {
        ZIRAN_TIMER();
        if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE == 1)
            gridToParticlesHelper<false, true, false>(dt);
        else if (transfer_scheme == FLIP_blend_PIC && ZIRAN_MPM_DEGREE != 1)
            gridToParticlesHelper<false, false, false>(dt);
        else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE == 1)
            gridToParticlesHelper<true, true, false>(dt);
        else if (transfer_scheme == APIC_blend_RPIC && ZIRAN_MPM_DEGREE != 1)
            gridToParticlesHelper<true, false, false>(dt);
    }
}

template <class T, int dim>
template <bool USE_APIC_BLEND_RPIC, bool USE_MPM_DEGREE_ONE, bool USE_MLS_MPM>
void MpmSimulationBase<T, dim>::gridToParticlesHelper(double dt)
{
    std::atomic<bool> faster_than_grid_cell(false);
    std::atomic<bool> faster_than_half_grid_cell(false);
    auto grid_array = grid.grid->Get_Array();
    auto* C = ((transfer_scheme == APIC_blend_RPIC) && interpolation_degree != 1)
        ? &particles.DataManager::get(C_name<TM>())
        : nullptr;

    tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
        for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
            bool local_faster_than_grid_cell = false;
            bool local_faster_than_half_grid_cell = false;

            int i = particle_order[idx];
            TV& Xp = particles.X[i];

            TV picV = TV::Zero();
            BSplineWeights<T, dim> spline(Xp, dx);
            if constexpr (!USE_APIC_BLEND_RPIC) {
                TV oldV = TV::Zero();
                TM& gradVp = scratch_gradV[i];
                gradVp = TM::Zero();
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    picV += w * g.new_v;
                    oldV += w * g.v;
                    gradVp.noalias() += g.new_v * dw.transpose();
                });
                particles.V[i] *= flip_pic_ratio;
                particles.V[i] += picV - flip_pic_ratio * oldV;
            }
            else if constexpr (!USE_MPM_DEGREE_ONE) {
                TM Bp = TM::Zero();
                TM& gradVp = scratch_gradV[i];
                gradVp = TM::Zero();
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    picV += w * g.new_v;
                    TV xi_minus_xp = node.template cast<T>() * dx - Xp;
                    Bp.noalias() += w * g.new_v * xi_minus_xp.transpose();
                    // !mls_mpm
                    if constexpr (!USE_MLS_MPM) {
                        gradVp.noalias() += g.new_v * dw.transpose();
                    }
                });
                particles.V[i] = picV;
                TM CC = Bp * D_inverse;
                (*C)[i] = ((apic_rpic_ratio + 1) * (T)0.5) * CC + ((apic_rpic_ratio - 1) * (T)0.5) * CC.transpose();
                // no check for NAN/unitialized data
                if constexpr (USE_MLS_MPM) {
                    gradVp = (*C)[i];
                }
            }
            else {
                TM& gradVp = scratch_gradV[i];
                gradVp = TM::Zero();
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    picV += w * g.new_v;
                    gradVp.noalias() += g.new_v * dw.transpose();
                });
                particles.V[i] = picV;
            }

            TV increment = dt * picV;
            particles.X[i] += increment; // particle position update is the same for all transfer schemes
            T inc = increment.squaredNorm();
            T dx2 = dx * dx;
            local_faster_than_grid_cell = local_faster_than_grid_cell + (inc > dx2);
            local_faster_than_half_grid_cell = local_faster_than_half_grid_cell + (inc > dx2 * (T)0.25 * (cfl * cfl));
            if (local_faster_than_half_grid_cell)
                faster_than_half_grid_cell.store(true, std::memory_order_relaxed);
            if (local_faster_than_grid_cell)
                faster_than_grid_cell.store(true, std::memory_order_relaxed);
        }
    });

    static int latest_frame_restarted = 0;
    if (faster_than_grid_cell && autorestart) {
        ZIRAN_WARN("Particle traveling more than a grid cell detected");
        int frame_to_restart = frame - 1;
        if (step.max_dt <= step.min_dt) {
            ZIRAN_WARN("Unable to shrink dt further");
            frame_to_restart--;
            ZIRAN_ASSERT(frame_to_restart >= start_frame, "Unstable intial conditions detected, shrink min_dt or change other parameters");
        }
        else {
            double new_max_dt = std::max(step.max_dt / 2, step.min_dt);
            ZIRAN_WARN("Shrinking max dt to ", new_max_dt);

            restart_callbacks.emplace_back([new_max_dt, this](int frame) {
                step.max_dt = new_max_dt;
            });
        }
        latest_frame_restarted = std::max(frame_to_restart, latest_frame_restarted);
        throw RestartException(frame_to_restart);
    }

    if (!faster_than_half_grid_cell && frame > latest_frame_restarted + 1) {
        double new_max_dt = std::min(step.max_dt_original, step.max_dt * (double)1);
        if (new_max_dt != step.max_dt) {
            ZIRAN_WARN("All particles traveled less than half a grid cell");
            step.max_dt = new_max_dt;
            ZIRAN_WARN("Increasing max dt to ", step.max_dt);
        }
    }

    force->evolveStrain(dt);

    applyPlasticity();
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::applyPlasticity()
{
    // Split all particles to subranges.
    // In parallel, each subrange goes through all plasticity appliers in serial.
    ZIRAN_QUIET_TIMER();
    auto ranges = particles.X.ranges;
    tbb::parallel_for(ranges,
        [&](DisjointRanges& subrange) {
            for (auto& p : plasticity_appliers)
                p->applyPlasticity(subrange, particles);
        });
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::sortParticlesAndPolluteGrid()
{
    auto& Xarray = particles.X.array;

    constexpr int index_bits = (32 - SparseMask::block_bits);
    ZIRAN_ASSERT(particles.count < (1 << index_bits));

    if (particles.count != (int)particle_sorter.size()) {
        particle_base_offset.resize(particles.count);
        particle_sorter.resize(particles.count);
        particle_order.resize(particles.count);
    }

    T one_over_dx = (T)1 / dx;
    tbb::parallel_for(0, particles.count, [&](int i) {
        uint64_t offset = SparseMask::Linear_Offset(to_std_array(
            baseNode<interpolation_degree, T, dim>(Xarray[i] * one_over_dx)));
        particle_sorter[i] = ((offset >> SparseMask::data_bits) << index_bits) + i;
    });

    tbb::parallel_sort(particle_sorter.begin(), particle_sorter.end());

    particle_group.clear();
    block_offset.clear();
    int last_index = 0;
    for (int i = 0; i < particles.count; ++i)
        if (i == particles.count - 1 || (particle_sorter[i] >> 32) != (particle_sorter[i + 1] >> 32)) {
            particle_group.push_back(std::make_pair(last_index, i));
            block_offset.push_back(particle_sorter[i] >> 32);
            last_index = i + 1;
        }

    grid.page_map->Clear();
    for (int i = 0; i < particles.count; ++i) {
        particle_order[i] = (int)(particle_sorter[i] & ((1ll << index_bits) - 1));
        uint64_t offset = (particle_sorter[i] >> index_bits) << SparseMask::data_bits;
        particle_base_offset[particle_order[i]] = offset;
        if (i == particles.count - 1 || (particle_sorter[i] >> 32) != (particle_sorter[i + 1] >> 32)) {
            grid.page_map->Set_Page(offset);
            if constexpr (dim == 2) {
                auto x = 1 << SparseMask::block_xbits;
                auto y = 1 << SparseMask::block_ybits;
                for (int i = 0; i < 2; ++i)
                    for (int j = 0; j < 2; ++j)
                        grid.page_map->Set_Page(SparseMask::Packed_Add(
                            offset, SparseMask::Linear_Offset(x * i, y * j)));
            }
            else {
                auto x = 1 << SparseMask::block_xbits;
                auto y = 1 << SparseMask::block_ybits;
                auto z = 1 << SparseMask::block_zbits;
                for (int i = 0; i < 2; ++i)
                    for (int j = 0; j < 2; ++j)
                        for (int k = 0; k < 2; ++k)
                            grid.page_map->Set_Page(SparseMask::Packed_Add(
                                offset, SparseMask::Linear_Offset(x * i, y * j, z * k)));
            }
        }
    }
    grid.page_map->Update_Block_Offsets();

    auto grid_array = grid.grid->Get_Array();
    auto blocks = grid.page_map->Get_Blocks();
    for (int b = 0; b < (int)blocks.second; ++b) {
        auto base_offset = blocks.first[b];
        std::memset(&grid_array(base_offset), 0, (size_t)(1 << MpmGrid<T, dim>::log2_page));
        GridState<T, dim>* g = reinterpret_cast<GridState<T, dim>*>(&grid_array(base_offset));
        for (int i = 0; i < (int)SparseMask::elements_per_block; ++i)
            g[i].idx = -1;
    }
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::buildInitialDvAndVnForNewton()
{
    ZIRAN_QUIET_TIMER();
    assert(!symplectic);

    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        TV old_v = g.v;
        TV vi = old_v;

        const TV xi = node.template cast<T>() * dx;
        TM normal_basis;

        bool any_collision = AnalyticCollisionObject<T, dim>::
            multiObjectCollision(collision_objects, xi, vi, normal_basis);

        int node_id = g.idx;
        assert(node_id < num_nodes);
        if (any_collision) {
            CollisionNode<T, dim> Z{ node_id, TM::Identity() - normal_basis * normal_basis.transpose() };
            collision_nodes.push_back(Z);

            // Newton initial guess for collided nodes.
            // Setting it this way so that the new velocity constructed after the solve agrees with collision.
            dv.col(node_id) = vi - old_v;
        }
        else
            dv.col(node_id) = gravity * dt; // Newton initial guess for non-collided nodes

        vn.col(node_id) = old_v;
    });
}

/**
       result(0) is purely particle speed maximum.
       result(1) is APIC enhanced 'speed' maximum.
    */
template <class T, int dim>
Eigen::Array<T, 2, 1> MpmSimulationBase<T, dim>::evalMaxParticleSpeed(TV& min_corner, TV& max_corner)
{
    ZIRAN_TIMER();
    using TArray = Eigen::Array<T, 2 * dim + 2, 1>;
    TArray init = TArray::Zero();
    T Tmin = std::numeric_limits<float>::lowest();
    TV TVmin = TV::Zero();
    TVmin.array() = Tmin;
    init << (T)0, (T)0, TVmin, TVmin;
    TArray result;
    result = particles.map_reduce(
        init,
        [](const TV& position, const TV& velocity) {
            T velocity_norm = velocity.norm();
            TArray r;
            r << velocity_norm, velocity_norm, position, -position;
            return r;
        },
        [](const TArray& a, const TArray& b) {
            return a.max(b);
        },
        particles.X_name(), particles.V_name());

    max_corner = result.template segment<dim>(2);
    min_corner = -result.template tail<dim>();
    return result.template head<2>();
}

template <class T, int dim>
void MpmSimulationBase<T, dim>::writeSimulationInformation()
{
    ZIRAN_INFO("Particle count:      ", particles.count);
    ZIRAN_INFO("Grid dx:             ", dx);
    ZIRAN_INFO("Interp degree:       ", interpolation_degree);
    ZIRAN_INFO("Gravity:             ", gravity.transpose());
    ZIRAN_INFO("");
    ZIRAN_INFO("CFL:                 ", cfl);
    ZIRAN_INFO("Transfer scheme:     ", transfer_scheme);
    ZIRAN_INFO("Print stats:         ", print_stats);
    ZIRAN_INFO("Using mls_mpm:       ", mls_mpm);
    writeParticlePerCellHistogram();
}

// Called by writeSimulationInformation
template <class T, int dim>
void MpmSimulationBase<T, dim>::writeParticlePerCellHistogram()
{
    typedef tbb::concurrent_unordered_map<uint64_t, size_t> ConcurrentHash;
    ConcurrentHash cellParticleCount; // mapping IV to size_t
    AttributeName<Vector<T, dim>> x_name = Particles<T, dim>::X_name();
    for (auto iter = particles.iter(x_name); iter; ++iter) {
        Vector<T, dim>& X = iter.template get<0>();
        Vector<int, dim> IX = (X / dx).template cast<int>();
        uint64_t index = SparseMask::Linear_Offset(to_std_array(IX));
        auto search = cellParticleCount.find(index);
        if (search == cellParticleCount.end())
            cellParticleCount[index] = 1;
        else
            cellParticleCount[index]++;
    }
    StdVector<size_t> number_of_cells_with_particles(101, 0);
    for (auto& iter : cellParticleCount)
        number_of_cells_with_particles[std::min(iter.second, (size_t)100)]++;

    for (int i = 0; i < 101; i++) {
        int count = number_of_cells_with_particles[i];
        if (count > 0 && i != 100)
            ZIRAN_INFO("Cells with ", std::setw(2), i, " particles：", count);
        else if (count > 0 && i == 100)
            ZIRAN_INFO("Cells with 100 or more particles: ", count);
    }
}

template <class T, int dim>
constexpr int MpmSimulationBase<T, dim>::interpolation_degree;
} // namespace ZIRAN
