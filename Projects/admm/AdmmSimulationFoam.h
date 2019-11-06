#ifndef ADMM_SIMULATION_FOAM_H
#define ADMM_SIMULATION_FOAM_H

#include <Ziran/CS/Util/AttributeNamesForward.h>
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/Physics/ConstitutiveModel/StvkWithHencky.h>
#include <Ziran/Physics/PlasticityApplier.h>
#include <Ziran/Math/Linear/Minres.h>
#include <MPM/MpmSimulationBase.h>
#include <Eigen/LU>
#include <Partio.h>
#include <map>
#include <mutex>

#include "AdmmFormula.h"
#include "AdmmSketch.h"
#include "GlobalSystem.h"
//#include "ConjugateGradient.h"
#include <Ziran/Math/Linear/Minres.h>
#include <Ziran/Math/Geometry/CollisionObject.h>

namespace ZIRAN {

template <class T, int dim>
class AdmmSimulation : public MpmSimulationBase<T, dim> {
public:
    using Base = MpmSimulationBase<T, dim>;

    typedef Vector<T, dim> TV;
    typedef Vector<int, dim> IV;
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, Eigen::Dynamic> TStack;
    typedef Matrix<T, dim, Eigen::Dynamic> TVStack;

    using Base::block_offset;
    using Base::collision_objects;
    using Base::dt;
    using Base::dv;
    using Base::dx;
    using Base::grid;
    using Base::mass_matrix;
    using Base::particle_base_offset;
    using Base::particle_group;
    using Base::particle_order;
    using Base::particles;
    using Base::step;
    using Base::vn;
    using typename Base::SparseMask;

    StdVector<TM>* Fn_pointer;
    StdVector<TM>* FFn_pointer;
    StdVector<TM> F;
    StdVector<TM> y;
    StdVector<T> omega;
    T rho_scale = 0.5;
    T over_relaxation_alpha = 1;

    using Simulation = AdmmSimulation<T, dim>;
    using Objective = GlobalSystem<Simulation, dim>;

    Objective cg_objective;
    Minres<T, Objective, TVStack> cg;

    T out_a = 0;
    int total_cg_iterations = 0;
    std::mutex mtx;

    bool use_ruiz = false;
    bool use_elasticity_plasticity = true;
    int admm_max_iterations = 100;
    T global_tolerance = 1e-6;
    T local_tolerance = 1e-8;
    int local_max_iteration;
    int global_max_iteration = 10000;

    static const bool build_matrix = true;

    TVStack preconditioner;
    std::vector<T> multiplier_value;
    std::vector<int> multiplier_idx;

    void visualize_particles(StdVector<TM>* pointer, std::string filename)
    {
        Partio::ParticlesDataMutable* parts = Partio::create();

        // visualize particles info
        Partio::ParticleAttribute posH, F00H, F01H, F10H, F11H;
        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        F00H = parts->addAttribute("F00", Partio::VECTOR, 1);
        F01H = parts->addAttribute("F01", Partio::VECTOR, 1);
        F10H = parts->addAttribute("F10", Partio::VECTOR, 1);
        F11H = parts->addAttribute("F11", Partio::VECTOR, 1);

        for (int k = 0; k < particles.count; k++) {
            int idx = parts->addParticle();
            float* posP = parts->dataWrite<float>(posH, idx);
            float* F00P = parts->dataWrite<float>(F00H, idx);
            float* F01P = parts->dataWrite<float>(F01H, idx);
            float* F10P = parts->dataWrite<float>(F10H, idx);
            float* F11P = parts->dataWrite<float>(F11H, idx);
            for (int d = 0; d < 3; ++d) posP[d] = 0;
            for (int d = 0; d < dim; ++d) posP[d] = particles.X.array[k](d);
            F00P[0] = (*pointer)[k](0, 0);
            F01P[0] = (*pointer)[k](0, 1);
            F10P[0] = (*pointer)[k](1, 0);
            F11P[0] = (*pointer)[k](1, 1);
        }

        Partio::write(filename.c_str(), *parts);
        parts->release();
    }

    void writeState(std::ostream& out)
    {
        Base::writeState(out);
        return;
        {
            std::string filename = SimulationBase::output_dir.absolutePath(SimulationBase::outputFileName("F_e", ".bgeo"));
            visualize_particles(Fn_pointer, filename);
        }
        {
            std::string filename = SimulationBase::output_dir.absolutePath(SimulationBase::outputFileName("F_ne", ".bgeo"));
            visualize_particles(FFn_pointer, filename);
        }
    }

    AdmmSimulation()
        : Base(), cg(10000)
    {
        auto& Xarray = particles.X.array;
        auto ranges = particles.X.ranges;

        if constexpr (build_matrix) {
            cg_objective.setMultiplier([&](const TVStack& x, TVStack& b) {
                b.setZero();
                int block_size = (dim == 2 ? 25 : 125);
                grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                    T* _value = &multiplier_value[g.idx * block_size];
                    int* _idx = &multiplier_idx[g.idx * block_size];
                    TV ret = g.m * x.col(g.idx);
                    for (int i = 0; i < block_size; ++i)
                        ret += x.col(_idx[i]) * _value[i];
                    b.col(g.idx) = ret;
                });
            });
            cg_objective.setPreconditioner([&](const TVStack& in, TVStack& out) {
                out = in.cwiseQuotient(preconditioner);
            });
        }
        else {
            ZIRAN_ASSERT(false, "moved into AdmmSketch.h");
        }
        cg_objective.setProjection([&](TVStack& r) {
            for (auto iter = Base::collision_nodes.begin(); iter != Base::collision_nodes.end(); ++iter) {
                int node_id = (*iter).node_id;
                TV tmp = r.col(node_id);
                (*iter).project(tmp);
                r.col(node_id) = tmp;
            }
        });
    }

    void dvStep(T tolerance, T relative_tolerance)
    {
        ZIRAN_TIMER();
        auto& Xarray = particles.X.array;

        TVStack rhs = TVStack::Zero(dim, Base::num_nodes);
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            rhs.col(g.idx) = dt * g.m * Base::gravity;
        });

        // G2P
        StdVector<TM> bFu(particles.count * 2, TM::Zero());
        tbb::parallel_for(0, (int)particles.count, [&](int i) {
            TV& Xp = Xarray[i];
            TM& Fn = (*Fn_pointer)[i];
            bFu[i] = F[i] - Fn + y[i] / (omega[i] * rho_scale);
            BSplineWeights<T, dim> spline(Xp, dx);
            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                bFu[i] -= dt * g.v * dw.transpose() * Fn;
            });
            if (enable_visco) {
                int pc = particles.count;
                TM& FFn = (*FFn_pointer)[i];
                bFu[i + pc] = F[i + pc] - FFn + y[i + pc] / (omega[i + pc] * rho_scale);
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    bFu[i + pc] -= dt * g.v * dw.transpose() * FFn;
                });
            }
        });
        // P2G
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = TV::Zero();
        });
        Base::parallel_for_updating_grid([&](int i) {
            TV& Xp = Xarray[i];
            TM& Fn = (*Fn_pointer)[i];
            BSplineWeights<T, dim> spline(Xp, dx);
            grid.iterateKernel(spline, particle_base_offset[i], [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                g.new_v += dt * rho_scale * omega[i] * omega[i] * bFu[i] * Fn.transpose() * dw;
            });
            if (enable_visco) {
                int pc = particles.count;
                TM& FFn = (*FFn_pointer)[i];
                grid.iterateKernel(spline, particle_base_offset[i], [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                    g.new_v += dt * rho_scale * omega[i + pc] * omega[i + pc] * bFu[i + pc] * FFn.transpose() * dw;
                });
            }
        });
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            rhs.col(g.idx) += g.new_v;
        });

        cg.setTolerance(tolerance);
        cg.setRelativeTolerance(relative_tolerance);
        cg.setMaxIteration(global_max_iteration);
        total_cg_iterations += cg.solve(cg_objective, dv, rhs, false);
        // outputFile(out_a += step.max_dt / 100.0, total_cg_iterations);
    }

    template <class TCONST, class TPCONST>
    void initOmega()
    {
        ZIRAN_TIMER();
        auto ranges = particles.X.ranges;
        auto constitutive_model_name = AttributeName<TCONST>(TCONST::name());
        auto plastic_name = AttributeName<TPCONST>(TPCONST::name());
        if (!particles.exist(constitutive_model_name) || !particles.exist(plastic_name))
            return;
        tbb::parallel_for(ranges, [&](DisjointRanges& subrange) {
            DisjointRanges subset(subrange, particles.commonRanges(constitutive_model_name, element_measure_name<T>(), F_name<T, dim>(), plastic_name));
            for (auto iter = particles.subsetIter(subset, constitutive_model_name, element_measure_name<T>(), F_name<T, dim>(), plastic_name); iter; ++iter) {
                auto& constitutive_model = iter.template get<0>();
                auto& vol = iter.template get<1>();
                auto& plasticity_model = iter.template get<3>();
                int i = iter.entryId();
                TM Fn = F[i];
                TM U, V;
                TV sigma;
                singularValueDecomposition(Fn, U, sigma, V);
                TM A = constitutive_model.firstPiolaDerivativeDiagonal(sigma);
                makePD(A);
                if constexpr (!std::is_same<TPCONST, StvkWithHencky<T, dim>>::value) {
                    sigma = plasticity_model.projectSigma(constitutive_model, sigma);
                }
                T total = A.squaredNorm();
                T clamp_value = 1e-6;
                total += (makePD2D(constitutive_model.Bij(sigma, 0, 1, clamp_value))).squaredNorm();
                if constexpr (dim == 3) {
                    total += (makePD2D(constitutive_model.Bij(sigma, 0, 2, clamp_value))).squaredNorm();
                    total += (makePD2D(constitutive_model.Bij(sigma, 1, 2, clamp_value))).squaredNorm();
                }
                // omega[i] = std::sqrt(vol * ((T)2*constitutive_model.mu/3+constitutive_model.lambda));
                omega[i] = std::sqrt(vol * std::sqrt(total));
                if constexpr (std::is_same<TPCONST, StvkWithHencky<T, dim>>::value) {
                    int pc = particles.count;
                    singularValueDecomposition(F[i + pc], U, sigma, V);
                    //  TODO
                    A = plasticity_model.firstPiolaDerivativeDiagonal(sigma);
                    makePD(A);
                    total = A.squaredNorm();
                    clamp_value = 1e-6;
                    total += (makePD2D(plasticity_model.Bij(sigma, 0, 1, clamp_value))).squaredNorm();
                    if constexpr (dim == 3) {
                        total += (makePD2D(plasticity_model.Bij(sigma, 0, 2, clamp_value))).squaredNorm();
                        total += (makePD2D(plasticity_model.Bij(sigma, 1, 2, clamp_value))).squaredNorm();
                    }
                    // omega[i + pc] = std::sqrt(vol * ((T)2*plasticity_model.mu/3+plasticity_model.lambda));
                    omega[i + pc] = std::sqrt(vol * std::sqrt(total));
                }
            }
        });
    }

    void updateCollisionNodes()
    {
        Base::collision_nodes.clear();
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            TV vi = dv.col(g.idx) + vn.col(g.idx);
            const TV xi = node.template cast<T>() * dx;
            TM normal_basis;
            bool any_collision = AnalyticCollisionObject<T, dim>::
                multiObjectCollisionWithFakeSeparate(collision_objects, xi, vi, normal_basis);
            if (any_collision) {
                CollisionNode<T, dim> Z{ (int)g.idx, TM::Identity() - normal_basis * normal_basis.transpose() };
                mtx.lock();
                Base::collision_nodes.push_back(Z);
                mtx.unlock();
                dv.col(g.idx) = vi - vn.col(g.idx);
            }
        });
    }

    StdVector<TV> last_dv;

    void initStep()
    {
        ZIRAN_TIMER();

        bool flag = (dv.cols() == 0);

        dv.resize(dim, Base::num_nodes);
        vn.resize(dim, Base::num_nodes);
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            vn.col(g.idx) = g.v;
            dv.col(g.idx) = TV::Zero();
            // dv.col(g.idx) = dt * Base::gravity;
        });
        auto& Xarray = particles.X.array;
        last_dv.resize((int)particles.count, TV::Zero());
        std::vector<T> dv_w(Base::num_nodes, 0);
        Base::parallel_for_updating_grid([&](int i) {
            TV& Xp = Xarray[i];
            BSplineWeights<T, dim> spline(Xp, dx);
            grid.iterateKernel(spline, particle_base_offset[i], [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                if (g.idx >= 0) {
                    dv.col(g.idx) += last_dv[i] * w;
                    dv_w[g.idx] += w;
                }
            });
        });
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            TV tmp = dv.col(g.idx) * dt / dv_w[g.idx];
            if ((tmp / dt).norm() < 0) {
                dv.col(g.idx) = TV::Zero();
            }
            else {
                dv.col(g.idx) = tmp;
            }
        });
        if (flag) {
            grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                dv.col(g.idx) = Base::gravity * dt * 0;
            });
        }
        updateCollisionNodes();
        // init F
        F.resize((int)particles.count * 2);
        grid.iterateTouchedGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = TV::Zero();
        });
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = g.v + dv.col(g.idx);
        });
        if (particles.count > 0) {
            Fn_pointer = &(particles.DataManager::get(F_name<T, dim>()).array);
            if (enable_visco) {
                FFn_pointer = &(particles.DataManager::get(FEN_name<T, dim>()).array);
            }
        }
        tbb::parallel_for(0, (int)particles.count, [&](int i) {
            TV& Xp = Xarray[i];
            TM& Fn = (*Fn_pointer)[i];
            F[i] = Fn;
            BSplineWeights<T, dim> spline(Xp, dx);
            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                F[i] += dt * g.new_v * dw.transpose() * Fn;
            });
            if (enable_visco) {
                int pc = particles.count;
                TM& FFn = (*FFn_pointer)[i];
                F[i + pc] = FFn;
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    F[i + pc] += dt * g.new_v * dw.transpose() * FFn;
                });
            }
        });
        // init omega
        StdVector<T> omega_old = omega;
        omega.resize((int)particles.count * 2, (T)1);
        if (use_elasticity_plasticity) {
            initOmega<StvkWithHencky<T, dim>, DummyPlasticity<T>>();
            initOmega<StvkWithHencky<T, dim>, DruckerPragerStvkHencky<T>>();
            initOmega<StvkWithHencky<T, dim>, VonMisesStvkHencky<T, dim>>();
        }
        else {
            initOmega<StvkWithHenckyIsotropic<T, dim>, StvkWithHencky<T, dim>>();
        }
        // init y
        StdVector<TM> y_old = y;
        int y_old_size = y.size();
        y.resize((int)particles.count * 2, TM::Zero());
        if (!enable_visco) {
            tbb::parallel_for(0, y_old_size / 2, [&](int i) {
                y[i] = y_old[i] * omega_old[i] / omega[i];
            });
        }
        else {
            int pc = particles.count;
            tbb::parallel_for(0, y_old_size / 2, [&](int i) {
                y[i] = y_old[i] * omega_old[i] / omega[i];
                y[i + pc] = y_old[i + y_old_size / 2] * omega_old[i + y_old_size / 2] / omega[i + pc];
            });
        }
        // init mass_matrix
        mass_matrix.resize((int)Base::num_nodes, 1);
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            mass_matrix(g.idx) = g.m;
        });
    }

    bool enable_visco = false;
    T viscosity_d;
    T viscosity_v;

    void initMultiplierAndPreconditioner()
    {
        ZIRAN_TIMER();
        if constexpr (build_matrix) {
            auto& Xarray = particles.X.array;
            // use new_v to store x
            // G2P
            int block_size = (dim == 2 ? 25 : 125);
            multiplier_value.resize((int)Base::num_nodes * block_size);
            multiplier_idx.resize((int)Base::num_nodes * block_size);
            std::fill(multiplier_value.begin(), multiplier_value.end(), 0);
            std::fill(multiplier_idx.begin(), multiplier_idx.end(), 0);
            Base::parallel_for_updating_grid([&](int i) {
                TV& Xp = Xarray[i];
                TM& Fn = (*Fn_pointer)[i];
                BSplineWeights<T, dim> spline(Xp, dx);
                std::vector<IV> cached_nodes;
                std::vector<int> cached_idx;
                std::vector<TV> cached_FTdw;
                grid.iterateKernel(spline, particle_base_offset[i], [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                    if (g.idx >= 0) {
                        cached_nodes.push_back(node);
                        cached_idx.push_back((int)g.idx);
                        cached_FTdw.push_back(Fn.transpose() * dw);
                    }
                });
                T coeff = rho_scale * omega[i] * omega[i] * dt * dt;
                int cached_size = cached_nodes.size();
                for (int p = 0; p < cached_size; ++p) {
                    T* _value = &multiplier_value[cached_idx[p] * block_size];
                    int* _idx = &multiplier_idx[cached_idx[p] * block_size];
                    for (int q = 0; q < cached_size; ++q) {
                        int offset = 0;
                        if constexpr (dim == 2) {
                            offset = 5 * (cached_nodes[q](0) - cached_nodes[p](0) + 2) + cached_nodes[q](1) - cached_nodes[p](1) + 2;
                        }
                        else {
                            offset = 25 * (cached_nodes[q](0) - cached_nodes[p](0) + 2) + 5 * (cached_nodes[q](1) - cached_nodes[p](1) + 2) + cached_nodes[q](2) - cached_nodes[p](2) + 2;
                        }
                        _idx[offset] = cached_idx[q];
                        _value[offset] += coeff * cached_FTdw[p].dot(cached_FTdw[q]);
                    }
                }
                if (enable_visco) {
                    int pc = particles.count;
                    TM& FFn = (*FFn_pointer)[i];
                    std::vector<TV> cached_FFFF;
                    grid.iterateKernel(spline, particle_base_offset[i], [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                        if (g.idx >= 0) {
                            cached_FFFF.push_back(FFn.transpose() * dw);
                        }
                    });
                    T ccccc = rho_scale * omega[i + pc] * omega[i + pc] * dt * dt;
                    for (int p = 0; p < cached_size; ++p) {
                        T* _value = &multiplier_value[cached_idx[p] * block_size];
                        for (int q = 0; q < cached_size; ++q) {
                            int offset = 0;
                            if constexpr (dim == 2) {
                                offset = 5 * (cached_nodes[q](0) - cached_nodes[p](0) + 2) + cached_nodes[q](1) - cached_nodes[p](1) + 2;
                            }
                            else {
                                offset = 25 * (cached_nodes[q](0) - cached_nodes[p](0) + 2) + 5 * (cached_nodes[q](1) - cached_nodes[p](1) + 2) + cached_nodes[q](2) - cached_nodes[p](2) + 2;
                            }
                            _value[offset] += ccccc * cached_FFFF[p].dot(cached_FFFF[q]);
                        }
                    }
                }
            });

            preconditioner.resize(dim, Base::num_nodes);
            grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                g.new_v = g.m * TV::Ones();
            });
            Base::parallel_for_updating_grid([&](int i) {
                TV& Xp = Xarray[i];
                TM& Fn = (*Fn_pointer)[i];
                BSplineWeights<T, dim> spline(Xp, dx);
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    T coeff = rho_scale * omega[i] * omega[i] * dt * dt;
                    TM tmp = coeff * (TV::Ones() * dw.transpose() * Fn);
                    TV value = (tmp * Fn.transpose() * dw);
                    g.new_v += value;
                });
                if (enable_visco) {
                    int pc = particles.count;
                    TM& FFn = (*FFn_pointer)[i];
                    grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                        T ccccc = rho_scale * omega[i + pc] * omega[i + pc] * dt * dt;
                        TM ttt = ccccc * (TV::Ones() * dw.transpose() * FFn);
                        TV value = (ttt * FFn.transpose() * dw);
                        g.new_v += value;
                    });
                }
            });
            grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                preconditioner.col(g.idx) = g.new_v;
            });
        }
    }

    template <class TCONST, class TPCONST>
    void FStepSolve(const T& localTol)
    {
        auto& Xarray = particles.X.array;
        auto constitutive_model_name = AttributeName<TCONST>(TCONST::name());
        auto plastic_name = AttributeName<TPCONST>(TPCONST::name());
        auto ranges = particles.X.ranges;
        if (!particles.exist(constitutive_model_name) || !particles.exist(plastic_name))
            return;
        tbb::parallel_for(ranges, [&](DisjointRanges& subrange) {
            DisjointRanges subset(subrange, particles.commonRanges(constitutive_model_name, plastic_name, element_measure_name<T>(), F_name<T, dim>()));
            for (auto iter = particles.subsetIter(subset, constitutive_model_name, plastic_name, element_measure_name<T>(), F_name<T, dim>()); iter; ++iter) {
                auto& c = iter.template get<0>();
                auto& p = iter.template get<1>();
                auto& vol = iter.template get<2>();
                int i = iter.entryId();
                TV& Xp = Xarray[i];
                TM& Fn = (*Fn_pointer)[i];
                // build rhs
                TM u_bar = -y[i] / (omega[i] * rho_scale);
                TM FDb = Fn + u_bar;
                BSplineWeights<T, dim> spline(Xp, dx);
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    FDb += dt * g.new_v * dw.transpose() * Fn;
                });
                // SVD rhs_raw
                TM U, V;
                TV sigma;
                singularValueDecomposition(FDb, U, sigma, V);
                TV zi = (U.transpose() * F[i] * V).diagonal();
                T coeff = rho_scale * omega[i] * omega[i];
                TV step;
                for (int j = 0; j < 5; j++) {
                    if constexpr (!std::is_same<TPCONST, StvkWithHencky<T, dim>>::value) {
                        TV zi_p;
                        TM zi_pd;
                        p.projectSigmaAndDerivative(c, zi, zi_p, zi_pd);
                        TV dE = c.firstPiolaDiagonal(zi_p);
                        TM ddE = c.firstPiolaDerivativeDiagonal(zi_p);
                        TV g = vol * dE + coeff * (zi - sigma);
                        TM P = vol * ddE * zi_pd + coeff * TM::Identity();
                        step = P.inverse() * (-g);
                        zi += step;
                        if (step.norm() < localTol)
                            break;
                    }
                    else {
                        TV dE = c.firstPiolaDiagonal(zi);
                        TM ddE = c.firstPiolaDerivativeDiagonal(zi);
                        TV g = vol * dE + coeff * (zi - sigma);
                        TM P = vol * ddE + coeff * TM::Identity();
                        step = P.inverse() * (-g);
                        zi += step;
                        if (step.norm() < localTol)
                            break;
                    }
                }
                // ZIRAN_ASSERT(step.norm() < localTol, "need larger rho_scale");
                F[i] = U * zi.asDiagonal() * V.transpose();
                if constexpr (std::is_same<TPCONST, StvkWithHencky<T, dim>>::value) {
                    int pc = particles.count;
                    TM& FFn = (*FFn_pointer)[i];
                    // build rhs
                    TM u_bar = -y[i + pc] / (omega[i + pc] * rho_scale);
                    TM FDb = FFn + u_bar;
                    BSplineWeights<T, dim> spline(Xp, dx);
                    grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                        FDb += dt * g.new_v * dw.transpose() * FFn;
                    });
                    // SVD rhs_raw
                    TM U, V;
                    TV sigma;
                    singularValueDecomposition(FDb, U, sigma, V);
                    TV zi = (U.transpose() * F[i + pc] * V).diagonal();
                    T coeff = rho_scale * omega[i + pc] * omega[i + pc];
                    TV step;
                    for (int j = 0; j < 20; j++) {
                        TV zi_p;
                        TM zi_pd;
                        {
                            T lambda = p.lambda;
                            T mu = p.mu;
                            T alpha = (T)2.0 * mu / viscosity_d;
                            T beta = (T)2.0 * ((T)2.0 * mu + lambda * dim) / ((T)9.0 * viscosity_v) - (T)2.0 * mu / (viscosity_d * dim);

                            TV epsilon_trial = zi.array().abs().log();
                            TV epsilon = (T)1.0 / ((T)1.0 + dt * alpha) * (epsilon_trial - dt * beta / ((T)1.0 + dt * (alpha + dim * beta)) * epsilon_trial.trace() * TV::Ones());
                            zi_p = epsilon.array().exp();

                            TM zi_p_m = zi_p.asDiagonal();
                            TV zi_inv = zi.array().inverse();
                            TM zi_inv_m = zi_inv.asDiagonal();
                            zi_pd = zi_p_m * ((T)1.0 / ((T)1.0 + dt * alpha)) * (zi_inv_m - dt * beta / ((T)1.0 + dt * (alpha + dim * beta)) * TV::Ones() * zi_inv.transpose());
                        }
                        TV dE = p.firstPiolaDiagonal(zi_p);
                        TM ddE = p.firstPiolaDerivativeDiagonal(zi_p);
                        TV g = vol * dE + coeff * (zi - sigma);
                        TM P = vol * ddE * zi_pd + coeff * TM::Identity();
                        step = P.inverse() * (-g);
                        zi += step;
                        if (step.norm() < localTol)
                            break;
                    }
                    // ZIRAN_ASSERT(step.norm() < localTol, "need larger rho_scale");
                    F[i + pc] = U * zi.asDiagonal() * V.transpose();
                }
            }
        });
    };

    // only work for non/ruiz
    void FStep(T localTol)
    {
        ZIRAN_TIMER();
        grid.iterateTouchedGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = TV::Zero();
        });
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = g.v + dv.col(g.idx);
        });
        local_max_iteration = 0;
        if (use_elasticity_plasticity) {
            FStepSolve<StvkWithHencky<T, dim>, DummyPlasticity<T>>(localTol);
            FStepSolve<StvkWithHencky<T, dim>, DruckerPragerStvkHencky<T>>(localTol);
            FStepSolve<StvkWithHencky<T, dim>, VonMisesStvkHencky<T, dim>>(localTol);
        }
        else {
            FStepSolve<StvkWithHenckyIsotropic<T, dim>, StvkWithHencky<T, dim>>(localTol);
        }
    }

    void uStep()
    {
        ZIRAN_TIMER();
        auto& Xarray = particles.X.array;

        // use new_v to store dv
        // G2P
        grid.iterateTouchedGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = TV::Zero();
        });
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = g.v + dv.col(g.idx);
        });
        tbb::parallel_for(0, (int)particles.count, [&](int i) {
            TV& Xp = Xarray[i];
            TM& Fn = (*Fn_pointer)[i];

            TM tmp = F[i] - Fn;
            BSplineWeights<T, dim> spline(Xp, dx);
            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                tmp -= dt * g.new_v * dw.transpose() * Fn;
            });
            y[i] += rho_scale * omega[i] * tmp;
            if (enable_visco) {
                int pc = particles.count;
                TM& FFn = (*FFn_pointer)[i];
                TM ttt = F[i + pc] - FFn;
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    ttt -= dt * g.new_v * dw.transpose() * FFn;
                });
                y[i + pc] += rho_scale * omega[i + pc] * ttt;
            }
        });
    }

    T admmResidual()
    {
        Base::force->updatePositionBasedState();
        Base::inertia->updatePositionBasedState();
        TVStack residual;
        residual.resize(dim, Base::num_nodes);
        tbb::parallel_for(tbb::blocked_range<int>(0, residual.cols(), 256),
            [&](const tbb::blocked_range<int>& range) {
                int start = range.begin();
                int length = (range.end() - range.begin());
                TV dtg = dt * Base::gravity;
                residual.middleCols(start, length) = dtg * mass_matrix.segment(start, length).transpose();
            });
        Base::addScaledForces(dt, residual);
        Base::inertia->addScaledForces(dt, residual);
        for (auto iter = Base::collision_nodes.begin(); iter != Base::collision_nodes.end(); ++iter) {
            int node_id = (*iter).node_id;
            (*iter).project(residual.col(node_id));
        }
        ZIRAN_ASSERT(residual.cols() == Base::num_nodes);
        T norm_sq = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, residual.cols(), 256), (T)0,
            [&](const tbb::blocked_range<int>& range, T ns) -> T {
                int start = range.begin();
                int length = (range.end() - range.begin());
                const auto& r_block = residual.middleCols(start, length);
                const auto& mass_block = mass_matrix.segment(start, length).transpose();
                return ns + ((r_block.colwise().squaredNorm()).array() / mass_block.array()).sum();
            },
            [](T a, T b) -> T { return a + b; });
        Base::force->restoreStrain();

        T residual_current = std::sqrt(norm_sq);
        outputResidual(out_a += step.max_dt / 200.0, residual_current);
        ZIRAN_INFO("admm residual ", residual_current, ", total cg iterations = ", total_cg_iterations, ", rho_scale = ", rho_scale);

        return residual_current;
    }

    TVStack old_dv;
    void updateRhoScaleAbsolute()
    {
        auto& Xarray = particles.X.array;
        // use new_v to store dv
        // G2P
        T prime_residual = 0;
        grid.iterateTouchedGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = TV::Zero();
        });
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = dv.col(g.idx);
        });
        tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
            for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                int i = particle_order[idx];
                TV& Xp = Xarray[i];
                TM& Fn = (*Fn_pointer)[i];
                BSplineWeights<T, dim> spline(Xp, dx);
                TM tmp = Fn - F[i];
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    tmp += dt * g.v * dw.transpose() * Fn;
                    tmp += dt * g.new_v * dw.transpose() * Fn;
                });
                mtx.lock();
                prime_residual += (tmp * omega[i]).squaredNorm();
                mtx.unlock();
            }
        });
        prime_residual = std::sqrt(prime_residual);

        T dual_residual = 0;
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v -= old_dv.col(g.idx);
        });
        tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
            for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                int i = particle_order[idx];
                TV& Xp = Xarray[i];
                TM& Fn = (*Fn_pointer)[i];
                BSplineWeights<T, dim> spline(Xp, dx);
                TM tmp = TM::Zero();
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    tmp += dt * g.new_v * dw.transpose() * Fn;
                });
                mtx.lock();
                dual_residual += (rho_scale * tmp * omega[i] * omega[i]).squaredNorm();
                mtx.unlock();
            }
        });
        dual_residual = std::sqrt(dual_residual);
        ZIRAN_INFO("prime residual : ", prime_residual, "\t\tdual residual : ", dual_residual, "\t\trho_scale :", rho_scale);
    }

    void admmStep()
    {
        ZIRAN_ASSERT(!use_ruiz, "please use AdmmSimulationFull.h");
        ZIRAN_ASSERT(over_relaxation_alpha == (T)1, "please use AdmmSimulationFull.h");

        Base::force->backupStrain();
        initStep();
        initMultiplierAndPreconditioner();

        out_a = step.time;
        // outputFile(out_a, total_cg_iterations);
        for (int tim = 0; tim < admm_max_iterations; ++tim) {
            updateCollisionNodes();
            // T residual_current = admmResidual();
            T residual_current = 1;
            FStep(local_tolerance);
            dvStep(global_tolerance, std::min((T)0.5, (T)1 * std::sqrt(std::max(residual_current, (T)1e-6))));
            uStep();
        }
        // same with BE method
        Base::constructNewVelocityFromNewtonResult();
        auto& Xarray = particles.X.array;
        tbb::parallel_for(0, (int)particles.count, [&](int i) {
            TV& Xp = Xarray[i];
            BSplineWeights<T, dim> spline(Xp, dx);
            T last_dv_w = 0;
            last_dv[i] = TV::Zero();
            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                if (g.idx >= 0) {
                    last_dv[i] += dv.col(g.idx) * w;
                    last_dv_w += w;
                }
            });
            last_dv[i] /= (last_dv_w * dt);
        });
    }

    virtual void advanceOneTimeStep(double dt)
    {
        ZIRAN_TIMER();
        ZIRAN_INFO("Advance one time step from time ", std::setw(7), step.time, " with                     dt = ", dt);

        Base::reinitialize();

        ZIRAN_ASSERT(!Base::symplectic, "It should be implicit P2G transfer.");

        Base::particlesToGrid();
        admmStep();
        for (size_t k = 0; k < collision_objects.size(); ++k)
            if (collision_objects[k]->updateState)
                collision_objects[k]->updateState(step.time + dt, *collision_objects[k]);
        Base::gridToParticles(dt);
        /*
        tbb::parallel_for(0, (int)particles.count, [&](int i) {
            TM& Fn = (*Fn_pointer)[i];
            TM U, V;
            TV sigma;
            singularValueDecomposition(Fn, U, sigma, V);
            MATH_TOOLS::clamp(sigma, 0.2, 5.0);
            Fn = U * sigma.asDiagonal() * V.transpose();
        });

        tbb::parallel_for(0, (int)particles.count, [&](int i) {
            TM& FFn = (*FFn_pointer)[i];
            TM U, V;
            TV sigma;
            singularValueDecomposition(FFn, U, sigma, V);
            MATH_TOOLS::clamp(sigma, 0.2, 5.0);
            FFn = U * sigma.asDiagonal() * V.transpose();
        });
*/
    }

    const char* name() override { return "admm"; }
};
} // namespace ZIRAN

#endif