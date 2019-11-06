#ifndef FRACTURE_SIMULATION_H
#define FRACTURE_SIMULATION_H

#include <MPM/MpmSimulationBase.h>
#include <Ziran/Math/MathTools.h>
#include <Ziran/Physics/ConstitutiveModel/NeoHookeanBorden.h>
#include <Partio.h>
#include "PhaseFieldSystem.h"
#include "ConjugateGradient.h"

#include "PhaseField.h"

#undef B2

namespace ZIRAN {

template <class T, int dim>
class FractureSimulation : public MpmSimulationBase<T, dim> {
public:
    using Base = MpmSimulationBase<T, dim>;
    using Base::autorestart;
    using Base::begin_time_step_callbacks;
    using Base::block_offset;
    using Base::cfl;
    using Base::dt;
    using Base::dx;
    using Base::flip_pic_ratio;
    using Base::force;
    using Base::frame;
    using Base::grid;
    using Base::mls_mpm;
    using Base::num_nodes;
    using Base::particle_base_offset;
    using Base::particle_group;
    using Base::particle_order;
    using Base::particles;
    using Base::restart_callbacks;
    using Base::scratch_gradV;
    using Base::start_frame;
    using Base::step;
    using Base::symplectic;
    using Base::transfer_scheme;

    typedef Vector<T, dim> TV;
    typedef Vector<int, dim> IV;
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, Eigen::Dynamic> Vec;
    typedef Matrix<T, dim, Eigen::Dynamic> Mat;

    bool use_phase_field = false;
    bool useNaiveDamage = false; //turn on so we can use the naive damage and compare with PFF!
    T sigmaF = 0;
    T parabolic_M = 0;
    T delete_particle_threshold = 0;
    bool lumping = true;

    using Simulation = FractureSimulation<T, dim>;
    using Objective = PhaseFieldSystem<Simulation>;

    Objective cg_objective;
    ConjugateGradient<T, Objective, Vec> cg;

    FractureSimulation()
        : cg_objective(*this), cg(100)
    {
    }

    inline static AttributeName<PhaseField<T, dim>> phase_field_range()
    {
        return AttributeName<PhaseField<T, dim>>("phase field");
    }

    inline static AttributeName<T> element_measure_range()
    {
        return AttributeName<T>("element measure");
    }

    inline static AttributeName<TM> F_range()
    {
        return AttributeName<TM>("F");
    }

    inline static AttributeName<NeoHookeanBorden<T, dim>> neohookean_borden_range()
    {
        return AttributeName<NeoHookeanBorden<T, dim>>(NeoHookeanBorden<T, dim>::name());
    }

    virtual void initialize()
    {
        Base::initialize();

        if (!use_phase_field) return;

        //        T unique_one_over_sigma_c = 0;
        //        tbb::parallel_for(particles.X.ranges,
        //            [&](DisjointRanges& subrange) {
        //                DisjointRanges subset(subrange,
        //                    particles.commonRanges(phase_field_range()));
        //                for (auto iter = particles.subsetIter(subset, phase_field_range()); iter; ++iter) {
        //                    auto& phase_field = iter.template get<0>();
        //                    if (unique_one_over_sigma_c == 0) {
        //                        unique_one_over_sigma_c = phase_field.one_over_sigma_c;
        //                    }
        //                    else {
        //                        ZIRAN_ASSERT(unique_one_over_sigma_c == phase_field.one_over_sigma_c, "all phase field particle should own same one_over_sigma_c");
        //                    }
        //                }
        //            });
        //        ZIRAN_INFO("Particles with phase field, one_over_sigma_c = ", unique_one_over_sigma_c);

        SimulationBase::end_frame_callbacks.emplace_back(
            [&](int frame) {
                std::string filename = SimulationBase::output_dir.absolutePath(SimulationBase::outputFileName("phaseField", ".bgeo"));

                Partio::ParticlesDataMutable* parts = Partio::create();
                Partio::ParticleAttribute posH, pfH, sigmaCH;
                posH = parts->addAttribute("position", Partio::VECTOR, 3);
                pfH = parts->addAttribute("pf", Partio::VECTOR, 1);
                sigmaCH = parts->addAttribute("sigmaC", Partio::VECTOR, 1); //NOTE: this is NOT the inverse of sigmaC!!

                // gather data as a flat array
                StdVector<T> pfData(particles.count, T(1));
                StdVector<T> sigmaCData(particles.count, T(1));
                for (auto iter = particles.iter(phase_field_range()); iter; ++iter) {
                    pfData[iter.entryId()] = iter.template get<0>().c;
                    sigmaCData[iter.entryId()] = iter.template get<0>().one_over_sigma_c;
                }

                // write to partio structure
                for (int k = 0; k < particles.count; k++) {
                    int idx = parts->addParticle();
                    float* posP = parts->dataWrite<float>(posH, idx);
                    float* pfP = parts->dataWrite<float>(pfH, idx);
                    float* sigmaCP = parts->dataWrite<float>(sigmaCH, idx);

                    for (int d = 0; d < 3; ++d)
                        posP[d] = 0;
                    for (int d = 0; d < dim; ++d)
                        posP[d] = (float)particles.X.array[k](d);
                    pfP[0] = (float)pfData[k];
                    sigmaCP[0] = (T)1 / (float)sigmaCData[k]; //store the inverse of the inverse so we have the actual sigmaC!
                }

                Partio::write(filename.c_str(), *parts);
                parts->release();
            });
    }

    virtual void particlesToGrid()
    {
        ZIRAN_TIMER();
        Base::particlesToGrid();

        if (useNaiveDamage) {
            solvePhaseFieldSystemNaive();
        }
        else {
            solvePhaseFieldSystem();
        }
    }

    virtual void reinitialize()
    {
        deleteParticles();

        Base::reinitialize();
    }

    void deleteParticles()
    {
        if (!use_phase_field) return;

        auto& Xarray = particles.X.array;
        auto& Varray = particles.V.array;
        auto* pf_pointer = &particles.DataManager::get(phase_field_range());
        auto* F_pointer = &particles.DataManager::get(F_name<T, dim>());
        auto* vol_pointer = &particles.DataManager::get(element_measure_name<T>());

        for (int i = 0; i < particles.count; ++i) {
            if ((*pf_pointer)[i].c < delete_particle_threshold) {
                if constexpr (dim == 2) {
                    Xarray[i] = TV(4, 4);
                    Varray[i] = TV(0, 0);
                    (*F_pointer)[i] = TM::Identity();
                    (*pf_pointer)[i].allow_damage = false;
                }
                else {
                    Xarray[i] = TV(4, 4, 4);
                    Varray[i] = TV(0, 0, 0);
                    (*F_pointer)[i] = TM::Identity();
                    (*pf_pointer)[i].allow_damage = false;
                }
            }
        }

        tbb::parallel_for(particles.X.ranges,
            [&](DisjointRanges& subrange) {
                DisjointRanges subset(subrange,
                    particles.commonRanges(phase_field_range(),
                        F_range()));
                for (auto iter = particles.subsetIter(subset, phase_field_range(), F_range()); iter; ++iter) {
                    auto& phase_field = iter.template get<0>();
                    auto& F = iter.template get<1>();
                    ZIRAN_ASSERT(!((phase_field.c < delete_particle_threshold) && (F.determinant() != (T)1)), "set F as identity wrongly");
                }
            });
    }

    void phaseP2G()
    {

        ZIRAN_TIMER();

        auto& Xarray = particles.X.array;
        auto& Varray = particles.V.array;
        auto& marray = particles.mass.array;
        auto* pf_pointer = &particles.DataManager::get(phase_field_range());

        if (parabolic_M > 0) {
            for (uint64_t color = 0; color < (1 << dim); ++color) {
                tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
                    if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                        return;
                    for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                        int i = particle_order[idx];
                        TV& Xp = Xarray[i];
                        BSplineWeights<T, dim> spline(Xp, dx);
                        grid.iterateKernel(spline, particle_base_offset[i],
                            [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                                g.phase_field_multiplier += w;
                                g.phase_field += (*pf_pointer)[i].c * w;
                            });
                    }
                });
            }

            grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                g.phase_field /= g.phase_field_multiplier;
            });
        }

        tbb::parallel_for(particles.X.ranges,
            [&](DisjointRanges& subrange) {
                DisjointRanges subset(subrange,
                    particles.commonRanges(phase_field_range(),
                        element_measure_range(),
                        F_range()));
                for (auto iter = particles.subsetIter(subset, phase_field_range(),
                         element_measure_range(), F_range());
                     iter; ++iter) {
                    auto& phase_field = iter.template get<0>();
                    auto& vol0 = iter.template get<1>();
                    auto& F = iter.template get<2>();
                    phase_field.vol = vol0 * F.determinant();
                }
            });
    }

    void solvePhaseFieldSystem()
    {
        ZIRAN_TIMER();

        if (!use_phase_field) return;

        auto& Xarray = particles.X.array;
        auto& Varray = particles.V.array;
        auto& marray = particles.mass.array;
        auto* pf_pointer = &particles.DataManager::get(phase_field_range());

        phaseP2G();

        /*if (parabolic_M > 0) {
            for (uint64_t color = 0; color < (1 << dim); ++color) {
                tbb::parallel_for(0, (int) particle_group.size(), [&](int group_idx) {
                    if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                        return;
                    for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                        int i = particle_order[idx];
                        TV &Xp = Xarray[i];
                        BSplineWeights<T, dim> spline(Xp, dx);
                        grid.iterateKernel(spline, particle_base_offset[i],
                                           [&](const IV &node, T w, const TV &dw, GridState<T, dim> &g) {
                                               g.phase_field_multiplier += w;
                                               g.phase_field += (*pf_pointer)[i].c * w;
                                           });
                    }
                });
            }

            grid.iterateGrid([&](IV node, GridState<T, dim> &g) {
                g.phase_field /= g.phase_field_multiplier;
            });
        }

        tbb::parallel_for(particles.X.ranges,
                          [&](DisjointRanges &subrange) {
                              DisjointRanges subset(subrange,
                                                    particles.commonRanges(phase_field_range(),
                                                                           element_measure_range(),
                                                                           F_range()));
                              for (auto iter = particles.subsetIter(subset, phase_field_range(),
                                                                    element_measure_range(), F_range()); iter; ++iter) {
                                  auto &phase_field = iter.template get<0>();
                                  auto &vol0 = iter.template get<1>();
                                  auto &F = iter.template get<2>();
                                  phase_field.vol = vol0 * F.determinant();
                              }
                          });*/

        Vec x = Vec::Zero(num_nodes, 1);
        Vec rhs = Vec::Zero(num_nodes, 1);

        // build rhs
        for (uint64_t color = 0; color < (1 << dim); ++color) {
            tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
                if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                    return;
                for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                    int i = particle_order[idx];
                    T vol = (*pf_pointer)[i].vol;
                    TV& Xp = particles.X[i];
                    BSplineWeights<T, dim> spline(Xp, dx);
                    grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                        int node_id = g.idx;
                        if (node_id < 0)
                            return;
                        if (parabolic_M > 0) {
                            rhs(node_id) += (parabolic_M + g.phase_field / dt) * vol * w;
                        }
                        else {
                            rhs(node_id) += vol * w;
                        }
                    });
                }
            });
        }

        cg_objective.setMultiplier([&](const Vec& x, Vec& b) {
            Vec c_scp = Vec::Zero(particles.count, 1);
            Mat gradc_scp = Mat::Zero(dim, particles.count);

            tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
                for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                    int i = particle_order[idx];
                    TV& Xp = particles.X[i];
                    T vol = (*pf_pointer)[i].vol;
                    BSplineWeights<T, dim> spline(Xp, dx);
                    grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                        int node_id = g.idx;
                        if (node_id < 0)
                            return;
                        if (!lumping) {
                            c_scp(i) += x(node_id) * w;
                        }
                        gradc_scp.col(i) += x(node_id) * dw;
                    });
                    if (!lumping) {
                        c_scp(i) *= vol * (*pf_pointer)[i].pf_Fp;
                    }
                    T l0 = (*pf_pointer)[i].l0;
                    gradc_scp.col(i) *= vol * (T)4 * l0 * l0;
                }
            });

            b.setZero();

            for (uint64_t color = 0; color < (1 << dim); ++color) {
                tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
                    if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                        return;
                    for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                        int i = particle_order[idx];
                        T vol = (*pf_pointer)[i].vol;
                        TV& Xp = particles.X[i];
                        BSplineWeights<T, dim> spline(Xp, dx);
                        grid.iterateKernel(spline, particle_base_offset[i],
                            [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                                int node_id = g.idx;
                                if (node_id < 0)
                                    return;
                                const PhaseField<T, dim>& phase = (*pf_pointer)[i];
                                if (!lumping) {
                                    b(node_id) += c_scp(i) * w;
                                }
                                else {
                                    if (parabolic_M > 0) {
                                        b(node_id) += vol * (phase.pf_Fp * parabolic_M + (T)1 / dt) * w * x(node_id);
                                    }
                                    else {
                                        b(node_id) += vol * phase.pf_Fp * w * x(node_id);
                                    }
                                }
                                if (parabolic_M > 0) {
                                    b(node_id) += gradc_scp.col(i).dot(dw) * parabolic_M;
                                }
                                else {
                                    b(node_id) += gradc_scp.col(i).dot(dw);
                                }
                            });
                    }
                });
            }
        });

        cg_objective.setPreconditioner([&](const Vec& in, Vec& out) {
            if (!lumping) {
                out = in;
            }
            out.setZero();
            for (uint64_t color = 0; color < (1 << dim); ++color) {
                tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
                    if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                        return;
                    for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                        int i = particle_order[idx];
                        T vol = (*pf_pointer)[i].vol;
                        TV& Xp = particles.X[i];
                        BSplineWeights<T, dim> spline(Xp, dx);
                        grid.iterateKernel(spline, particle_base_offset[i],
                            [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                                int node_id = g.idx;
                                if (node_id < 0)
                                    return;
                                const PhaseField<T, dim>& phase = (*pf_pointer)[i];
                                if (parabolic_M > 0) {
                                    out(node_id) += vol * (parabolic_M * phase.pf_Fp + (T)1 / dt) * w;
                                }
                                else {
                                    out(node_id) += vol * phase.pf_Fp * w;
                                }
                            });
                    }
                });
            }
            for (int i = 0; i < num_nodes; ++i)
                out(i) = in(i) / out(i);
        });

        T tolerance = 1e-6;
        cg.setTolerance(tolerance);
        cg.solve(cg_objective, x, rhs, false);

        if (parabolic_M > 0) {
            grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                g.phase_field = x(g.idx) - g.phase_field;
            });
        }
        else {
            grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                g.phase_field = x(g.idx);
            });
        }

        phaseG2P(); //transfer to particle view to update particle phase
    }

    void phaseG2P()
    {

        ZIRAN_TIMER();

        auto* pf_pointer = &particles.DataManager::get(phase_field_range());

        tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
            for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                int i = particle_order[idx];
                if (!(*pf_pointer)[i].allow_damage)
                    continue;
                TV& Xp = particles.X[i];
                T pf = (T)0;
                BSplineWeights<T, dim> spline(Xp, dx);
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    if (g.idx >= 0)
                        pf += w * g.phase_field;
                });
                // FLIP
                T& c = (*pf_pointer)[i].c;
                T new_c;
                if (parabolic_M > 0) {
                    new_c = c + pf;
                }
                else {
                    new_c = pf;
                }
                c = std::max(std::min(c, new_c), (T)0);
            }
        });

        auto ranges = particles.X.ranges;
        tbb::parallel_for(ranges,
            [&](DisjointRanges& subrange) {
                DisjointRanges subset(subrange,
                    particles.commonRanges(phase_field_range(),
                        neohookean_borden_range(),
                        F_range()));
                for (auto iter = particles.subsetIter(subset, phase_field_range(), neohookean_borden_range(), F_range()); iter; ++iter) {
                    auto& phase_field = iter.template get<0>();
                    auto& model = iter.template get<1>();
                    auto& F = iter.template get<2>();
                    T tmp = model.get_psi_pos(F);
                    phase_field.Update_Phase_Field_Fp(tmp);
                    T c = phase_field.c;
                    T k = phase_field.residual_phase;
                    model.g = c * c * (1 - k) + k;
                }
            });
    }

    void solvePhaseFieldSystemNaive()
    {

        auto ranges = particles.X.ranges;
        tbb::parallel_for(ranges,
            [&](DisjointRanges& subrange) {
                DisjointRanges subset(subrange,
                    particles.commonRanges(phase_field_range(),
                        neohookean_borden_range(),
                        F_range()));
                for (auto iter = particles.subsetIter(subset, phase_field_range(), neohookean_borden_range(), F_range()); iter; ++iter) {
                    auto& phase_field = iter.template get<0>();
                    auto& model = iter.template get<1>();
                    auto& F = iter.template get<2>();
                    //T tmp = model.get_psi_pos(F);
                    //phase_field.Update_Phase_Field_Fp(tmp);
                    T c = phase_field.c;
                    T k = phase_field.residual_phase;

                    //Now, compute c based on naive damage assumption (using sigmaF as our defined threshold)
                    T J = F.determinant();
                    TM tau;
                    typename NeoHookeanBorden<T, dim>::Scratch s;
                    s.F = F;
                    s.J = J;

                    //set model.g to 1 so we compute based on the undegraded principal stress
                    model.g = 1;

                    model.kirchhoff(s, tau); //get tau (kirchoff stress)
                    TM cauchy = ((T)1 / J) * tau;

                    //Now that we have cauchy stress take eigen value decomposition of it
                    Eigen::VectorXcd eigenVals = cauchy.eigenvalues(); //TODO: NEED TO CHANGE THIS WHEN SWITCHING BETWEEN FLOAT AND DOUBLE!!!!
                    T sigmaMax = (T)eigenVals[0].real(); //get the real part of the first eigen value (which shouldnt have be complex anyway lol)
                    T newC = c;
                    if (sigmaMax > sigmaF) {
                        newC = std::min(c, (sigmaF / sigmaMax));
                    }

                    //Update c and g
                    phase_field.c = newC;
                    model.g = newC * newC * (1 - k) + k;
                }
            });
    }
};

} // namespace ZIRAN

#endif
