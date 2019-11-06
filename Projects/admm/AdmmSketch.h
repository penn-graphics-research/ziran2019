#ifndef ZIRAN_ADMMSKETCH_H
#define ZIRAN_ADMMSKETCH_H

#include <Ziran/CS/Util/Forward.h>
#include <Eigen/Eigenvalues>

namespace ZIRAN {
/*

            cg_objective.setMultiplier([&](const TVStack& x, TVStack& b) {
                b.setZero();

                grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                    b.col(g.idx) = g.m * x.col(g.idx).cwiseProduct(Q[g.idx]).cwiseProduct(Q[g.idx]);
                });

                // use new_v to store x
                // G2P
                StdVector<TM> WTWD;
                WTWD.resize((int)particles.count);
                grid.iterateTouchedGrid([&](IV node, GridState<T, dim>& g) {
                    g.new_v = TV::Zero();
                });
                grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                    g.new_v = x.col(g.idx);
                });
                tbb::parallel_for(0, (int)particles.count, [&](int i) {
                    TV& Xp = Xarray[i];
                    TM& Fn = (*Fn_pointer)[i];
                    BSplineWeights<T, dim> spline(Xp, dx);
                    WTWD[i] = TM::Zero();
                    grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                        TM tmp = dt * rho_scale * (g.new_v.cwiseProduct(Q[g.idx])) * dw.transpose() * Fn;
                        WTWD[i] += tmp.cwiseProduct(omega[i]).cwiseProduct(omega[i]).cwiseProduct(R[i]).cwiseProduct(R[i]);
                    });
                });
                // P2G
                grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                    g.new_v = TV::Zero();
                });
                for (uint64_t color = 0; color < (1 << dim); ++color) {
                    tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
                        if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                            return;
                        for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                            int i = particle_order[idx];
                            TV& Xp = Xarray[i];
                            TM& Fn = (*Fn_pointer)[i];
                            BSplineWeights<T, dim> spline(Xp, dx);
                            grid.iterateKernel(spline, particle_base_offset[i], [&](const IV& node, T w, const TV& dw, GridState<T, dim>& g) {
                                TV tmp = dt * WTWD[i] * Fn.transpose() * dw;
                                g.new_v += tmp.cwiseProduct(Q[g.idx]);
                            });
                        }
                    });
                }
                grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                    b.col(g.idx) += g.new_v;
                });
            });

            cg_objective.setPreconditioner([&](const TVStack& in, TVStack& out) {
                grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                    g.new_v = g.m * (Q[g.idx].cwiseProduct(Q[g.idx]));
                });
                for (uint64_t color = 0; color < (1 << dim); ++color) {
                    tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
                        if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                            return;
                        for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                            int i = particle_order[idx];
                            TV& Xp = Xarray[i];
                            TM& Fn = (*Fn_pointer)[i];
                            BSplineWeights<T, dim> spline(Xp, dx);
                            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                                TM tmp = dt * dt * rho_scale * Q[g.idx] * dw.transpose() * Fn;
                                TM tmpooRR = tmp.cwiseProduct(omega[i]).cwiseProduct(R[i]).cwiseProduct(R[i]).cwiseProduct(omega[i]);
                                TV opt = tmpooRR * Fn.transpose() * dw;
                                g.new_v += opt.cwiseProduct(Q[g.idx]);
                            });
                        }
                    });
                }
                grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
                    out.col(g.idx) = (in.col(g.idx)).cwiseQuotient(g.new_v);
                });
            });
 */
} // namespace ZIRAN

#endif
