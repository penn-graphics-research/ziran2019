#ifndef ADMM_FORMULA_H
#define ADMM_FORMULA_H

#include <Ziran/CS/Util/Forward.h>
#include <Eigen/Eigenvalues>

namespace ZIRAN {

#define TV Vector<T, dim>
#define TM Matrix<T, dim, dim>

template <class T, int dim>
TM makePD(const TM& symMtr)
{
    Eigen::SelfAdjointEigenSolver<TM> eigenSolver(symMtr);
    TV D(eigenSolver.eigenvalues());
    for (int i = 0; i < dim; i++) {
        if (D[i] < 0.0) {
            D[i] = 0.0;
        }
    }
    return eigenSolver.eigenvectors() * D.asDiagonal() * eigenSolver.eigenvectors().transpose();
};

template <class T>
Matrix<T, 2, 2> makePD2D(const Matrix<T, 2, 2>& symMtr)
{
    Eigen::SelfAdjointEigenSolver<Matrix<T, 2, 2>> eigenSolver(symMtr);
    Vector<T, 2> D(eigenSolver.eigenvalues());
    for (int i = 0; i < 2; i++) {
        if (D[i] < 0.0)
            D[i] = 0.0;
    }
    return eigenSolver.eigenvectors() * D.asDiagonal() * eigenSolver.eigenvectors().transpose();
};

void outputResidual(double a, double b)
{
    static bool first_line_output = true;
    if (first_line_output) {
        FILE* f = fopen("output/residual.txt", "w");
        fclose(f);
        first_line_output = false;
    }
    FILE* f = fopen("output/residual.txt", "a");
    fprintf(f, "%.10lf %.10lf\n", a, b);
    fclose(f);
}

#undef TV
#undef TM

} // namespace ZIRAN

#endif

/*
void outputGlobal()
{
    FILE* f = fopen("output/result_original.txt", "w");
    TVStack x = dv, y = dv;
    for (int i = 0; i < Base::num_nodes; ++i)
        for (int j = 0; j < dim; ++j) {
            printf("%d %d %d %d\n", i, j, Base::num_nodes, dim);
            x.setZero();
            x.col(i)(j) = 1;
            cg_objective.multiply(x, y);
            for (int p = 0; p < Base::num_nodes; ++p)
                for (int q = 0; q < dim; ++q)
                    fprintf(f, "%.30f ", y.col(p)(q));
            fprintf(f, "\n");
        }
    fclose(f);
}

void outputGlobalOOT()
{
    auto multiply = [&](const TStack& x, TStack& b) {
        b.setZero();
        // P2G
        auto grid_array = grid.grid->Get_Array();
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            g.new_v = TV::Zero();
        });
        for (int i = 0; i < (int)div_nodes.size(); ++i) {
            uint64_t offset00 = div_nodes[i];
            uint64_t offset10 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 0));
            uint64_t offset01 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(0, 1));
            uint64_t offset11 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 1));
            auto& g00 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset00));
            auto& g10 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset10));
            auto& g01 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset01));
            auto& g11 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset11));
            T tmp = x(i) / (T)2 / dx;
            g00.new_v += TV(-tmp, -tmp).cwiseProduct(Q[g00.idx]);
            g10.new_v += TV(tmp, -tmp).cwiseProduct(Q[g10.idx]);
            g01.new_v += TV(-tmp, tmp).cwiseProduct(Q[g01.idx]);
            g11.new_v += TV(tmp, tmp).cwiseProduct(Q[g11.idx]);
        }

        tbb::parallel_for(0, (int)div_nodes.size(), [&](int i) {
            uint64_t offset00 = div_nodes[i];
            uint64_t offset10 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 0));
            uint64_t offset01 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(0, 1));
            uint64_t offset11 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 1));
            auto& g00 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset00));
            auto& g10 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset10));
            auto& g01 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset01));
            auto& g11 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset11));
            TV v00 = g00.new_v.cwiseProduct(Q[g00.idx]);
            TV v10 = g10.new_v.cwiseProduct(Q[g10.idx]);
            TV v01 = g01.new_v.cwiseProduct(Q[g01.idx]);
            TV v11 = g11.new_v.cwiseProduct(Q[g11.idx]);
            b(i) = rho_scale * omegain[i] * Rin[i] * Rin[i] * omegain[i] * get_div(v00, v10, v01, v11);
        });
    };

    FILE* f = fopen("output/result_oot.txt", "w");
    TStack x = TStack::Zero(div_nodes.size());
    TStack y = TStack::Zero(div_nodes.size());
    for (int i = 0; i < (int)div_nodes.size(); ++i) {
        x.setZero();
        x(i) = 1;
        multiply(x, y);
        for (int p = 0; p < (int)div_nodes.size(); ++p)
            fprintf(f, "%.30f ", y(p));
        fprintf(f, "\n");
    }
    fclose(f);
}

void outputGlobalMass()
{
    auto& Xarray = particles.X.array;
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = g.m * (Q[g.idx].cwiseProduct(Q[g.idx]));
    });
    TVStack pre = dv;
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        pre.col(g.idx) = g.new_v;
    });
    FILE* f = fopen("output/result_mass.txt", "w");
    TVStack x = dv, y = dv;
    for (int i = 0; i < Base::num_nodes; ++i)
        for (int j = 0; j < dim; ++j)
            fprintf(f, "%.30f ", pre.col(i)(j));
    fclose(f);
}

void outputGlobalPreconditioner()
{
    auto& Xarray = particles.X.array;
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = TV::Zero();
    });
//        for (uint64_t color = 0; color < (1 << dim); ++color) {
//            tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
//                if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
//                    return;
//                for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
//                    int i = particle_order[idx];
//                    TV& Xp = Xarray[i];
//                    TM& Fn = (*Fn_pointer)[i];
//                    BSplineWeights<T, dim> spline(Xp, dx);
//                    grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
//                        TM tmp = dt * dt * rho_scale * Q[g.idx] * dw.transpose() * Fn;
//                        TM tmpooRR = tmp.cwiseProduct(omega[i]).cwiseProduct(R[i]).cwiseProduct(R[i]).cwiseProduct(omega[i]);
//                        TV opt = tmpooRR * Fn.transpose() * dw;
//                        g.new_v += opt.cwiseProduct(Q[g.idx]);
//                    });
//                }
//            });
//        }
    if constexpr (use_incompressibility && dim == 2) {
        auto grid_array = grid.grid->Get_Array();
        tbb::parallel_for(0, (int)div_nodes.size(), [&](int i) {
            uint64_t offset00 = div_nodes[i];
            uint64_t offset10 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 0));
            uint64_t offset01 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(0, 1));
            uint64_t offset11 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 1));
            auto& g00 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset00));
            auto& g10 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset10));
            auto& g01 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset01));
            auto& g11 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset11));
            mtx.lock();
            g00.new_v += (rho_scale * Q[g00.idx] * omegain[i] * Rin[i] * Rin[i] * omegain[i]).cwiseProduct(Q[g00.idx]) / (T)4 / dx / dx;
            g10.new_v += (rho_scale * Q[g10.idx] * omegain[i] * Rin[i] * Rin[i] * omegain[i]).cwiseProduct(Q[g10.idx]) / (T)4 / dx / dx;
            g01.new_v += (rho_scale * Q[g01.idx] * omegain[i] * Rin[i] * Rin[i] * omegain[i]).cwiseProduct(Q[g01.idx]) / (T)4 / dx / dx;
            g11.new_v += (rho_scale * Q[g11.idx] * omegain[i] * Rin[i] * Rin[i] * omegain[i]).cwiseProduct(Q[g11.idx]) / (T)4 / dx / dx;
            mtx.unlock();
        });
    }
    TVStack pre = dv;
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        pre.col(g.idx) = g.new_v;
    });
    FILE* f = fopen("output/result_diagonal.txt", "w");
    TVStack x = dv, y = dv;
    for (int i = 0; i < Base::num_nodes; ++i)
        for (int j = 0; j < dim; ++j)
            fprintf(f, "%.30f ", pre.col(i)(j));
    fclose(f);
}

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
            prime_residual += (tmp.cwiseProduct(omega[i]).cwiseProduct(R[i])).squaredNorm();
            mtx.unlock();
        }
    });
    if constexpr (use_incompressibility && dim == 2) {
        auto grid_array = grid.grid->Get_Array();
        tbb::parallel_for(0, (int)div_nodes.size(), [&](int i) {
            uint64_t offset00 = div_nodes[i];
            uint64_t offset10 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 0));
            uint64_t offset01 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(0, 1));
            uint64_t offset11 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 1));
            auto& g00 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset00));
            auto& g10 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset10));
            auto& g01 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset01));
            auto& g11 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset11));
            TV v00 = g00.v + dv.col(g00.idx);
            TV v10 = g10.v + dv.col(g10.idx);
            TV v01 = g01.v + dv.col(g01.idx);
            TV v11 = g11.v + dv.col(g11.idx);
            mtx.lock();
            T tmp = Rin[i] * omegain[i] * get_div(v00, v10, v01, v11);
            prime_residual += tmp * tmp;
            mtx.unlock();
        });
    }
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
            dual_residual += (rho_scale * tmp.cwiseProduct(omega[i]).cwiseProduct(R[i]).cwiseProduct(R[i]).cwiseProduct(omega[i])).squaredNorm();
            mtx.unlock();
        }
    });
    dual_residual = std::sqrt(dual_residual);

    if (prime_residual > 10 * dual_residual) {
        rho_scale *= 2;
    }
    else if (10 * prime_residual < dual_residual) {
        rho_scale *= 0.5;
    }
}

void updateRhoScaleRelative()
{
    auto& Xarray = particles.X.array;
    // use new_v to store dv
    // G2P
    T prime_residual = 0;
    T prime_invidual = 0;
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
            TM tmp1 = Fn, tmp2 = F[i];
            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                tmp2 += dt * g.v * dw.transpose() * Fn;
                tmp2 += dt * g.new_v * dw.transpose() * Fn;
            });
            TM tmp = tmp1 - tmp2;
            mtx.lock();
            prime_residual = std::max(prime_residual, (tmp.cwiseProduct(omega[i]).cwiseProduct(R[i])).cwiseAbs().maxCoeff());
            prime_invidual = std::max(prime_invidual, (tmp1.cwiseProduct(omega[i]).cwiseProduct(R[i])).cwiseAbs().maxCoeff());
            prime_invidual = std::max(prime_invidual, (tmp2.cwiseProduct(omega[i]).cwiseProduct(R[i])).cwiseAbs().maxCoeff());
            //prime_invidual = std::max(prime_invidual, (tmp3.cwiseProduct(omega[i]).cwiseProduct(R[i])).cwiseAbs().maxCoeff());
            mtx.unlock();
        }
    });
    prime_residual = std::sqrt(prime_residual / prime_invidual);
    T dual_residual = 0;
    T dual_invidual = 0;

    auto ce_name = constitutive_model_name<CorotatedElasticity<T, dim>>;
    auto ranges = particles.X.ranges;
    tbb::parallel_for(ranges, [&](DisjointRanges& subrange) {
        DisjointRanges subset(subrange, particles.commonRanges(ce_name(), element_measure_name<T>(), F_name()));
        for (auto iter = particles.subsetIter(subset, ce_name(), element_measure_name<T>(), F_name()); iter; ++iter) {
            auto& constitutive_model = iter.template get<0>();
            auto& vol = iter.template get<1>();
            auto& Fn = iter.template get<2>();
            int i = iter.entryId();
            TV& Xp = Xarray[i];
            BSplineWeights<T, dim> spline(Xp, dx);
            TM tmp1 = TM::Zero(), tmp2 = TM::Zero();
            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                if (g.idx >= 0) {
                    tmp1 += dt * g.new_v * dw.transpose() * Fn;
                    tmp2 += dt * old_dv.col(g.idx) * dw.transpose() * Fn;
                }
            });
            TM tmp = tmp1 - tmp2;
            mtx.lock();
            TM u_bar = -y[i].cwiseQuotient(omega[i]) / rho_scale;
            dual_residual = std::max(dual_residual, (rho_scale * tmp.cwiseProduct(omega[i]).cwiseProduct(R[i]).cwiseProduct(R[i]).cwiseProduct(omega[i])).cwiseAbs().maxCoeff());
            dual_invidual = std::max(dual_invidual, (rho_scale * (u_bar.cwiseProduct(omega[i]).cwiseProduct(omega[i]))).cwiseAbs().maxCoeff());

            typename CorotatedElasticity<T, dim>::Scratch s;
            constitutive_model.updateScratch(F[i], s);
            TM firstPiola;
            constitutive_model.firstPiola(s, firstPiola);
            TM VP = vol * firstPiola;

            dual_invidual = std::max(dual_invidual, VP.cwiseAbs().maxCoeff());
            mtx.unlock();
        }
    });

    dual_residual = std::sqrt(dual_residual / dual_invidual);
    if (prime_residual > 10 * dual_residual) {
        rho_scale *= 2;
    }
    else if (10 * prime_residual < dual_residual) {
        rho_scale *= 0.5;
    }
}

void calcDivergence()
{
    T result = 0;
    auto grid_array = grid.grid->Get_Array();
    tbb::parallel_for(0, (int)div_nodes.size(), [&](int i) {
        uint64_t offset00 = div_nodes[i];
        uint64_t offset10 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 0));
        uint64_t offset01 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(0, 1));
        uint64_t offset11 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 1));
        auto& g00 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset00));
        auto& g10 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset10));
        auto& g01 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset01));
        auto& g11 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset11));
        TV v00 = g00.idx < 0 ? TV::Zero() : (g00.v + dv.col(g00.idx)).eval();
        TV v10 = g10.idx < 0 ? TV::Zero() : (g10.v + dv.col(g10.idx)).eval();
        TV v01 = g01.idx < 0 ? TV::Zero() : (g01.v + dv.col(g01.idx)).eval();
        TV v11 = g11.idx < 0 ? TV::Zero() : (g11.v + dv.col(g11.idx)).eval();
        mtx.lock();
        T tmp = get_div(v00, v10, v01, v11);
        result += tmp * tmp;
        mtx.unlock();
    });
    printf("=========== %.20lf\n", result);
}

void FStep(T localTol)
{
    ZIRAN_TIMER();
    auto& Xarray = particles.X.array;
    grid.iterateTouchedGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = TV::Zero();
    });
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        g.new_v = dv.col(g.idx);
    });

    auto ce_name = constitutive_model_name<CorotatedElasticity<T, dim>>;
    auto ranges = particles.X.ranges;
    tbb::parallel_for(ranges, [&](DisjointRanges& subrange) {
        DisjointRanges subset(subrange, particles.commonRanges(ce_name(), element_measure_name<T>(), F_name()));
        for (auto iter = particles.subsetIter(subset, ce_name(), element_measure_name<T>(), F_name()); iter; ++iter) {
            auto& constitutive_model = iter.template get<0>();
            auto& vol = iter.template get<1>();
            auto& Fn = iter.template get<2>();
            int i = iter.entryId();
            TV& Xp = Xarray[i];
            // build rhs
            TM u_bar = -y[i].cwiseQuotient(omega[i]) / rho_scale;
            TM FDb = Fn + u_bar.cwiseQuotient(R[i]).cwiseQuotient(R[i]);
            BSplineWeights<T, dim> spline(Xp, dx);
            grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                FDb += dt * g.v * dw.transpose() * Fn;
                FDb += dt * g.new_v * dw.transpose() * Fn;
            });
            typedef typename CorotatedElasticity<T, dim>::Hessian Hessian;
            // SVD rhs_raw
            auto computeEnergyVal_zUpdate = [&](const FlattenTM& flatten_F, T& E) {
                TM F = TM(flatten_F.data());
                typename CorotatedElasticity<T, dim>::Scratch s;
                constitutive_model.updateScratch(F, s);
                E = vol * constitutive_model.psi(s);
                for (int col = 0; col < dim; ++col)
                    for (int row = 0; row < dim; ++row) {
                        T ooRR = omega[i](row, col) * omega[i](row, col) * R[i](row, col) * R[i](row, col);
                        E += 0.5 * rho_scale * ooRR * (F(row, col) - FDb(row, col)) * (F(row, col) - FDb(row, col));
                    }
                // E = vol * constitutive_model.psi(s) + 0.5 * rho_scale * ((F - FDb).cwiseProduct(omega[i]).cwiseProduct(R[i])).squaredNorm();
            };
            auto computeGradient_zUpdate = [&](const FlattenTM& flatten_F, FlattenTM& g) {
                TM F = TM(flatten_F.data());
                typename CorotatedElasticity<T, dim>::Scratch s;
                constitutive_model.updateScratch(F, s);
                TM firstPiola;
                constitutive_model.firstPiola(s, firstPiola);
                g = vol * FlattenTM(firstPiola.data());
                for (int col = 0; col < dim; ++col)
                    for (int row = 0; row < dim; ++row) {
                        int idx = col * dim + row;
                        T ooRR = omega[i](row, col) * omega[i](row, col) * R[i](row, col) * R[i](row, col);
                        g(idx) += rho_scale * ooRR * (F(row, col) - FDb(row, col));
                    }
                // g = FlattenTM((vol * firstPiola + rho_scale * (F - FDb).cwiseProduct(omega[i]).cwiseProduct(omega[i]).cwiseProduct(R[i]).cwiseProduct(R[i])).eval().data());
            };
            auto computeHessianProxy_zUpdate = [&](const FlattenTM& flatten_F, Hessian& P) {
                TM F = TM(flatten_F.data());
                typename CorotatedElasticity<T, dim>::Scratch s;
                constitutive_model.updateScratch(F, s);
                Hessian firstPiolaDerivative;
                constitutive_model.firstPiolaDerivative(s, firstPiolaDerivative);
                // TODO: diagonal for this
                makePD(firstPiolaDerivative);
                P = vol * firstPiolaDerivative;
                for (int col = 0; col < dim; ++col)
                    for (int row = 0; row < dim; ++row) {
                        int idx = col * dim + row;
                        T ooRR = omega[i](row, col) * omega[i](row, col) * R[i](row, col) * R[i](row, col);
                        P(idx, idx) += rho_scale * ooRR;
                    }
                // P = vol * firstPiolaDerivative + rho_scale * Hessian(FlattenTM((omega[i].cwiseProduct(omega[i]).cwiseProduct(R[i]).cwiseProduct(R[i])).eval().data()).asDiagonal());
            };
            FlattenTM Fi(F[i].data());
            FlattenTM g;
            for (int j = 0;; j++) {
                computeGradient_zUpdate(Fi, g);
                if (g.norm() < localTol) {
                    break;
                }
                Hessian P;
                computeHessianProxy_zUpdate(Fi, P);
                FlattenTM p = P.ldlt().solve(-g);
                T alpha = 1.0;
                FlattenTM Fi0 = Fi;
                T E0;
                computeEnergyVal_zUpdate(Fi0, E0);
                Fi = Fi0 + alpha * p;
                T E;
                computeEnergyVal_zUpdate(Fi, E);
                while (E > E0) {
                    alpha /= 2.0;
                    Fi = Fi0 + alpha * p;
                    computeEnergyVal_zUpdate(Fi, E);
                }
                if (E0 - E < 1e-6 * E0) {
                    break;
                }
            }
            F[i] = TM(Fi.data());
            //                F[i] = TM::Identity();
        }
    });
}

void ruizEqulibration()
{
    ZIRAN_TIMER();
    int Np = particles.count;
    int Nn = Base::num_nodes;
    rGeps = TStack::Zero(Np * dim * dim);
    rQeps = TStack::Zero(Nn * dim);
    rReps = TStack::Zero(Np * dim * dim);
    rRineps = TStack::Zero((int)div_nodes.size());
    T eps = 1e-12;
    for (int tim = 0;; ++tim) {
        T Gresidual = (TStack::Ones(Np * dim * dim) - rGeps).cwiseAbs().maxCoeff();
        T Qresidual = (TStack::Ones(Nn * dim) - rQeps).cwiseAbs().maxCoeff();
        T Rresidual = (TStack::Ones(Np * dim * dim) - rReps).cwiseAbs().maxCoeff();
        T Rinresidual = div_nodes.size() == 0 ? (T)0 : (TStack::Ones((int)div_nodes.size()) - rRineps).cwiseAbs().maxCoeff();
        if (std::max(std::max(Gresidual, Qresidual), std::max(Rresidual, Rinresidual)) < eps) {
            ZIRAN_WARN("Residuals : ", Gresidual, " ", Qresidual, " ", Rresidual, " ", Rinresidual, ". Convergence number : ", tim);
            break;
        }
        // rGeps
        DEFINE_CE_NAME
        //auto ce_name = constitutive_model_name<CorotatedElasticity<T, dim>>;
        auto ranges = particles.X.ranges;
        tbb::parallel_for(ranges, [&](DisjointRanges& subrange) {
            DisjointRanges subset(subrange, particles.commonRanges(ce_name(), element_measure_name<T>(), F_name()));
            for (auto iter = particles.subsetIter(subset, ce_name(), element_measure_name<T>(), F_name()); iter; ++iter) {
                auto& constitutive_model = iter.template get<0>();
                auto& vol = iter.template get<1>();
                auto& Fn = iter.template get<2>();
                int i = iter.entryId();
                typename CorotatedElasticity<T, dim>::Scratch s;
                constitutive_model.updateScratch(Fn, s);
                typename CorotatedElasticity<T, dim>::Hessian dPdF;
                constitutive_model.firstPiolaDerivative(s, dPdF);
                for (int col = 0; col < dim; ++col)
                    for (int row = 0; row < dim; ++row) {
                        int idx = col * dim + row;
                        FlattenTM cp = (vol * dPdF.row(idx).transpose()).cwiseProduct(FlattenTM(G[i].data()));
                        T Zcom = std::abs(cp.cwiseAbs().maxCoeff() * G[i](row, col));
                        T WTcom = std::abs(omega[i](row, col) * R[i](row, col) * G[i](row, col));
                        rGeps(i * dim * dim + idx) = std::max(Zcom, WTcom);
                    }
            }
        });
        // rQeps
        grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
            for (int d = 0; d < dim; ++d)
                rQeps(g.idx * dim + d) = std::abs(g.m * Q[g.idx](d) * Q[g.idx](d));
        });
        auto& Xarray = particles.X.array;
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
                        for (int d = 0; d < dim; ++d) {
                            TV ts = ((R[i].row(d)).cwiseProduct(omega[i].row(d))).transpose();
                            T DTWTcom = std::abs(ts.cwiseProduct(dt * Fn.transpose() * dw).cwiseAbs().maxCoeff() * Q[g.idx](d));
                            rQeps(g.idx * dim + d) = std::max(rQeps(g.idx * dim + d), DTWTcom);
                        }
                    });
                }
            });
        }
        if constexpr (use_incompressibility && dim == 2) {
            auto grid_array = grid.grid->Get_Array();
            tbb::parallel_for(0, (int)div_nodes.size(), [&](int i) {
                uint64_t offset00 = div_nodes[i];
                uint64_t offset10 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 0));
                uint64_t offset01 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(0, 1));
                uint64_t offset11 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 1));
                auto& g00 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset00));
                auto& g10 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset10));
                auto& g01 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset01));
                auto& g11 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset11));
                rQeps(g00.idx * dim + 0) = std::max(rQeps(g00.idx * dim + 0), std::abs(Q[g00.idx](0) * omegain[i] * Rin[i] / (T)2 / dx));
                rQeps(g00.idx * dim + 1) = std::max(rQeps(g00.idx * dim + 1), std::abs(Q[g00.idx](1) * omegain[i] * Rin[i] / (T)2 / dx));
                rQeps(g10.idx * dim + 0) = std::max(rQeps(g10.idx * dim + 0), std::abs(Q[g10.idx](0) * omegain[i] * Rin[i] / (T)2 / dx));
                rQeps(g10.idx * dim + 1) = std::max(rQeps(g10.idx * dim + 1), std::abs(Q[g10.idx](1) * omegain[i] * Rin[i] / (T)2 / dx));
                rQeps(g01.idx * dim + 0) = std::max(rQeps(g01.idx * dim + 0), std::abs(Q[g01.idx](0) * omegain[i] * Rin[i] / (T)2 / dx));
                rQeps(g01.idx * dim + 1) = std::max(rQeps(g01.idx * dim + 1), std::abs(Q[g01.idx](1) * omegain[i] * Rin[i] / (T)2 / dx));
                rQeps(g11.idx * dim + 0) = std::max(rQeps(g11.idx * dim + 0), std::abs(Q[g11.idx](0) * omegain[i] * Rin[i] / (T)2 / dx));
                rQeps(g11.idx * dim + 1) = std::max(rQeps(g11.idx * dim + 1), std::abs(Q[g11.idx](1) * omegain[i] * Rin[i] / (T)2 / dx));
            });
        }
        // rReps
        tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
            for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                int i = particle_order[idx];
                TV& Xp = Xarray[i];
                TM& Fn = (*Fn_pointer)[i];
                for (int col = 0; col < dim; ++col)
                    for (int row = 0; row < dim; ++row) {
                        int idx = col * dim + row;
                        rReps(i * dim * dim + idx) = std::abs(omega[i](row, col) * R[i](row, col) * G[i](row, col));
                    }
                BSplineWeights<T, dim> spline(Xp, dx);
                grid.iterateKernel(spline, particle_base_offset[i], [&](IV node, T w, TV dw, GridState<T, dim>& g) {
                    for (int col = 0; col < dim; ++col)
                        for (int row = 0; row < dim; ++row) {
                            int idx = col * dim + row;
                            T WDcom = std::abs(Q[g.idx](row) * (dt * dw.transpose() * Fn).transpose()(col) * omega[i](row, col) * R[i](row, col));
                            rReps(i * dim * dim + idx) = std::max(rReps(i * dim * dim + idx), WDcom);
                        }
                });
            }
        });
        // rRineps
        if constexpr (use_incompressibility && dim == 2) {
            auto grid_array = grid.grid->Get_Array();
            tbb::parallel_for(0, (int)div_nodes.size(), [&](int i) {
                uint64_t offset00 = div_nodes[i];
                uint64_t offset10 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 0));
                uint64_t offset01 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(0, 1));
                uint64_t offset11 = SparseMask::Packed_Add(offset00, SparseMask::Linear_Offset(1, 1));
                auto& g00 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset00));
                auto& g10 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset10));
                auto& g01 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset01));
                auto& g11 = reinterpret_cast<GridState<T, dim>&>(grid_array(offset11));
                T Q00 = Q[g00.idx].cwiseAbs().maxCoeff();
                T Q10 = Q[g10.idx].cwiseAbs().maxCoeff();
                T Q01 = Q[g01.idx].cwiseAbs().maxCoeff();
                T Q11 = Q[g11.idx].cwiseAbs().maxCoeff();
                rRineps(i) = std::abs(Rin[i] * omegain[i] * std::max(std::max(Q00, Q10), std::max(Q01, Q11)) / (T)2 / dx);
            });
        }
        // update
        StdVector<TM> GG(G);
        StdVector<TV> QQ(Q);
        StdVector<TM> RR(R);
        StdVector<T> RRin(Rin);
        for (int i = 0; i < Np; ++i)
            G[i] = GG[i].cwiseQuotient(TM(rGeps.template block<dim * dim, 1>(i * dim * dim, 0).data()).cwiseSqrt());
        for (int i = 0; i < Nn; ++i)
            Q[i] = QQ[i].cwiseQuotient(TV(rQeps.template block<dim, 1>(i * dim, 0).data()).cwiseSqrt());
        for (int i = 0; i < Np; ++i)
            R[i] = RR[i].cwiseQuotient(TM(rReps.template block<dim * dim, 1>(i * dim * dim, 0).data()).cwiseSqrt());
        for (int i = 0; i < (int)div_nodes.size(); ++i)
            Rin[i] = RRin[i] / std::sqrt(rRineps(i));
    }
}

void writeState(std::ostream& out)
{
    Base::writeState(out);

    std::string filename = SimulationBase::output_dir.absolutePath(SimulationBase::outputFileName("residual", ".bgeo"));

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, typeH, residualH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    typeH = parts->addAttribute("type", Partio::VECTOR, 1);
    residualH = parts->addAttribute("residual", Partio::VECTOR, 1);

    // write to partio structure
    for (int k = 0; k < particles.count; k++) {
        int idx = parts->addParticle();
        float* posP = parts->dataWrite<float>(posH, idx);
        float* typeP = parts->dataWrite<float>(typeH, idx);
        float* residualP = parts->dataWrite<float>(residualH, idx);
        for (int d = 0; d < 3; ++d)
            posP[d] = 0;
        for (int d = 0; d < dim; ++d)
            posP[d] = (float)particles.X.array[k](d);
        typeP[0] = 0;
        residualP[0] = 0;
    }

    StdVector<TV> grid_pos;
    grid_pos.resize(Base::num_nodes);
    grid.iterateGrid([&](IV node, GridState<T, dim>& g) {
        grid_pos[g.idx] = node.template cast<T>() * dx;
    });
    for (int k = 0; k < Base::num_nodes; k++) {
        int idx = parts->addParticle();
        float* posP = parts->dataWrite<float>(posH, idx);
        float* typeP = parts->dataWrite<float>(typeH, idx);
        float* residualP = parts->dataWrite<float>(residualH, idx);
        for (int d = 0; d < 3; ++d)
            posP[d] = 0;
        for (int d = 0; d < dim; ++d)
            posP[d] = (float)grid_pos[k](d);
        typeP[0] = 1;
        residualP[0] = _residual.col(k).norm();
    }

    Partio::write(filename.c_str(), *parts);
    parts->release();
}

*/