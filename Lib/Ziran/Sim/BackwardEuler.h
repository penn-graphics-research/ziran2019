#ifndef BACKWARD_EULER_H
#define BACKWARD_EULER_H

#include <Ziran/Math/Linear/DenseExt.h>
#include <Ziran/Math/Linear/KrylovSolvers.h>
#include <Ziran/Math/Linear/Minres.h>
#include <Ziran/Math/Linear/GeneralizedMinimalResidual.h>
#include <Ziran/Math/Linear/ConjugateGradient.h>
#include <Ziran/Math/MathTools.h>
#include <Ziran/Sim/SimulationBase.h>
#include <Ziran/Math/Linear/EigenSparseLU.h>

#include <tbb/tbb.h>

namespace ZIRAN {

/**
   Simulation owns Objective and NewtonsMethod.
   This will allow options to do Newton.
   It stores the parameters for them. But it is a little weird.
   This is supposed to be used by both MPM and Elasticity.
 */
template <class Simulation>
class BackwardEulerLagrangianForceObjective {
public:
    using T = typename Simulation::Scalar;
    static const int dim = Simulation::dim;
    using Objective = BackwardEulerLagrangianForceObjective<Simulation>;
    using TV = Vector<T, dim>;
    using TVStack = Matrix<T, dim, Eigen::Dynamic>;
    using Vec = Vector<T, Eigen::Dynamic>;
    using SpMat = Eigen::SparseMatrix<T, Eigen::RowMajor>;

    using Scalar = T;
    using NewtonVector = TVStack;

    Simulation& simulation;

    SpMat implicit_system_matrix;

    std::function<void(TVStack&)> project;
    std::function<void(const TVStack&, TVStack&)> precondition;

    GeneralizedMinimalResidual<T, Objective, TVStack> gmres;
    std::unique_ptr<LinearSolver<T, Objective, TVStack>> linear_solver;
    Minres<T, Objective, TVStack>& minres;

    bool matrix_free;
    bool lu;

    BackwardEulerLagrangianForceObjective(Simulation& simulation)
        : simulation(simulation)
        , project([](TVStack&) {})
        , precondition([](const TVStack& x, TVStack& b) { b = x; })
        , gmres(20)
        , linear_solver(std::make_unique<Minres<T, Objective, TVStack>>(20))
        , minres(*reinterpret_cast<Minres<T, Objective, TVStack>*>(linear_solver.get()))
        , matrix_free(false)
        , lu(false)
    {
        gmres.setTolerance(1);
        linear_solver->setTolerance(1);
    }

    BackwardEulerLagrangianForceObjective(const BackwardEulerLagrangianForceObjective& objective) = delete;

    ~BackwardEulerLagrangianForceObjective() {}

    template <class Projection>
    void initialize(const Projection& projection_function)
    {
        ZIRAN_INFO("Krylov matrix free : ", matrix_free);
        project = projection_function;
    }

    void reinitialize()
    {
        if (!matrix_free)
            initializeImplicitTimeSteppingMatrix(); //set up BE/QS matrix
    }

    void switchLinearSolver(SimulationBase::LinearSolverType p_linear_solver_type)
    {
        switch (p_linear_solver_type) {
        case SimulationBase::MINRES:
            linear_solver = std::make_unique<Minres<T, Objective, TVStack>>(20);
            break;

        case SimulationBase::CG:
            linear_solver = std::make_unique<ConjugateGradient<T, Objective, TVStack>>(20);
            break;

        default:
            ZIRAN_ASSERT(0 && "Invalid linear solver type!");
            break;
        }
        simulation.linear_solver_type = p_linear_solver_type;
    }

    // called by Newton
    void computeResidual(TVStack& residual)
    {
        ZIRAN_ASSERT(residual.cols() == simulation.num_nodes);
        tbb::parallel_for(tbb::blocked_range<int>(0, residual.cols(), 256),
            [&](const tbb::blocked_range<int>& range) {
                int start = range.begin();
                int length = (range.end() - range.begin());
                TV dtg = simulation.dt * simulation.gravity;
                residual.middleCols(start, length) = dtg * simulation.mass_matrix.segment(start, length).transpose();
            });

        simulation.addScaledForces(simulation.dt, residual);

        if (!simulation.quasistatic) {
            simulation.inertia->addScaledForces(simulation.dt, residual);
        }
        project(residual);

        if (simulation.full_implicit) {
            for (auto& lf : simulation.forces) {
                lf->updateStrainWithFullImplicit();
            }
        }
    }

    // called by Newton
    T computeNorm(const TVStack& residual) const
    {
        ZIRAN_ASSERT(residual.cols() == simulation.num_nodes);
        T norm_sq = tbb::parallel_reduce(tbb::blocked_range<int>(0, residual.cols(), 256), (T)0,
            [&](const tbb::blocked_range<int>& range, T ns) -> T {
                int start = range.begin();
                int length = (range.end() - range.begin());
                const auto& r_block = residual.middleCols(start, length);
                const auto& mass_block = simulation.mass_matrix.segment(start, length).transpose();
                return ns + ((r_block.colwise().squaredNorm()).array() / mass_block.array()).sum();
            },
            [](T a, T b) -> T { return a + b; });
        return std::sqrt(norm_sq);
    }

    T innerProduct(const TVStack& ddv, const TVStack& residual) const
    {
        ZIRAN_ASSERT(residual.cols() == ddv.cols());
        T result = tbb::parallel_reduce(tbb::blocked_range<int>(0, residual.cols(), 256), (T)0,
            [&](const tbb::blocked_range<int>& range, T ns) -> T {
                int start = range.begin();
                int length = (range.end() - range.begin());
                const auto& ddv_block = ddv.middleCols(start, length);
                const auto& r_block = residual.middleCols(start, length);
                const auto& mass_block = simulation.mass_matrix.segment(start, length);
                // ZIRAN_DEBUG(ddv_block);
                // ZIRAN_DEBUG(r_block);
                // ZIRAN_DEBUG(mass_block);

                T to_add = ((ddv_block.transpose() * r_block).diagonal().array() / mass_block.array()).sum();
                // ZIRAN_DEBUG(to_add);
                return ns += to_add;
            },
            [](T a, T b) -> T { return a + b; });
        return result;
    }

    // called by Newton
    void updateState(const TVStack& dv)
    {
        ZIRAN_ASSERT(dv.cols() == simulation.num_nodes);
        simulation.moveNodes(dv);
        for (auto& lf : simulation.forces) {
            lf->updatePositionBasedState();
        }
        if (!simulation.quasistatic) {
            simulation.inertia->updatePositionBasedState();
        }
    }

    T totalEnergy()
    {
        T result = 0;
        for (auto& lf : simulation.forces) {
            result += lf->totalEnergy();
        }
        if (!simulation.quasistatic) {
            result += simulation.inertia->totalEnergy();
        }

        result += tbb::parallel_reduce(tbb::blocked_range<int>(0, simulation.dv.cols(), 256), (T)0,
            [&](const tbb::blocked_range<int>& range, T ns) -> T {
                int start = range.begin();
                int length = (range.end() - range.begin());
                const auto& dv_block = simulation.dv.middleCols(start, length);
                const auto& mass_block = simulation.mass_matrix.segment(start, length);
                return ns - simulation.dt * ((simulation.gravity.transpose() * dv_block).dot(mass_block));
            },
            [](T a, T b) -> T { return a + b; });
        return result;
    }

    void computeStep(TVStack& ddv, const TVStack& residual, const T linear_solve_relative_tolerance)
    {
        ZIRAN_ASSERT(ddv.cols() == simulation.num_nodes);
        ddv.setZero();
        if (!matrix_free) {
            buildMatrix();
            if (lu) {
                // integrate boundary condition to the matrix
                /*
                TVStack proj_finder;
                proj_finder.resizeLike(ddv);
                proj_finder.fill(1);
                project(proj_finder);
                for (int j = 0; j < proj_finder.cols(); ++j) {
                    for (int i = 0; i < dim; ++i) {
                        if (proj_finder(i, j) == 0) {
                            // note residual is already projected
                            // TODO: sparse matrix doesn't have a convenient `fill'
                            // implicit_system_matrix.col(j*dim + i).fill(0);
                            // implicit_system_matrix.row(j*dim + i).fill(0);
                            implicit_system_matrix(j*dim + i, j*dim + i) = 1;
                        }
                    }
                }
                */
                eigenSparseLUSolve(implicit_system_matrix, EIGEN_EXT::vec(residual), EIGEN_EXT::vec(ddv));
                TVStack linear_solve_res;
                linear_solve_res.resizeLike(ddv);
                multiply(ddv, linear_solve_res);
                linear_solve_res -= residual;
                ZIRAN_DEBUG(computeNorm(linear_solve_res));
                return;
            }
        }

        if (simulation.full_implicit) {
            gmres.setRelativeTolerance(linear_solve_relative_tolerance);
            gmres.solve(*this, ddv, residual, simulation.verbose);
        }
        else {
            linear_solver->setRelativeTolerance(linear_solve_relative_tolerance);
            linear_solver->solve(*this, ddv, residual, simulation.verbose);
        }
    }

    T computeDescentStep(TVStack& ddv, const TVStack& residual, const T linear_solve_relative_tolerance)
    {
        computeStep(ddv, residual, linear_solve_relative_tolerance);

        T inner_prod = innerProduct(ddv, residual);
        T tol = 1e-5 * computeNorm(ddv) * computeNorm(residual);
        if (inner_prod < -tol) {
            ZIRAN_INFO("Newton step direction is of same direction as gradient. Look backwards!!");
            inner_prod = -inner_prod;
            ddv = -ddv;
        }
        else if (std::abs(inner_prod) < tol) {
            // ZIRAN_DEBUG(simulation.mass_matrix);
            // ZIRAN_DEBUG(ddv);
            // ZIRAN_DEBUG(residual);
            // ZIRAN_DEBUG(inner_prod);

            ZIRAN_INFO("Newton step direciton is almost orthogonal to gradient. Look towards negative gradient direction!");
            T res_norm = computeNorm(residual);
            T scale = computeNorm(ddv) / res_norm;
            ddv = scale * residual;
            inner_prod = scale * res_norm * res_norm;
        }
        else if (!std::isfinite(inner_prod)) {
            ddv = residual;
            T res_norm = computeNorm(residual);
            inner_prod = res_norm * res_norm;
        }
        return inner_prod;
    }

    /**
       This sets up the sparsity pattern in the implicit time stepping system matrix. This sparsity pattern will not change over time, although it's entries will.
    */
    void initializeImplicitTimeSteppingMatrix()
    {
        assert(!matrix_free);

        typedef Eigen::Triplet<T> Triplet;
        StdVector<Triplet> tripletList;
        implicit_system_matrix.resize(simulation.num_nodes * dim, simulation.num_nodes * dim);

        if (!simulation.quasistatic)
            simulation.inertia->initializeStiffnessSparsityPattern(tripletList);

        //add in the stiffness entries
        for (auto& lf : simulation.forces)
            lf->initializeStiffnessSparsityPattern(tripletList);

        //initialize from triplet list
        implicit_system_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
        implicit_system_matrix.makeCompressed();

        for (auto& lf : simulation.forces)
            lf->updateStiffnessSparsityPatternBasedState(implicit_system_matrix);
    }

    void buildMatrix()
    {
        assert(!matrix_free);
        using MATH_TOOLS::sqr;

        //update the matrix
        for (int k = 0; k < implicit_system_matrix.outerSize(); ++k)
            for (typename SpMat::InnerIterator it(implicit_system_matrix, k); it; ++it) {
                it.valueRef() = (T)0;
            }
        if (!simulation.quasistatic)
            simulation.inertia->addScaledStiffnessEntries(sqr(simulation.dt), implicit_system_matrix);
        //stiffness matrix
        for (auto& lf : simulation.forces) {
            lf->addScaledStiffnessEntries(sqr(simulation.dt), implicit_system_matrix);
        }

        precondition = JacobiPreconditioner<TVStack>(implicit_system_matrix);
        // ZIRAN_DEBUG(implicit_system_matrix);
    }

    void multiply(const TVStack& x, TVStack& b) const
    {
        ZIRAN_ASSERT(x.cols() == simulation.num_nodes);
        ZIRAN_ASSERT(b.cols() == simulation.num_nodes);
        using EIGEN_EXT::vec;
        using MATH_TOOLS::sqr;
        if (matrix_free) {
            b.setZero();
            // Note the relationship H dx = - df, where H is the stiffness matrix
            if (!simulation.quasistatic)
                simulation.inertia->addScaledForceDifferential(-sqr(simulation.dt), x, b);

            simulation.addScaledForceDifferentials(-sqr(simulation.dt), x, b);

            // for (auto& lf : simulation.forces) {
            //     lf->addScaledForceDifferential(-sqr(simulation.dt), x, b);
            // }
        }
        else {
            ZIRAN_ASSERT(implicit_system_matrix.cols() == simulation.num_nodes * dim);
            ZIRAN_ASSERT(implicit_system_matrix.rows() == simulation.num_nodes * dim);
            tbb::parallel_for(tbb::blocked_range<int>(0, b.cols(), (b.cols() * dim) / 32),
                [&](const tbb::blocked_range<int>& range) {
                    int start = range.begin() * dim;
                    int length = (range.end() - range.begin()) * dim;
                    auto vec_b = vec(b);
                    vec_b.segment(start, length) = implicit_system_matrix.middleRows(start, length) * vec(x);
                });
        }
    }

    template <class Func>
    void setProjection(Func project_func)
    {
        project = project_func;
    }

    template <class Func>
    void setPreconditioner(Func preconditioner_func)
    {
        precondition = preconditioner_func;
    }
};
} // namespace ZIRAN

#endif
