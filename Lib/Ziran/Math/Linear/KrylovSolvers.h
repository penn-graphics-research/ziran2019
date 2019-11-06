#ifndef KRYLOV_SOLVERS_H
#define KRYLOV_SOLVERS_H
#include <Ziran/Math/Linear/DenseExt.h>
#include <Ziran/Math/Linear/Givens.h>
#include <Ziran/CS/Util/Timer.h>
#include <Ziran/CS/Util/Logging.h>
#include <Eigen/Core>

namespace ZIRAN {

template <class TV>
struct DirichletDofProjection {

    const StdVector<int>& dirichlet_dofs;
    DirichletDofProjection(const StdVector<int>& dirichlet_dofs)
        : dirichlet_dofs(dirichlet_dofs)
    {
    }

    void operator()(TV& v)
    {
        using T = ScalarType<TV>;
        Eigen::Map<Vector<T, Eigen::Dynamic>> v_vec = EIGEN_EXT::vec(v);
        for (unsigned int i = 0; i < dirichlet_dofs.size(); i++) {
            int index = dirichlet_dofs[i];
            v_vec(index) = (T)0;
        }
    }
};

template <class TV>
struct JacobiPreconditioner {

    TV inverse_diagonal;

    template <class TM>
    JacobiPreconditioner(const TM& mat, typename TV::Scalar tol = 16 * std::numeric_limits<typename TV::Scalar>::epsilon())
    {
        const int rows = TV::RowsAtCompileTime;
        const int cols = TV::ColsAtCompileTime;
        using T = typename TV::Scalar;
        static_assert(rows != Eigen::Dynamic || cols != Eigen::Dynamic, "rows and columns can't both be dynamic");
        if (rows == Eigen::Dynamic)
            inverse_diagonal.resize(mat.cols() / cols, cols);
        if (cols == Eigen::Dynamic)
            inverse_diagonal.resize(rows, mat.rows() / rows);
        EIGEN_EXT::vec(inverse_diagonal) = mat.diagonal().unaryExpr(
            [tol](T x) -> T {
                return (T)1 / std::max(tol, std::abs(x));
            });
    }

    void operator()(const TV& in, TV& out)
    {
        out = in.array() * inverse_diagonal.array();
    }
};

template <class T, class TV>
class EigenSparseKrylovMatrix {
public:
    using TM = Eigen::SparseMatrix<T, Eigen::RowMajor>;
    TM& A;
    std::function<void(TV&)> project;
    std::function<void(const TV&, TV&)> precondition;

    // StdVector<int>* dirichlet_dofs;

    EigenSparseKrylovMatrix(TM& A)
        : A(A)
        , project([](TV&) {})
        , precondition([](const TV& x, TV& b) { b = x; })

    {
    }

    template <class Func1, class Func2>
    EigenSparseKrylovMatrix(TM& A, Func1 project, Func2 precondition)
        : A(A)
        , project(project)
        , precondition(precondition)
    {
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

    void multiply(const TV& x, TV& b) const
    {
        using EIGEN_EXT::vec;
        vec(b).noalias() = A * vec(x);
    }
};
} // namespace ZIRAN

#endif
