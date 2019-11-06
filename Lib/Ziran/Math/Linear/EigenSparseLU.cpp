#include <Ziran/CS/Util/Debug.h>
#include <Eigen/SparseLU>
#include "EigenSparseLU.h"

namespace ZIRAN {

template <class T>
void eigenSparseLUSolve(const Eigen::SparseMatrix<T, Eigen::RowMajor>& mat, Eigen::Map<const Vector<T, Eigen::Dynamic>> rhs, Eigen::Map<Vector<T, Eigen::Dynamic>> to_solve)
{
    using SpMatCol = Eigen::SparseMatrix<T, Eigen::ColMajor>;
    SpMatCol copy = mat;
    Eigen::SparseLU<SpMatCol> solver;
    solver.analyzePattern(copy);
    solver.factorize(copy);
    to_solve = solver.solve(rhs);
    ZIRAN_VERB_IF(solver.info() == Eigen::NumericalIssue, "NumericalIssue");
    ZIRAN_VERB_IF(solver.info() == Eigen::InvalidInput, "InvalidInput");
    ZIRAN_VERB_IF(solver.info() == Eigen::NoConvergence, "NoConvergence");
    ZIRAN_ASSERT(solver.info() == Eigen::Success);
}

template void eigenSparseLUSolve<double>(const Eigen::SparseMatrix<double, Eigen::RowMajor>&, Eigen::Map<const Vector<double, Eigen::Dynamic>>, Eigen::Map<Vector<double, Eigen::Dynamic>>);
template void eigenSparseLUSolve<float>(const Eigen::SparseMatrix<float, Eigen::RowMajor>&, Eigen::Map<const Vector<float, Eigen::Dynamic>>, Eigen::Map<Vector<float, Eigen::Dynamic>>);
} // namespace ZIRAN
