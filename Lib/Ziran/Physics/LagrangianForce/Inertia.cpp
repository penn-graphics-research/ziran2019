#include <Ziran/CS/Util/Logging.h>
#include <Ziran/Physics/LagrangianForce/Inertia.h>

namespace ZIRAN {

template <class T, int dim>
MassLumpedInertia<T, dim>::MassLumpedInertia(const Vec& mass_matrix, const TVStack& dv, const T& dt)
    : mass_matrix(mass_matrix)
    , dv(dv)
    , dt(dt)
{
}

// E
template <class T, int dim>
T MassLumpedInertia<T, dim>::totalEnergy() const
{
    T ke = tbb::parallel_reduce(tbb::blocked_range<int>(0, dv.cols(), 256), (T)0,
        [&](const tbb::blocked_range<int>& range, T ns) -> T {
            int start = range.begin();
            int length = (range.end() - range.begin());
            const auto& dv_block = dv.middleCols(start, length);
            const auto& mass_block = mass_matrix.segment(start, length).transpose();
            return ns + ((dv_block.colwise().squaredNorm()).array() * mass_block.array()).sum();
        },
        [](T a, T b) -> T { return a + b; });
    // ZIRAN_INFO("Inertia Called ", ke);
    return ke / 2;
}

// f= -dE/dx
template <class T, int dim>
void MassLumpedInertia<T, dim>::addScaledForces(T scale, TVStack& forces) const
{
    scale /= dt; // dE/dx = ddv/dx dE/ddv = 1/dt dE/ddv
    tbb::parallel_for(tbb::blocked_range<int>(0, forces.cols(), 256),
        [&](const tbb::blocked_range<int>& range) {
            for (int p = range.begin(), s = range.end(); p < s; p++)
                forces.col(p) -= scale * mass_matrix(p) * dv.col(p);
        });
}

// df
template <class T, int dim>
void MassLumpedInertia<T, dim>::addScaledForceDifferential(T scale, const TVStack& dx, TVStack& df) const
{
    scale /= dt * dt;
    tbb::parallel_for(tbb::blocked_range<int>(0, df.cols(), 256),
        [&](const tbb::blocked_range<int>& range) {
            for (int p = range.begin(), s = range.end(); p < s; p++)
                df.col(p) -= scale * mass_matrix(p) * dx.col(p);
        });
}

// -df/dx
template <class T, int dim>
void MassLumpedInertia<T, dim>::addScaledStiffnessEntries(T scale, Eigen::SparseMatrix<T, Eigen::RowMajor>& newton_matrix) const
{
    //mass matrix
    scale /= dt * dt;
    for (int p = 0; p < mass_matrix.rows(); p++) {
        for (int d = 0; d < dim; d++) {
            int index = p * dim + d;
            newton_matrix.coeffRef(index, index) += scale * mass_matrix(p);
        }
    }
}

template <class T, int dim>
void MassLumpedInertia<T, dim>::initializeStiffnessSparsityPattern(StdVector<Eigen::Triplet<T>>& tripletList) const
{
    tripletList.reserve(tripletList.size() + dim * mass_matrix.rows());
    for (int p = 0; p < dim * mass_matrix.rows(); p++)
        tripletList.emplace_back(p, p, (T)0);
}
template class MassLumpedInertia<double, 2>;
template class MassLumpedInertia<double, 3>;
template class MassLumpedInertia<float, 2>;
template class MassLumpedInertia<float, 3>;
} // namespace ZIRAN
