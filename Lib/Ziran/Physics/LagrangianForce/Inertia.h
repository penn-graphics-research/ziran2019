#ifndef INERTIA_H
#define INERTIA_H
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/Physics/LagrangianForce/LagrangianForce.h>
#include <tbb/tbb.h>

namespace ZIRAN {

template <class T, int dim>
class LagrangianForce;

/**
This is the class for computing the mass lumped discretization of inertia
**/
template <class T, int dim>
class MassLumpedInertia : public LagrangianForce<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    using TVStack = Matrix<T, dim, Eigen::Dynamic>;
    using Vec = Vector<T, Eigen::Dynamic>;

    const Vec& mass_matrix;
    const TVStack& dv;
    const T& dt;

    MassLumpedInertia(const Vec& mass_matrix, const TVStack& dv, const T& dt);

    // E
    T totalEnergy() const;

    // f= -dE/dx
    void addScaledForces(T scale, TVStack& forces) const;

    // df
    void addScaledForceDifferential(T scale, const TVStack& dx, TVStack& df) const;

    // -df/dx
    void addScaledStiffnessEntries(T scale, Eigen::SparseMatrix<T, Eigen::RowMajor>& newton_matrix) const;

    void initializeStiffnessSparsityPattern(StdVector<Eigen::Triplet<T>>& tripletList) const;
};
} // namespace ZIRAN
#endif
