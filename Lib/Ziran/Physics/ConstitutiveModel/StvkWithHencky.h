#ifndef STVK_WITH_HENCKY_H
#define STVK_WITH_HENCKY_H
#include <iostream>
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/Util/BinaryIO.h>
#include <Ziran/Physics/ConstitutiveModel/HyperelasticConstitutiveModel.h>

namespace ZIRAN {

template <typename Derived>
struct HyperelasticTraits;

// scratch (non-state) variables for the consitutive model
template <class T, int dim>
struct StvkWithHenckyScratch {
    using TM = Matrix<T, dim, dim>;
    using TV = Vector<T, dim>;
    TM F, U, V;
    TV sigma;
    TV log_sigma;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    StvkWithHenckyScratch()
        : F(TM::Identity())
    {
    }

    static const char* name()
    {
        return "StvkWithHenckyScratch";
    }
};

template <class T, int _dim>
class StvkWithHencky : public HyperelasticConstitutiveModel<StvkWithHencky<T, _dim>> {
public:
    static const int dim = _dim;
    using Base = HyperelasticConstitutiveModel<StvkWithHencky<T, dim>>;
    using TM = typename Base::TM;
    using TV = typename Base::TV;

    using Strain = TM;
    using Hessian = typename Base::Hessian;
    using Scalar = typename Base::Scalar;
    using Scratch = typename HyperelasticTraits<StvkWithHencky<T, dim>>::ScratchType;
    using Base::hessianImplemented;
    using Vec = Vector<T, Eigen::Dynamic>;
    using VecBlock = Eigen::VectorBlock<Vec>;

    T mu, lambda;

    StvkWithHencky(const T E = (T)1, const T nu = (T)0.3);

    void setLameParameters(const T E, const T nu);

    Matrix<T, 2, 2> Bij(const TV& sigma, int i, int j, T clamp_value) const;

    Matrix<T, 2, 2> BijFull(const TV& zi, const TV& sigma, const T& second_term, int i, int j, T clamp_value) const;

    void updateScratch(const TM& new_F, Scratch& scratch) const;

    static constexpr bool diagonalDifferentiable()
    {
        return true;
    }

    /**
       psi = mu tr((log S)^2) + 1/2 lambda (tr(log S))^2
     */
    T psi(const Scratch& s) const;

    /**
       P = U (2 mu S^{-1} (log S) + lambda tr(log S) S^{-1}) V^T
     */
    void firstPiola(const Scratch& s, TM& P) const;

    void firstPiolaDerivative(const Scratch& s, Hessian& dPdF) const;

    T psiDiagonal(const TV& sigma) const;

    TV firstPiolaDiagonal(const TV& sigma) const;

    TM firstPiolaDerivativeDiagonal(const TV& sigma) const;

    void firstPiolaDifferential(const Scratch& s, const TM& dF, TM& dP) const;

    bool isC2(const Scratch& s, T tolerance) const;

    /**
       Returns whether dP (or dPdF) is implemented
    */
    bool hessianImplemented() const;

    void write(std::ostream& out) const;

    static StvkWithHencky<T, _dim> read(std::istream& in);

    static const char* name() { return "StvkWithHencky"; }

    static const char* scratch_name() { return Scratch::name(); }
};

template <class T, int dim>
struct HyperelasticTraits<StvkWithHencky<T, dim>> {
    using ScratchType = StvkWithHenckyScratch<T, dim>;
};

template <class T, int dim>
struct RW<StvkWithHenckyScratch<T, dim>> {
    using Tag = NoWriteTag<StvkWithHenckyScratch<T, dim>>;
};
} // namespace ZIRAN

#endif
