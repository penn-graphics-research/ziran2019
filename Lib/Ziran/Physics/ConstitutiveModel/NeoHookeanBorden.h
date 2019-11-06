#ifndef NEO_HOOKEAN_BORDEN_H
#define NEO_HOOKEAN_BORDEN_H
#include <iostream>
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/Util/BinaryIO.h>
#include <Ziran/Physics/ConstitutiveModel/HyperelasticConstitutiveModel.h>

namespace ZIRAN {

template <class T, int dim>
struct NeoHookeanBordenScratch {
    using TM = Matrix<T, dim, dim>;

    T J, logJ;
    TM F, FinvT;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    NeoHookeanBordenScratch()
        : F(TM::Identity())
    {
    }

    static const char* name()
    {
        return "NeoHookeanBordenScratch";
    }
};

template <class T, int _dim>
class NeoHookeanBorden : public HyperelasticConstitutiveModel<NeoHookeanBorden<T, _dim>> {
public:
    static const int dim = _dim;
    using Base = HyperelasticConstitutiveModel<NeoHookeanBorden<T, dim>>;
    using TV = typename Base::TV;
    using TM = typename Base::TM;
    using Strain = TM;
    using Hessian = typename Base::Hessian;
    using Scalar = typename Base::Scalar;
    using Scratch = typename HyperelasticTraits<NeoHookeanBorden<T, dim>>::ScratchType;

    T mu, kappa;

    T psi_pos;

    T g = 1;

    NeoHookeanBorden(const T E = (T)1, const T nu = (T)0.3);

    void setLameParameters(const T E, const T nu);

    void updateScratch(const TM& new_F, Scratch& scratch);

    static constexpr bool diagonalDifferentiable() { return true; }

    T psi(const Scratch& s) const;

    T get_psi_pos(const TM& s);

    void kirchhoff(const Scratch& s, TM& tau) const;

    //void firstPiolaDerivative(const Scratch& s, Hessian& dPdF) const;

    void firstPiolaDifferential(const Scratch& s, const TM& dF, TM& dP) const;

    bool isC2(const Scratch& s, T tolerance) const;

    bool hessianImplemented() const;

    void write(std::ostream& out) const;

    static NeoHookeanBorden<T, _dim> read(std::istream& in);

    static const char* name() { return "NeoHookeanBorden"; }

    static const char* scratch_name() { return Scratch::name(); }
};

template <class T, int dim>
struct HyperelasticTraits<NeoHookeanBorden<T, dim>> {
    using ScratchType = NeoHookeanBordenScratch<T, dim>;
};

template <class T, int dim>
struct RW<NeoHookeanBordenScratch<T, dim>> {
    using Tag = NoWriteTag<NeoHookeanBordenScratch<T, dim>>;
};
} // namespace ZIRAN

#endif
