#ifndef HYPERELASTIC_CONSTITUTIVE_MODEL_H
#define HYPERELASTIC_CONSTITUTIVE_MODEL_H
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/Math/Linear/DenseExt.h>

namespace ZIRAN {

template <class Derived>
class HyperelasticConstitutiveModel;

template <typename Derived>
struct HyperelasticTraits;

template <template <class, int> class Derived, class T, int dim>
class HyperelasticConstitutiveModel<Derived<T, dim>> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    using Hessian = Eigen::Matrix<T, dim * dim, dim * dim>;
    using Scratch = typename HyperelasticTraits<Derived<T, dim>>::ScratchType;
    using Scalar = T;

    static constexpr bool diagonalDifferentiable() { return false; }

    T psi(const TM& F) const
    {
        Scratch s;
        static_cast<const Derived<T, dim>*>(this)->updateScratch(F, s);
        return static_cast<const Derived<T, dim>*>(this)->psi(s);
    }

    /**
      Overload me to compute the elastic potential
    */
    T psi(const Scratch& s) const
    {
        return 0;
    }

    void kirchhoff(const TM& F, TM& tau) const
    {
        Scratch s;
        static_cast<const Derived<T, dim>*>(this)->updateScratch(F, s);
        static_cast<const Derived<T, dim>*>(this)->kirchhoff(s, tau);
    }

    void firstPiola(const TM& F, TM& P) const
    {
        Scratch s;
        static_cast<const Derived<T, dim>*>(this)->updateScratch(F, s);
        static_cast<const Derived<T, dim>*>(this)->firstPiola(s, P);
    }

    void firstPiolaDifferential(const TM& F, const TM& dF, TM& dP) const
    {
        Scratch s;
        static_cast<const Derived<T, dim>*>(this)->updateScratch(F, s);
        static_cast<const Derived<T, dim>*>(this)->firstPiolaDifferential(s, dF, dP);
    }

    void firstPiolaDerivative(const TM& F, Hessian& dPdF) const
    {
        Scratch s;
        static_cast<const Derived<T, dim>*>(this)->updateScratch(F, s);
        static_cast<const Derived<T, dim>*>(this)->firstPiolaDerivative(s, dPdF);
    }

    /**
      Overload at least one of
      void kirchhoff(const Scratch& s, TM& tau) const
      void firstPiola(const Scratch& s, TM& P) const
      */

    void kirchhoff(const Scratch& s, TM& tau) const
    {
        TM P;
        static_cast<const Derived<T, dim>*>(this)->firstPiola(s, P);
        tau.noalias() = P * s.F.transpose();
    }

    void firstPiola(const Scratch& s, TM& P) const
    {
        TM tau;
        static_cast<const Derived<T, dim>*>(this)->kirchhoff(s, tau);
        P.noalias() = tau * s.FinvT;
    }

    /**
       Overload at least one of
       firstPiolaDifferential(TM& dP,const TM& dF, const Scratch& s)
       firstPiolaDerivative(Hessian& dPdF, const Scratch& s)
       */
    void firstPiolaDifferential(const Scratch& s, const TM& dF, TM& dP) const
    {
        Hessian dPdF;
        static_cast<const Derived<T, dim>*>(this)->firstPiolaDerivative(s, dPdF);
        EIGEN_EXT::contract(dP, dPdF, dF);
    }

    /**
       Overload at least one of
       firstPiolaDifferential(TM& dP,const TM& dF, const Scratch& s)
       firstPiolaDerivative(Hessian& dPdF, const Scratch& s)
       */
    void firstPiolaDerivative(const Scratch& s, Hessian& dPdF) const
    {
        TM dP;
        TM dF = TM::Zero();
        for (int i = 0; i < TM::SizeAtCompileTime; i++) {
            dF(i) = 1;
            static_cast<const Derived<T, dim>*>(this)->firstPiolaDifferential(s, dF, dP);
            dPdF.col(i) = EIGEN_EXT::vec(dP);
            dF(i) = 0;
        }
    }

    /**
      Returns whether the constitutive model is C2 in a neighborhood of radius tol
      about the current state
      */
    bool isC2(const TM& F, T tol) const
    {
        Scratch s;
        static_cast<const Derived<T, dim>*>(this)->updateScratch(F, s);
        return static_cast<const Derived<T, dim>*>(this)->isC2(tol, s);
    }

    /**
      Returns whether dP (or dPdF) is implemented
      */
    bool hessianImplemented() const
    {
        return true;
    }
};
} // namespace ZIRAN
#endif
