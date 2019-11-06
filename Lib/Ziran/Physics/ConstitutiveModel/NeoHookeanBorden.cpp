#include <Ziran/Physics/ConstitutiveModel/NeoHookeanBorden.h>
#include <Ziran/Physics/ConstitutiveModel/HyperelasticConstitutiveModel.h>
#include <Ziran/Math/Linear/DenseExt.h>
#include <Ziran/CS/Util/BinaryIO.h>
#include <cmath>

namespace ZIRAN {

template <class T, int _dim>
NeoHookeanBorden<T, _dim>::NeoHookeanBorden(const T E, const T nu)
{
    setLameParameters(E, nu);
}

template <class T, int _dim>
void NeoHookeanBorden<T, _dim>::setLameParameters(const T E, const T nu)
{
    T lambda = E * nu / (((T)1 + nu) * ((T)1 - (T)2 * nu));
    mu = E / ((T)2 * ((T)1 + nu));
    kappa = (T).666666 * mu + lambda;
}

template <class T, int _dim>
void NeoHookeanBorden<T, _dim>::updateScratch(const TM& new_F, Scratch& scratch)
{
    using std::log;
    using namespace EIGEN_EXT;
    scratch.F = new_F;
    scratch.J = scratch.F.determinant();

    TM JFinvT;
    EIGEN_EXT::cofactorMatrix(scratch.F, JFinvT);
    scratch.FinvT = ((T)1 / scratch.J) * JFinvT;
    scratch.logJ = log(scratch.J);
}

template <class T, int _dim>
T NeoHookeanBorden<T, _dim>::psi(const Scratch& s) const
{
    // step 1 compute psi pos
    TM JaF = std::pow(s.J, -(T)1 / (T)dim) * s.F;
    T psi_dev = mu / 2 * ((JaF.transpose() * JaF).trace() - dim);
    T psi_vol = kappa / 2 * ((s.J * s.J - 1) / 2 - log(s.J));
    T psi_pos = (s.J >= 1) ? (psi_vol + psi_dev) : psi_dev;

    // step 2 compute psi neg
    T psi_neg = (s.J >= 1) ? (T)0 : psi_vol;

    // step 3
    return g * psi_pos + psi_neg;
}

template <class T, int _dim>
T NeoHookeanBorden<T, _dim>::get_psi_pos(const TM& F)
{
    Scratch s;
    updateScratch(F, s);

    // step 1 compute psi pos
    TM JaF = std::pow(s.J, -(T)1 / (T)dim) * s.F;
    T psi_dev = mu / 2 * ((JaF.transpose() * JaF).trace() - dim);
    T psi_vol = kappa / 2 * ((s.J * s.J - 1) / 2 - log(s.J));
    T psi_pos = (s.J >= 1) ? (psi_vol + psi_dev) : psi_dev;

    return psi_pos;
}

template <class T, int _dim>
void NeoHookeanBorden<T, _dim>::kirchhoff(const Scratch& s, TM& tau) const
{
    // T scale = lambda * s.logJ - mu;
    // tau = (mu * s.F * s.F.transpose() + scale * TM::Identity());

    TM B = s.F * s.F.transpose();
    TM devB = B - TM::Identity() * (T)1 / dim * B.trace();
    TM tau_dev = mu * std::pow(s.J, -(T)2 / (T)dim) * devB;

    T prime = kappa / 2 * (s.J - 1 / s.J);
    TM tau_vol = s.J * prime * TM::Identity();

    if (s.J >= 1)
        tau = g * (tau_dev + tau_vol);
    else
        tau = g * tau_dev + tau_vol;
}

template <class T, int _dim>
void NeoHookeanBorden<T, _dim>::firstPiolaDifferential(const Scratch& s, const TM& dF, TM& dP) const
{
    using namespace EIGEN_EXT;

    T deltaJ = s.J * contractMatrices(s.FinvT, dF);
    T FcontractF = contractMatrices(s.F, s.F);
    TM deltaFinvT = -s.FinvT * dF.transpose() * s.FinvT;

    T a = -(T)1 / dim;
    T muJ2a = mu * std::pow(s.J, (T)2 * a);

    T prime = kappa / 2 * (s.J - (T)1 / s.J);
    T prime2 = kappa / 2 * ((T)1 + (T)1 / (s.J * s.J));

    TM dP_dev = 2 * a * muJ2a * contractMatrices(s.FinvT, dF) * (a * FcontractF * s.FinvT + s.F);
    dP_dev += muJ2a * (2 * a * contractMatrices(s.F, dF) * s.FinvT + a * FcontractF * deltaFinvT + dF);

    TM dP_vol = deltaJ * prime * s.FinvT + s.J * deltaJ * prime2 * s.FinvT + s.J * prime * deltaFinvT;

    if (s.J >= 1)
        dP = g * (dP_dev + dP_vol);
    else
        dP = g * dP_dev + dP_vol;
}

template <class T, int _dim>
bool NeoHookeanBorden<T, _dim>::isC2(const Scratch& s, T tolerance) const
{
    return s.J > tolerance;
}

template <class T, int _dim>
void NeoHookeanBorden<T, _dim>::write(std::ostream& out) const
{
    writeEntry(out, mu);
    writeEntry(out, kappa);
    writeEntry(out, psi_pos);
    writeEntry(out, g);
}

template <class T, int _dim>
NeoHookeanBorden<T, _dim> NeoHookeanBorden<T, _dim>::read(std::istream& in)
{
    NeoHookeanBorden<T, _dim> model;
    model.mu = readEntry<T>(in);
    model.kappa = readEntry<T>(in);
    model.psi_pos = readEntry<T>(in);
    model.g = readEntry<T>(in);
    return model;
}

template <class T, int _dim>
bool NeoHookeanBorden<T, _dim>::hessianImplemented() const
{
    return true;
}

template class NeoHookeanBorden<double, 1>;
template class NeoHookeanBorden<double, 2>;
template class NeoHookeanBorden<double, 3>;
template class NeoHookeanBorden<float, 1>;
template class NeoHookeanBorden<float, 2>;
template class NeoHookeanBorden<float, 3>;
} // namespace ZIRAN
