#ifndef PHASE_FIELD_H
#define PHASE_FIELD_H

#include <Ziran/Math/MathTools.h>
#include <Ziran/Physics/ConstitutiveModel/NeoHookeanBorden.h>
#include <Ziran/CS/Util/BinaryIO.h>

namespace ZIRAN {

template <class T, int dim>
class PhaseField {
public:
    T residual_phase;
    T c; // 1 is undamaged, 0 is fully damaged
    T H;
    T l0;
    T one_over_sigma_c;
    T pf_Fp;
    T H_max;
    T vol;
    bool allow_damage;

    PhaseField() {}

    PhaseField(const T residual_phase, const T c, const T H, const T l0, const T one_over_sigma_c, const T pf_Fp, const T H_max, const T vol, const bool allow_damage)
        : residual_phase(residual_phase), c(c), H(H), l0(l0), one_over_sigma_c(one_over_sigma_c), pf_Fp(pf_Fp), H_max(H_max), vol(vol), allow_damage(allow_damage)
    {
    }

    // Maybe I don't need read/write because the data are trivial.
    void write(std::ostream& out) const
    {
        writeEntry(out, residual_phase);
        writeEntry(out, c);
        writeEntry(out, H);
        writeEntry(out, l0);
        writeEntry(out, one_over_sigma_c);
        writeEntry(out, pf_Fp);
        writeEntry(out, H_max);
        writeEntry(out, vol);
        writeEntry(out, allow_damage);
    }
    static PhaseField<T, dim> read(std::istream& in)
    {
        PhaseField<T, dim> target;
        target.residual_phase = readEntry<T>(in);
        target.c = readEntry<T>(in);
        target.H = readEntry<T>(in);
        target.l0 = readEntry<T>(in);
        target.one_over_sigma_c = readEntry<T>(in);
        target.pf_Fp = readEntry<T>(in);
        target.H_max = readEntry<T>(in);
        target.vol = readEntry<T>(in);
        target.allow_damage = readEntry<bool>(in);
        return target;
    }

    template <class TCONST>
    static T Get_Sigma_C_Based_On_Max_Deformation(const T percentage, const TCONST& model)
    {
        typename NeoHookeanBorden<T, dim>::Scratch s;
        s.F = Matrix<T, dim, dim>::Identity() * ((T)1 + percentage);
        s.J = s.F.determinant();
        T G = model.psi(s);
        T one_over_sigma_c = (T)1 / G;
        ZIRAN_INFO("G Value: ", G);
        return one_over_sigma_c;
    }

    void Update_Phase_Field_Fp(const T psi_pos)
    {
        if (psi_pos > H) {
            H = std::min(psi_pos, H_max);
            pf_Fp = 4 * l0 * (1 - residual_phase) * H * one_over_sigma_c + 1;
        }
    }

    static const char* name()
    {
        return "PhaseField";
    }
};

} // namespace ZIRAN

#endif