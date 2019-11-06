#include "PlasticityApplier.h"

namespace ZIRAN {
///////////////////////////////////////////////////////////////////////////////
/**
   This is NonAssociativeCamClay
 */
///////////////////////////////////////////////////////////////////////////////

template <class T>
NonAssociativeCamClay<T>::NonAssociativeCamClay(T logJp, T friction_angle, T beta, T xi, int dim, bool hardeningOn)
    : logJp(logJp)
    , beta(beta)
    , xi(xi)
    , hardeningOn(hardeningOn)
{
    T sin_phi = std::sin(friction_angle / (T)180 * (T)3.141592653);
    T mohr_columb_friction = std::sqrt((T)2 / (T)3) * (T)2 * sin_phi / ((T)3 - sin_phi);
    M = mohr_columb_friction * (T)dim / std::sqrt((T)2 / ((T)6 - dim));
}

template <class T, int dim>
void Compare_With_Phybam_Numerical_Check(
    const Vector<T, dim>& Se1,
    const Vector<T, dim>& Se2,
    const T logJp1,
    const T logJp2)
{

    ZIRAN_ASSERT(std::abs(logJp1 - logJp2) < 1e-4, logJp1, logJp2);
    for (int i = 0; i < dim; ++i)
        ZIRAN_ASSERT(std::abs(Se1(i) - Se2(i)) < 1e-4, Se1(i), Se2(i));
}

template <class T, int dim>
void Compare_With_Physbam(
    const T mu,
    const T kappa,
    const T cam_clay_M,
    const T cam_clay_beta,
    const T p0,
    Vector<T, dim>& Se,
    T& cam_clay_logJp)
{
    using namespace EIGEN_EXT;
    using namespace MATH_TOOLS;
    //typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, dim> TV;

    T a = ((T)1 + (T)2 * cam_clay_beta) * ((T)6 - (T)dim) / (T)2;
    T b = cam_clay_beta * p0;
    T c = p0;
    T M2 = sqr(cam_clay_M);

    T Je = 1.;
    for (int i = 0; i < dim; ++i) Je *= Se(i);

    TV Se2;
    for (int i = 0; i < dim; ++i)
        Se2(i) = sqr(Se(i));

    TV s_hat_trial = mu * std::pow(Je, -(T)2 / (T)dim) * deviatoric(Se2);

    TV tau_hat;
    T Uprime = kappa / (T)2 * (Je - 1 / Je);
    T p_trial = -Uprime * Je;

    // Projecting to the tips
    if (p_trial > c) {
        T Je_new = sqrt(-2 * c / kappa + 1);
        Se = TV::Ones() * pow(Je_new, (T)1 / dim);
        cam_clay_logJp += log(Je / Je_new);
        return;
    }
    else if (p_trial < -b) {
        T Je_new = sqrt(2 * b / kappa + 1);
        Se = TV::Ones() * pow(Je_new, (T)1 / dim);
        cam_clay_logJp += log(Je / Je_new);
        return;
    }
#if 1
    T k = sqrt(-M2 * (p_trial + b) * (p_trial - c) / a);

    T s_hat_trial_norm = s_hat_trial.norm();
    T y = a * sqr(s_hat_trial_norm) + M2 * (p_trial + b) * (p_trial - c);

    if (y < 1e-4) return; // inside the yield surface

    // Fake hardening by computing intersection to center
    T pc = ((T)1 - cam_clay_beta) * p0 / 2;
    if (p0 > 1e-4 && p_trial < p0 - (1e-4) && p_trial > -cam_clay_beta * p0 + (1e-4)) {
        T aa = M2 * sqr(p_trial - pc) / (a * sqr(s_hat_trial_norm));
        T dd = 1 + aa;
        T ff = aa * cam_clay_beta * p0 - aa * p0 - 2 * pc;
        T gg = sqr(pc) - aa * cam_clay_beta * sqr(p0);
        T zz = std::sqrt(std::abs(sqr(ff) - 4 * dd * gg));
        T p1 = (-ff + zz) / (2 * dd);
        T p2 = (-ff - zz) / (2 * dd);
        T p_fake = (p_trial - pc) * (p1 - pc) > 0 ? p1 : p2;
        T Je_new_fake = sqrt(std::abs(-2 * p_fake / kappa + 1));
        if (Je_new_fake > 1e-4)
            cam_clay_logJp += log(Je / Je_new_fake);
    }

    TV be_new = k / mu * std::pow(Je, (T)2 / (T)dim) * s_hat_trial / s_hat_trial_norm + (T)1 / dim * Se2.sum() * TV::Ones();

    for (int i = 0; i < dim; ++i)
        Se(i) = std::sqrt(be_new(i));
#endif
}

template <class T>
template <class TConst>
bool NonAssociativeCamClay<T>::projectStrain(TConst& c, Matrix<T, TConst::dim, TConst::dim>& strain)
{
    using namespace EIGEN_EXT;
    static const int dim = TConst::dim;
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, dim> TV;
    TM U, V;
    TV sigma;

    // TODO: this is inefficient because next time step updateState will do the svd again!
    singularValueDecomposition(strain, U, sigma, V);

    T p0 = c.kappa * (T(0.00001) + std::sinh(xi * std::max(-logJp, (T)0)));

    // debug with physbam only

    TV XXXSe = sigma;
    T XXXcam_clay_logJp = logJp;

    Compare_With_Physbam(c.mu, c.kappa, M, beta, p0, XXXSe, XXXcam_clay_logJp);
    logJp = XXXcam_clay_logJp;
    strain = U * XXXSe.asDiagonal() * V.transpose();
    return true;

    T J = 1.;
    for (int i = 0; i < dim; ++i) J *= sigma(i);

    // step 0 compute the value of y
    TV B_hat_trial;
    for (int i = 0; i < dim; ++i)
        B_hat_trial(i) = sigma(i) * sigma(i);
    TV s_hat_trial = c.mu * std::pow(J, -(T)2 / (T)dim) * deviatoric(B_hat_trial);

    T prime = c.kappa / (T)2 * (J - 1 / J);
    T p_trial = -prime * J;

    T y_s_half_coeff = ((T)6 - dim) / (T)2 * ((T)1 + (T)2 * beta);
    T y_p_half = M * M * (p_trial + beta * p0) * (p_trial - p0);
    T y = y_s_half_coeff * s_hat_trial.squaredNorm() + y_p_half;

    // project sigma's onto the tips
    T p_min = beta * p0;
    if (p_trial > p0) {
        T Je_new = std::sqrt(-2 * p0 / c.kappa + 1);
        sigma = TV::Ones() * std::pow(Je_new, (T)1 / dim);
        Eigen::DiagonalMatrix<T, dim, dim> sigma_m(sigma);
        TM Fe = U * sigma_m * V.transpose();
        strain = Fe;
        if (hardeningOn) {
            logJp += log(J / Je_new);
        }
        if (hardeningOn) { Compare_With_Phybam_Numerical_Check(XXXSe, sigma, XXXcam_clay_logJp, logJp); }
        return false;
    }
    else if (p_trial < -p_min) {
        T Je_new = std::sqrt(2 * p_min / c.kappa + 1);
        sigma = TV::Ones() * std::pow(Je_new, (T)1 / dim);
        Eigen::DiagonalMatrix<T, dim, dim> sigma_m(sigma);
        TM Fe = U * sigma_m * V.transpose();
        strain = Fe;
        if (hardeningOn) {
            logJp += log(J / Je_new);
        }

        if (hardeningOn) { Compare_With_Phybam_Numerical_Check(XXXSe, sigma, XXXcam_clay_logJp, logJp); }
        return false;
    }

    if (y < 1e-4) return false;

    // project to yield surface
    TV B_hat_new = std::pow(J, (T)2 / (T)dim) / c.mu * std::sqrt(-y_p_half / y_s_half_coeff) * s_hat_trial / s_hat_trial.norm();
    B_hat_new += (T)1 / dim * B_hat_trial.sum() * TV::Ones();

    for (int i = 0; i < dim; ++i)
        sigma(i) = std::sqrt(B_hat_new(i));
    Eigen::DiagonalMatrix<T, dim, dim> sigma_m(sigma);
    TM Fe = U * sigma_m * V.transpose();
    strain = Fe;

    // step 2 hack the hardening by computing a fake delta_p

    if (p0 > 1e-4 && p_trial < p0 - 1e-4 && p_trial > 1e-4 - p_min) {
        T p_center = (p0 - p_min) * .5;
        T q_trial = std::sqrt(((T)6 - (T)dim) / (T)2) * s_hat_trial.norm();
        Vector<T, 2> direction;
        direction(0) = p_center - p_trial;
        direction(1) = 0 - q_trial;
        direction = direction / direction.norm();

        T C = M * M * (p_center + beta * p0) * (p_center - p0);
        T B = M * M * direction(0) * (2 * p_center - p0 + beta * p0);
        T A = M * M * direction(0) * direction(0) + (1 + 2 * beta) * direction(1) * direction(1);

        T l1 = (-B + std::sqrt(B * B - 4 * A * C)) / (2 * A);
        T l2 = (-B - std::sqrt(B * B - 4 * A * C)) / (2 * A);

        T p1 = p_center + l1 * direction(0);
        T p2 = p_center + l2 * direction(0);

        T p_fake = (p_trial - p_center) * (p1 - p_center) > 0 ? p1 : p2;
        T Je_new_fake = sqrt(std::abs(-2 * p_fake / c.kappa + 1));
        if (Je_new_fake > 1e-4 && hardeningOn) //only update logJp if hardening is turned on!
            logJp += log(J / Je_new_fake);
    }

    if (hardeningOn) { Compare_With_Phybam_Numerical_Check(XXXSe, sigma, XXXcam_clay_logJp, logJp); }
    return false;
}

//#########################################################################
// Function: fillAttributesToVec3
//
// This is for PartioIO purposes so we can dump out attributes out from the plasticity model.
//#########################################################################
template <class T>
void NonAssociativeCamClay<T>::fillAttributesToVec3(Vector<T, 3>& data)
{
    for (int d = 0; d < 3; d++) { data(d) = logJp; }
}

template <class T>
const char* NonAssociativeCamClay<T>::name()
{
    return "NonAssociativeCamClay";
}

///////////////////////////////////////////////////////////////////////////////
/**
   This is the Drucker Prager plasticity model from

   Drucker-Prager Elastoplasticity for Sand Animation,
   G. Klar, T. Gast, A. Pradhana, C. Fu, C. Schroeder, C. Jiang, J. Teran,
   ACM Transactions on Graphics (SIGGRAPH 2016).

   It assumes the StvkHencky elasticity consitutive model.
 */
///////////////////////////////////////////////////////////////////////////////
template <class T>
DruckerPragerStvkHencky<T>::DruckerPragerStvkHencky(const T friction_angle, const T beta, const T cohesion, const bool volume_correction)
    : beta(beta)
    , logJp(0)
    , cohesion(cohesion)
    , volume_correction(volume_correction)
{
    T sin_phi = std::sin(friction_angle / (T)180 * (T)3.141592653);
    alpha = std::sqrt((T)2 / (T)3) * (T)2 * sin_phi / ((T)3 - sin_phi);
}

template <class T>
void DruckerPragerStvkHencky<T>::setParameters(const T friction_angle_in, const T beta_in, const T cohesion_in)
{
    beta = beta_in;
    cohesion = cohesion_in;
    T sin_phi = std::sin(friction_angle_in / (T)180 * (T)3.141592653);
    alpha = std::sqrt((T)2 / (T)3) * (T)2 * sin_phi / ((T)3 - sin_phi);
}

// bool is for fracture
template <class T>
template <class TConst>
bool DruckerPragerStvkHencky<T>::projectStrain(TConst& c, Matrix<T, TConst::dim, TConst::dim>& strain)
{
    static const int dim = TConst::dim;
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, dim> TV;
    TM U, V;
    TV sigma;

    // TODO: this is inefficient because next time step updateState will do the svd again!
    singularValueDecomposition(strain, U, sigma, V);

    TV epsilon = sigma.array().abs().max(1e-4).log() - cohesion;
    T trace_epsilon = epsilon.sum() + logJp;
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();

    if (trace_epsilon >= (T)0) // case II: project to tip
    {
        strain = U * std::exp(cohesion) * V.transpose();
        if (volume_correction) {
            logJp = beta * epsilon.sum() + logJp;
        }
        return false;
    }
    else if (c.mu == 0) {
        return false;
    }
    logJp = 0;
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm + (dim * c.lambda + 2 * c.mu) / (2 * c.mu) * trace_epsilon * alpha;
    TV H;
    if (delta_gamma <= 0) // case I: inside yield surface
    {
        H = epsilon + TV::Constant(cohesion);
    }
    else {
        H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat + TV::Constant(cohesion); // case III: projection
    }
    TV exp_H = H.array().exp();
    strain = U * exp_H.asDiagonal() * V.transpose();
    return false;
}

// bool is for fracture
template <class T>
template <class TConst>
Vector<T, TConst::dim> DruckerPragerStvkHencky<T>::projectSigma(TConst& c, const Vector<T, TConst::dim>& sigma)
{
    static const int dim = TConst::dim;
    typedef Vector<T, dim> TV;
    TV epsilon = sigma.array().abs().max(1e-4).log() - cohesion;
    T trace_epsilon = epsilon.sum();
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
    if (trace_epsilon >= 0) // case II: project to tip
    {
        TV ret = std::exp(cohesion) * TV::Ones();
        return ret;
    }
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm + (dim * c.lambda + 2 * c.mu) / (2 * c.mu) * trace_epsilon * alpha;
    TV H;
    if (delta_gamma <= 0) // case I: inside yield surface
    {
        H = epsilon + TV::Constant(cohesion);
    }
    else {
        H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat + TV::Constant(cohesion); // case III: projection
    }
    TV ret = H.array().exp();
    return ret;
}

// bool is for fracture
template <class T>
template <class TConst>
Matrix<T, TConst::dim, TConst::dim> DruckerPragerStvkHencky<T>::projectSigmaDerivative(TConst& c, const Vector<T, TConst::dim>& sigma)
{
    // const T eps = (T)1e-6;
    static const int dim = TConst::dim;
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, dim> TV;
    TV epsilon = sigma.array().abs().max(1e-4).log() - cohesion;
    T trace_epsilon = epsilon.sum();
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
    if (trace_epsilon >= (T)0) // case II: project to tip
    {
        TM ret = TM::Zero();
        return ret;
    }
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm + (dim * c.lambda + 2 * c.mu) / (2 * c.mu) * trace_epsilon * alpha;
    TV H;
    if (delta_gamma <= 0) // case I: inside yield surface
    {
        TM ret = TM::Identity();
        return ret;
    }
    else {
        TV w = sigma.array().inverse();
        T k = trace_epsilon;
        TV s = epsilon - k / dim * TV::Ones();
        TV s_hat = s / s.norm();
        T p = alpha * k * (dim * c.lambda + 2 * c.mu) / (2 * c.mu * s.norm());
        H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat + TV::Constant(cohesion); // case III: projection
        TV Z_hat = H.array().exp();
        TM ret = Z_hat.asDiagonal() * ((((T)1 + 2 * p) / dim * TV::Ones() - p / k * epsilon) * w.transpose() - p * (TM::Identity() - s_hat * s_hat.transpose()) * w.asDiagonal());
        return ret;
    }
}

template <class T>
template <class TConst>
void DruckerPragerStvkHencky<T>::projectSigmaAndDerivative(TConst& c, const Vector<T, TConst::dim>& sigma, Vector<T, TConst::dim>& projectedSigma, Matrix<T, TConst::dim, TConst::dim>& projectedSigmaDerivative)
{
    // const T eps = (T)1e-6;
    static const int dim = TConst::dim;
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, dim> TV;
    TV epsilon = sigma.array().abs().max(1e-4).log() - cohesion;
    T trace_epsilon = epsilon.sum();
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm + (dim * c.lambda + 2 * c.mu) / (2 * c.mu) * trace_epsilon * alpha;
    TV H;
    if (delta_gamma <= 0) // case I: inside yield surface
    {
        H = epsilon + TV::Constant(cohesion);
        projectedSigma = H.array().exp();
        projectedSigmaDerivative = TM::Identity();
    }
    else if (trace_epsilon > (T)0 || epsilon_hat_norm == 0) // case II: project to tip
    {
        projectedSigma = std::exp(cohesion) * TV::Ones();
        projectedSigmaDerivative = TM::Zero();
    }
    else {
        TV w = sigma.array().inverse();
        T k = trace_epsilon;
        TV s = epsilon - k / dim * TV::Ones();
        TV s_hat = s / s.norm();
        T p = alpha * k * (dim * c.lambda + 2 * c.mu) / (2 * c.mu * s.norm());
        H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat + TV::Constant(cohesion); // case III: projection
        projectedSigma = H.array().exp();
        projectedSigmaDerivative = projectedSigma.asDiagonal() * ((((T)1 + 2 * p) / dim * TV::Ones() - p / k * epsilon) * w.transpose() - p * (TM::Identity() - s_hat * s_hat.transpose()) * w.asDiagonal());
    }
}

template <class T>
template <class TConst>
void DruckerPragerStvkHencky<T>::computeSigmaPInverse(TConst& c, const Vector<T, TConst::dim>& sigma_e, Vector<T, TConst::dim>& sigma_p_inv)
{
    using TV = typename TConst::TV;
    TV sigma_proj = sigma_e;
    projectStrainDiagonal(c, sigma_proj);
    sigma_p_inv.array() = sigma_proj.array() / sigma_e.array();

    ZIRAN_WARN("Drucker Prager lambda step not fully implemented yet. ");
}

template <class T>
const char* DruckerPragerStvkHencky<T>::name()
{
    return "DruckerPragerStvkHencky";
}

template <class T, int dim>
VonMisesStvkHencky<T, dim>::VonMisesStvkHencky(const T yield_stress, const T fail_stress, const T xi)
    : yield_stress(yield_stress)
    , xi(xi)
    , fail_stress(fail_stress)
{
}

template <class T, int dim>
void VonMisesStvkHencky<T, dim>::setParameters(const T yield_stress_in, const T xi_in, const T fail_stress_in)
{
    yield_stress = yield_stress_in;
    xi = xi_in;
    if (fail_stress_in > -1)
        fail_stress = fail_stress_in;
}

// strain s is deformation F
//TODO no support for secondary cuts yet.
template <class T, int dim>
template <class TConst>
bool VonMisesStvkHencky<T, dim>::projectStrain(TConst& c, Matrix<T, TConst::dim, TConst::dim>& strain)
{
    static_assert(dim == TConst::dim, "Plasticity model has a different dimension as the Constitutive model!");
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, dim> TV;
    TM U, V;
    TV sigma;

    // TODO: this is inefficient because next time step updateState will do the svd again!
    singularValueDecomposition(strain, U, sigma, V);

    //TV epsilon = sigma.array().log();
    TV epsilon = sigma.array().max(1e-4).log(); //TODO: need the max part?
    T trace_epsilon = epsilon.sum();
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm - yield_stress / (2 * c.mu);
    if (delta_gamma <= 0) // case I
    {
        return false;
    }
    //hardening
    yield_stress -= xi * delta_gamma; //supposed to only increase yield_stress
    //yield_stress = std::max((T)0, yield_stress);

    TV H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat; // case II
    T tau_0 = 2 * c.mu * H(0) + c.lambda * H.sum();
    TV exp_H = H.array().exp();
    strain = U * exp_H.asDiagonal() * V.transpose();

    if (tau_0 >= fail_stress) {
        broken = true;
        crack_normal = V.col(0);
    }

    return false;
}

//TODO this is not up to date
template <class T, int dim>
template <class TConst>
void VonMisesStvkHencky<T, dim>::projectStrainDiagonal(TConst& c, Vector<T, TConst::dim>& sigma)
{
    static_assert(dim == TConst::dim, "Plasticity model has a different dimensiona s the Constitutive model!");
    typedef Vector<T, dim> TV;

    //TV epsilon = sigma.array().log();
    TV epsilon = sigma.array().max(1e-4).log();
    T trace_epsilon = epsilon.sum();
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm - yield_stress / (2 * c.mu);
    if (delta_gamma <= 0) // case I
    {
        return;
    }
    TV H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat; // case II
    sigma = H.array().exp();
}

template <class T, int dim>
template <class TConst>
Vector<T, TConst::dim> VonMisesStvkHencky<T, dim>::projectSigma(TConst& c, const Vector<T, TConst::dim>& sigma)
{
    static_assert(dim == TConst::dim, "Plasticity model has a different dimensiona s the Constitutive model!");
    typedef Vector<T, dim> TV;

    //TV epsilon = sigma.array().log();
    TV epsilon = sigma.array().max(1e-4).log();
    T trace_epsilon = epsilon.sum();
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm - yield_stress / (2 * c.mu);
    if (delta_gamma <= 0) // case I
    {
        return sigma;
    }
    TV H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat; // case II
    TV ret = H.array().exp();
    return ret;
}

template <class T, int dim>
template <class TConst>
Matrix<T, TConst::dim, TConst::dim> VonMisesStvkHencky<T, dim>::projectSigmaDerivative(TConst& c, const Vector<T, TConst::dim>& sigma)
{
    static_assert(dim == TConst::dim, "Plasticity model has a different dimensiona s the Constitutive model!");
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, dim> TV;

    //TV epsilon = sigma.array().log();
    TV epsilon = sigma.array().max(1e-4).log();
    T trace_epsilon = epsilon.sum();
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm - yield_stress / (2 * c.mu);
    if (delta_gamma <= 0) // case I
    {
        return TM::Identity();
    }
    TV w = sigma.array().inverse();
    T k = trace_epsilon;
    TV s = epsilon - k / dim * TV::Ones();
    TV s_hat = s / s.norm();
    TV H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat; // case II
    TV Z_hat = H.array().exp();
    TM ret = Z_hat.asDiagonal() * (TM::Identity() * w.asDiagonal() - ((T)1 - yield_stress / (2 * c.mu * s.norm())) * (TM::Identity() * w.asDiagonal() - TV::Ones() * w.transpose() / (T)dim) - yield_stress / (2 * c.mu * s.norm()) * s_hat * s_hat.transpose() * w.asDiagonal());
    return ret;
}

template <class T, int dim>
template <class TConst>
void VonMisesStvkHencky<T, dim>::projectSigmaAndDerivative(TConst& c, const Vector<T, TConst::dim>& sigma, Vector<T, TConst::dim>& projectedSigma, Matrix<T, TConst::dim, TConst::dim>& projectedSigmaDerivative)
{
    static_assert(dim == TConst::dim, "Plasticity model has a different dimensiona s the Constitutive model!");
    typedef Matrix<T, dim, dim> TM;
    typedef Vector<T, dim> TV;

    //TV epsilon = sigma.array().log();
    TV epsilon = sigma.array().max(1e-4).log();
    T trace_epsilon = epsilon.sum();
    TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
    T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
    T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
    T delta_gamma = epsilon_hat_norm - yield_stress / (2 * c.mu);
    if (delta_gamma <= 0) // case I
    {
        projectedSigma = sigma;
        projectedSigmaDerivative = TM::Identity();
        return;
    }
    TV w = sigma.array().inverse();
    T k = trace_epsilon;
    TV s = epsilon - k / dim * TV::Ones();
    TV s_hat = s / s.norm();
    TV H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat; // case II
    projectedSigma = H.array().exp();
    projectedSigmaDerivative = projectedSigma.asDiagonal() * (TM::Identity() * w.asDiagonal() - ((T)1 - yield_stress / (2 * c.mu * s.norm())) * (TM::Identity() * w.asDiagonal() - TV::Ones() * w.transpose() / (T)dim) - yield_stress / (2 * c.mu * s.norm()) * s_hat * s_hat.transpose() * w.asDiagonal());
}

template <class T, int dim>
template <class TConst>
void VonMisesStvkHencky<T, dim>::computeSigmaPInverse(TConst& c, const Vector<T, TConst::dim>& sigma_e, Vector<T, TConst::dim>& sigma_p_inv)
{
    using TV = typename TConst::TV;
    TV eps = sigma_e.array().log();
    TV tau = c.lambda * eps.sum() + 2 * c.mu * eps.array();
    TV tau_trfree;
    tau_trfree.array() = tau.array() - ((T)1 / TConst::dim) * tau.sum();

    T tau_trfree_norm = tau_trfree.norm();
    if (tau_trfree_norm > yield_stress) {
        // Lambda = 1 / (2*dt*mu*(1+ k/(tau_trefree.norm()-k)))*tau_trfree
        sigma_p_inv = tau_trfree / tau_trfree_norm;
        sigma_p_inv *= (yield_stress - tau_trfree_norm) / (2.0 * c.mu);
        sigma_p_inv = sigma_p_inv.array().exp();
    }
    else
        sigma_p_inv = TV::Ones();
}

template <class T, int dim>
const char* VonMisesStvkHencky<T, dim>::name()
{
    return "VonMisesStvkHencky";
}

template class NonAssociativeCamClay<double>;
template class NonAssociativeCamClay<float>;

template class DruckerPragerStvkHencky<double>;
template class DruckerPragerStvkHencky<float>;

template class VonMisesStvkHencky<double, 2>;
template class VonMisesStvkHencky<float, 2>;
template class VonMisesStvkHencky<double, 3>;
template class VonMisesStvkHencky<float, 3>;

template class PlasticityApplier<NeoHookeanBorden<double, 2>, NonAssociativeCamClay<double>, Eigen::Matrix<double, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<NeoHookeanBorden<double, 3>, NonAssociativeCamClay<double>, Eigen::Matrix<double, 3, 3, 0, 3, 3>>;
template class PlasticityApplier<NeoHookeanBorden<float, 2>, NonAssociativeCamClay<float>, Eigen::Matrix<float, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<NeoHookeanBorden<float, 3>, NonAssociativeCamClay<float>, Eigen::Matrix<float, 3, 3, 0, 3, 3>>;

template class PlasticityApplier<StvkWithHencky<double, 2>, DruckerPragerStvkHencky<double>, Eigen::Matrix<double, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<StvkWithHencky<double, 2>, VonMisesStvkHencky<double, 2>, Eigen::Matrix<double, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<StvkWithHencky<double, 3>, DruckerPragerStvkHencky<double>, Eigen::Matrix<double, 3, 3, 0, 3, 3>>;
template class PlasticityApplier<StvkWithHencky<double, 3>, VonMisesStvkHencky<double, 3>, Eigen::Matrix<double, 3, 3, 0, 3, 3>>;
template class PlasticityApplier<StvkWithHencky<float, 2>, DruckerPragerStvkHencky<float>, Eigen::Matrix<float, 2, 2, 0, 2, 2>>;

template class PlasticityApplier<StvkWithHencky<float, 2>, VonMisesStvkHencky<float, 2>, Eigen::Matrix<float, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<StvkWithHencky<float, 3>, DruckerPragerStvkHencky<float>, Eigen::Matrix<float, 3, 3, 0, 3, 3>>;
template class PlasticityApplier<StvkWithHencky<float, 3>, VonMisesStvkHencky<float, 3>, Eigen::Matrix<float, 3, 3, 0, 3, 3>>;
template class PlasticityApplier<StvkWithHenckyIsotropic<double, 2>, DruckerPragerStvkHencky<double>, Eigen::Matrix<double, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<StvkWithHenckyIsotropic<double, 2>, VonMisesStvkHencky<double, 2>, Eigen::Matrix<double, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<StvkWithHenckyIsotropic<double, 3>, DruckerPragerStvkHencky<double>, Eigen::Matrix<double, 3, 3, 0, 3, 3>>;
template class PlasticityApplier<StvkWithHenckyIsotropic<double, 3>, VonMisesStvkHencky<double, 3>, Eigen::Matrix<double, 3, 3, 0, 3, 3>>;
template class PlasticityApplier<StvkWithHenckyIsotropic<float, 2>, DruckerPragerStvkHencky<float>, Eigen::Matrix<float, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<StvkWithHenckyIsotropic<float, 2>, VonMisesStvkHencky<float, 2>, Eigen::Matrix<float, 2, 2, 0, 2, 2>>;
template class PlasticityApplier<StvkWithHenckyIsotropic<float, 3>, DruckerPragerStvkHencky<float>, Eigen::Matrix<float, 3, 3, 0, 3, 3>>;
template class PlasticityApplier<StvkWithHenckyIsotropic<float, 3>, VonMisesStvkHencky<float, 3>, Eigen::Matrix<float, 3, 3, 0, 3, 3>>;

template void VonMisesStvkHencky<double, 2>::projectStrainDiagonal<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, Eigen::Matrix<double, StvkWithHencky<double, 2>::dim, 1, 0, StvkWithHencky<double, 2>::dim, 1>&);
template void VonMisesStvkHencky<double, 3>::projectStrainDiagonal<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, Eigen::Matrix<double, StvkWithHencky<double, 3>::dim, 1, 0, StvkWithHencky<double, 3>::dim, 1>&);
template void VonMisesStvkHencky<float, 2>::projectStrainDiagonal<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, Eigen::Matrix<float, StvkWithHencky<float, 2>::dim, 1, 0, StvkWithHencky<float, 2>::dim, 1>&);
template void VonMisesStvkHencky<float, 3>::projectStrainDiagonal<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, Eigen::Matrix<float, StvkWithHencky<float, 3>::dim, 1, 0, StvkWithHencky<float, 3>::dim, 1>&);
template void VonMisesStvkHencky<double, 2>::computeSigmaPInverse<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, Eigen::Matrix<double, StvkWithHencky<double, 2>::dim, 1, 0, StvkWithHencky<double, 2>::dim, 1> const&, Eigen::Matrix<double, StvkWithHencky<double, 2>::dim, 1, 0, StvkWithHencky<double, 2>::dim, 1>&);
template void VonMisesStvkHencky<double, 3>::computeSigmaPInverse<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, Eigen::Matrix<double, StvkWithHencky<double, 3>::dim, 1, 0, StvkWithHencky<double, 3>::dim, 1> const&, Eigen::Matrix<double, StvkWithHencky<double, 3>::dim, 1, 0, StvkWithHencky<double, 3>::dim, 1>&);
template void VonMisesStvkHencky<float, 2>::computeSigmaPInverse<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, Eigen::Matrix<float, StvkWithHencky<float, 2>::dim, 1, 0, StvkWithHencky<float, 2>::dim, 1> const&, Eigen::Matrix<float, StvkWithHencky<float, 2>::dim, 1, 0, StvkWithHencky<float, 2>::dim, 1>&);
template void VonMisesStvkHencky<float, 3>::computeSigmaPInverse<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, Eigen::Matrix<float, StvkWithHencky<float, 3>::dim, 1, 0, StvkWithHencky<float, 3>::dim, 1> const&, Eigen::Matrix<float, StvkWithHencky<float, 3>::dim, 1, 0, StvkWithHencky<float, 3>::dim, 1>&);

template Eigen::Matrix<double, 3, 1, 0, 3, 1> DruckerPragerStvkHencky<double>::projectSigma<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> DruckerPragerStvkHencky<float>::projectSigma<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<double, 2, 1, 0, 2, 1> DruckerPragerStvkHencky<double>::projectSigma<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<float, 2, 1, 0, 2, 1> DruckerPragerStvkHencky<float>::projectSigma<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<double, 3, 3, 0, 3, 3> DruckerPragerStvkHencky<double>::projectSigmaDerivative<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<float, 3, 3, 0, 3, 3> DruckerPragerStvkHencky<float>::projectSigmaDerivative<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<double, 2, 2, 0, 2, 2> DruckerPragerStvkHencky<double>::projectSigmaDerivative<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<float, 2, 2, 0, 2, 2> DruckerPragerStvkHencky<float>::projectSigmaDerivative<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&);
template void DruckerPragerStvkHencky<double>::projectSigmaAndDerivative<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&, Eigen::Matrix<double, 3, 1, 0, 3, 1>&, Eigen::Matrix<double, 3, 3, 0, 3, 3>&);
template void DruckerPragerStvkHencky<float>::projectSigmaAndDerivative<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&, Eigen::Matrix<float, 3, 1, 0, 3, 1>&, Eigen::Matrix<float, 3, 3, 0, 3, 3>&);
template void DruckerPragerStvkHencky<double>::projectSigmaAndDerivative<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&, Eigen::Matrix<double, 2, 1, 0, 2, 1>&, Eigen::Matrix<double, 2, 2, 0, 2, 2>&);
template void DruckerPragerStvkHencky<float>::projectSigmaAndDerivative<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&, Eigen::Matrix<float, 2, 1, 0, 2, 1>&, Eigen::Matrix<float, 2, 2, 0, 2, 2>&);

template Eigen::Matrix<double, 3, 1, 0, 3, 1> VonMisesStvkHencky<double, 3>::projectSigma<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> VonMisesStvkHencky<float, 3>::projectSigma<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<double, 2, 1, 0, 2, 1> VonMisesStvkHencky<double, 2>::projectSigma<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<float, 2, 1, 0, 2, 1> VonMisesStvkHencky<float, 2>::projectSigma<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<double, 3, 3, 0, 3, 3> VonMisesStvkHencky<double, 3>::projectSigmaDerivative<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<float, 3, 3, 0, 3, 3> VonMisesStvkHencky<float, 3>::projectSigmaDerivative<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<double, 2, 2, 0, 2, 2> VonMisesStvkHencky<double, 2>::projectSigmaDerivative<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<float, 2, 2, 0, 2, 2> VonMisesStvkHencky<float, 2>::projectSigmaDerivative<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&);
template void VonMisesStvkHencky<double, 3>::projectSigmaAndDerivative<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&, Eigen::Matrix<double, 3, 1, 0, 3, 1>&, Eigen::Matrix<double, 3, 3, 0, 3, 3>&);
template void VonMisesStvkHencky<float, 3>::projectSigmaAndDerivative<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&, Eigen::Matrix<float, 3, 1, 0, 3, 1>&, Eigen::Matrix<float, 3, 3, 0, 3, 3>&);
template void VonMisesStvkHencky<double, 2>::projectSigmaAndDerivative<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&, Eigen::Matrix<double, 2, 1, 0, 2, 1>&, Eigen::Matrix<double, 2, 2, 0, 2, 2>&);
template void VonMisesStvkHencky<float, 2>::projectSigmaAndDerivative<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&, Eigen::Matrix<float, 2, 1, 0, 2, 1>&, Eigen::Matrix<float, 2, 2, 0, 2, 2>&);

template Eigen::Matrix<double, 3, 1, 0, 3, 1> DummyPlasticity<double>::projectSigma<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> DummyPlasticity<float>::projectSigma<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<double, 2, 1, 0, 2, 1> DummyPlasticity<double>::projectSigma<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<float, 2, 1, 0, 2, 1> DummyPlasticity<float>::projectSigma<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<double, 3, 3, 0, 3, 3> DummyPlasticity<double>::projectSigmaDerivative<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<float, 3, 3, 0, 3, 3> DummyPlasticity<float>::projectSigmaDerivative<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&);
template Eigen::Matrix<double, 2, 2, 0, 2, 2> DummyPlasticity<double>::projectSigmaDerivative<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&);
template Eigen::Matrix<float, 2, 2, 0, 2, 2> DummyPlasticity<float>::projectSigmaDerivative<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&);
template void DummyPlasticity<double>::projectSigmaAndDerivative<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3>&, const Eigen::Matrix<double, 3, 1, 0, 3, 1>&, Eigen::Matrix<double, 3, 1, 0, 3, 1>&, Eigen::Matrix<double, 3, 3, 0, 3, 3>&);
template void DummyPlasticity<float>::projectSigmaAndDerivative<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3>&, const Eigen::Matrix<float, 3, 1, 0, 3, 1>&, Eigen::Matrix<float, 3, 1, 0, 3, 1>&, Eigen::Matrix<float, 3, 3, 0, 3, 3>&);
template void DummyPlasticity<double>::projectSigmaAndDerivative<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2>&, const Eigen::Matrix<double, 2, 1, 0, 2, 1>&, Eigen::Matrix<double, 2, 1, 0, 2, 1>&, Eigen::Matrix<double, 2, 2, 0, 2, 2>&);
template void DummyPlasticity<float>::projectSigmaAndDerivative<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2>&, const Eigen::Matrix<float, 2, 1, 0, 2, 1>&, Eigen::Matrix<float, 2, 1, 0, 2, 1>&, Eigen::Matrix<float, 2, 2, 0, 2, 2>&);

template bool NonAssociativeCamClay<double>::projectStrain<NeoHookeanBorden<double, 2>>(NeoHookeanBorden<double, 2>&, Eigen::Matrix<double, NeoHookeanBorden<double, 2>::dim, NeoHookeanBorden<double, 2>::dim, 0, NeoHookeanBorden<double, 2>::dim, NeoHookeanBorden<double, 2>::dim>&);
template bool NonAssociativeCamClay<float>::projectStrain<NeoHookeanBorden<float, 2>>(NeoHookeanBorden<float, 2>&, Eigen::Matrix<float, NeoHookeanBorden<float, 2>::dim, NeoHookeanBorden<float, 2>::dim, 0, NeoHookeanBorden<float, 2>::dim, NeoHookeanBorden<float, 2>::dim>&);

} // namespace ZIRAN
