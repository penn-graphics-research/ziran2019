#ifndef SOUND_SPEED_CFL_H
#define SOUND_SPEED_CFL_H
#include <cmath>
namespace ZIRAN {

//#########################################################################
// Function: evaluateSoundSpeedLinearElasticityAnalysis
//#########################################################################
template <class T>
inline T evaluateSoundSpeedLinearElasticityAnalysis(const T E, const T nu, const T rho)
{
    return std::sqrt(E * (1 - nu) / ((1 + nu) * (1 - 2 * nu) * rho));
}

//#########################################################################
// Function: evaluateTimestepLinearElasticityAnalysis
//#########################################################################
template <class T>
inline T evaluateTimestepLinearElasticityAnalysis(const T E, const T nu, const T rho, const T dx, const T cfl)
{
    return cfl * dx / evaluateSoundSpeedLinearElasticityAnalysis(E, nu, rho);
}

} // namespace ZIRAN
#endif
