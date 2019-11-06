#include <MPM/MpmParticleHandleBase.cpp>
namespace ZIRAN {
template class MpmParticleHandleBase<double, 3>;
extern template class MpmForceBase<double, 3>;
template void MpmParticleHandleBase<double, 3>::addFBasedMpmForceWithPhaseField<NeoHookeanBorden<double, 3>>(const double&, const double&, NeoHookeanBorden<double, 3> const&, bool, const double);
template void MpmParticleHandleBase<double, 3>::addFBasedMpmForce<NeoHookeanBorden<double, 3>>(NeoHookeanBorden<double, 3> const&);
template void MpmParticleHandleBase<double, 3>::addFBasedMpmForce<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3> const&);
template void MpmParticleHandleBase<double, 3>::addFBasedMpmForce<StvkWithHenckyIsotropic<double, 3>>(StvkWithHenckyIsotropic<double, 3> const&);
template void MpmParticleHandleBase<double, 3>::addFElasticNonequilibratedBasedMpmForce<StvkWithHencky<double, 3>>(StvkWithHencky<double, 3> const&, double, double);
template void MpmParticleHandleBase<double, 3>::addPlasticity<NeoHookeanBorden<double, 3>, NonAssociativeCamClay<double>>(NeoHookeanBorden<double, 3> const&, NonAssociativeCamClay<double> const&, std::string);
template void MpmParticleHandleBase<double, 3>::addPlasticity<StvkWithHencky<double, 3>, DruckerPragerStvkHencky<double>>(StvkWithHencky<double, 3> const&, DruckerPragerStvkHencky<double> const&, std::string);
template void MpmParticleHandleBase<double, 3>::addPlasticity<StvkWithHencky<double, 3>, VonMisesStvkHencky<double, 3>>(StvkWithHencky<double, 3> const&, VonMisesStvkHencky<double, 3> const&, std::string);
template void MpmParticleHandleBase<double, 3>::addPlasticity<StvkWithHenckyIsotropic<double, 3>, VonMisesStvkHencky<double, 3>>(StvkWithHenckyIsotropic<double, 3> const&, VonMisesStvkHencky<double, 3> const&, std::string);
} // namespace ZIRAN
