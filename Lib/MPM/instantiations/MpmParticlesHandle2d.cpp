#include <MPM/MpmParticleHandleBase.cpp>
namespace ZIRAN {
template class MpmParticleHandleBase<double, 2>;
extern template class MpmForceBase<double, 2>;
template void MpmParticleHandleBase<double, 2>::addFBasedMpmForceWithPhaseField<NeoHookeanBorden<double, 2>>(const double&, const double&, NeoHookeanBorden<double, 2> const&, bool, const double);
template void MpmParticleHandleBase<double, 2>::addFBasedMpmForce<NeoHookeanBorden<double, 2>>(NeoHookeanBorden<double, 2> const&);
template void MpmParticleHandleBase<double, 2>::addFBasedMpmForce<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2> const&);
template void MpmParticleHandleBase<double, 2>::addFBasedMpmForce<StvkWithHenckyIsotropic<double, 2>>(StvkWithHenckyIsotropic<double, 2> const&);
template void MpmParticleHandleBase<double, 2>::addFElasticNonequilibratedBasedMpmForce<StvkWithHencky<double, 2>>(StvkWithHencky<double, 2> const&, double, double);
template void MpmParticleHandleBase<double, 2>::addPlasticity<NeoHookeanBorden<double, 2>, NonAssociativeCamClay<double>>(NeoHookeanBorden<double, 2> const&, NonAssociativeCamClay<double> const&, std::string);
template void MpmParticleHandleBase<double, 2>::addPlasticity<StvkWithHencky<double, 2>, DruckerPragerStvkHencky<double>>(StvkWithHencky<double, 2> const&, DruckerPragerStvkHencky<double> const&, std::string);
template void MpmParticleHandleBase<double, 2>::addPlasticity<StvkWithHencky<double, 2>, VonMisesStvkHencky<double, 2>>(StvkWithHencky<double, 2> const&, VonMisesStvkHencky<double, 2> const&, std::string);
template void MpmParticleHandleBase<double, 2>::addPlasticity<StvkWithHenckyIsotropic<double, 2>, VonMisesStvkHencky<double, 2>>(StvkWithHenckyIsotropic<double, 2> const&, VonMisesStvkHencky<double, 2> const&, std::string);
} // namespace ZIRAN
