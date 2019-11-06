#include <MPM/MpmParticleHandleBase.cpp>
namespace ZIRAN {
template class MpmParticleHandleBase<float, 3>;
extern template class MpmForceBase<float, 3>;
template void MpmParticleHandleBase<float, 3>::addFBasedMpmForceWithPhaseField<NeoHookeanBorden<float, 3>>(const float&, const float&, NeoHookeanBorden<float, 3> const&, bool, const float);
template void MpmParticleHandleBase<float, 3>::addFBasedMpmForce<NeoHookeanBorden<float, 3>>(NeoHookeanBorden<float, 3> const&);
template void MpmParticleHandleBase<float, 3>::addFBasedMpmForce<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3> const&);
template void MpmParticleHandleBase<float, 3>::addFBasedMpmForce<StvkWithHenckyIsotropic<float, 3>>(StvkWithHenckyIsotropic<float, 3> const&);
template void MpmParticleHandleBase<float, 3>::addFElasticNonequilibratedBasedMpmForce<StvkWithHencky<float, 3>>(StvkWithHencky<float, 3> const&, float, float);
template void MpmParticleHandleBase<float, 3>::addPlasticity<NeoHookeanBorden<float, 3>, NonAssociativeCamClay<float>>(NeoHookeanBorden<float, 3> const&, NonAssociativeCamClay<float> const&, std::string);
template void MpmParticleHandleBase<float, 3>::addPlasticity<StvkWithHencky<float, 3>, DruckerPragerStvkHencky<float>>(StvkWithHencky<float, 3> const&, DruckerPragerStvkHencky<float> const&, std::string);
template void MpmParticleHandleBase<float, 3>::addPlasticity<StvkWithHencky<float, 3>, VonMisesStvkHencky<float, 3>>(StvkWithHencky<float, 3> const&, VonMisesStvkHencky<float, 3> const&, std::string);
template void MpmParticleHandleBase<float, 3>::addPlasticity<StvkWithHenckyIsotropic<float, 3>, VonMisesStvkHencky<float, 3>>(StvkWithHenckyIsotropic<float, 3> const&, VonMisesStvkHencky<float, 3> const&, std::string);
} // namespace ZIRAN
