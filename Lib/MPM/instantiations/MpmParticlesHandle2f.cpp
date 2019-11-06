#include <MPM/MpmParticleHandleBase.cpp>
namespace ZIRAN {
template class MpmParticleHandleBase<float, 2>;
extern template class MpmForceBase<float, 2>;
template void MpmParticleHandleBase<float, 2>::addFBasedMpmForceWithPhaseField<NeoHookeanBorden<float, 2>>(const float&, const float&, NeoHookeanBorden<float, 2> const&, bool, const float);
template void MpmParticleHandleBase<float, 2>::addFBasedMpmForce<NeoHookeanBorden<float, 2>>(NeoHookeanBorden<float, 2> const&);
template void MpmParticleHandleBase<float, 2>::addFBasedMpmForce<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2> const&);
template void MpmParticleHandleBase<float, 2>::addFBasedMpmForce<StvkWithHenckyIsotropic<float, 2>>(StvkWithHenckyIsotropic<float, 2> const&);
template void MpmParticleHandleBase<float, 2>::addFElasticNonequilibratedBasedMpmForce<StvkWithHencky<float, 2>>(StvkWithHencky<float, 2> const&, float, float);
template void MpmParticleHandleBase<float, 2>::addPlasticity<NeoHookeanBorden<float, 2>, NonAssociativeCamClay<float>>(NeoHookeanBorden<float, 2> const&, NonAssociativeCamClay<float> const&, std::string);
template void MpmParticleHandleBase<float, 2>::addPlasticity<StvkWithHencky<float, 2>, DruckerPragerStvkHencky<float>>(StvkWithHencky<float, 2> const&, DruckerPragerStvkHencky<float> const&, std::string);
template void MpmParticleHandleBase<float, 2>::addPlasticity<StvkWithHencky<float, 2>, VonMisesStvkHencky<float, 2>>(StvkWithHencky<float, 2> const&, VonMisesStvkHencky<float, 2> const&, std::string);
template void MpmParticleHandleBase<float, 2>::addPlasticity<StvkWithHenckyIsotropic<float, 2>, VonMisesStvkHencky<float, 2>>(StvkWithHenckyIsotropic<float, 2> const&, VonMisesStvkHencky<float, 2> const&, std::string);
} // namespace ZIRAN
