#include "../Force/FBasedMpmForceHelper.cpp"
namespace ZIRAN {
template class FBasedMpmForceHelper<NeoHookeanBorden<float, 2>>;
template class FBasedMpmForceHelper<StvkWithHencky<float, 2>>;
template class FBasedMpmForceHelper<StvkWithHenckyIsotropic<float, 2>>;
} // namespace ZIRAN
