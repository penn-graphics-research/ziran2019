#include "../Force/FBasedMpmForceHelper.cpp"
namespace ZIRAN {
template class FBasedMpmForceHelper<NeoHookeanBorden<float, 3>>;
template class FBasedMpmForceHelper<StvkWithHencky<float, 3>>;
template class FBasedMpmForceHelper<StvkWithHenckyIsotropic<float, 3>>;
} // namespace ZIRAN
