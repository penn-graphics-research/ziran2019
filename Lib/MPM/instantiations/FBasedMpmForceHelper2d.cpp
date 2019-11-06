#include "../Force/FBasedMpmForceHelper.cpp"
namespace ZIRAN {
template class FBasedMpmForceHelper<NeoHookeanBorden<double, 2>>;
template class FBasedMpmForceHelper<StvkWithHencky<double, 2>>;
template class FBasedMpmForceHelper<StvkWithHenckyIsotropic<double, 2>>;
} // namespace ZIRAN
