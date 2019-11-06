#include "../Force/FBasedMpmForceHelper.cpp"
namespace ZIRAN {
template class FBasedMpmForceHelper<NeoHookeanBorden<double, 3>>;
template class FBasedMpmForceHelper<StvkWithHencky<double, 3>>;
template class FBasedMpmForceHelper<StvkWithHenckyIsotropic<double, 3>>;
} // namespace ZIRAN
