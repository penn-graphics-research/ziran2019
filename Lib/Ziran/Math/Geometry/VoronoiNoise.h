#ifndef VORONOI_NOISE_H
#define VORONOI_NOISE_H
#include <Ziran/CS/Util/Forward.h>
namespace ZIRAN {

// Generate voronoi noise
float voronoiDistance(const Vector<float, 3>& x, const Vector<int, 3>& seed = Vector<int, 3>::Zero());

// Generate voronoi noise
float voronoiDistance(const Vector<float, 2>& x, const Vector<int, 2>& seed = Vector<int, 2>::Zero());
} // namespace ZIRAN
#endif
