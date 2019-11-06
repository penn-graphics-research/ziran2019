#ifndef ELEMENT_MANAGER_FOR_H
#define ELEMENT_MANAGER_FOR_H
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/Math/Geometry/SimplexMesh.h>
#include <Ziran/Math/Geometry/SimplexElements.h>
namespace ZIRAN {
namespace INTERNAL {
//*****************************************************************************
// ElementManagerFor is a type that encodes SimplexElements
//*****************************************************************************

template <class T, class TMesh, int dim>
struct ElementManagerForHelper;
template <class T, int manifold_dim, int dim>
struct ElementManagerForHelper<T, SimplexMesh<manifold_dim>, dim> {
    using type = SimplexElements<T, manifold_dim, dim>;
};

} // namespace INTERNAL
template <class T, class TMesh, int dim>
using ElementManagerFor = typename INTERNAL::ElementManagerForHelper<T, TMesh, dim>::type;
} // namespace ZIRAN
#endif
