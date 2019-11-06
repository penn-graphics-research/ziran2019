#ifndef MESH_CONSTRUCTION_H
#define MESH_CONSTRUCTION_H
#include <Ziran/CS/Util/Logging.h>
#include <Ziran/Math/Geometry/Grid.h>
#include <Ziran/Math/Geometry/SimplexMesh.h>
#include <Ziran/CS/Util/Forward.h>

/**
This is a collectiong of helpers for building various types of Lagrangian meshes
**/

namespace ZIRAN {
namespace MESHING {

template <class T, int dim>
void constructUnitSimplex(SimplexMesh<dim>& mesh, StdVector<Vector<T, dim>>& X);

template <class T, int dim>
void appendStrandToSegmesh(const Vector<T, dim>& start, const Vector<T, dim>& end,
    const int N_cells,
    SimplexMesh<1>& mesh, StdVector<Vector<T, dim>>& X);

template <class T, int dim>
void appendCircleToSegmesh(const Vector<T, dim>& center, const T radius,
    const int N_cells,
    SimplexMesh<1>& mesh, StdVector<Vector<T, dim>>& X);

template <class T, int dim>
void constructMattressMesh(const Grid<T, 2>& grid, SimplexMesh<2>& mesh, StdVector<Vector<T, dim>>& X);

template <class T, int dim>
void constructMattressMeshWithSurfaceMarker(const Grid<T, 2>& grid, SimplexMesh<2>& mesh, StdVector<Vector<T, dim>>& X, StdVector<Vector<T, dim>>& normals);

template <class T, int dim>
void constructHairMeshWithSurfaceMarker(const Grid<T, 2>& grid, SimplexMesh<2>& mesh, StdVector<Vector<T, dim>>& X, StdVector<Vector<T, dim>>& normals);

template <class T>
void constructMattressMesh(const Grid<T, 3>& grid, SimplexMesh<3>& mesh, StdVector<Vector<T, 3>>& X);
}
} // namespace ZIRAN::MESHING
#endif
