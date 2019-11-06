#include "MeshConstruction.h"
#include <Ziran/CS/Util/Logging.h>
#include <Ziran/Math/Geometry/Grid.h>
#include <Ziran/Math/Geometry/SimplexMesh.h>

namespace ZIRAN {
namespace MESHING {

template <class T, int dim>
void constructUnitSimplex(SimplexMesh<dim>& mesh, StdVector<Vector<T, dim>>& X)
{
    using TV = Vector<T, dim>;
    using IV = Eigen::Matrix<int, dim + 1, 1>;

    X.reserve(dim + 1);
    TV z = TV::Zero();
    X.push_back(z);
    for (int i = 0; i < dim; i++)
        X.push_back(TV::Unit(i));

    IV indices;
    for (int i = 0; i < dim + 1; i++)
        indices(i) = i;
    mesh.indices.push_back(indices);
}

template <class T, int dim>
void appendStrandToSegmesh(const Vector<T, dim>& start, const Vector<T, dim>& end,
    const int N_cells,
    SimplexMesh<1>& mesh, StdVector<Vector<T, dim>>& X)
{
    using IVE = Vector<int, 2>;
    using TV = Vector<T, dim>;

    TV dx = (end - start) / N_cells;

    //first create the mesh
    ZIRAN_INFO("Creating mesh ...");
    size_t starting_index = 0;
    for (int i = 0; i < N_cells; i++) {
        IVE indices;
        indices(0) = starting_index + i;
        indices(1) = indices(0) + 1;
        mesh.indices.push_back(indices);
    }

    // now create the particles
    ZIRAN_INFO("constructing particles ...");
    X.reserve(N_cells + 1);
    for (int i = 0; i < N_cells + 1; i++) {
        TV x = start + i * dx;
        X.push_back(x);
    }
}

template <class T, int dim>
void appendCircleToSegmesh(const Vector<T, dim>& center, const T radius,
    const int N_cells,
    SimplexMesh<1>& mesh, StdVector<Vector<T, dim>>& X)
{
    using IVE = Vector<int, 2>;
    using TV = Vector<T, dim>;

    T dtheta = 2 * M_PI / N_cells;

    //first create the mesh
    ZIRAN_INFO("Creating mesh ...");
    size_t starting_index = 0;
    for (int i = 0; i < N_cells; i++) {
        IVE indices;
        indices(0) = starting_index + i;
        indices(1) = indices(0) + 1;
        mesh.indices.push_back(indices);
    }
    IVE indices;
    indices(0) = starting_index;
    indices(1) = starting_index + N_cells;
    mesh.indices.push_back(indices);

    // now create the particles
    ZIRAN_INFO("constructing particles ...");
    X.reserve(N_cells);
    for (int i = 0; i < N_cells + 1; i++) {
        T theta = i * dtheta;
        TV radial;
        radial << std::cos(theta) * radius, std::sin(theta) * radius;
        TV x = center + radial;
        X.push_back(x);
    }
}

template <class T, int dim>
void constructMattressMesh(const Grid<T, 2>& grid, SimplexMesh<2>& mesh, StdVector<Vector<T, dim>>& X)
{
    using IVE = Vector<int, 3>;
    using IV = Vector<int, 2>;
    using TV = Vector<T, dim>;

    //first create the mesh
    IV index;
    IVE indices;
    ZIRAN_INFO("Creating mesh ...");

    for (int i = 0; i < 2; i++)
        ZIRAN_INFO("Grid number cells ", i, ": ", grid.numberCellsPerDimension(i));

    for (int i = 0; i < grid.numberCellsPerDimension(0); i++) {
        for (int j = 0; j < grid.numberCellsPerDimension(1); j++) {
            index << i, j;
            indices(0) = grid.nodeIndex(index);
            index << i + 1, j;
            indices(1) = grid.nodeIndex(index);
            index << i + 1, j + 1;
            indices(2) = grid.nodeIndex(index);
            mesh.indices.push_back(indices);

            index << i, j;
            indices(0) = grid.nodeIndex(index);
            index << i + 1, j + 1;
            indices(1) = grid.nodeIndex(index);
            index << i, j + 1;
            indices(2) = grid.nodeIndex(index);
            mesh.indices.push_back(indices);
        }
    }

    // now create the particles
    ZIRAN_INFO("constructing particles ...");
    X.reserve(grid.numberNodes());
    TV x = TV::Zero();
    for (int i = 0; i < grid.numberNodesPerDimension(0); i++) {
        for (int j = 0; j < grid.numberNodesPerDimension(1); j++) {
            index << i, j;
            x.template head<2>() = grid.node(index);
            X.push_back(x);
        }
    }
}

template <class T, int dim>
void constructMattressMeshWithSurfaceMarker(const Grid<T, 2>& grid, SimplexMesh<2>& mesh, StdVector<Vector<T, dim>>& X, StdVector<Vector<T, dim>>& normals)
{
    using IVE = Vector<int, 3>;
    using IV = Vector<int, 2>;
    using TV = Vector<T, dim>;

    //first create the mesh
    IV index;
    IVE indices;
    ZIRAN_INFO("Creating mesh ...");

    for (int i = 0; i < 2; i++)
        ZIRAN_INFO("Grid number cells ", i, ": ", grid.numberCellsPerDimension(i));

    for (int i = 0; i < grid.numberCellsPerDimension(0); i++) {
        for (int j = 0; j < grid.numberCellsPerDimension(1); j++) {
            index << i, j;
            indices(0) = grid.nodeIndex(index);
            index << i + 1, j;
            indices(1) = grid.nodeIndex(index);
            index << i + 1, j + 1;
            indices(2) = grid.nodeIndex(index);
            mesh.indices.push_back(indices);

            index << i, j;
            indices(0) = grid.nodeIndex(index);
            index << i + 1, j + 1;
            indices(1) = grid.nodeIndex(index);
            index << i, j + 1;
            indices(2) = grid.nodeIndex(index);
            mesh.indices.push_back(indices);
        }
    }

    // now create the particles
    ZIRAN_INFO("constructing particles ...");
    X.reserve(grid.numberNodes());
    TV x = TV::Zero();
    for (int i = 0; i < grid.numberNodesPerDimension(0); i++) {
        for (int j = 0; j < grid.numberNodesPerDimension(1); j++) {
            index << i, j;
            x.template head<2>() = grid.node(index);
            X.push_back(x);
            if (i == 0 || i == grid.numberNodesPerDimension(0) - 1 || j == 0 || j == grid.numberNodesPerDimension(1) - 1) {
                normals.push_back(Vector<T, dim>::Unit(0));
            }
            else {
                normals.push_back(Vector<T, dim>::Zero());
            }
        }
    }
}

template <class T, int dim>
void constructHairMeshWithSurfaceMarker(const Grid<T, 2>& grid, SimplexMesh<2>& mesh, StdVector<Vector<T, dim>>& X, StdVector<Vector<T, dim>>& normals)
{
    using IVE = Vector<int, 3>;
    using IV = Vector<int, 2>;
    using TV = Vector<T, dim>;

    //first create the mesh
    IV index;
    IVE indices;
    ZIRAN_INFO("Creating mesh ...");

    for (int i = 0; i < 2; i++)
        ZIRAN_INFO("Grid number cells ", i, ": ", grid.numberCellsPerDimension(i));

    for (int i = 0; i < grid.numberCellsPerDimension(0); i++) {
        for (int j = 0; j < grid.numberCellsPerDimension(1); j++) {
            index << i, j;
            indices(0) = grid.nodeIndex(index);
            index << i + 1, j;
            indices(1) = grid.nodeIndex(index);
            index << i + 1, j + 1;
            indices(2) = grid.nodeIndex(index);
            mesh.indices.push_back(indices);

            index << i, j;
            indices(0) = grid.nodeIndex(index);
            index << i + 1, j + 1;
            indices(1) = grid.nodeIndex(index);
            index << i, j + 1;
            indices(2) = grid.nodeIndex(index);
            mesh.indices.push_back(indices);
        }
    }

    // now create the particles
    ZIRAN_INFO("constructing particles ...");
    X.reserve(grid.numberNodes());
    TV x = TV::Zero();
    for (int i = 0; i < grid.numberNodesPerDimension(0); i++) {
        for (int j = 0; j < grid.numberNodesPerDimension(1); j++) {
            index << i, j;
            x.template head<2>() = grid.node(index);
            X.push_back(x);
            if (j == 0) {
                if (i == 0 || i == grid.numberNodesPerDimension(0) - 1) {
                    normals.push_back(Vector<T, dim>::Unit(0) * (-1));
                }
                else {
                    normals.push_back(Vector<T, dim>::Unit(0));
                }
            }
            else {
                normals.push_back(Vector<T, dim>::Zero());
            }
        }
    }
}

template <class T>
void constructMattressMesh(const Grid<T, 3>& grid, SimplexMesh<3>& mesh, StdVector<Vector<T, 3>>& X)
{
    using IVE = Vector<int, 4>;
    using IV = Vector<int, 3>;

    //first create the mesh
    IV index;
    IVE indices;
    ZIRAN_INFO("creating mesh ...");

    for (int i = 0; i < 3; i++)
        ZIRAN_INFO("grid number cells ", i, ": ", grid.numberCellsPerDimension(i));

    for (int i = 0; i < grid.numberCellsPerDimension(0); i++)
        for (int j = 0; j < grid.numberCellsPerDimension(1); j++)
            for (int k = 0; k < grid.numberCellsPerDimension(2); k++) {

                if ((i + j + k) % 2 != 0) {
                    index << i, j, k;
                    indices(0) = grid.nodeIndex(index);
                    index << i + 1, j, k;
                    indices(1) = grid.nodeIndex(index);
                    index << i, j + 1, k;
                    indices(2) = grid.nodeIndex(index);
                    index << i, j, k + 1;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);

                    index << i + 1, j, k;
                    indices(0) = grid.nodeIndex(index);
                    index << i + 1, j, k + 1;
                    indices(1) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k + 1;
                    indices(2) = grid.nodeIndex(index);
                    index << i, j, k + 1;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);

                    index << i, j + 1, k;
                    indices(0) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k;
                    indices(1) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k + 1;
                    indices(2) = grid.nodeIndex(index);
                    index << i + 1, j, k;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);

                    index << i, j + 1, k + 1;
                    indices(0) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k + 1;
                    indices(1) = grid.nodeIndex(index);
                    index << i, j, k + 1;
                    indices(2) = grid.nodeIndex(index);
                    index << i, j + 1, k;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);

                    index << i + 1, j, k;
                    indices(0) = grid.nodeIndex(index);
                    index << i, j, k + 1;
                    indices(1) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k + 1;
                    indices(2) = grid.nodeIndex(index);
                    index << i, j + 1, k;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);
                }
                else {
                    index << i, j, k;
                    indices(0) = grid.nodeIndex(index);
                    index << i + 1, j, k;
                    indices(1) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k;
                    indices(2) = grid.nodeIndex(index);
                    index << i + 1, j, k + 1;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);

                    index << i, j, k;
                    indices(0) = grid.nodeIndex(index);
                    index << i, j + 1, k;
                    indices(1) = grid.nodeIndex(index);
                    index << i, j + 1, k + 1;
                    indices(2) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);

                    index << i, j + 1, k + 1;
                    indices(0) = grid.nodeIndex(index);
                    index << i + 1, j, k + 1;
                    indices(1) = grid.nodeIndex(index);
                    index << i, j, k + 1;
                    indices(2) = grid.nodeIndex(index);
                    index << i, j, k;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);

                    index << i, j + 1, k + 1;
                    indices(0) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k + 1;
                    indices(1) = grid.nodeIndex(index);
                    index << i + 1, j, k + 1;
                    indices(2) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);

                    index << i, j + 1, k + 1;
                    indices(0) = grid.nodeIndex(index);
                    index << i, j, k;
                    indices(1) = grid.nodeIndex(index);
                    index << i + 1, j + 1, k;
                    indices(2) = grid.nodeIndex(index);
                    index << i + 1, j, k + 1;
                    indices(3) = grid.nodeIndex(index);
                    mesh.indices.push_back(indices);
                }
            }
    // now create the particles
    ZIRAN_INFO("constructing particles ...");

    X.reserve(grid.numberNodes());
    for (int i = 0; i < grid.numberNodesPerDimension(0); i++)
        for (int j = 0; j < grid.numberNodesPerDimension(1); j++)
            for (int k = 0; k < grid.numberNodesPerDimension(2); k++) {
                index << i, j, k;
                X.push_back(grid.node(index));
            }
}
} // namespace MESHING
template void MESHING::appendStrandToSegmesh<double, 2>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, int, SimplexMesh<1>&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&);
template void MESHING::appendStrandToSegmesh<double, 3>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, int, SimplexMesh<1>&, std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1, 0, 3, 1>>>&);
template void MESHING::appendStrandToSegmesh<float, 2>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, int, SimplexMesh<1>&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&);
template void MESHING::appendStrandToSegmesh<float, 3>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, int, SimplexMesh<1>&, std::vector<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 3, 1, 0, 3, 1>>>&);
template void MESHING::appendCircleToSegmesh<double, 2>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, double, int, SimplexMesh<1>&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&);
template void MESHING::appendCircleToSegmesh<double, 3>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, double, int, SimplexMesh<1>&, std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1, 0, 3, 1>>>&);
template void MESHING::appendCircleToSegmesh<float, 2>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, float, int, SimplexMesh<1>&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&);
template void MESHING::appendCircleToSegmesh<float, 3>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, float, int, SimplexMesh<1>&, std::vector<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 3, 1, 0, 3, 1>>>&);

template void MESHING::constructMattressMesh<double, 2>(Grid<double, 2> const&, SimplexMesh<2>&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&);
template void MESHING::constructMattressMesh<double, 3>(Grid<double, 2> const&, SimplexMesh<2>&, std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1, 0, 3, 1>>>&);
template void MESHING::constructMattressMesh<double>(Grid<double, 3> const&, SimplexMesh<3>&, std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1, 0, 3, 1>>>&);
template void MESHING::constructMattressMesh<float, 2>(Grid<float, 2> const&, SimplexMesh<2>&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&);
template void MESHING::constructMattressMesh<float, 3>(Grid<float, 2> const&, SimplexMesh<2>&, std::vector<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 3, 1, 0, 3, 1>>>&);
template void MESHING::constructMattressMesh<float>(Grid<float, 3> const&, SimplexMesh<3>&, std::vector<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 3, 1, 0, 3, 1>>>&);
template void MESHING::constructUnitSimplex<double, 2>(SimplexMesh<2>&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&);
template void MESHING::constructUnitSimplex<double, 3>(SimplexMesh<3>&, std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1, 0, 3, 1>>>&);
template void MESHING::constructUnitSimplex<float, 2>(SimplexMesh<2>&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&);
template void MESHING::constructUnitSimplex<float, 3>(SimplexMesh<3>&, std::vector<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 3, 1, 0, 3, 1>>>&);

template void MESHING::constructMattressMeshWithSurfaceMarker<double, 2>(Grid<double, 2> const&, SimplexMesh<2>&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&);
template void MESHING::constructMattressMeshWithSurfaceMarker<float, 2>(Grid<float, 2> const&, SimplexMesh<2>&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&);

template void MESHING::constructHairMeshWithSurfaceMarker<double, 2>(Grid<double, 2> const&, SimplexMesh<2>&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&);
template void MESHING::constructHairMeshWithSurfaceMarker<float, 2>(Grid<float, 2> const&, SimplexMesh<2>&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&);
} // namespace ZIRAN
