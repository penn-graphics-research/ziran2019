#pragma once

#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/DataStructure/HashTable.h>
#include <vector>
#include <Ziran/CS/Util/BinaryIO.h>
#include "SimplexMesh.h"

namespace ZIRAN {

// Meshes with real vertices, instead of particles as vertices
template <typename T, int dim>
class VertexMesh {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    using TMAffine = Matrix<T, dim + 1, dim + 1>;
    using Mesh = SimplexMesh<dim - 1>;
    StdVector<TV> vertices;
    Mesh mesh;

    VertexMesh()
    {
    }

    // TODO: fix rvalue hack here
    VertexMesh(const StdVector<TV>& vertices, StdVector<typename Mesh::IV> indices)
        : vertices(vertices)
        , mesh(std::move(indices))
    {
    }

    VertexMesh(StdVector<TV>&& vertices, StdVector<typename Mesh::IV>&& indices)
        : vertices(std::move(vertices))
        , mesh(std::move(indices))
    {
    }
};

template <class T, int dim>
struct VertexMeshTag {
};

template <class T, int dim>
struct RW<VertexMesh<T, dim>> {
    using Tag = VertexMeshTag<T, dim>;
};

template <class T, int dim>
void writeHelper(std::ostream& out, const VertexMesh<T, dim>& mesh, VertexMeshTag<T, dim>)
{
    writeSTDVector(out, mesh.vertices);
    //TODO: IO for a whole SimplexMesh (needs HashTable IO first.)
    writeSTDVector(out, mesh.mesh.indices);
}
template <class T, int dim>
VertexMesh<T, dim> readHelper(std::istream& in, VertexMeshTag<T, dim>)
{
    VertexMesh<T, dim> mesh;
    readSTDVector(in, mesh.vertices);
    readSTDVector(in, mesh.mesh.indices);
    return mesh;
}
} // namespace ZIRAN