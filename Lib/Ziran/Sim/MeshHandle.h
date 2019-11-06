#ifndef MESH_HANDLE_H
#define MESH_HANDLE_H

#include <Ziran/Math/Geometry/SimplexMesh.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/CS/DataStructure/DisjointRanges.h>
#include <Ziran/CS/DataStructure/Ref.h>
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {

template <class T, int dim, class TMesh>
class MeshReader;

template <class T, int dim>
class MeshReader<T, dim, SimplexMesh<1>> {
public:
    using TV = Vector<T, dim>;

    std::string filename;
    T percentage_to_keep;

    MeshReader(const std::string& filename, const T percentage_to_keep_in = 1);

    void operator()(SimplexMesh<1>& mesh, StdVector<TV>& X0);
};

template <class T, int dim>
class MeshReader<T, dim, SimplexMesh<2>> {
public:
    using TV = Vector<T, dim>;

    std::string filename;
    MeshReader(const std::string& filename);

    void operator()(SimplexMesh<2>& mesh, StdVector<TV>& X0);
};

template <class T, int dim>
class MeshReader<T, dim, SimplexMesh<3>> {
public:
    using TV = Vector<T, dim>;

    std::string filename;
    MeshReader(const std::string& filename);

    void operator()(SimplexMesh<3>& mesh, StdVector<TV>& X0);
};

template <class T, int dim, class TMesh>
class MeshHandle {
public:
    Particles<T, dim>& particles;
    std::shared_ptr<TMesh> mesh;
    Range particle_range;
    int particle_index_offset;

    using TV = Vector<T, dim>;

    MeshHandle(Particles<T, dim>& particles, std::shared_ptr<TMesh> mesh, Range particle_range);

    MeshHandle(Particles<T, dim>& particles, std::shared_ptr<TMesh> mesh_in, Range particle_range, int particle_index_offset);

    MeshHandle(Particles<T, dim>& particles, std::shared_ptr<TMesh> mesh, const std::function<void(TMesh&, StdVector<TV>&)> mesh_constructor);

    MeshHandle(Particles<T, dim>& particles, std::shared_ptr<TMesh> mesh, const std::function<void(TMesh&, StdVector<TV>&, StdVector<TV>&)> mesh_constructor);

    // Creates a copy with new particles
    MeshHandle copy();

    void transform(const std::function<void(int, Ref<T>, TV&, TV&)>& mapping);
};

/**
   This constructs a segmesh handle from a trimesh handle using its boundary. Works in both 2D and 3D.
   This will append new particles to the same 'particles' reference of the trimeshhandle.
   The segmesh index into these new particles.
*/
template <class T, int dim>
MeshHandle<T, dim, SimplexMesh<1>> constructTrimeshBoundaryWithParticles(MeshHandle<T, dim, SimplexMesh<2>>& trimesh);

/*
 * The particle range is only for the newly added control points
 */

template <class T, int dim>
void correctInverted(MeshHandle<T, dim, SimplexMesh<3>>& tet_handle);
} // namespace ZIRAN
#endif
