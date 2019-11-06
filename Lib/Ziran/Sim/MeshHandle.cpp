#include <Ziran/Sim/MeshHandle.h>
#include <Ziran/CS/Util/DataDir.h>
#include <Ziran/Math/Geometry/ObjIO.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Math/Geometry/UscHairDataIO.h>
#include <Ziran/Math/Geometry/VtkIO.h>
#include <Ziran/Math/Geometry/PolyIO.h>

#include <Eigen/LU>
namespace ZIRAN {

template <class T, int dim>
MeshReader<T, dim, SimplexMesh<1>>::MeshReader(const std::string& filename, const T percentage_to_keep_in)
    : filename(filename)
    , percentage_to_keep(percentage_to_keep_in)
{
}

// readSegMesh
template <class T, int dim>
void MeshReader<T, dim, SimplexMesh<1>>::operator()(SimplexMesh<1>& mesh, StdVector<TV>& X0)
{
    std::string ext = filename.substr(filename.find_last_of('.') + 1);
    std::string absolute_path = DataDir().absolutePath(filename);
    ZIRAN_INFO("Reading mesh ", filename);
    if (ext == "data") {
        readUscHairData(absolute_path, X0, mesh.indices, percentage_to_keep);
    }
    else if (ext == "poly") {
        ZIRAN_ASSERT(percentage_to_keep == (T)1, "reading poly segmesh doesn't support downsample");
        readSegmeshPoly(absolute_path, X0, mesh.indices);
    }
    else {
        ZIRAN_ASSERT(false, "Unsupported segmesh file ", filename);
    }
}

template <class T, int dim>
MeshReader<T, dim, SimplexMesh<2>>::MeshReader(const std::string& filename)
    : filename(filename)
{
}

template <class T, int dim>
void MeshReader<T, dim, SimplexMesh<2>>::operator()(SimplexMesh<2>& mesh, StdVector<TV>& X0)
{
    std::string ext = filename.substr(filename.find_last_of('.') + 1);
    std::string absolute_path = DataDir().absolutePath(filename);
    ZIRAN_INFO("Reading mesh ", filename);
    if (ext == "obj") {
        if (dim == 2) {
            StdVector<Vector<T, 3>> X3d;
            readTrimeshObj(absolute_path, X3d, mesh.indices);
            for (auto& p : X3d) {
                TV X2d;
                X2d << p(0), p(1);
                X0.emplace_back(X2d);
            }
        }
        else
            readTrimeshObj(absolute_path, X0, mesh.indices);
    }
    else {
        ZIRAN_ASSERT(false, "Unsupported trimesh file ", filename);
    }
}

template <class T, int dim>
MeshReader<T, dim, SimplexMesh<3>>::MeshReader(const std::string& filename)
    : filename(filename)
{
}

template <class T, int dim>
void MeshReader<T, dim, SimplexMesh<3>>::operator()(SimplexMesh<3>& mesh, StdVector<TV>& X0)
{
    std::string ext = filename.substr(filename.find_last_of('.') + 1);
    std::string absolute_path = DataDir().absolutePath(filename);
    ZIRAN_INFO("Reading mesh ", filename);
    if (ext == "vtk") {
        readTetmeshVtk(absolute_path, X0, mesh.indices);
    }
    else {
        ZIRAN_ASSERT(false, "Unsupported mesh file ", filename);
    }
}

template <class T, int dim, class TMesh>
MeshHandle<T, dim, TMesh>::MeshHandle(Particles<T, dim>& particles,
    std::shared_ptr<TMesh> mesh_in,
    Range particle_range)
    : particles(particles)
    , mesh(std::move(mesh_in))
    , particle_range(particle_range)
    , particle_index_offset(particle_range.lower)
{
}

template <class T, int dim, class TMesh>
MeshHandle<T, dim, TMesh>::MeshHandle(Particles<T, dim>& particles,
    std::shared_ptr<TMesh> mesh_in,
    Range particle_range,
    int particle_index_offset)
    : particles(particles)
    , mesh(std::move(mesh_in))
    , particle_range(particle_range)
    , particle_index_offset(particle_index_offset)
{
}

template <class T, int dim, class TMesh>
MeshHandle<T, dim, TMesh>::MeshHandle(Particles<T, dim>& particles,
    std::shared_ptr<TMesh> mesh_in,
    const std::function<void(TMesh&, StdVector<TV>&)> mesh_constructor)
    : particles(particles)
    , mesh(std::move(mesh_in))
    , particle_index_offset(particles.count)
{
    StdVector<TV>& X0 = particles.X.array;
    particle_range.lower = X0.size();
    mesh_constructor(*mesh, X0);
    particle_range.upper = X0.size();
    particles.X.ranges.append(particle_range);
    particles.add(particles.V_name(), particle_range, TV::Zero());
    particles.add(particles.mass_name(), particle_range, (T)0);
}

template <class T, int dim, class TMesh>
MeshHandle<T, dim, TMesh>::MeshHandle(Particles<T, dim>& particles,
    std::shared_ptr<TMesh> mesh_in,
    const std::function<void(TMesh&, StdVector<TV>&, StdVector<TV>&)> mesh_constructor)
    : particles(particles)
    , mesh(std::move(mesh_in))
    , particle_index_offset(particles.count)
{
    StdVector<TV>& X0 = particles.X.array;
    particle_range.lower = X0.size();
    StdVector<TV> normals;
    mesh_constructor(*mesh, X0, normals);
    particle_range.upper = X0.size();
    particles.X.ranges.append(particle_range);
    particles.add(particles.V_name(), particle_range, TV::Zero());
    particles.add(particles.mass_name(), particle_range, (T)0);
    particles.add(AttributeName<TV>("normal"), particle_range, normals);
}

template <class T, int dim, class TMesh>
MeshHandle<T, dim, TMesh> MeshHandle<T, dim, TMesh>::copy()
{
    return MeshHandle(particles, mesh, [=](TMesh&, StdVector<TV>& X0) {
        X0.insert(X0.end(),
            X0.begin() + particle_range.lower,
            X0.begin() + particle_range.upper);
    });
}

template <class T, int dim, class TMesh>
void MeshHandle<T, dim, TMesh>::transform(const std::function<void(int, Ref<T>, TV&, TV&)>& mapping)
{
    for (int i = particle_range.lower; i < particle_range.upper; ++i) {
        mapping(i - particle_range.lower, particles.mass[i], particles.X[i], particles.V[i]);
    }
}

template <class T, int dim>
MeshHandle<T, dim, SimplexMesh<1>> constructTrimeshBoundaryWithParticles(MeshHandle<T, dim, SimplexMesh<2>>& trimesh)
{
    int trimesh_offset = trimesh.particle_index_offset;

    std::shared_ptr<SimplexMesh<1>> segmesh = std::make_shared<SimplexMesh<1>>();

    trimesh.mesh->initializeBoundaryElements();

    int old_count = trimesh.particles.count;
    int offset = 0;

    HashTable<int, int> vertex_hash; // map tri mesh point to seg mesh point (both global indices)

    {
        auto ap = trimesh.particles.appender();
        for (auto segment : trimesh.mesh->boundary_indices) {
            int A = segment(0), B = segment(1);
            if (vertex_hash.get(A) == nullptr) {
                vertex_hash[A] = offset++;
                ap.append(/*mass*/ (T)0, trimesh.particles.X[A + trimesh_offset], trimesh.particles.V[A + trimesh_offset]);
            }
            if (vertex_hash.get(B) == nullptr) {
                vertex_hash[B] = offset++;
                ap.append(/*mass*/ (T)0, trimesh.particles.X[B + trimesh_offset], trimesh.particles.V[B + trimesh_offset]);
            }
            Vector<int, 2> seg;
            seg << vertex_hash[A], vertex_hash[B];
            segmesh->indices.emplace_back(seg);
        }
    }

    MeshHandle<T, dim, SimplexMesh<1>> handle(trimesh.particles, segmesh, Range{ old_count, old_count + offset });
    return handle;
}

//invert the tets with wrong orientation in initial tetmesh
template <class T, int dim>
void correctInverted(MeshHandle<T, dim, SimplexMesh<3>>& tet_handle)
{
    using TMesh = SimplexMesh<3>;
    Matrix<T, 4, 4> tet;
    int inverted_count = 0;
    std::shared_ptr<TMesh> mesh = tet_handle.mesh;
    for (size_t i = 0; i < mesh->indices.size(); i++) {
        Vector<int, 4> element = mesh->indices[i];
        tet.col(0) = Vector<T, dim + 1>::Ones();
        tet.block(0, 1, 1, 3) = tet_handle.particles.X[element(0)].transpose();
        tet.block(1, 1, 1, 3) = tet_handle.particles.X[element(1)].transpose();
        tet.block(2, 1, 1, 3) = tet_handle.particles.X[element(2)].transpose();
        tet.block(3, 1, 1, 3) = tet_handle.particles.X[element(3)].transpose();

        T volume = tet.determinant();
        if (volume < 0) {
            inverted_count++;
            std::swap(element(2), element(3));
            ZIRAN_INFO("corrected tet number ", i, "with vertices ", element);
        }
        mesh->indices[i] = element;
    }
    ZIRAN_INFO("total inverted tets count: ", inverted_count);
}

template class MeshReader<double, 2, SimplexMesh<1>>;
template class MeshReader<double, 3, SimplexMesh<1>>;
template class MeshReader<double, 2, SimplexMesh<2>>;
template class MeshReader<double, 3, SimplexMesh<2>>;
template class MeshReader<double, 3, SimplexMesh<3>>;
template class MeshReader<float, 2, SimplexMesh<1>>;
template class MeshReader<float, 3, SimplexMesh<1>>;
template class MeshReader<float, 2, SimplexMesh<2>>;
template class MeshReader<float, 3, SimplexMesh<2>>;
template class MeshReader<float, 3, SimplexMesh<3>>;

template class MeshHandle<double, 2, SimplexMesh<1>>;
template class MeshHandle<double, 2, SimplexMesh<2>>;
template class MeshHandle<double, 3, SimplexMesh<1>>;
template class MeshHandle<double, 3, SimplexMesh<2>>;
template class MeshHandle<double, 3, SimplexMesh<3>>;
template class MeshHandle<float, 2, SimplexMesh<1>>;
template class MeshHandle<float, 2, SimplexMesh<2>>;
template class MeshHandle<float, 3, SimplexMesh<1>>;
template class MeshHandle<float, 3, SimplexMesh<2>>;
template class MeshHandle<float, 3, SimplexMesh<3>>;
template MeshHandle<double, 2, SimplexMesh<1>> constructTrimeshBoundaryWithParticles<double, 2>(MeshHandle<double, 2, SimplexMesh<2>>&);
template MeshHandle<double, 3, SimplexMesh<1>> constructTrimeshBoundaryWithParticles<double, 3>(MeshHandle<double, 3, SimplexMesh<2>>&);
template MeshHandle<float, 2, SimplexMesh<1>> constructTrimeshBoundaryWithParticles<float, 2>(MeshHandle<float, 2, SimplexMesh<2>>&);
template MeshHandle<float, 3, SimplexMesh<1>> constructTrimeshBoundaryWithParticles<float, 3>(MeshHandle<float, 3, SimplexMesh<2>>&);

template void correctInverted<double, 3>(MeshHandle<double, 3, SimplexMesh<3>>&);
template void correctInverted<float, 3>(MeshHandle<float, 3, SimplexMesh<3>>&);
} // namespace ZIRAN
