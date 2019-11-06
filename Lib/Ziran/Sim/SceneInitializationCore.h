#ifndef SCENE_INITIALIZATION_CORE_H
#define SCENE_INITIALIZATION_CORE_H
#include <Ziran/Sim/Scene.h>
#include <Ziran/CS/Util/Forward.h>
#include <functional>

namespace ZIRAN {

template <class T, int dim, class TMesh>
class MeshReader;

template <class T, int dim>
class SceneInitialization {
public:
    using TV = Vector<T, dim>;
    using TVI = Vector<int, dim>;

    Scene<T, dim>& scene;
    /*
       scene has:

       particles
       dirichlet_position_dofs
       segmesh_to_write
       trimesh_to_write
       forces
       getElementManager()
     */

    SceneInitialization(Scene<T, dim>& scene_in)
        : scene(scene_in)
    {
    }

    void setBoundaryConditions(const std::function<bool(int, TV&)>& is_dirichlet)
    {
        for (int i = 0; i < scene.particles.count; i++)
            for (int d = 0; d < dim; d++) {
                int dof = i * dim + d;
                if (is_dirichlet(dof, scene.particles.X[i])) {
                    scene.dirichlet_position_dofs.push_back(dof);
                }
            }
    }

    void setDirichletBoundaryFromFile(const std::string& file_name)
    {
        std::ifstream fs;
        std::string absolute_path = DataDir().absolutePath(file_name);
        fs.open(absolute_path);
        ZIRAN_ASSERT(fs, "could not open ", file_name);
        int index;
        while (fs >> index) {
            for (int d = 0; d < dim; ++d) {
                int dof_id = index * dim + d;
                scene.dirichlet_position_dofs.push_back(dof_id);
            }
        }
    }

    void addMeshToWrite(SimplexMesh<1>& mesh, int offset)
    {
        auto& indices = scene.segmesh_to_write.indices;
        for (const auto& i : mesh.indices)
            indices.emplace_back(i.array() + offset);
    }

    void addMeshToWrite(SimplexMesh<2>& mesh, int offset)
    {
        auto& indices = scene.trimesh_to_write.indices;
        for (const auto& i : mesh.indices)
            indices.emplace_back(i.array() + offset);
    }

    void addMeshToWrite(SimplexMesh<3>& mesh, int offset)
    {
        mesh.initializeBoundaryElements();
        auto& indices = scene.trimesh_to_write.indices;
        for (const auto& i : mesh.boundary_indices)
            indices.emplace_back(i.array() + offset);
    }

    template <class TMesh>
    MeshHandle<T, dim, TMesh> createMesh(std::function<void(TMesh&, StdVector<TV>&)> mesh_constructor)
    {
        return MeshHandle<T, dim, TMesh>(scene.particles, std::make_shared<TMesh>(), mesh_constructor);
    }

    template <class TMesh>
    MeshHandle<T, dim, TMesh> createMesh(std::function<void(TMesh&, StdVector<TV>&, StdVector<TV>&)> mesh_constructor)
    {
        return MeshHandle<T, dim, TMesh>(scene.particles, std::make_shared<TMesh>(), mesh_constructor);
    }
};
} // namespace ZIRAN
#endif
