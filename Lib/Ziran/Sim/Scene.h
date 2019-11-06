#ifndef SCENE_H
#define SCENE_H

#include <Ziran/CS/Util/DataDir.h>
#include <Ziran/Math/Geometry/Elements.h>
#include <Ziran/Math/Geometry/ObjIO.h>
#include <Ziran/Math/Geometry/SimplexMesh.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Math/Geometry/PolyIO.h>
#include <Ziran/Physics/LagrangianForce/LagrangianForce.h>

namespace ZIRAN {
/**
       scene needs:

       particles
       dirichlet_position_dofs
       segmesh_to_write
       trimesh_to_write
       forces
       element_managers
       getElementManager()
     */

/**
   Helper struct for dynamic lookup of the element manager
*/
template <class T, int dim>
struct ElementManagerProxy {
    StdVector<std::unique_ptr<ElementManager<T, dim>>>& element_managers;
    ElementManagerProxy(StdVector<std::unique_ptr<ElementManager<T, dim>>>& element_managers)
        : element_managers(element_managers)
    {
    }

    template <int manifold_dim>
    operator SimplexElements<T, manifold_dim, dim>&()
    {
        return get<SimplexElements<T, manifold_dim, dim>>();
    }

    /**
       Get the element manager with type Type or insert it if it's not there
       Not threadsafe
    */
    template <class Type>
    Type& get()
    {
        for (auto& emp : element_managers)
            if (Type* em = dynamic_cast<Type*>(emp.get()))
                return *em;
        // We didn't find it
        std::unique_ptr<Type> emp = std::make_unique<Type>();
        element_managers.emplace_back(std::move(emp));
        return *dynamic_cast<Type*>(element_managers.back().get());
    }
};

template <class T, int dim>
class Scene {
public:
    using TVStack = Matrix<T, dim, Eigen::Dynamic>;
    Particles<T, dim> particles;
    StdVector<int> dirichlet_position_dofs;
    StdVector<std::unique_ptr<ElementManager<T, dim>>> element_managers;
    SimplexMesh<2> trimesh_to_write;
    SimplexMesh<1> segmesh_to_write;

    Scene()
    {
        particles.X.registerReorderingCallback(
            [&](const DataArrayBase::ValueIDToNewValueID& old_to_new) {
                for (auto& element : trimesh_to_write.indices)
                    for (int d = 0; d < 3; d++)
                        element(d) = old_to_new(element(d));

                for (auto& element : segmesh_to_write.indices)
                    for (int d = 0; d < 2; d++)
                        element(d) = old_to_new(element(d));

                for (auto& dirichlet_dof : dirichlet_position_dofs) {
                    int p = dirichlet_dof / dim;
                    int off = dirichlet_dof % dim;
                    dirichlet_dof = dim * old_to_new(p) + off;
                }
            });
    }

    // These are meshed forces.
    // SceneInitialization will create it with mesh etc.
    StdVector<std::unique_ptr<LagrangianForce<T, dim>>> forces;

    T totalEnergy()
    {
        return tbb::parallel_reduce(tbb::blocked_range<int>(0, forces.size()), (double)0,
            [&](const tbb::blocked_range<int>& range, const double& e) -> double {
                double psi = e;
                for (int i = range.begin(); i < range.end(); ++i)
                    psi += forces[i]->totalEnergy();
                return psi;
            },
            [&](const double& x, const double& y) -> double {
                return x + y;
            });
    }

    void addScaledForces(const T scale, TVStack& f)
    {
        if (element_managers.size() == 0)
            return;
        int element_partitions = element_managers[0]->pads.size();

        if (element_partitions) { // parallel
            // splat forces to pads
            for (auto& em : element_managers)
                tbb::parallel_for(tbb::blocked_range<int>(0, em->pads.size()),
                    [&](const tbb::blocked_range<int>& r) {
                        for (int k = r.begin(); k < r.end(); ++k) {
                            em->pads[k].setZero(); // zero out pads
                            for (auto& lf : forces)
                                if (lf->isThisMyElementManager(*em))
                                    lf->splatToPad(1, k); // k is pad_id
                        }
                    });

            // gather forces
            for (auto& em : element_managers)
                em->gatherFromPadsToTVStack(scale, f); // this is internally a parallel_for
        }

        else { // serial
            for (auto& lf : forces)
                lf->addScaledForces(scale, f);
        }
    }

    void addScaledForceDifferentials(const T scale, const TVStack& dx, TVStack& df)
    {
        if (element_managers.size() == 0)
            return;
        int element_partitions = element_managers[0]->pads.size();

        if (element_partitions) { // parallel

            // splat force differentials to pads
            for (auto& em : element_managers)
                tbb::parallel_for(tbb::blocked_range<int>(0, em->pads.size()),
                    [&](const tbb::blocked_range<int>& r) {
                        for (int k = r.begin(); k < r.end(); ++k) { // k is pad id
                            em->pads[k].setZero(); // zero out pads

                            // build pad_dx
                            TVStack& pad_dx = em->pads_dx[k];
                            for (int t = 0; t < pad_dx.cols(); t++) // t is pad k's local particle index
                                pad_dx.col(t) = dx.col(em->particle_l2g[k][t]);

                            for (auto& lf : forces)
                                if (lf->isThisMyElementManager(*em))
                                    lf->splatDifferentialToPad(1, k, pad_dx);
                        }
                    });

            // gather force differentials
            for (auto& em : element_managers)
                em->gatherFromPadsToTVStack(scale, df); // this is internally a parallel_for
        }

        else { // serial
            for (auto& lf : forces)
                lf->addScaledForceDifferential(scale, dx, df);
        }
    }

    ElementManagerProxy<T, dim> getElementManager()
    {
        return ElementManagerProxy<T, dim>(element_managers);
    }

    template <class Type>
    bool exist()
    {
        for (auto& emp : element_managers)
            if (Type* em = dynamic_cast<Type*>(emp.get()))
                return true;
        return false;
    }

    /* binary ver 1 */
    void writeState(std::ostream& out, const std::string& boundary_filename, const std::string& segment_filename, DataDir& output_dir, int binary_ver, bool write_meshes = true)
    {
        if (write_meshes && (trimesh_to_write.indices.size() != 0)) {
            std::ofstream boundary_file = output_dir.openTextOutput(boundary_filename);
            writePositionObj(boundary_file, particles.X.array);
            writeFacesObj(boundary_file, trimesh_to_write.indices);
        }
        if (write_meshes && segmesh_to_write.indices.size() != 0) {
            std::ofstream segmesh_file = output_dir.openTextOutput(segment_filename);
            writeSegmeshPoly(segmesh_file, particles.X.array, segmesh_to_write.indices);
        }
        particles.writeData(out);
        for (const auto& em : element_managers)
            em->writeData(out);

        writeSTDVector(out, trimesh_to_write.indices);
        writeSTDVector(out, segmesh_to_write.indices);
    }

    /* binary ver 1 */
    void readState(std::istream& in, int binary_ver)
    {
        particles.readData(in);
        for (auto& em : element_managers)
            em->readData(in);

        readSTDVector(in, trimesh_to_write.indices);
        readSTDVector(in, segmesh_to_write.indices);
    }
};
} // namespace ZIRAN

#endif
