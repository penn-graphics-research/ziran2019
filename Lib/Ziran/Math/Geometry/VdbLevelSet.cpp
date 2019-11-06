#ifdef ZIRAN_WITH_VDB
#include <Ziran/Math/Geometry/VdbLevelSet.h>
#include <Ziran/Math/Geometry/OutputPolyMesh.h>
#include <Ziran/CS/Util/DataDir.h>

#include <openvdb/io/File.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/VolumeToMesh.h>
#undef B2

namespace ZIRAN {

template <class T, int dim>
VdbLevelSet<T, dim>::VdbLevelSet(const std::string file_in)
{
    constructFromVdbFile(file_in);
}

template <class T, int dim>
VdbLevelSet<T, dim>::VdbLevelSet()
{
}

template <class T, int dim>
VdbLevelSet<T, dim>::~VdbLevelSet()
{
}

template <class T, int dim>
void VdbLevelSet<T, dim>::constructFromVdbFile(const std::string& file_in)
{
    std::string filename = DataDir().absolutePath(file_in);

    openvdb::io::File file(filename);
    file.open();
    openvdb::GridPtrVecPtr my_grids = file.getGrids();
    file.close();
    int count = 0;
    for (openvdb::GridPtrVec::iterator iter = my_grids->begin(); iter != my_grids->end(); ++iter) {
        grid = openvdb::gridPtrCast<GridT>(*iter);
        count++;
    }
    ZIRAN_ASSERT(count == 1, "Vdb file to load should only contain one levelset.");

    recomputeGradient();
}

template <class T, int dim>
void VdbLevelSet<T, dim>::recomputeGradient()
{
    // construct gradient field
    openvdb::tools::Gradient<GridT> mg(*grid);
    grad_phi = mg.process();
}

// This is only correct at narrow band.
// Otherwise, it only gets the sign right.
template <class T, int dim>
T VdbLevelSet<T, dim>::signedDistance(const TV& X_input) const
{
    TV3 X;
    X.setZero();
    for (int d = 0; d < dim; d++)
        X(d) = X_input(d);
    openvdb::tools::GridSampler<TreeT, openvdb::tools::BoxSampler> interpolator(grid->constTree(), grid->transform());
    openvdb::math::Vec3<T> P(X(0), X(1), X(2));
    float phi = interpolator.wsSample(P); //ws denotes world space
    return (T)phi;
}

template <class T, int dim>
bool VdbLevelSet<T, dim>::inside(const TV& X_input) const
{
    return signedDistance(X_input) <= 0;
}

// Normals are evaluated as interpolated gradients.
// It is only safe to use when the vdb narrow band is wide enough.
// A safe way to achieve this is to fill the interor when creating the levelset.
template <class T, int dim>
auto VdbLevelSet<T, dim>::normal(const TV& X_input) const -> TV
{
    TV3 X;
    X.setZero();
    for (int d = 0; d < dim; d++)
        X(d) = X_input(d);
    openvdb::tools::GridSampler<GradientTreeT, openvdb::tools::BoxSampler> interpolator(grad_phi->constTree(), grad_phi->transform());
    openvdb::math::Vec3<T> P(X(0), X(1), X(2));
    auto grad_phi = interpolator.wsSample(P); //ws denotes world space
    TV result;
    for (int d = 0; d < dim; d++)
        result(d) = grad_phi(d);
    T norm = result.norm();
    if (norm != 0)
        return result / norm;
    else
        return TV::Zero();
}

template <class T, int dim>
bool VdbLevelSet<T, dim>::query(const TV& X, T& signed_distance, TV& n, TM* hessian) const
{
    signed_distance = signedDistance(X);
    bool inside = (signed_distance <= 0);
    n = normal(X);
    ZIRAN_ASSERT(!hessian, "not implemented");
    return inside;
}

template <class T, int dim>
void VdbLevelSet<T, dim>::getBounds(TV& min_corner, TV& max_corner) const
{
    min_corner.setZero();
    max_corner.setZero();
    // TODO: double check is this world space or material space
    openvdb::CoordBBox box = grid->evalActiveVoxelBoundingBox();
    auto world_min = grid->indexToWorld(box.min());
    auto world_max = grid->indexToWorld(box.max());

    for (size_t d = 0; d < dim; d++) {
        min_corner(d) = world_min(d);
        max_corner(d) = world_max(d);
    }
}

template <class T, int dim>
std::unique_ptr<AnalyticLevelSet<T, dim>> VdbLevelSet<T, dim>::createPtr() const
{
    return std::make_unique<VdbLevelSet<T, dim>>(*this);
}

template <class T, int dim>
std::function<OutputPolyMesh()> VdbLevelSet<T, dim>::createMesh(T dx) const
{
    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    openvdb::tools::volumeToMesh(*grid, points, quads, 0.0);
    std::vector<int32_t> face_vertex_count(quads.size(), 4);
    return [points(std::move(points)),
               quads(std::move(quads)),
               face_vertex_count(std::move(face_vertex_count))]() {
        return OutputPolyMesh{
            points.data()->asPointer(), points.size(),
            (int32_t*)quads.data()->asPointer(),
            quads.size() * 4,
            face_vertex_count.data(),
            face_vertex_count.size()
        };
    };
}
template class VdbLevelSet<float, 2>;
template class VdbLevelSet<float, 3>;
template class VdbLevelSet<double, 2>;
template class VdbLevelSet<double, 3>;
} // namespace ZIRAN
#endif
