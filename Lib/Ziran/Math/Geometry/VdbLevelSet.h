#ifdef ZIRAN_WITH_VDB
#ifndef VDB_LEVEL_SET_H
#define VDB_LEVEL_SET_H

#include <Ziran/Math/Geometry/AnalyticLevelSet.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/GridOperators.h>
#undef B2

namespace ZIRAN {

/**
   This is a level set wrapper class that internally uses an openvdb level set to represent a 3D sdf.
   phi and normals are only correct near narrow band.
   Otherwise, only the sign of phi is meaningful.
*/
template <class T, int dim>
class VdbLevelSet : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    using TV3 = Vector<T, 3>;
    typedef typename openvdb::Grid<typename openvdb::tree::Tree4<float, 5, 4, 3>::Type> GridT;
    typedef typename GridT::TreeType TreeT;
    typedef typename openvdb::tools::ScalarToVectorConverter<GridT>::Type GradientGridT;
    typedef typename GradientGridT::TreeType GradientTreeT;

    typename GridT::Ptr grid;
    typename GradientGridT::Ptr grad_phi;

    VdbLevelSet(const std::string file_in);

    VdbLevelSet();

    ~VdbLevelSet() override;

    void constructFromVdbFile(const std::string& file_in);

    void recomputeGradient();

    // This is only correct at narrow band.
    // Otherwise, it only gets the sign right.
    T signedDistance(const TV& X_input) const override;

    bool inside(const TV& X_input) const override;

    // Normals are evaluated as interpolated gradients.
    // It is only safe to use when the vdb narrow band is wide enough.
    // A safe way to achieve this is to fill the interor when creating the levelset.
    TV normal(const TV& X_input) const override;

    bool query(const TV& X, T& signed_distance, TV& n, TM* hessian = 0) const override;

    void getBounds(TV& min_corner, TV& max_corner) const override;

    std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override;
    std::function<OutputPolyMesh()> createMesh(T dx) const override;
};
} // namespace ZIRAN

#endif
#endif
