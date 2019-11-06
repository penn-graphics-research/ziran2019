#ifndef ANALYTIC_LEVEL_SET_H
#define ANALYTIC_LEVEL_SET_H
#include <Ziran/CS/Util/Forward.h>
#include <Eigen/Core>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/Math/MathTools.h>
#include <Ziran/Math/Nonlinear/AutoDiff.h>
#include <Ziran/Math/Geometry/OutputPolyMesh.h>
#include <math.h>
#include <memory>

namespace ZIRAN {

template <class T, int dim>
class AnalyticLevelSet {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual ~AnalyticLevelSet() {}

    // Inclusive on points on the surface
    virtual bool inside(const TV& X) const;

    /**
       Signed distance function. Negative if X is inside the object.
     */
    virtual T signedDistance(const TV& X) const = 0;

    virtual TV closestPointOnSurface(const TV& X) const;

    /**
       Derivative of signed distance function at X.
     */
    virtual TV normal(const TV& X) const = 0;

    /**
       Hessian of signed distance function at X.
     */
    virtual TM hessian(const TV& X) const;

    virtual bool query(const TV& X, T& signed_distance, TV& normal, TM* hessian = 0) const;

    virtual bool queryInside(const TV& X, T& signed_distance, TV& normal, const T phi_critical = 0) const;

    virtual void getBounds(TV& min_corner, TV& max_corner) const;

    virtual T volume() const;

    virtual std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const = 0;

    // By default this uses marching cubes to construct a mesh
    virtual std::function<OutputPolyMesh()> createMesh(T dx) const;
};

template <class T, int dim>
class DisjointUnionLevelSet : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>> level_sets;

    DisjointUnionLevelSet() {}

    DisjointUnionLevelSet(StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>>&& level_sets)
        : level_sets(std::move(level_sets))
    {
    }

    DisjointUnionLevelSet(const StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>>& level_sets_input);

    virtual ~DisjointUnionLevelSet() {}

    void add(const AnalyticLevelSet<T, dim>& ls);

    T signedDistance(const TV& X) const override;

    virtual TV normal(const TV& X) const override;

    void getBounds(TV& min_corner, TV& max_corner) const override;

    T volume() const override;

    virtual std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override;
};

/**
   setA - setB
   This operation takes intersection of A and -B.
   setB should be completely inside of setA for this operation to be accurate
*/
template <class T, int dim>
class DifferenceLevelSet : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>> level_sets;

    DifferenceLevelSet() {}

    DifferenceLevelSet(StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>>&& level_sets);

    DifferenceLevelSet(const StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>>& level_sets_input);

    virtual ~DifferenceLevelSet() {}

    void add(const AnalyticLevelSet<T, dim>& ls1, const AnalyticLevelSet<T, dim>& ls2);

    T signedDistance(const TV& X) const override;

    virtual TV normal(const TV& X) const override;

    void getBounds(TV& min_corner, TV& max_corner) const override;

    T volume() const override;

    virtual std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override;
};

template <class T, int dim>
class HalfSpace : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;

    Eigen::Matrix<T, dim, 1, Eigen::DontAlign> origin;
    Eigen::Matrix<T, dim, 1, Eigen::DontAlign> outward_normal;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    HalfSpace(const TV& origin_in, const TV& outward_normal_in);

    ~HalfSpace() override {}

    // Inclusive on points on the surface
    bool inside(const TV& X) const override;

    T signedDistance(const TV& X) const override;

    TV closestPointOnSurface(const TV& X) const override;

    TV normal(const TV& X) const override;

    TM hessian(const TV& X) const override;

    bool query(const TV& X, T& signed_distance, TV& normal, TM* hessian = 0) const override;

    std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override;
    std::function<OutputPolyMesh()> createMesh(T dx) const override;
};

template <class T, int dim>
class AnalyticBox;

template <class T, int dim>
class AxisAlignedAnalyticBox : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;

    AnalyticBox<T, dim> box;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AxisAlignedAnalyticBox(const TV& min_corner, const TV& max_corner);

    ~AxisAlignedAnalyticBox() override {}

    T signedDistance(const TV& X) const override;

    TV normal(const TV& X) const override;

    void getBounds(TV& min_corner, TV& max_corner) const override;

    T volume() const override;

    std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override;
};

//
// TODO: Fix for (X == center) case.
// Most functions are broken in that case.
// Should just choose n = (1,0,0) in that case.
//
template <class T, int dim>
class Sphere : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;

    Eigen::Matrix<T, dim, 1, Eigen::DontAlign> center;
    T radius;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Sphere(const TV& center, const T& radius);

    ~Sphere() override {}

    // Inclusive on points on the surface
    bool inside(const TV& X) const override;

    T signedDistance(const TV& X) const override;

    TV closestPointOnSurface(const TV& X) const override;

    TV normal(const TV& X) const override;

    TM hessian(const TV& X) const override;

    bool queryInside(const TV& X, T& signed_distance, TV& norm, const T phi_critical = 0) const override;

    bool query(const TV& X, T& signed_distance, TV& norm, TM* hess = 0) const override;

    void getBounds(TV& min_corner, TV& max_corner) const override;

    T volume() const override;

    std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override;
};

template <class T, int dim>
class CappedCylinder : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TV2 = Vector<T, 2>;
    using TV3 = Vector<T, 3>;
    using TM = Matrix<T, dim, dim>;

    T radius;
    T height;
    Eigen::Matrix<T, dim, dim, Eigen::DontAlign> R;
    Eigen::Matrix<T, dim, dim, Eigen::DontAlign> R_inverse;
    Eigen::Matrix<T, dim, 1, Eigen::DontAlign> b;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
       The primitive cylinder is centered at origin, along y axis.
       The input transformation (rotation R, translation b)
       applies to the primitive as R*c+b
       Note that q_in is quartetion <w,x,y,z>
     */
    CappedCylinder(const T radius_in, const T height_in, const Vector<T, 4>& q, const TV& b_in)
        : radius(radius_in)
        , height(height_in)
    {
        ZIRAN_ASSERT(dim == 3);
        // the following is basically:
        //R = Eigen::Quaternion<T>(q(0), q(1), q(2), q(3)).normalized().toRotationMatrix().template topLeftCorner<dim, dim>();
        Vector<T, 4> vv = q;
        vv.normalize();
        Eigen::Quaternion<T> qv(vv);
        R = qv.toRotationMatrix().template topLeftCorner<dim, dim>();
        R_inverse = R.inverse();
        b = b_in;
    }

    ~CappedCylinder() override {}

    template <class S>
    S signedDistancePrimitive(const Vector<S, dim>& X) const
    {
        using std::max;
        using std::min;
        using SV2 = Vector<S, 2>;
        SV2 h;
        h << radius, (S)0.5 * height;
        SV2 xz;
        xz << X(0), X(2);
        SV2 qq;
        qq << xz.norm(), X(1);
        SV2 d = qq.array().abs() - h.array();
        qq << max(d(0), (S)0), max(d(1), (S)0);
        S result = min(max(d(0), d(1)), (S)0) + qq.norm();
        return result;
    }

    T signedDistance(const TV& X) const override
    {
        TV X_primitive = R_inverse * (X - b);
        return signedDistancePrimitive(X_primitive);
    }

    TV normal(const TV& X) const override
    {
        TV X_primitive = R_inverse * (X - b);
        ADVec<T, dim> x = vars(X_primitive);
        ADScalar<T, dim> s = signedDistancePrimitive(x);
        TV n = R * s.ds;

        return n;
    }

    void getBounds(TV& min_corner, TV& max_corner) const override
    {
        using namespace MATH_TOOLS;
        T bounding_sphere_radius = std::sqrt(sqr(radius) + sqr((T)0.5 * height));
        min_corner = TV::Constant(-bounding_sphere_radius) + b;
        max_corner = TV::Constant(bounding_sphere_radius) + b;
    }

    T volume() const override
    {
        return M_PI * MATH_TOOLS::sqr(radius) * height;
    }

    std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override
    {
        return std::make_unique<CappedCylinder<T, dim>>(*this);
    }
};

template <class T, int dim>
class AnalyticBox : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TV2 = Vector<T, 2>;
    using TV3 = Vector<T, 3>;
    using TM = Matrix<T, dim, dim>;

    TV half_edges;

    Eigen::Matrix<T, dim, dim, Eigen::DontAlign> R;
    Eigen::Matrix<T, dim, dim, Eigen::DontAlign> R_inverse;
    Eigen::Matrix<T, dim, 1, Eigen::DontAlign> b;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
       The primitive cylinder is centered at origin,
       with min corner -half_edges and max corner half_edges.
       The input transformation (rotation R, translation b)
       applies to the primitive as R*c+b
       Note that q_in is quartetion <w,x,y,z> in 3d, <theta,*,*,*> in 2d.
     */
    AnalyticBox(const TV& half_edges_in, const Vector<T, 4>& q, const TV& b_in);
    /*AnalyticBox(const TV& half_edges_in, const Vector<T, 4>& q, const TV& b_in)
    {
        using std::sin;
        using std::cos;
        half_edges = half_edges_in;
        if (dim == 3) {
	    R = Eigen::Quaternion<T>(q(0), q(1), q(2), q(3)).normalized().toRotationMatrix().template topLeftCorner<dim, dim>();
        }
        else if (dim == 2) {
            T theta = q(0);
            R << cos(theta), -sin(theta), sin(theta), cos(theta);
        }
        R_inverse = R.inverse();
        b = b_in;
    }*/

    ~AnalyticBox() override {}

    template <class S>
    S signedDistancePrimitive(const Vector<S, dim>& X) const;

    T signedDistance(const TV& X) const override;

    TV normal(const TV& X) const override;

    void getBounds(TV& min_corner, TV& max_corner) const override;

    T volume() const override;

    std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override;
};

template <class T, int dim>
class Torus : public AnalyticLevelSet<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TV2 = Vector<T, 2>;
    using TV3 = Vector<T, 3>;
    using TM = Matrix<T, dim, dim>;

    T r0; // major radius (radius of the center of the ring)
    T r1; // minor radius
    // use -0.5 and 1 to make a football

    TM R;
    TM R_inverse;
    TV b;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
       The primitive torus is centered at origin, hole is through y axis.
       The input transformation (rotation R, translation b)
       applies to the primitive as R*c+b
       Note that q_in is quartetion <w,x,y,z>
     */
    Torus(const T r0_in, const T r1_in, const Vector<T, 4>& q, const TV& b_in);

    ~Torus() override {}

    template <class S>
    S signedDistancePrimitive(const Vector<S, dim>& X) const;

    T signedDistance(const TV& X) const override;

    TV normal(const TV& X) const override;

    void getBounds(TV& min_corner, TV& max_corner) const override;

    T volume() const override;

    std::unique_ptr<AnalyticLevelSet<T, dim>> createPtr() const override;
};
} // namespace ZIRAN

#endif
