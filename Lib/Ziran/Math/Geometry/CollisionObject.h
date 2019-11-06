#ifndef COLLISION_OBJECT_H
#define COLLISION_OBJECT_H
#include <Ziran/Math/Geometry/Rotation.h>
#include <Ziran/Math/Geometry/OutputPolyMesh.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <memory>

namespace ZIRAN {

struct OutputPolyMesh;
template <class T, int dim>
class AnalyticLevelSet;

template <class T, int dim>
struct CollisionNode {
    int node_id;
    Matrix<T, dim, dim> P;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
       P projects v to the feasible set.

       Assume there are d basis vectors bi of slip collision normals,
       they are the columns of a dim by dim matrix K.
       (For STICKY collision, K is identity)
       (For less then d normals, the rest columns are zero)
       The operation to project v is:
       v -= \sum_i (bi.dot(v)) bi
       This computation simplifies to
       v = P*v, where P = I - KK'
     */
    template <class TV>
    inline void project(const TV& cv)
    {
        TV& v = const_cast<TV&>(cv);
        v = P * v;
    }
};

template <class T, int dim>
class AnalyticCollisionObject {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    enum COLLISION_OBJECT_TYPE {
        STICKY = 1,
        SLIP = 2,
        SEPARATE = 3,
        GHOST = 4
    };

    std::unique_ptr<AnalyticLevelSet<T, dim>> ls;
    COLLISION_OBJECT_TYPE type;
    T friction;

    Rotation<T, dim> R; // Rotation
    AngularVelocity<T, dim> omega; // angular velocity
    Eigen::Matrix<T, dim, 1, Eigen::DontAlign> b; // translation
    Eigen::Matrix<T, dim, 1, Eigen::DontAlign> dbdt; // translation velocity
    T s; // uniform scaling
    T dsdt; // uniform scaling velocity
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     function for updating transformations and their derivatives (veloties) wrt time
    */
    std::function<void(T, AnalyticCollisionObject<T, dim>&)> updateState;

    static void trialCollision(const StdVector<std::unique_ptr<AnalyticCollisionObject<T, dim>>>& objects,
        const TV& xi, const T dt, TV& vi);

    static bool multiObjectCollision(const StdVector<std::unique_ptr<AnalyticCollisionObject<T, dim>>>& objects,
        const TV& xi, TV& vi, TM& normal_basis);

    static bool multiObjectCollisionWithFakeSeparate(const StdVector<std::unique_ptr<AnalyticCollisionObject<T, dim>>>& objects,
        const TV& xi, TV& vi, TM& normal_basis);

    /**
       A moving collision object supporting x = R S x0 + b, where R is rotation, S is uniform scaling, b is translation
    */
    AnalyticCollisionObject(std::function<void(T, AnalyticCollisionObject<T, dim>&)> updateState, std::unique_ptr<AnalyticLevelSet<T, dim>>&& ls_in, COLLISION_OBJECT_TYPE type_in);

    /**
       if no updateState is passed in, the collision object is static
     */
    AnalyticCollisionObject(std::unique_ptr<AnalyticLevelSet<T, dim>>&& ls_in, COLLISION_OBJECT_TYPE type_in);

    /**
       A moving collision object supporting x = R S x0 + b, where R is rotation, S is uniform scaling, b is translation
    */
    AnalyticCollisionObject(std::function<void(T, AnalyticCollisionObject<T, dim>&)> updateState, AnalyticLevelSet<T, dim>& ls_in, COLLISION_OBJECT_TYPE type_in);

    /**
       if no updateState is passed in, the collision object is static
     */
    AnalyticCollisionObject(AnalyticLevelSet<T, dim>& ls_in, COLLISION_OBJECT_TYPE type_in);

    AnalyticCollisionObject(AnalyticCollisionObject&& other) = default;
    AnalyticCollisionObject(const AnalyticCollisionObject& other) = delete;

    virtual ~AnalyticCollisionObject()
    {
    }

    virtual TV getMaterialVelocity(const TV& X) const
    {
        return TV::Zero();
    }

    virtual T evalMaxSpeed(const TV& p_min_corner, const TV& p_max_corner) const;

    virtual std::function<OutputPolyMesh()> createMesh(T dx) const
    {
        return ls->createMesh(dx);
    }

    /**
      return true if collision object transforms only
      by scaling, rotation and translation
    */
    virtual bool rigid() const
    {
        return true;
    }

    void
    setStatic()
    {
        R.setIdentityRotation();
        omega.setZero();
        b = TV::Zero();
        dbdt = TV::Zero();
        s = 1;
        dsdt = 0;
    }

    void setFriction(const T& friction_in)
    {
        friction = friction_in;
    }

    void setTranslation(const TV& b_in, const TV& dbdt_in)
    {
        b = b_in;
        dbdt = dbdt_in;
    }

    // For 3d, q represents a quaternion <w,x,y,z>
    // For 2d, only the first number of q matters, representing the angle
    void setRotation(const Vector<T, 4>& q)
    {
        R = Rotation<T, dim>(q);
    }

    void setAngularVelocity(const Vector<T, 3>& q)
    {
        omega.set(q);
    }

    void setScaling(const T& s_in, const T& dsdt_in)
    {
        s = s_in;
        dsdt = dsdt_in;
    }

    void getTransform(Matrix<T, 4, 4>& t) const
    {
        t.setIdentity();
        t.template topLeftCorner<dim, dim>() = s * R.rotation;
        t.template topRightCorner<dim, 1>() = b;
    }

    // find the world space position and velocity for any material space point based on the transformation
    void transformedPoint(const TV& X, TV& x, TV& v) const;

    void transformedVector(const TV& N, TV& n) const;

    void untransformedPoint(const TV& x, TV& X) const;

    T signedDistance(const TV& x) const;

    bool queryInside(const TV& x, T& phi, TV& n, const T phi_critical = 0) const;

    void particleCollision(const T phi_critical, const TV& x, TV& v) const;

    bool getBoundaryVelocityAndNormal(const TV& x, TV& v, TV&) const;

    /**
       This takes a position and its velocity,
       project the velocity, and returns a normal if the collision happened as a SLIP collsion.
    */
    bool detectAndResolveCollision(const TV& x, TV& v, TV& n) const;
};
} // namespace ZIRAN
#endif
