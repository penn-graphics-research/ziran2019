#include <Ziran/Math/Geometry/CollisionObject.h>
#include <Ziran/CS/Util/Forward.h>
#include <openvdb/openvdb.h>
#undef B2
namespace ZIRAN {

template <class T, int dim>
class SourceCollisionObject : public AnalyticCollisionObject<T, dim> {
public:
    using Base = AnalyticCollisionObject<T, dim>;
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;

    using Base::friction;
    using Base::ls;
    using Base::type;

    using Base::b; //translation
    using Base::dbdt; // translation velocity
    using Base::dsdt; // uniform scaling velocity
    using Base::omega; // angular velocity
    using Base::R; // Rotation
    using Base::s; // uniform scaling
    using Base::setStatic;
    using Base::updateState;

    std::function<void(const TV&, TV&)> materialVelocity;
    TV const_material_velocity;
    StdVector<TV> samples;
    StdVector<TV> add_samples;
    StdVector<int> inside;

    T uniform_mass;

    /**
       A moving collision object supporting x = R S x0 + b, where R is rotation, S is uniform scaling, b is translation
    */
    SourceCollisionObject(std::function<void(T, AnalyticCollisionObject<T, dim>&)> updateState, std::unique_ptr<AnalyticLevelSet<T, dim>>&& ls_in, const TV& const_vel)
        : AnalyticCollisionObject<T, dim>(updateState, std::move(ls_in), Base::COLLISION_OBJECT_TYPE::STICKY)
    {
        setStatic();
        materialVelocity = [const_vel](const TV& x, TV& v) { v = const_vel; };
        const_material_velocity = const_vel;
    }

    /**
       if no updateState is passed in, the collision object is static
     */
    SourceCollisionObject(std::unique_ptr<AnalyticLevelSet<T, dim>>&& ls_in, const TV& const_vel)
        : AnalyticCollisionObject<T, dim>(std::move(ls_in), Base::COLLISION_OBJECT_TYPE::STICKY)
    {
        setStatic();
        materialVelocity = [const_vel](const TV& x, TV& v) { v = const_vel; };
        const_material_velocity = const_vel;
    }

    /**
       A moving collision object supporting x = R S x0 + b, where R is rotation, S is uniform scaling, b is translation
    */
    SourceCollisionObject(std::function<void(T, AnalyticCollisionObject<T, dim>&)> updateState, AnalyticLevelSet<T, dim>& ls_in, const TV& const_vel)
        : AnalyticCollisionObject<T, dim>(updateState, ls_in, Base::COLLISION_OBJECT_TYPE::STICKY)
    {
        setStatic();
        materialVelocity = [const_vel](const TV& x, TV& v) { v = const_vel; };
        const_material_velocity = const_vel;
    }

    /**
       if no updateState is passed in, the collision object is static
     */
    SourceCollisionObject(AnalyticLevelSet<T, dim>& ls_in, const TV& const_vel)
        : AnalyticCollisionObject<T, dim>(ls_in, Base::COLLISION_OBJECT_TYPE::STICKY)
    {
        setStatic();
        materialVelocity = [const_vel](const TV& x, TV& v) { v = const_vel; };
        const_material_velocity = const_vel;
    }

    void setMaterialVelocity(std::function<void(const TV&, TV&)> new_material_velocity)
    {
        materialVelocity = new_material_velocity;
    }

    void setMaterialVelocity(const TV& new_material_velocity)
    {
        materialVelocity = [new_material_velocity](const TV& x, TV& v) { v = new_material_velocity; };
        const_material_velocity = new_material_velocity;
    }

    SourceCollisionObject(SourceCollisionObject&& other) = default;
    SourceCollisionObject(const SourceCollisionObject& other) = delete;

    virtual TV getMaterialVelocity(const TV& X) const override
    {
        TV v;
        materialVelocity(X, v);
        return v;
    }

    virtual T evalMaxSpeed(const TV& p_min_corner, const TV& p_max_corner) const override;

    void uniformSample(T grid_dx, int total_particles_number);

    void poissonSample(T grid_dx, T particle_per_cell, T density);

    void sampleAndPrune(T dt, T grid_dx);

    void addParticlesFromSample(T dt, T grid_dx);

    void addParticlesFromSource(T density, T total_volume, T particles_per_cell);

    using Base::detectAndResolveCollision;
    using Base::getTransform;
    using Base::multiObjectCollision;
    using Base::particleCollision;
    using Base::queryInside;
    using Base::setAngularVelocity;
    using Base::setFriction;
    using Base::setRotation;
    using Base::setScaling;
    using Base::setTranslation;
    using Base::signedDistance;
    using Base::trialCollision;
};
} // namespace ZIRAN
