#include <Ziran/Math/Geometry/CollisionObject.h>
#include <Ziran/Math/Geometry/AnalyticLevelSet.h>

#include <memory>

namespace ZIRAN {

template <class T, int dim>
AnalyticCollisionObject<T, dim>::AnalyticCollisionObject(std::function<void(T, AnalyticCollisionObject<T, dim>&)> updateState, std::unique_ptr<AnalyticLevelSet<T, dim>>&& ls_in, COLLISION_OBJECT_TYPE type_in)
    : ls(std::move(ls_in))
    , type(type_in)
    , friction(0)
    , updateState(updateState)
{
    setStatic();
}

template <class T, int dim>
AnalyticCollisionObject<T, dim>::AnalyticCollisionObject(std::unique_ptr<AnalyticLevelSet<T, dim>>&& ls_in, COLLISION_OBJECT_TYPE type_in)
    : ls(std::move(ls_in))
    , type(type_in)
    , friction(0)
    , updateState(nullptr)
{
    setStatic();
}

template <class T, int dim>
AnalyticCollisionObject<T, dim>::AnalyticCollisionObject(std::function<void(T, AnalyticCollisionObject<T, dim>&)> updateState, AnalyticLevelSet<T, dim>& ls_in, COLLISION_OBJECT_TYPE type_in)
    : ls(ls_in.createPtr())
    , type(type_in)
    , friction(0)
    , updateState(updateState)
{
    setStatic();
}

template <class T, int dim>
AnalyticCollisionObject<T, dim>::AnalyticCollisionObject(AnalyticLevelSet<T, dim>& ls_in, COLLISION_OBJECT_TYPE type_in)
    : ls(ls_in.createPtr())
    , type(type_in)
    , friction(0)
    , updateState(nullptr)
{
    setStatic();
}

template <class T, int dim>
void AnalyticCollisionObject<T, dim>::trialCollision(const StdVector<std::unique_ptr<AnalyticCollisionObject<T, dim>>>& objects,
    const TV& xi, const T dt, TV& vi)
{
    for (size_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->type == GHOST)
            continue;

        assert(objects[k]->ls != nullptr);
        T phi;
        TV n;
        TV x_trial = xi + dt * vi;
        bool colliding = objects[k]->queryInside(x_trial, phi, n);
        if (colliding)
            vi -= (phi / dt) * n;
    }
}

template <class T, int dim>
bool AnalyticCollisionObject<T, dim>::multiObjectCollision(const StdVector<std::unique_ptr<AnalyticCollisionObject<T, dim>>>& objects,
    const TV& xi, TV& vi, TM& normal_basis)
{
    bool any_collision = false;
    int slip_collision_count = 0;
    normal_basis = TM::Zero();
    for (size_t k = 0; k < objects.size(); ++k) {
        assert(objects[k]->ls != nullptr);
        if (objects[k]->type == GHOST)
            continue;

        TV n;
        bool collide = objects[k]->detectAndResolveCollision(xi, vi, n);
        any_collision = (collide || any_collision);
        if (collide) {
            if (objects[k]->type == STICKY) {
                normal_basis = TM::Identity();
                break;
            }
            else { // SLIP
                /*
                   Gram-Schmit the normals
                   */
                for (int c = 0; c < slip_collision_count; c++) {
                    TV n_old = normal_basis.col(c);
                    T dot = n_old.dot(n);
                    n -= dot * n_old;
                }
                T length = n.norm();
                if (length) {
                    normal_basis.col(slip_collision_count) = n / length;
                    if (++slip_collision_count == dim)
                        break;
                }
            }
        }
    }
    return any_collision;
}

template <class T, int dim>
bool AnalyticCollisionObject<T, dim>::multiObjectCollisionWithFakeSeparate(const StdVector<std::unique_ptr<AnalyticCollisionObject<T, dim>>>& objects,
    const TV& xi, TV& vi, TM& normal_basis)
{
    bool any_collision = false;
    int slip_collision_count = 0;
    normal_basis = TM::Zero();
    for (size_t k = 0; k < objects.size(); ++k) {
        assert(objects[k]->ls != nullptr);
        if (objects[k]->type == GHOST)
            continue;

        TV n;
        TV old_vi = vi;
        bool collide = objects[k]->detectAndResolveCollision(xi, vi, n);
        if (collide && objects[k]->type == SEPARATE) {
            TV boundary_velocity;
            TV boundary_normal;
            objects[k]->getBoundaryVelocityAndNormal(xi, boundary_velocity, boundary_normal);
            if ((old_vi - boundary_velocity).dot(boundary_normal) > 0) {
                vi = old_vi;
                collide = false;
            }
        }
        any_collision = (collide || any_collision);
        if (collide) {
            if (objects[k]->type == STICKY) {
                normal_basis = TM::Identity();
                break;
            }
            else { // SLIP
                /*
                   Gram-Schmit the normals
                   */
                for (int c = 0; c < slip_collision_count; c++) {
                    TV n_old = normal_basis.col(c);
                    T dot = n_old.dot(n);
                    n -= dot * n_old;
                }
                T length = n.norm();
                if (length) {
                    normal_basis.col(slip_collision_count) = n / length;
                    if (++slip_collision_count == dim)
                        break;
                }
            }
        }
    }
    return any_collision;
}

template <class T, int dim>
T AnalyticCollisionObject<T, dim>::evalMaxSpeed(const TV& p_min_corner, const TV& p_max_corner) const
{
    if (dsdt != 0 || omega.norm() != 0) {
        // Find intersection of bounding boxes
        // Max speed occurs at corner
        TV min_corner, max_corner;
        ls->getBounds(min_corner, max_corner);
        StdVector<TV> intersection_corners;
        assert(s);
        T one_over_s = 1 / s;
        for (int i = 0; i < 1 << dim; i++) {
            TV x;
            for (int d = 0; d < dim; d++)
                x(d) = (i & 1 << d) ? p_min_corner(d) : p_max_corner(d);
            TV X = R.rotation.transpose() * (x - b) * one_over_s;
            if ((min_corner.array() < X.array()).any() && (X.array() < max_corner.array()).any())
                intersection_corners.emplace_back(x);
        }

        for (int i = 0; i < 1 << dim; i++) {
            TV X;
            for (int d = 0; d < dim; d++)
                X(d) = (i & 1 << d) ? min_corner(d) : max_corner(d);
            TV x = R.rotation * s * X + b;
            if ((p_min_corner.array() < x.array()).any() && (x.array() < p_max_corner.array()).any())
                intersection_corners.emplace_back(x);
        }
        T max_speed = 0;
        for (int i = 0; i < (int)intersection_corners.size(); i++) {
            TV x = intersection_corners[i];
            TV x_minus_b = x - b;
            TV v = omega.cross(x_minus_b) + (dsdt * one_over_s) * x_minus_b + dbdt;
            max_speed = std::max(max_speed, v.norm());
        }
        return max_speed;
    }
    return dbdt.norm();
}

template <class T, int dim>
void AnalyticCollisionObject<T, dim>::transformedPoint(const TV& X, TV& x, TV& v) const
{
    TV x_minus_b = R.rotation * s * X;
    x = x_minus_b + b;
    v = omega.cross(x_minus_b) + (dsdt / s) * x_minus_b + R.rotation * s * getMaterialVelocity(X) + dbdt;
}

template <class T, int dim>
void AnalyticCollisionObject<T, dim>::transformedVector(const TV& N, TV& n) const
{
    n = R.rotation * s * N;
}

template <class T, int dim>
void AnalyticCollisionObject<T, dim>::untransformedPoint(const TV& x, TV& X) const
{
    TV x_minus_b = x - b;
    assert(s);
    T one_over_s = 1 / s;
    X = R.rotation.transpose() * x_minus_b * one_over_s;
}

template <class T, int dim>
T AnalyticCollisionObject<T, dim>::signedDistance(const TV& x) const
{
    TV x_minus_b = x - b;
    assert(s);
    T one_over_s = 1 / s;
    TV X = R.rotation.transpose() * x_minus_b * one_over_s; // material space
    return ls->signedDistance(X) * s;
}

template <class T, int dim>
bool AnalyticCollisionObject<T, dim>::queryInside(const TV& x, T& phi, TV& n, const T phi_critical) const
{
    TV x_minus_b = x - b;
    assert(s);
    T one_over_s = 1 / s;
    TV X = R.rotation.transpose() * x_minus_b * one_over_s; // material space
    TV N;
    bool inside = ls->queryInside(X, phi, N, phi_critical);
    if (inside) {
        n = R.rotation * N; // we care about world space normal, otherwise it is nan
        phi *= s;
        return true;
    }
    return false;
}

template <class T, int dim>
void AnalyticCollisionObject<T, dim>::particleCollision(const T phi_critical, const TV& x, TV& v) const
{
    if (type == GHOST)
        return;

    TV n;
    n.setConstant(NAN);

    TV x_minus_b = x - b;
    assert(s);
    T one_over_s = 1 / s;
    TV X = R.rotation.transpose() * x_minus_b * one_over_s; // material space
    T phi_dummy;
    TV N; // material space normal
    bool colliding = ls->queryInside(X, phi_dummy, N, phi_critical);

    if (colliding) {
        TV v_object = omega.cross(x_minus_b) + (dsdt * one_over_s) * x_minus_b + R.rotation * s * getMaterialVelocity(X) + dbdt;
        v -= v_object;
        if (type == STICKY) {
            v.setZero();
        }
        else if (type == SLIP) {
            n = R.rotation * N; // we care about world space normal, otherwise it is nan
            T dot = v.dot(n);
            TV v_normal = n * dot;
            v -= v_normal;

            // kinematics friction
            if (friction != 0) {
                if (dot < 0) {
                    if (-dot * friction < v.norm())
                        v += v.normalized() * dot * friction;
                    else
                        v.setZero();
                }
            }

            // v += v_normal;
        }
        else if (type == SEPARATE) {
            n = R.rotation * N; // we care about world space normal, otherwise it is nan
            T dot = v.dot(n);
            if (dot < 0) {

                TV v_normal = n * dot;
                v -= v_normal;

                // kinematics friction
                if (friction != 0) {
                    if (-dot * friction < v.norm())
                        v += v.normalized() * dot * friction;
                    else
                        v.setZero();
                }

                // v += v_normal;
            }
        }
        v += v_object;
    }
}

template <class T, int dim>
bool AnalyticCollisionObject<T, dim>::getBoundaryVelocityAndNormal(const TV& x, TV& v, TV& n) const
{
    if (type == GHOST)
        return false;

    n.setConstant(NAN);

    TV x_minus_b = x - b;
    assert(s);
    T one_over_s = 1 / s;
    TV X = R.rotation.transpose() * x_minus_b * one_over_s; // material space
    T phi_dummy;
    TV N; // material space normal
    bool colliding = ls->queryInside(X, phi_dummy, N);

    if (colliding) {
        v = omega.cross(x_minus_b) + (dsdt * one_over_s) * x_minus_b + R.rotation * s * getMaterialVelocity(X) + dbdt;
        n = R.rotation * N;
        return true;
    }
    return false;
}

/**
   This takes a position and its velocity,
   project the velocity, and returns a normal if the collision happened as a SLIP collsion.
*/
template <class T, int dim>
bool AnalyticCollisionObject<T, dim>::detectAndResolveCollision(const TV& x, TV& v, TV& n) const
{
    if (type == GHOST)
        return false;

    /** derivation:
     
      x = \phi(X,t) = R(t)s(t)X+b(t)
      X = \phi^{-1}(x,t) = (1/s) R^{-1} (x-b)
      V(X,t) = \frac{\partial \phi}{\partial t}
             = R'sX + Rs'X + RsX' + b'
      v(x,t) = V(\phi^{-1}(x,t),t)
             = R'R^{-1}(x-b) + (s'/s)(x-b) + RsX' + b'
             = omega \cross (x-b) + (s'/s)(x-b) +b'
    */
    n.setConstant(NAN);

    TV x_minus_b = x - b;
    assert(s);
    T one_over_s = 1 / s;
    TV X = R.rotation.transpose() * x_minus_b * one_over_s; // material space
    T phi_dummy;
    TV N; // material space normal
    bool colliding = ls->queryInside(X, phi_dummy, N);

    if (colliding) {
        TV v_object = omega.cross(x_minus_b) + (dsdt * one_over_s) * x_minus_b + R.rotation * s * getMaterialVelocity(X) + dbdt;
        v -= v_object;
        if (type == STICKY) {
            v.setZero();
        }
        else if (type == SLIP) {
            n = R.rotation * N; // we care about world space normal, otherwise it is nan
            T dot = v.dot(n);
            v -= n * dot;

            // kinematics friction
            if (friction != 0) {
                if (dot < 0) {
                    if (-dot * friction < v.norm())
                        v += v.normalized() * dot * friction;
                    else
                        v.setZero();
                }
            }
        }
        else if (type == SEPARATE) {
            n = R.rotation * N; // we care about world space normal, otherwise it is nan
            T dot = v.dot(n);
            if (dot < 0) {
                v -= n * dot;

                // kinematics friction
                if (friction != 0) {
                    if (-dot * friction < v.norm())
                        v += v.normalized() * dot * friction;
                    else
                        v.setZero();
                }
            }
        }
        v += v_object;
    }
    return colliding;
}

template class AnalyticCollisionObject<double, 2>;
template class AnalyticCollisionObject<double, 3>;
template class AnalyticCollisionObject<float, 2>;
template class AnalyticCollisionObject<float, 3>;
} // namespace ZIRAN
