#ifndef ROTATION_H
#define ROTATION_H
#include <Ziran/CS/Util/Forward.h>
#include <Eigen/Geometry>
#include <tick/requires.h>
#include <type_traits>

namespace ZIRAN {

// template <class T, int dim>
// struct RotationTypeHelper;

// template <class T>
// struct RotationTypeHelper<T, 2> {
//     using type = Eigen::Rotation2D<T>;
// };

// template <class T>
// struct RotationTypeHelper<T, 3> {
//     using type = Eigen::Quaternion<T>;
// };

// template <class T, int dim>
// using RotationType = typename RotationTypeHelper<T, dim>::type;

/**
   This is a ZIRAN rotation class that wraps Eigen's Rotation2D (in 2D) and Quaternion (in 3D).
 */
template <class T, int dim>
class Rotation {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;

    // RotationType<T, dim> rotation;
    Eigen::Matrix<T, dim, dim, Eigen::DontAlign> rotation;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Rotation()
    {
        setIdentityRotation();
    }

    TICK_MEMBER_REQUIRES(dim == 2)
    Rotation(const Vector<T, 4>& q)
    {
        Eigen::Rotation2D<T> r(q(0)); // set from angle that is stored in q(0)
        rotation = r.toRotationMatrix();
    }

    TICK_MEMBER_REQUIRES(dim == 3)
    Rotation(const Vector<T, 4>& q)
    {
        Eigen::Quaternion<T> r(q(0), q(1), q(2), q(3)); // q is quaternion <w,x,y,z>
        rotation = r.toRotationMatrix();
    }

    // The resulting rotation rotates a to b. a or b don't need to be normalized and can have different length (not zero through).
    TICK_MEMBER_REQUIRES(dim == 3)
    Rotation(const TV& a, const TV& b)
    {
        Eigen::Quaternion<T> r = Eigen::Quaternion<T>::FromTwoVectors(a, b);
        rotation = r.toRotationMatrix();
    }

    // The resulting rotation rotates a to b. a or b don't need to be normalized and can have different length (not zero through).
    TICK_MEMBER_REQUIRES(dim == 2)
    Rotation(const TV& a, const TV& b)
    {
        TV aa = a.normalized();
        TV bb = b.normalized();
        TM R;
        R(0, 0) = aa(0) * bb(0) + aa(1) * bb(1);
        R(0, 1) = -(aa(0) * bb(1) - bb(0) * aa(1));
        R(1, 0) = aa(0) * bb(1) - bb(0) * aa(1);
        R(1, 1) = aa(0) * bb(0) + aa(1) * bb(1);
        rotation = R;
    }

    ~Rotation() {}

    void setIdentityRotation()
    {
        rotation = TM::Identity();
    }
};

template <class T, int dim, class enable = void>
class AngularVelocity;

/**
   2D angular velocity is a number
 */
template <class T, int dim>
class AngularVelocity<T, dim, std::enable_if_t<dim == 2>> {
public:
    typedef Vector<T, dim> TV;

    T angular_velocity;

    AngularVelocity() { setZero(); }
    AngularVelocity(T w) { set(w); }
    ~AngularVelocity() {}

    void setZero() { angular_velocity = (T)0; }
    void set(T w) { angular_velocity = w; }
    void set(const Vector<T, 3>& w) { angular_velocity = w(0); }

    void operator=(const AngularVelocity& another) { set(another.angular_velocity); }
    void operator+=(const AngularVelocity& another) { angular_velocity += another.angular_velocity; }
    void operator-=(const AngularVelocity& another) { angular_velocity -= another.angular_velocity; }
    void operator*=(T alpha) { angular_velocity *= alpha; }
    void operator/=(T alpha) { angular_velocity /= alpha; }

    AngularVelocity operator+(const AngularVelocity& another) const { return AngularVelocity(angular_velocity + another.angular_velocity); }
    AngularVelocity operator*(T alpha) const { return AngularVelocity(angular_velocity * alpha); }

    Vector<T, dim> cross(const Vector<T, dim>& x) const
    {
        Vector<T, dim> result;
        result << -angular_velocity * x(1), angular_velocity * x(0);
        return result;
    }

    T norm() const
    {
        using std::abs;
        return abs(angular_velocity);
    }

    template <class Type>
    AngularVelocity<Type, dim> cast()
    {
        return AngularVelocity<Type, dim>((Type)angular_velocity);
    }
};

/**
   3D angular velocity is a 3D vector
 */
template <class T, int dim>
class AngularVelocity<T, dim, std::enable_if_t<dim == 3>> {
public:
    typedef Vector<T, dim> TV;

    TV angular_velocity;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AngularVelocity() { setZero(); }
    AngularVelocity(const Vector<T, dim>& w) { set(w); }
    ~AngularVelocity() {}

    void setZero() { angular_velocity = TV::Zero(); }
    void set(const Vector<T, dim>& w) { angular_velocity = w; }

    void operator=(const AngularVelocity& another) { set(another.angular_velocity); }
    void operator+=(const AngularVelocity& another) { angular_velocity += another.angular_velocity; }
    void operator-=(const AngularVelocity& another) { angular_velocity -= another.angular_velocity; }
    void operator*=(T alpha) { angular_velocity *= alpha; }
    void operator/=(T alpha) { angular_velocity /= alpha; }

    AngularVelocity operator+(const AngularVelocity& another) const { return AngularVelocity(angular_velocity + another.angular_velocity); }
    AngularVelocity operator*(T alpha) const { return AngularVelocity(angular_velocity * alpha); }

    Vector<T, dim> cross(const Vector<T, dim>& x) const { return angular_velocity.cross(x); }

    T norm() const
    {
        return angular_velocity.norm();
    }

    template <class Type>
    AngularVelocity<Type, dim> cast()
    {
        return AngularVelocity<Type, dim>(angular_velocity.template cast<Type>());
    }
};

// B : epsilon
template <class T>
inline AngularVelocity<T, 2> LeviCivitaContract(const Matrix<T, 2, 2>& B)
{
    return AngularVelocity<T, 2>(B(0, 1) - B(1, 0));
}

// B : epsilon
template <class T>
inline AngularVelocity<T, 3> LeviCivitaContract(const Matrix<T, 3, 3>& B)
{
    Vector<T, 3> result;
    result(0) = B(1, 2) - B(2, 1);
    result(1) = B(2, 0) - B(0, 2);
    result(2) = B(0, 1) - B(1, 0);
    return AngularVelocity<T, 3>(result);
}
} // namespace ZIRAN
#endif
