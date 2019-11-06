#include "AnalyticLevelSet.h"
#include <Ziran/Math/Linear/GivensQR.h>
#include <Ziran/Math/Geometry/OutputPolyMesh.h>
#ifdef ZIRAN_WITH_VDB
#include <openvdb/tools/SignedFloodFill.h>
#include <openvdb/tools/VolumeToMesh.h>
#undef B2
#endif

namespace ZIRAN {

template <class T, int dim>
std::function<OutputPolyMesh()> AnalyticLevelSet<T, dim>::createMesh(T dx) const
{
#ifdef ZIRAN_WITH_VDB
    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    using GridT = typename openvdb::Grid<typename openvdb::tree::Tree4<float, 5, 4, 3>::Type>;
    using IV = Vector<int, 3>;
    using std::abs;
    typename GridT::Ptr grid;
    T max_distance = 3 * dx;
    TV min_corner = TV::Constant(-max_distance);
    TV max_corner = TV::Constant(max_distance);
    Vector<T, dim> min;
    Vector<T, dim> max;
    getBounds(min, max);
    min_corner.template head<dim>() = min;
    max_corner.template head<dim>() = max;

    grid = GridT::create(max_distance);
    grid->setTransform(openvdb::math::Transform::createLinearTransform(dx));

    typename GridT::Accessor acc = grid->getAccessor();
    IV min_index;
    IV max_index;
    for (size_t d = 0; d < 3; ++d) {
        min_index(d) = std::floor(min_corner(d) / dx) - 4;
        max_index(d) = std::ceil(max_corner(d) / dx) + 4;
    }
    for (int i = min_index(0); i <= max_index(0); i++)
        for (int j = min_index(1); j <= max_index(1); j++)
            for (int k = min_index(2); k <= max_index(2); k++) {
                openvdb::Coord ijk(i, j, k);
                Vector<T, dim> x;
                for (int d = 0; d < dim; d++)
                    x(d) = ijk[d] * dx;
                T distance = signedDistance(x);
                if (dim == 2) {
                    T distance_to_wall = (abs(ijk[2]) - (T)0.5) * dx;
                    distance = std::max(distance_to_wall, distance);
                }
                if (abs(distance) < max_distance) {
                    acc.setValue(ijk, distance);
                }
            }
    openvdb::tools::signedFloodFill(grid->tree());
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
#else
    ZIRAN_ASSERT(false, "Analytic level set create mesh not implemented");
#endif
}

template <class T, int dim>
bool AnalyticLevelSet<T, dim>::inside(const TV& X) const
{
    return signedDistance(X) <= (T)0;
}

template <class T, int dim>
typename AnalyticLevelSet<T, dim>::TV AnalyticLevelSet<T, dim>::closestPointOnSurface(const TV& X) const
{
    TV normal;
    T signed_distance;
    query(X, signed_distance, normal);
    return X - signed_distance * normal;
}

/**
Hessian of signed distance function at X.
*/
template <class T, int dim>
typename AnalyticLevelSet<T, dim>::TM AnalyticLevelSet<T, dim>::hessian(const TV& X) const
{
    ZIRAN_ASSERT(false, "Analytic level set hessian not implemented.");
}

template <class T, int dim>
bool AnalyticLevelSet<T, dim>::query(const TV& X, T& signed_distance, TV& normal, TM* hessian) const
{
    signed_distance = this->signedDistance(X);
    normal = this->normal(X);
    if (hessian)
        *hessian = this->hessian(X);
    return signed_distance <= (T)0;
}

template <class T, int dim>
bool AnalyticLevelSet<T, dim>::queryInside(const TV& X, T& signed_distance, TV& normal, const T phi_critical) const
{
    signed_distance = this->signedDistance(X);
    if (signed_distance <= phi_critical) {
        normal = this->normal(X);
        return true;
    }
    return false;
}

template <class T, int dim>
void AnalyticLevelSet<T, dim>::getBounds(TV& min_corner, TV& max_corner) const
{
    ZIRAN_ASSERT(false, "No bounds available.");
}

template <class T, int dim>
T AnalyticLevelSet<T, dim>::volume() const
{
    ZIRAN_ASSERT(false, "No volume available.");
}

// DisjointUnioniLevelSet
template <class T, int dim>
DisjointUnionLevelSet<T, dim>::DisjointUnionLevelSet(const StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>>& level_sets_input)
{
    for (const auto& ls : level_sets_input)
        level_sets.emplace_back(ls->createPtr());
}

template <class T, int dim>
void DisjointUnionLevelSet<T, dim>::add(const AnalyticLevelSet<T, dim>& ls)
{
    level_sets.emplace_back(ls.createPtr());
}

template <class T, int dim>
T DisjointUnionLevelSet<T, dim>::signedDistance(const TV& X) const
{
    T min = std::numeric_limits<T>::max();
    for (const auto& ls : level_sets)
        min = std::min(min, ls->signedDistance(X));
    return min;
}

template <class T, int dim>
void DisjointUnionLevelSet<T, dim>::getBounds(TV& min_corner, TV& max_corner) const
{
    min_corner.setConstant(std::numeric_limits<T>::max());
    max_corner.setConstant(std::numeric_limits<T>::min());
    for (const auto& tls : level_sets) {
        TV local_min, local_max;
        tls->getBounds(local_min, local_max);
        min_corner = min_corner.array().min(local_min.array());
        max_corner = max_corner.array().max(local_max.array());
    }
}
template <class T, int dim>
T DisjointUnionLevelSet<T, dim>::volume() const
{
    T result = 0;
    for (const auto& tls : level_sets)
        result += tls->volume();
    return result;
}

template <class T, int dim>
typename DisjointUnionLevelSet<T, dim>::TV DisjointUnionLevelSet<T, dim>::normal(const TV& X) const
{
    T min = std::numeric_limits<T>::max();
    TV n = TV::Zero();
    for (const auto& ls : level_sets) {
        T distance = ls->signedDistance(X);
        if (distance < min) {
            min = distance;
            n = ls->normal(X);
        }
    }
    return n;
}

template <class T, int dim>
std::unique_ptr<AnalyticLevelSet<T, dim>> DisjointUnionLevelSet<T, dim>::createPtr() const
{
    return std::make_unique<DisjointUnionLevelSet<T, dim>>(level_sets);
}

// DifferenceLevelSet
template <class T, int dim>
DifferenceLevelSet<T, dim>::DifferenceLevelSet(StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>>&& level_sets)
    : level_sets(std::move(level_sets))
{
}

template <class T, int dim>
DifferenceLevelSet<T, dim>::DifferenceLevelSet(const StdVector<std::unique_ptr<AnalyticLevelSet<T, dim>>>& level_sets_input)
{
    for (const auto& ls : level_sets_input)
        level_sets.emplace_back(ls->createPtr());
}

template <class T, int dim>
void DifferenceLevelSet<T, dim>::add(const AnalyticLevelSet<T, dim>& ls1, const AnalyticLevelSet<T, dim>& ls2)
{
    level_sets.emplace_back(ls1.createPtr());
    level_sets.emplace_back(ls2.createPtr());
}

template <class T, int dim>
T DifferenceLevelSet<T, dim>::signedDistance(const TV& X) const
{
    auto& setA = level_sets[0];
    auto& setB = level_sets[1];
    return std::max(setA->signedDistance(X), -(setB->signedDistance(X)));
}

template <class T, int dim>
typename DifferenceLevelSet<T, dim>::TV DifferenceLevelSet<T, dim>::normal(const TV& X) const
{
    auto& setA = level_sets[0];
    auto& setB = level_sets[1];
    if (-(setB->signedDistance(X)) > setA->signedDistance(X))
        return -(setB->normal(X));
    else
        return setA->normal(X);
}

template <class T, int dim>
void DifferenceLevelSet<T, dim>::getBounds(TV& min_corner, TV& max_corner) const
{
    level_sets[0]->getBounds(min_corner, max_corner);
}

template <class T, int dim>
T DifferenceLevelSet<T, dim>::volume() const
{
    return level_sets[0]->volume() - level_sets[1]->volume();
}

template <class T, int dim>
std::unique_ptr<AnalyticLevelSet<T, dim>> DifferenceLevelSet<T, dim>::createPtr() const
{
    return std::make_unique<DifferenceLevelSet<T, dim>>(level_sets);
}

// HalfSpace

template <class T, int dim>
HalfSpace<T, dim>::HalfSpace(const TV& origin_in, const TV& outward_normal_in)
    : origin(origin_in)
    , outward_normal(outward_normal_in.normalized())
{
}

template <class T, int dim>
bool HalfSpace<T, dim>::inside(const TV& X) const
{
    return (signedDistance(X) <= (T)0);
}

template <class T, int dim>
T HalfSpace<T, dim>::signedDistance(const TV& X) const
{
    return outward_normal.dot(X - origin);
}

template <class T, int dim>
typename HalfSpace<T, dim>::TV HalfSpace<T, dim>::closestPointOnSurface(const TV& X) const
{
    return X - outward_normal.dot(X - origin) * outward_normal;
}

template <class T, int dim>
typename HalfSpace<T, dim>::TV HalfSpace<T, dim>::normal(const TV& X) const
{
    return outward_normal;
}

template <class T, int dim>
typename HalfSpace<T, dim>::TM HalfSpace<T, dim>::hessian(const TV& X) const
{
    return TM::Zero();
}

template <class T, int dim>
bool HalfSpace<T, dim>::query(const TV& X, T& signed_distance, TV& normal, TM* hessian) const
{
    signed_distance = signedDistance(X);
    normal = outward_normal;
    if (hessian) {
        hessian->setZero();
    }
    return (signed_distance <= (T)0);
}

template <class T, int dim>
std::unique_ptr<AnalyticLevelSet<T, dim>> HalfSpace<T, dim>::createPtr() const
{
    return std::make_unique<HalfSpace<T, dim>>(*this);
}

template <class T, int dim>
std::function<OutputPolyMesh()> HalfSpace<T, dim>::createMesh(T dx) const
{
    StdVector<Eigen::Vector3f> points;
    StdVector<Eigen::Vector4i> quads;
    Matrix<T, 3, 3> Q;
    Matrix<T, 3, 1> R;
    Vector<T, 3> n = Vector<T, 3>::Zero();
    n.template head<dim>() = outward_normal;
    GivensQR(n, Q, R);
    for (int i = -10; i <= 10; i++)
        for (int j = -10; j <= 10; j++) {
            Vector<T, 3> x = Q * Vector<T, 3>(0, dx * i, dx * j);
            x.template head<dim>() += origin;
            points.emplace_back(x(0), x(1), x(2));
            if (i < 10 && j < 10) {
                int q0 = 21 * (j + 10) + (i + 10);
                int q1 = q0 + 1;
                int q2 = q1 + 21;
                int q3 = q0 + 21;
                quads.emplace_back(q0, q1, q2, q3);
            }
        }
    std::vector<int32_t> face_vertex_count(quads.size(), 4);
    return [points(std::move(points)),
               quads(std::move(quads)),
               face_vertex_count(std::move(face_vertex_count))]() {
        return OutputPolyMesh{
            (float*)points.data(), points.size(),
            (int32_t*)quads.data(),
            quads.size() * 4,
            face_vertex_count.data(),
            face_vertex_count.size()
        };
    };
}
//template <class T, int dim>
//class AnalyticBox;

// AxisAlignedAnalyticBox
template <class T, int dim>
AxisAlignedAnalyticBox<T, dim>::AxisAlignedAnalyticBox(const TV& min_corner, const TV& max_corner)
    : box((max_corner - min_corner) / (T)2, Vector<T, 4>(dim == 3, 0, 0, 0), (min_corner + max_corner) / (T)2)
{
}

template <class T, int dim>
T AxisAlignedAnalyticBox<T, dim>::signedDistance(const TV& X) const
{
    return box.signedDistance(X);
}

template <class T, int dim>
typename AxisAlignedAnalyticBox<T, dim>::TV AxisAlignedAnalyticBox<T, dim>::normal(const TV& X) const
{
    return box.normal(X);
}

template <class T, int dim>
void AxisAlignedAnalyticBox<T, dim>::getBounds(TV& min_corner, TV& max_corner) const
{
    box.getBounds(min_corner, max_corner);
}

template <class T, int dim>
T AxisAlignedAnalyticBox<T, dim>::volume() const
{
    return box.volume();
}

template <class T, int dim>
std::unique_ptr<AnalyticLevelSet<T, dim>> AxisAlignedAnalyticBox<T, dim>::createPtr() const
{
    return std::make_unique<AxisAlignedAnalyticBox<T, dim>>(*this);
}

// Sphere
template <class T, int dim>
Sphere<T, dim>::Sphere(const TV& center, const T& radius)
    : center(center)
    , radius(radius)
{
}

template <class T, int dim>
bool Sphere<T, dim>::inside(const TV& X) const
{
    return (X - center).squaredNorm() <= radius * radius;
}

template <class T, int dim>
T Sphere<T, dim>::signedDistance(const TV& X) const
{
    return (X - center).norm() - radius;
}

template <class T, int dim>
typename Sphere<T, dim>::TV Sphere<T, dim>::closestPointOnSurface(const TV& X) const
{
    TV outward_normal = X - center;
    return radius * (outward_normal).normalized();
}

template <class T, int dim>
typename Sphere<T, dim>::TV Sphere<T, dim>::normal(const TV& X) const
{
    TV outward_normal = X - center;
    if (outward_normal.squaredNorm() < (T)1e-7)
        return TV::Unit(0);
    return (outward_normal).normalized();
}

template <class T, int dim>
typename Sphere<T, dim>::TM Sphere<T, dim>::hessian(const TV& X) const
{
    TV outward_normal = X - center;
    if (outward_normal.squaredNorm() < (T)1e-7)
        return TM::Zero();
    T inv_norm = (T)1 / outward_normal.norm();
    return inv_norm * TM::Identity() - std::pow(inv_norm, 3) * outward_normal * outward_normal.transpose();
}

template <class T, int dim>
bool Sphere<T, dim>::queryInside(const TV& X, T& signed_distance, TV& norm, const T phi_critical) const
{
    using MATH_TOOLS::sqr;
    TV to_center = X - center;
    T to_center_dist2 = to_center.squaredNorm();
    T r2 = sqr(radius + phi_critical);
    if (to_center_dist2 < r2) {
        T to_center_dist = std::sqrt(to_center_dist2);
        signed_distance = to_center_dist - radius;
        if (to_center_dist < (T)1e-7)
            norm = TV::Unit(0);
        else
            norm = (1 / to_center_dist) * to_center;
        return true;
    }
    return false;
}

template <class T, int dim>
bool Sphere<T, dim>::query(const TV& X, T& signed_distance, TV& norm, TM* hess) const
{
    signed_distance = signedDistance(X);
    norm = normal(X);
    if (hess) {
        *hess = hessian(X);
    }
    return (signed_distance <= (T)0);
}

template <class T, int dim>
void Sphere<T, dim>::getBounds(TV& min_corner, TV& max_corner) const
{
    min_corner = center - TV::Constant(radius);
    max_corner = center + TV::Constant(radius);
}

template <class T, int dim>
T Sphere<T, dim>::volume() const
{
    using std::pow;
    constexpr T factor = (dim == 2) ? M_PI : 4 * M_PI / 3;
    return factor * pow(radius, dim);
}

template <class T, int dim>
std::unique_ptr<AnalyticLevelSet<T, dim>> Sphere<T, dim>::createPtr() const
{
    return std::make_unique<Sphere<T, dim>>(*this);
}

template <class T, int dim>
AnalyticBox<T, dim>::AnalyticBox(const TV& half_edges_in, const Vector<T, 4>& q, const TV& b_in)
{
    using std::cos;
    using std::sin;
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
}

template <class T, int dim>
template <class S>
S AnalyticBox<T, dim>::signedDistancePrimitive(const Vector<S, dim>& X) const
{
    using std::abs;
    using std::max;
    using std::min;
    using SV = Vector<S, dim>;
    SV d;
    if (dim == 2)
        d << abs(X(0)) - half_edges(0), abs(X(1)) - half_edges(1);
    else if (dim == 3)
        d << abs(X(0)) - half_edges(0), abs(X(1)) - half_edges(1), abs(X(2)) - half_edges(2);
    S dd = d.array().maxCoeff();
    SV qq = d;
    for (int i = 0; i < dim; i++)
        if (qq(i) < (S)0)
            qq(i) = (S)0;
    S result = min(dd, (S)0) + qq.norm();
    return result;
}

template <class T, int dim>
T AnalyticBox<T, dim>::signedDistance(const TV& X) const
{
    TV X_primitive = R_inverse * (X - b);
    return signedDistancePrimitive(X_primitive);
}

template <class T, int dim>
typename AnalyticBox<T, dim>::TV AnalyticBox<T, dim>::normal(const TV& X) const
{
    TV X_primitive = R_inverse * (X - b);
    ADVec<T, dim> x = vars(X_primitive);
    ADScalar<T, dim> s = signedDistancePrimitive(x);
    TV n = R * s.ds;
    return n;
}

template <class T, int dim>
void AnalyticBox<T, dim>::getBounds(TV& min_corner, TV& max_corner) const
{
    using namespace MATH_TOOLS;
    T bounding_sphere_radius = half_edges.norm();
    min_corner = TV::Constant(-bounding_sphere_radius) + b;
    max_corner = TV::Constant(bounding_sphere_radius) + b;
}

template <class T, int dim>
T AnalyticBox<T, dim>::volume() const
{
    T result = (T)1;
    for (int i = 0; i < dim; i++)
        result *= (2 * std::abs(half_edges(i)));
    return result;
}

template <class T, int dim>
std::unique_ptr<AnalyticLevelSet<T, dim>> AnalyticBox<T, dim>::createPtr() const
{
    return std::make_unique<AnalyticBox<T, dim>>(*this);
}

// Torus

template <class T, int dim>
Torus<T, dim>::Torus(const T r0_in, const T r1_in, const Vector<T, 4>& q, const TV& b_in)
    : r0(r0_in)
    , r1(r1_in)
    , b(b_in)
{
    ZIRAN_ASSERT(dim == 3);
    R = Eigen::Quaternion<T>(q(0), q(1), q(2), q(3)).normalized().toRotationMatrix().template topLeftCorner<dim, dim>();
    R_inverse = R.inverse();
}

template <class T, int dim>
template <class S>
S Torus<T, dim>::signedDistancePrimitive(const Vector<S, dim>& X) const
{
    using std::max;
    using std::min;
    using SV2 = Vector<S, 2>;
    SV2 xz;
    xz << X(0), X(2);
    SV2 q;
    q << xz.norm() - r0, X(1);
    S result = q.norm() - r1;
    return result;
}

template <class T, int dim>
T Torus<T, dim>::signedDistance(const TV& X) const
{
    TV X_primitive = R_inverse * (X - b);
    return signedDistancePrimitive(X_primitive);
}

template <class T, int dim>
typename Torus<T, dim>::TV Torus<T, dim>::normal(const TV& X) const
{
    TV X_primitive = R_inverse * (X - b);
    ADVec<T, dim> x = vars(X_primitive);
    ADScalar<T, dim> s = signedDistancePrimitive(x);
    TV n = R * s.ds;

    return n;
}

template <class T, int dim>
void Torus<T, dim>::getBounds(TV& min_corner, TV& max_corner) const
{
    using namespace MATH_TOOLS;
    T bounding_sphere_radius = r0 + r1;
    min_corner = TV::Constant(-bounding_sphere_radius) + b;
    max_corner = TV::Constant(bounding_sphere_radius) + b;
}

template <class T, int dim>
T Torus<T, dim>::volume() const
{
    return std::abs(M_PI * MATH_TOOLS::sqr(r1) * 2 * M_PI * r0); // the abs is for football (negative r0). TODO: fix.
}

template <class T, int dim>
std::unique_ptr<AnalyticLevelSet<T, dim>> Torus<T, dim>::createPtr() const
{
    return std::make_unique<Torus<T, dim>>(*this);
}

template class DisjointUnionLevelSet<double, 2>;
template class DisjointUnionLevelSet<double, 3>;
template class DisjointUnionLevelSet<float, 2>;
template class DisjointUnionLevelSet<float, 3>;

template class DifferenceLevelSet<double, 2>;
template class DifferenceLevelSet<double, 3>;
template class DifferenceLevelSet<float, 2>;
template class DifferenceLevelSet<float, 3>;

template class HalfSpace<double, 2>;
template class HalfSpace<double, 3>;
template class HalfSpace<float, 2>;
template class HalfSpace<float, 3>;

template class AxisAlignedAnalyticBox<double, 2>;
template class AxisAlignedAnalyticBox<double, 3>;
template class AxisAlignedAnalyticBox<float, 2>;
template class AxisAlignedAnalyticBox<float, 3>;

template class AnalyticBox<double, 2>;
template class AnalyticBox<double, 3>;
template class AnalyticBox<float, 2>;
template class AnalyticBox<float, 3>;
template class Sphere<double, 2>;
template class Sphere<double, 3>;
template class Sphere<float, 2>;
template class Sphere<float, 3>;
template class Torus<double, 2>;
template class Torus<double, 3>;
template class Torus<float, 2>;
template class Torus<float, 3>;
} // namespace ZIRAN
