#include <Ziran/Math/Geometry/PoissonDisk.h>
#include <Ziran/Math/Geometry/RandomSampling.h>
#include <Ziran/Math/Geometry/AnalyticLevelSet.h>
#include <Ziran/Math/Geometry/SourceCollisionObject.h>

namespace ZIRAN {
template <class T, int dim>
void SourceCollisionObject<T, dim>::uniformSample(T grid_dx, int total_particles_number)
{
    static int random_seed = 123;
    random_seed++;
    TV min_corner, max_corner;
    ZIRAN_ASSERT(ls != nullptr);
    ls->getBounds(min_corner, max_corner);

    RandomSampling<T, dim> rs(random_seed, min_corner, max_corner, total_particles_number);
    rs.sample(samples, [&](const TV& x) { TV v; materialVelocity(x, v); return (ls->inside(x) && !ls->inside(x)); });
    ZIRAN_ASSERT(samples.size() == (size_t)total_particles_number);
}

template <class T, int dim>
void SourceCollisionObject<T, dim>::poissonSample(T grid_dx, T particles_per_cell, T density)
{
    static int random_seed = 123;
    random_seed++;
    TV min_corner, max_corner;
    ZIRAN_ASSERT(ls != nullptr);
    ls->getBounds(min_corner, max_corner);
    AxisAlignedAnalyticBox<T, dim> bb(min_corner, max_corner); // bounding box

    PoissonDisk<T, dim> pd(random_seed, 0, min_corner, max_corner, 30, /*periodic=*/true);
    pd.setDistanceByParticlesPerCell(grid_dx, particles_per_cell);
    if (dim == 3)
        pd.sampleFromPeriodicData(samples, [&](const TV& x) { return bb.inside(x); });
    else
        pd.sample(samples, [&](const TV& x) { return bb.inside(x); });

    // compute total volume and set mass accordingly
    T total_volume = 1;
    for (int d = 0; d < dim; ++d) {
        total_volume *= (max_corner(d) - min_corner(d));
    }
    // set uniform mass
    int N = samples.size();
    uniform_mass = total_volume * density / N;
}

template <class T, int dim>
T SourceCollisionObject<T, dim>::evalMaxSpeed(const TV& p_min_corner, const TV& p_max_corner) const
{
    TV min_corner, max_corner, center;
    ls->getBounds(min_corner, max_corner);
    center = (T)0.5 * (min_corner + max_corner);

    // get the speed at the center
    TV v;
    materialVelocity(center, v);
    T max_speed = v.norm();
    materialVelocity(min_corner, v);
    if (v.norm() > max_speed) {
        max_speed = v.norm();
    }
    materialVelocity(max_corner, v);
    if (v.norm() > max_speed) {
        max_speed = v.norm();
    }
    return max_speed;

    /*for (int d = dim - 1; d >= 0; --d) {

    }

    for (int d = 0; d < dim; ++d) {
        for (int low = 0; low < 1; ++low) {
            x(d) = low ? max_corner(d) : min_corner(d);
            T speed = v.norm();
            if (speed > max_speed) {
                max_speed = speed;
            }
        }
    }*/
}

template <class T, int dim>
void SourceCollisionObject<T, dim>::addParticlesFromSample(T dt, T grid_dx)
{
    add_samples.clear();

    for (int i = samples.size() - 1; i >= 0; --i) {
        TV s = samples[i];
        TV v;
        materialVelocity(s, v);
        s += dt * v;
        if (!ls->inside(s)) {
            add_samples.emplace_back(s);
            samples[i] = samples.back();
            int kk = samples.size() - 1;
            if (kk > 0)
                samples.resize(kk);
        }
    }
}

template <class T, int dim>
void SourceCollisionObject<T, dim>::sampleAndPrune(T dt, T grid_dx)
{
    // move the particles in the sample modulo bounding box of the Analytic Level Set
    ZIRAN_ASSERT(ls != nullptr);
    TV min_corner, max_corner;
    ls->getBounds(min_corner, max_corner);

    // clear add_samples
    add_samples.clear();
    inside.clear();

    for (auto& s : samples) {
        // condition A: check if current position is in the sample
        bool originally_inside = ls->inside(s);
        // condition B: check if after the particles being advected (modulo the bounding box) that the particle is outside of the sample
        TV v;
        materialVelocity(s, v);
        s += dt * v;

        // check if the new position is outside of the level set
        bool endedup_outside = !(ls->inside(s));
        inside.push_back(endedup_outside);

        // if both A and B are true, then add particles to add_samples
        if (originally_inside && endedup_outside) {
            add_samples.emplace_back(s + b);
        }

        // move s modulo the bounding box
        for (int i = 0; i < dim; ++i) {
            if (s(i) > max_corner(i)) {
                s(i) = min_corner(i) - max_corner(i) + s(i);
            }
            else if (s(i) < min_corner(i)) {
                s(i) = max_corner(i) - min_corner(i) + s(i);
            }
        }
    }
    std::cout << "inside.size() inside sampleAndPrune = " << inside.size() << std::endl;
}

template class SourceCollisionObject<double, 2>;
template class SourceCollisionObject<double, 3>;
template class SourceCollisionObject<float, 2>;
template class SourceCollisionObject<float, 3>;
} // namespace ZIRAN
