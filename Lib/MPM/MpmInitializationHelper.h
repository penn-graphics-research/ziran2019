#ifndef MPM_INITIALIZATION_HELPER_H
#define MPM_INITIALIZATION_HELPER_H

#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/Math/Geometry/CollisionObject.h>
#include <Ziran/Math/Geometry/SourceCollisionObject.h>
#include <Ziran/Math/Geometry/PoissonDisk.h>
#include <Ziran/Math/Geometry/RandomSampling.h>
#include <Ziran/Math/Geometry/CurveFinder.h>
#include <Ziran/Math/Geometry/Grid.h>
#include <Ziran/Physics/PlasticityApplier.h>
#include <MPM/Force/FBasedMpmForceHelper.h>
#include <MPM/Forward/MpmForward.h>
#include <MPM/MpmParticleHandleBase.h>
#include <MPM/MpmSimulationBase.h>
#include <float.h>

namespace ZIRAN {

template <class T, int dim>
class SourceCollisionObject;

template <class T, int dim>
class Sphere;

template <class T, int dim>
class Grid;

template <class T, int dim>
class VdbLevelSet;

template <class T, int dim>
class PoissonDisk;

template <class T, int dim>
class MpmInitializationHelper {
public:
    static constexpr int interpolation_degree = ZIRAN_MPM_DEGREE;
    using TV = Vector<T, dim>;
    using IV = Vector<int, dim>;

    MpmSimulationBase<T, dim>& mpm;
    Scene<T, dim>& scene;

    MpmInitializationHelper(MpmSimulationBase<T, dim>& mpm)
        : mpm(mpm)
        , scene(mpm.scene)
    {
    }

    template <class BaryStack>
    void sampleBarycentricWeights(const int count, BaryStack& barys, RandomNumber<typename BaryStack::Scalar>& rand)
    {
        constexpr static int manifold_dim = BaryStack::RowsAtCompileTime - 1;
        static const T bary_table[4][3] = { { (T)1 / 3, (T)1 / 3, (T)1 / 3 }, { (T)1 / 6, (T)1 / 6, (T)4 / 6 }, { (T)4 / 6, (T)1 / 6, (T)1 / 6 }, { (T)1 / 6, (T)4 / 6, (T)1 / 6 } };

        for (int p = 0; p < count; p++) {
            if (manifold_dim == 1) { // hair
                barys(0, p) = (1 + p) * (T)1 / (count + 1);
                barys(1, p) = 1 - barys(0, p);
            }
            else if (manifold_dim == 2) { // 2d cloth in 3d
                if (p < 4) {
                    barys(0, p) = bary_table[p][0];
                    barys(1, p) = bary_table[p][1];
                    barys(2, p) = bary_table[p][2];
                }
                else {
                    barys.col(p) = rand.template randomBarycentricWeights<manifold_dim + 1>();
                }
            }
            else
                ZIRAN_ASSERT(false);
        }
    }

    MpmParticleHandleBase<T, dim> sampleOneParticle(const TV& position, const TV& velocity, T density = 1, T total_volume = 1)
    {
        ZIRAN_INFO("Sampling one particle");
        size_t N = 1;
        Range particle_range = mpm.particles.getNextRange(N);
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(position));
        mpm.particles.add(mpm.particles.V_name(), particle_range, std::move(velocity));
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> sampleFromVdbFile(std::string filename, T density, T particles_per_cell = (1 << dim))
    {
        T per_particle_volume = std::pow(mpm.dx, dim) / particles_per_cell;
        VdbLevelSet<float, dim> vdbls(filename); // houdini's vdb files are stored with floats
        Vector<float, dim> min_corner, max_corner;
        vdbls.getBounds(min_corner, max_corner);

        PoissonDisk<T, dim> pd(/*random seed*/ 123, 0, min_corner.template cast<T>(), max_corner.template cast<T>());
        pd.setDistanceByParticlesPerCell(mpm.dx, particles_per_cell);
        StdVector<TV> samples;
        if (dim == 3)
            pd.sampleFromPeriodicData(samples, [&](const TV& x) { return vdbls.inside(x.template cast<float>()); });
        else
            pd.sample(samples, [&](const TV& x) { return vdbls.inside(x.template cast<float>()); });
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, per_particle_volume * density);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, per_particle_volume * N);
    }

    MpmParticleHandleBase<T, dim> sampleFromVdbFileWithExistingPoints(StdVector<TV>& existing_samples, std::string filename, T density, T particles_per_cell = (1 << dim))
    {
        T per_particle_volume = std::pow(mpm.dx, dim) / particles_per_cell;
        VdbLevelSet<float, dim> vdbls(filename); // houdini's vdb files are stored with floats
        Vector<float, dim> min_corner, max_corner;
        vdbls.getBounds(min_corner, max_corner);

        PoissonDisk<T, dim> pd(/*random seed*/ 123, 0, min_corner.template cast<T>(), max_corner.template cast<T>());
        pd.setDistanceByParticlesPerCell(mpm.dx, particles_per_cell);
        StdVector<TV> samples;
        if (dim == 3)
            pd.sampleFromPeriodicData(samples, [&](const TV& x) { return vdbls.inside(x.template cast<float>()); });
        else
            pd.sample(samples, [&](const TV& x) { return vdbls.inside(x.template cast<float>()); });

        existing_samples.insert(existing_samples.end(), samples.begin(), samples.end());
        samples = existing_samples;

        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, per_particle_volume * density);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, per_particle_volume * N);
    }

    void sampleSourceAtTheBeginning(int pos, T density, T particles_per_cell = (1 << dim))
    {
        auto source = dynamic_cast<SourceCollisionObject<T, dim>*>(mpm.collision_objects[pos].get());
        ZIRAN_ASSERT(source != nullptr, "mpm.collision_objects[pos] is not a source collision object!");
        source->poissonSample(mpm.dx, particles_per_cell, density);
    }

    size_t sourceSampleAndPrune(int pos, T density, T particles_per_cell = (1 << dim))
    {
        auto source = dynamic_cast<SourceCollisionObject<T, dim>*>(mpm.collision_objects[pos].get());
        ZIRAN_ASSERT(source != nullptr, "mpm.collision_objects[pos] is not a source collision object!");
        source->sampleAndPrune(mpm.dt, mpm.dx);
        return source->add_samples.size();
    }

    MpmParticleHandleBase<T, dim> getZeroParticle()
    {
        Range particle_range = mpm.particles.getNextRange(0);
        mpm.particles.add(mpm.particles.X_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, 0);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, 0);
    }

    MpmParticleHandleBase<T, dim> getParticlesFromSource(int pos, T density, T particles_per_cell = (1 << dim))
    {
        auto source = dynamic_cast<SourceCollisionObject<T, dim>*>(mpm.collision_objects[pos].get());
        ZIRAN_ASSERT(source != nullptr, "mpm.collision_objects[pos] is not a source collision object!");
        T per_particle_volume = std::pow(mpm.dx, dim) / particles_per_cell;
        size_t N = source->add_samples.size();
        ZIRAN_ASSERT(N);
        Range particle_range = mpm.particles.getNextRange(N);
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(source->add_samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, source->const_material_velocity);
        mpm.particles.add(mpm.particles.mass_name(), particle_range, per_particle_volume * density);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, per_particle_volume * N);
    }

    MpmParticleHandleBase<T, dim> sampleFromObjPointCloudFile(std::string filename, const T mass_per_particle, const T volume_per_particle)
    {
        StdVector<TV> samples;
        std::string absolute_path = DataDir().absolutePath(filename);
        readPositionObj(absolute_path, samples);
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, mass_per_particle);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, N * volume_per_particle);
    }

    T sampleInAnalyticLevelSetHelper(AnalyticLevelSet<T, dim>& levelset, T particles_per_cell, StdVector<Vector<T, dim>>& samples)
    {
        using TV = Vector<T, dim>;
        ZIRAN_INFO("Sampling ", particles_per_cell, " particles per cell in the levelset");
        TV min_corner, max_corner;
        levelset.getBounds(min_corner, max_corner);
        T total_volume = levelset.volume();
        PoissonDisk<T, dim> pd(/*random seed*/ 123, 0, min_corner, max_corner);
        pd.setDistanceByParticlesPerCell(mpm.dx, particles_per_cell);

        if (dim == 3)
            pd.sampleFromPeriodicData(samples, [&](const TV& x) { return levelset.inside(x); });
        else
            pd.sample(samples, [&](const TV& x) { return levelset.inside(x); });
        return total_volume;
    }

    T sampleWaterInWaterLevelSetOutsideSandLevelSetHelper(AnalyticLevelSet<T, dim>& water_levelset, AnalyticLevelSet<T, dim>& sand_levelset, T particles_per_cell, StdVector<Vector<T, dim>>& samples)
    {
        using TV = Vector<T, dim>;
        ZIRAN_INFO("Sampling ", particles_per_cell, " particles per cell in the levelset");
        TV min_corner, max_corner;
        water_levelset.getBounds(min_corner, max_corner);
        T total_volume = water_levelset.volume();
        PoissonDisk<T, dim> pd(/*random seed*/ 123, 0, min_corner, max_corner);
        pd.setDistanceByParticlesPerCell(mpm.dx, particles_per_cell);

        if (dim == 3)
            pd.sampleFromPeriodicData(samples, [&](const TV& x) { return water_levelset.inside(x) && !sand_levelset.inside(x); });
        else
            pd.sample(samples, [&](const TV& x) { return water_levelset.inside(x) && !sand_levelset.inside(x); });
        return total_volume;
    }

    T sampleInAnalyticLevelSetHelperSpecial(AnalyticLevelSet<T, dim>& levelset, T particles_per_cell, StdVector<Vector<T, dim>>& samples, const Vector<T, dim>& min_bound, const Vector<T, dim>& max_bound)
    {
        using TV = Vector<T, dim>;
        ZIRAN_INFO("Sampling ", particles_per_cell, " particles per cell in the levelset");
        TV min_corner, max_corner;
        levelset.getBounds(min_corner, max_corner);
        T total_volume = levelset.volume();
        PoissonDisk<T, dim> pd(/*random seed*/ 123, 0, min_corner, max_corner);
        pd.setDistanceByParticlesPerCell(mpm.dx, particles_per_cell);
        if (dim == 3)
            pd.sampleFromPeriodicData(samples, [&](const TV& x) { return levelset.inside(x); });
        else
            pd.sample(samples, [&](const TV& x) { return (levelset.inside(x) && (x(0) >= min_bound(0) && x(0) <= max_bound(0)) && (x(1) >= min_bound(1) && x(1) <= max_bound(1))); });
        return total_volume;
    }

    T sampleInAnalyticLevelSetHelperSpecial(AnalyticLevelSet<T, dim>& levelset_in, AnalyticLevelSet<T, dim>& levelset_out, T particles_per_cell, StdVector<Vector<T, dim>>& samples, const Vector<T, dim>& min_bound, const Vector<T, dim>& max_bound)
    {
        using TV = Vector<T, dim>;
        ZIRAN_INFO("Sampling ", particles_per_cell, " particles per cell in the levelset");
        TV min_corner, max_corner;
        levelset_in.getBounds(min_corner, max_corner);
        T total_volume = levelset_in.volume();
        PoissonDisk<T, dim> pd(/*random seed*/ 123, 0, min_corner, max_corner);
        pd.setDistanceByParticlesPerCell(mpm.dx, particles_per_cell);
        if (dim == 3)
            pd.sampleFromPeriodicData(samples, [&](const TV& x) { return levelset_in.inside(x); });
        else
            pd.sample(samples, [&](const TV& x) { return levelset_in.inside(x) && !((levelset_out.inside(x) && (x(0) >= min_bound(0) && x(0) <= max_bound(0)) && (x(1) >= min_bound(1) && x(1) <= max_bound(1)))); });
        return total_volume;
    }

    MpmParticleHandleBase<T, dim> sampleWaterInWaterLevelSetOutsideSandLevelSet(AnalyticLevelSet<T, dim>& water_levelset, AnalyticLevelSet<T, dim>& sand_levelset, T density, T particles_per_cell)
    {
        StdVector<TV> samples;
        T total_volume = sampleWaterInWaterLevelSetOutsideSandLevelSetHelper(water_levelset, sand_levelset, particles_per_cell, samples);
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> readFromFile(const char* filename, int number, T density)
    {
        // comment another element_measure_name() function first

        StdVector<TV> samples;
        StdVector<TV> vs;

        TV v = TV::Zero();
        T mass, total_mass = 0, vol;
        StdVector<T> vols;
        StdVector<T> masss;
        FILE* f = fopen(filename, "r");
        for (int i = 0; i < number; ++i) {
            TV pos;
            if constexpr (std::is_same<T, float>::value) {
                for (int d = 0; d < dim; ++d)
                    fscanf(f, "%f", &pos(d));
                samples.push_back(pos);
                for (int d = 0; d < dim; ++d)
                    fscanf(f, "%f", &v(d));
                vs.push_back(v);
                fscanf(f, "%f", &mass);
                masss.push_back(mass);
                total_mass += mass;
                fscanf(f, "%f", &vol);
                vols.push_back(vol);
            }
            else {
                for (int d = 0; d < dim; ++d)
                    fscanf(f, "%lf", &pos(d));
                samples.push_back(pos);
                for (int d = 0; d < dim; ++d)
                    fscanf(f, "%lf", &v(d));
                vs.push_back(v);
                fscanf(f, "%lf", &mass);
                masss.push_back(mass);
                total_mass += mass;
                fscanf(f, "%lf", &vol);
                vols.push_back(vol);
            }
        }
        fclose(f);

        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, vs);
        mpm.particles.add(mpm.particles.mass_name(), particle_range, masss);
        mpm.particles.add(element_measure_name<T>(), particle_range, vols);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_mass / density);
    }

    MpmParticleHandleBase<T, dim> sampleInAnalyticLevelSet(AnalyticLevelSet<T, dim>& levelset, T density, T particles_per_cell)
    {
        StdVector<TV> samples;
        T total_volume = sampleInAnalyticLevelSetHelper(levelset, particles_per_cell, samples);
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> sampleInAnalyticLevelSetSpecial(AnalyticLevelSet<T, dim>& levelset, T density, T particles_per_cell, const Vector<T, dim>& min_bound, const Vector<T, dim>& max_bound)
    {
        StdVector<TV> samples;
        T total_volume = sampleInAnalyticLevelSetHelperSpecial(levelset, particles_per_cell, samples, min_bound, max_bound);
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> sampleInAnalyticLevelSetSpecial(AnalyticLevelSet<T, dim>& levelset_in, AnalyticLevelSet<T, dim>& levelset_out, T density, T particles_per_cell, const Vector<T, dim>& min_bound, const Vector<T, dim>& max_bound)
    {
        StdVector<TV> samples;
        T total_volume = sampleInAnalyticLevelSetHelperSpecial(levelset_in, levelset_out, particles_per_cell, samples, min_bound, max_bound);
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> sampleInAnalyticLevelSetWithExistingPoints(StdVector<TV>& existing_samples, AnalyticLevelSet<T, dim>& levelset, T density, T particles_per_cell)
    {
        StdVector<TV> samples;
        T total_volume = sampleInAnalyticLevelSetHelper(levelset, particles_per_cell, samples);
        existing_samples.insert(existing_samples.end(), samples.begin(), samples.end());
        samples = existing_samples;
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> sampleInAnalyticLevelSetPartial0(AnalyticLevelSet<T, dim>& levelset, T density, T particles_per_cell)
    {
        StdVector<TV> _samples;
        T total_volume = sampleInAnalyticLevelSetHelper(levelset, particles_per_cell, _samples);
        StdVector<TV> samples;
        for (const auto& pos : _samples) {
            T y = pos[1];
            if (y < 5)
                samples.push_back(pos);
        }
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> sampleInAnalyticLevelSetPartial(AnalyticLevelSet<T, dim>& levelset, T density, T particles_per_cell)
    {
        StdVector<TV> _samples;
        T total_volume = sampleInAnalyticLevelSetHelper(levelset, particles_per_cell, _samples);
        StdVector<TV> samples;
        for (const auto& pos : _samples) {
            T sqrt2 = std::sqrt((T)2.0);
            T x = pos[0];
            T y = pos[1];
            if (y > -x + 2.5 + 0.2 * sqrt2 && y > x - 0.1 + 0.2 * sqrt2)
                samples.push_back(pos);
        }
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> sampleInAnalyticLevelSetPartial2(AnalyticLevelSet<T, dim>& levelset, T density, T particles_per_cell)
    {
        StdVector<TV> _samples;
        T total_volume = sampleInAnalyticLevelSetHelper(levelset, particles_per_cell, _samples);
        StdVector<TV> samples;
        for (const auto& pos : _samples) {
            T y = pos[1];
            if (y < 3)
                samples.push_back(pos);
        }
        size_t N = samples.size();
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> sampleOnLine(const TV& A, const TV& B, T density, T h, T thickness)
    {
        int N = int((A - B).norm() / h) + 1;
        TV d = (B - A) / (N - 1);
        StdVector<TV> samples;
        for (int i = 0; i < N; i++)
            samples.push_back(A + i * d);
        T total_volume = (N - 1) * h * thickness;
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> samplePackedSpheresInAnalyticLevelSet(AnalyticLevelSet<T, dim>& levelset, T density, T particles_per_cell, T sphere_radius, T gap)
    {
        TV min_corner, max_corner;
        levelset.getBounds(min_corner, max_corner);
        T total_volume = levelset.volume();
        PoissonDisk<T, dim> pd(/*random seed*/ 123, sphere_radius * 2 + gap, min_corner, max_corner);
        StdVector<TV> center_samples;
        if (dim == 3)
            pd.sampleFromPeriodicData(center_samples, [&](const TV& x) { return levelset.inside(x); });
        else
            pd.sample(center_samples, [&](const TV& x) { return levelset.inside(x); });
        size_t N_spheres = center_samples.size();

        ZIRAN_INFO("N spheres:", N_spheres);

        StdVector<TV> samples;
        for (size_t k = 0; k < N_spheres; k++) {
            const TV center = center_samples[k];
            Sphere<T, dim> local_sphere(center, sphere_radius);
            PoissonDisk<T, dim> pd_local(/*random seed*/ 123, 0, center - sphere_radius * TV::Ones(), center + sphere_radius * TV::Ones());
            pd_local.setDistanceByParticlesPerCell(mpm.dx, particles_per_cell);

            StdVector<TV> local_samples;
            if (dim == 3)
                pd_local.sampleFromPeriodicData(local_samples, [&](const TV& x) { return local_sphere.inside(x); });
            else
                pd_local.sample(local_samples, [&](const TV& x) { return local_sphere.inside(x); });

            for (auto s : local_samples)
                samples.push_back(s);
        }
        size_t N = samples.size();
        ZIRAN_INFO("sampled particle: ", N);
        Range particle_range = mpm.particles.getNextRange(samples.size());
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / N);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> uniformGridSample(const Grid<T, dim>& grid, T density)
    {
        ZIRAN_INFO("Sampling totally ", grid.numberNodes(), " particles on a Cartesian grid");
        StdVector<TV> samples;
        grid.getAllPositions(samples);
        int point_number = samples.size();
        ZIRAN_ASSERT(point_number == grid.numberNodes());
        T total_volume = grid.gridVolume();

        Range particle_range = mpm.particles.getNextRange(point_number);
        ZIRAN_ASSERT(particle_range.upper - particle_range.lower == point_number);
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / point_number);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    //!!! point_number is total number of points sampled instead of particles per cell
    MpmParticleHandleBase<T, dim> uniformSampleInAnalyticLevelSet(AnalyticLevelSet<T, dim>& levelset, T density, int point_number)
    {
        ZIRAN_INFO("Sampling totally ", point_number, " particles in the levelset");
        TV min_corner, max_corner;
        levelset.getBounds(min_corner, max_corner);
        RandomSampling<T, dim> rs(123, min_corner, max_corner, point_number);
        StdVector<TV> samples;
        rs.sample(samples, [&](const TV& x) { return levelset.inside(x); });
        ZIRAN_ASSERT(samples.size() == (size_t)point_number);
        Range particle_range = mpm.particles.getNextRange(point_number);
        ZIRAN_ASSERT(particle_range.upper - particle_range.lower == point_number);
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        T total_volume = levelset.volume();
        ZIRAN_ASSERT(total_volume == total_volume);
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / point_number);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    MpmParticleHandleBase<T, dim> uniformSampleInAnalyticLevelSetSpecifyPPC(AnalyticLevelSet<T, dim>& levelset, T density, T ppc)
    {
        ZIRAN_INFO("Sampling on average ", ppc, " particles per cell randomly in the levelset");
        T volume = levelset.volume();
        T cells = volume / std::pow(mpm.dx, dim);
        int point_number = cells * ppc;
        TV min_corner, max_corner;
        levelset.getBounds(min_corner, max_corner);
        RandomSampling<T, dim> rs(123, min_corner, max_corner, point_number);
        StdVector<TV> samples;
        rs.sample(samples, [&](const TV& x) { return levelset.inside(x); });
        ZIRAN_ASSERT(samples.size() == (size_t)point_number);
        Range particle_range = mpm.particles.getNextRange(point_number);
        ZIRAN_ASSERT(particle_range.upper - particle_range.lower == point_number);
        mpm.particles.add(mpm.particles.X_name(), particle_range, std::move(samples));
        mpm.particles.add(mpm.particles.V_name(), particle_range, TV::Zero());
        T total_volume = levelset.volume();
        ZIRAN_ASSERT(total_volume == total_volume);
        mpm.particles.add(mpm.particles.mass_name(), particle_range, total_volume * density / point_number);
        return MpmParticleHandleBase<T, dim>(mpm.getParticles(), mpm.getScene(), mpm.getMpmForce(), mpm.getPlasticityAppliers(), mpm.getScratchXp(), mpm.getDt(), particle_range, total_volume);
    }

    void addExplicitVelocityField(const std::function<void(T&, const TV&, TV&)>& mapping)
    {
        mpm.explicit_velocity_field = mapping;
    }

    void addAnalyticCollisionObject(AnalyticCollisionObject<T, dim>& c)
    {
        mpm.addCollisionObject(std::move(c));
    }

    int addSourceCollisionObject(SourceCollisionObject<T, dim>& c)
    {
        return mpm.addSourceCollisionObject(std::move(c));
    }

    void addExternalBodyForce(T scale, T frequency, const TV& direction)
    {
        auto ff = [=](const TV& x, const T t) -> TV {
            TV f = TV::Zero();
            f += scale * (std::sin(frequency * t) + 1) * direction;
            return f;
        };
        mpm.fext.emplace_back(std::move(ff));
    }

    int getNumNodes()
    {
        return mpm.num_nodes;
    }

    void scaleCotangent(Matrix<T, dim, dim> scale)
    {
        using TM = Matrix<T, dim, dim>;
        Range particle_range(0, mpm.particles.count);

        DisjointRanges subset(DisjointRanges{ particle_range },
            mpm.particles.commonRanges(AttributeName<TM>("cotangent")));

        for (auto iter = mpm.particles.subsetIter(subset, AttributeName<TM>("cotangent")); iter; ++iter) {
            auto& F = iter.template get<0>();
            F.array() = F.array() * scale.array(); // coefficient wise scale
        }
    }

    void applyPlasticity()
    {
        mpm.applyPlasticity();
    }

    T signedDistanceToCollisionObject(const int pid, const int object_id)
    {
        return mpm.collision_objects[object_id]->signedDistance(mpm.particles.X(pid));
    }

    template <int splat_size>
    void addFS(const MpmParticleHandleBase<T, dim>& particles_handle)
    {
        using TMFS = Matrix<T, splat_size - 1, dim>;
        Matrix<T, splat_size - 1, dim> FS = Matrix<T, splat_size - 1, dim>::Zero();
        particles_handle.particles.add(AttributeName<TMFS>("fullS"), particles_handle.particle_range, FS);
    }

    template <int splat_size>
    void addFSToAllParticles()
    {
        using TMFS = Matrix<T, splat_size - 1, dim>;
        Range range;
        range.lower = 0;
        range.upper = mpm.particles.count;
        Matrix<T, splat_size - 1, dim> FS = Matrix<T, splat_size - 1, dim>::Zero();
        mpm.particles.add(AttributeName<TMFS>("fullS"), range, FS);
    }

    //#########################################################################
    // A helper function to add all walls in the [0,max_domain_boundary] domain with offset.
    //#########################################################################
    template <class BoundaryType>
    void addAllWallsInDomain(const T max_domain_boundary, const T offset, BoundaryType boundary_type) // boundary_type = AnalyticCollisionObject<T, dim>::SLIP/STICKY/SEPARATE
    {
        for (int d = 0; d < dim; d++) {
            for (int s = -1; s <= 1; s += 2) {
                TV O;
                if (s == 1)
                    O = offset * TV::Unit(d);
                else
                    O = (TV::Unit(d) * (max_domain_boundary - offset));
                TV N = TV::Unit(d) * s;
                HalfSpace<T, dim> ls(O, N);
                AnalyticCollisionObject<T, dim> ground([&](T, AnalyticCollisionObject<T, dim>&) {}, ls, boundary_type);
                addAnalyticCollisionObject(ground);
            }
        }
    }
};
} // namespace ZIRAN

#endif
