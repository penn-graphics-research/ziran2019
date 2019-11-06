#pragma once

#include "MpmParticleHandleBase.h"
#include <MPM/MpmSimulationBase.h>

#include <Ziran/CS/DataStructure/DataManager.h>
#include <Ziran/CS/DataStructure/KdTree.h>
#include <Ziran/CS/Util/DataDir.h>
#include <Ziran/CS/Util/AttributeNamesForward.h>
#include <Ziran/Math/Geometry/ObjIO.h>
#include <Ziran/Math/Geometry/PoissonDisk.h>
#include <Ziran/Math/Geometry/VoronoiNoise.h>
#include <Ziran/Math/Geometry/AnalyticLevelSet.h>
#include <Ziran/Math/Geometry/VdbLevelSet.h>
#include <Ziran/Math/Geometry/VtkIO.h>
#include <Ziran/Physics/LagrangianForce/LagrangianForce.h>
#include <Ziran/Physics/ConstitutiveModel/ConstitutiveModel.h>
#include <Ziran/Physics/PlasticityApplier.h>
#include <Ziran/Sim/Scene.h>

#include <MPM/Force/FBasedMpmForceHelper.h>
#include <MPM/Force/FElasticNonequilibratedBasedMpmForceHelper.h>

#include "../../Projects/fracture/PhaseField.h"

namespace ZIRAN {

template <class T, int dim>
MpmParticleHandleBase<T, dim>::
    MpmParticleHandleBase(Particles<T, dim>& particles, Scene<T, dim>& scene, MpmForceBase<T, dim>* mpmforce,
        StdVector<std::unique_ptr<PlasticityApplierBase>>& plasticity_appliers,
        StdVector<TV>& scratch_xp, T& dt, Range particle_range, T total_volume, int cotangent_manifold_dim)
    : particles(particles)
    , scene(scene)
    , mpmforce(mpmforce)
    , plasticity_appliers(plasticity_appliers)
    , scratch_xp(scratch_xp)
    , dt(dt)
    , particle_range(particle_range)
    , total_volume(total_volume)
    , cotangent_manifold_dim(cotangent_manifold_dim)
{
}

// Creates a copy with new particles
template <class T, int dim>
MpmParticleHandleBase<T, dim>
MpmParticleHandleBase<T, dim>::copy()
{
    Range new_particle_range;
    new_particle_range.lower = particles.count;
    {
        auto ap = particles.appender();
        for (int i = particle_range.lower; i < particle_range.upper; i++)
            ap.append(particles.mass[i], particles.X[i], particles.V[i]);
    }
    new_particle_range.upper = particles.count;
    return MpmParticleHandleBase(particles, scene, mpmforce, plasticity_appliers, scratch_xp, dt, new_particle_range, total_volume);
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::transform(const std::function<void(int, Ref<T>, Vector<T, dim>&, Vector<T, dim>&)>& mapping)
{
    for (int i = particle_range.lower; i < particle_range.upper; ++i) {
        // lua does not support passing scalars by reference. This is a work around to actually change mass.
        mapping(i - particle_range.lower, particles.mass[i], particles.X[i], particles.V[i]);
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    addVolumeFraction(const T b)
{
    particles.add(volume_fraction_name<T>(), particle_range, b);
}

template <class T, int dim>
template <class TCONST>
void MpmParticleHandleBase<T, dim>::
    addFBasedMpmForce(const TCONST& model)
{
    if (cotangent_manifold_dim == 1)
        ZIRAN_ASSERT(false, "not implemented");
    else if (cotangent_manifold_dim == 2)
        ZIRAN_ASSERT(false, "not implemented");
    else if (cotangent_manifold_dim == 0)
        addFBasedMpmForceWithMeasure(model, particle_range, total_volume);
    else
        ZIRAN_ASSERT(false);
}

template <class T, int dim>
template <class TCONST>
void MpmParticleHandleBase<T, dim>::
    addFBasedMpmForceWithPhaseField(const T& percentage, const T& l0, const TCONST& model, bool allow_damage, const T random_fluctuation_percentage)
{
    ZIRAN_ASSERT(cotangent_manifold_dim == 0);
    FBasedMpmForceHelper<TCONST>& helper = mpmforce->getHelper(); // this will let mpmforce create a consitutive model helper
    particles.add(helper.constitutive_model_name(), particle_range, model);
    if (total_volume != 0) {
        particles.add(element_measure_name<T>(), particle_range, total_volume / particle_range.length());
        particles.add(helper.F_name(), particle_range, TM::Identity());
    }
    particles.add(helper.constitutive_model_scratch_name(), particle_range, typename FBasedMpmForceHelper<TCONST>::Scratch());
    PhaseField<T, dim> pf;
    pf.residual_phase = (T)0.001;
    pf.c = (T)1;
    pf.H = (T)0;
    pf.l0 = l0;
    pf.one_over_sigma_c = PhaseField<T, dim>::Get_Sigma_C_Based_On_Max_Deformation(percentage, model);
    pf.pf_Fp = (T)1;
    pf.H_max = std::numeric_limits<T>::max();
    pf.vol = total_volume / particle_range.length();
    pf.allow_damage = allow_damage;
    particles.add(phase_field_name<PhaseField<T, dim>>(), particle_range, pf);

    if (random_fluctuation_percentage) randomizePhaseFieldSigmaC(random_fluctuation_percentage);
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::addOriginalPositionAsAttribute()
{
    // Create X0 as attribute
    particles.add(X0_name<T, dim>(), particle_range, TV::Zero());

    // Fill in X0 with current X
    for (auto iter = particles.iter(X0_name<T, dim>()); iter; ++iter) {
        Vector<T, dim>& X0 = iter.template get<0>();
        X0 = particles.X(iter.entryId());
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::randomizePhaseFieldSigmaC(const T random_fluctuation_percentage)
{
    RandomNumber<T> rn(123);
    for (auto iter = particles.subsetIter(DisjointRanges{ particle_range }, phase_field_name<PhaseField<T, dim>>()); iter; ++iter) {
        auto& pf = iter.template get<0>();
        T scale = 1 + rn.randReal(-random_fluctuation_percentage, random_fluctuation_percentage);
        pf.one_over_sigma_c *= scale;
    }
}

//Read in a file of voronoi points and use them to initialize sigmaC to be scaled based on point distance from voronoi surfaces
template <class T, int dim>
void MpmParticleHandleBase<T, dim>::voronoiSigmaC(std::string voronoiPointsFile, T radius)

//void MpmParticleHandleBase<T, dim>::voronoiSigmaC(std::string voronoiPointsFile, T radius, T minPercent, T maxPercent)
{
    //Read the file to grab our points from obj file
    StdVector<TV> vPoints;
    readPositionObj(voronoiPointsFile, vPoints);

    //Set up our KDTree and add each of our voronoi points to it!
    KdTree<dim> voronoiPoints;
    for (int i = 0; i < (int)vPoints.size(); i++) {
        voronoiPoints.addPoint(i, vPoints[i]);
    }
    voronoiPoints.optimize();

    //Now we can query this voronoi surface point KDTree for distances from our actual points!
    int i = 0;
    T zeta = -1; //set -1 until we calculate it for the first time
    for (auto iter = particles.subsetIter(DisjointRanges{ particle_range }, phase_field_name<PhaseField<T, dim>>()); iter; ++iter) {

        auto& pf = iter.template get<0>();

        int id;
        TV p;
        T dist;
        voronoiPoints.findNearest(particles.X[i], id, p, dist); //fill in "dist" with distance from nearest point

        //4th Way: use base to reverse engineer what value of zeta we need to fit the values to the curve y = e^(zeta * x) - 1
        if (zeta == -1) { //only compute once
            T baseSigmaC = (T)1 / pf.one_over_sigma_c;
            zeta = std::log(baseSigmaC + 1) / radius;
        }
        T newSigmaC = (T)1 / pf.one_over_sigma_c;
        if (dist < radius) {
            newSigmaC = std::exp(zeta * dist) - 1;
        }
        pf.one_over_sigma_c = (T)1 / newSigmaC;

        //3rd way: use 1 - log(dist/radius) and clamp all scales past (dist/radius)=1 to be 1
        //NOTE: this relies on scaling SPECIFICALLY one over sigma C!!!!!
        /*T x = dist / radius;
        T scale = 1;
        if (x < 1) {
            scale = (1 - std::log(x)) * magnitude; //scale by magnitude to allow for greater disparity between large and small sigmaC
        }
        pf.one_over_sigma_c *= scale; //NEED to scale one over sigma C here
        */

        //2nd Attempt: Exponential Way
        //T scale = std::exp(zeta * dist);

        //1st Way: linear scale based on dist, scale sigmaC!
        /*T proportion = (std::exp((dist / radius)) - 1) * (maxPercent - minPercent); //proportion through the interval min% to max% we are based on dist --> but using tanh instead of linear
        T scale = 1 + (minPercent + proportion);

        //Clamp all points outside the radius to be the max percent increase
        if (dist >= radius) {
            scale = 1 + maxPercent;
        }

        //clamp values to not be greater than 1 + maxPercent
        if(scale > (1+maxPercent)){
            scale = 1 + maxPercent;
        }

        //Clamp the scale to never be non-zero
        if (scale <= 0) {
            scale = 0.0000000001; //arbitrary small epsilon (to avoid divide by 0)
        }*/

        //want to scale the actual sigmaC, not its inverse!
        //T currSigmaC = (T)1 / pf.one_over_sigma_c;
        //T newSigmaC = currSigmaC * scale;
        //pf.one_over_sigma_c = (T)1 / newSigmaC;
        //pf.one_over_sigma_c *= scale; //for if we want to directly scale one over sigmaC

        i++; //increment particle index since our iterator is just for pf itself
    }
}

//Read in a file of voronoi points and use them to initialize sigmaC to be scaled based on point distance from voronoi surfaces
template <class T, int dim>
void MpmParticleHandleBase<T, dim>::yDimSigmaCInit(T yMin, T yMax, T maxScale)
{
    //Now we can query this voronoi surface point KDTree for distances from our actual points!
    int i = 0;
    for (auto iter = particles.subsetIter(DisjointRanges{ particle_range }, phase_field_name<PhaseField<T, dim>>()); iter; ++iter) {

        auto& pf = iter.template get<0>();
        T yVal = particles.X[i][1]; //grab y val

        //std::cout << "Y val: " << yVal << std::endl;

        if (yVal < yMin) {
            i++;
            continue;
        }
        else if (yVal > yMax) {
            i++;
            continue;
        }

        T dist = (yMax - yVal) / (yMax - yMin);

        T scale = dist * maxScale; //linear scale based on dist

        if (scale < 1) {
            scale = 1;
        }

        T newSigmaC = ((T)1 / pf.one_over_sigma_c) * scale;
        pf.one_over_sigma_c = (T)1 / newSigmaC;

        i++; //increment particle index since our iterator is just for pf itself
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    scaleF(const T scale)
{
    for (auto iter = particles.subsetIter(DisjointRanges{ particle_range }, F_name<T, dim>()); iter; ++iter) {
        auto& F = iter.template get<0>();
        F *= scale;
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    scaleF(const TV scale)
{
    for (auto iter = particles.subsetIter(DisjointRanges{ particle_range }, F_name<T, dim>()); iter; ++iter) {
        auto& F = iter.template get<0>();
        if constexpr (dim >= 1)
            F.row(0) *= scale[0];
        if constexpr (dim >= 2)
            F.row(1) *= scale[1];
        if constexpr (dim >= 3)
            F.row(2) *= scale[2];
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    scaleJ(const T scale)
{
    for (auto iter = particles.subsetIter(DisjointRanges{ particle_range }, J_name<T>()); iter; ++iter) {
        auto& J = iter.template get<0>();
        J *= scale;
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    scaleFCurve(int frame, const std::function<T(int)>& growCurve)
{
    DisjointRanges subset{ particle_range };
    for (auto iter = particles.subsetIter(subset, F_name<T, dim>()); iter; ++iter) {
        auto& F = iter.template get<0>();
        if (frame == 0)
            F /= growCurve(0);
        // else
        //      F *= (growCurve(frame - 1) / growCurve(frame));
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    resetDeformation()
{
    for (auto iter = particles.subsetIter(DisjointRanges{ particle_range }, F_name<T, dim>()); iter; ++iter) {
        auto& F = iter.template get<0>();
        F = Matrix<T, dim, dim>::Identity();
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    addElementMeasure()
{
    if (total_volume != 0) {
        particles.add(element_measure_name<T>(), particle_range, total_volume / particle_range.length());
    }
}

template <class T, int dim>
template <class TConst, class TPlastic>
void MpmParticleHandleBase<T, dim>::
    addPlasticity(const TConst& cons, const TPlastic& plasticity, std::string strain_name)
{
    using TStrain = typename TConst::Strain;
    PlasticityApplier<TConst, TPlastic, TStrain>* plasticity_model = nullptr;
    for (auto& p : plasticity_appliers) {
        plasticity_model = dynamic_cast<PlasticityApplier<TConst, TPlastic, TStrain>*>(p.get());
        if (plasticity_model && plasticity_model->strain_name.name == strain_name)
            break;
        else
            plasticity_model = nullptr;
    }
    if (plasticity_model == nullptr)
        plasticity_appliers.push_back(std::make_unique<PlasticityApplier<TConst, TPlastic, TStrain>>(strain_name));
    particles.add(TPlastic::name(), particle_range, plasticity);
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    addDummyPlasticity()
{
    auto p = DummyPlasticity<T>();
    particles.add(DummyPlasticity<T>::name(), particle_range, p);
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    prescorePhaseFieldSigmaC(AnalyticLevelSet<T, dim>& levelset, T grain_size, T scale_lower_clamp)
{
    ZIRAN_INFO("Prescoring phase field sigma_c");
    TV min_corner, max_corner;
    levelset.getBounds(min_corner, max_corner);
    PoissonDisk<T, dim> pd(/*random seed*/ 123, grain_size, min_corner, max_corner);
    StdVector<TV> grain_centroids;
    if (dim == 3)
        pd.sampleFromPeriodicData(grain_centroids, [&](const TV& x) { return levelset.inside(x); });
    else
        pd.sample(grain_centroids, [&](const TV& x) { return levelset.inside(x); });
    T phi = T(0.5) * (1 - std::sqrt((T)5));
    RandomNumber<T> rand;
    for (auto iter = particles.subsetIter({ particle_range }, Particles<T, dim>::X_name(), phase_field_name<PhaseField<T, dim>>()); iter; ++iter) {
        TV& x = iter.template get<0>();
        for (const TV& c : grain_centroids) {
            TV v = c - x;
            T alpha = std::min(v.norm() / ((T)1.2 * grain_size), (T)1);
            T scale = phi + 1 / (alpha - phi);
            scale *= rand.randReal(0.75, 1.0);
            //x = x + scale * v;
        }
        T min_distance2 = std::numeric_limits<T>::max();
        for (const TV& c : grain_centroids)
            min_distance2 = std::min(min_distance2, (x - c).squaredNorm());
        auto& pf = iter.template get<1>();
        T scale = 1 - std::min(std::sqrt(min_distance2) * rand.randReal(0.75, 1.25) / (2 * grain_size), (T)1);
        scale = scale * (1 - scale_lower_clamp) + scale_lower_clamp;
        T old_sigma_c = (T)1 / pf.one_over_sigma_c;
        T new_sigma_c = old_sigma_c * scale;
        pf.one_over_sigma_c = (T)1 / new_sigma_c;
    }
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::
    prescoreWetSand2D(AnalyticLevelSet<T, dim>& levelset, T grain_size, T density_min, T Jp_min, T scale_min)
{
    ZIRAN_INFO("Prescoring Wet Sand Particles");
    TV min_corner, max_corner;
    levelset.getBounds(min_corner, max_corner);
    PoissonDisk<T, dim> pd(/*random seed*/ 123, grain_size, min_corner, max_corner);
    StdVector<TV> grain_centroids;
    if (dim == 3)
        pd.sampleFromPeriodicData(grain_centroids, [&](const TV& x) { return levelset.inside(x); });
    else
        pd.sample(grain_centroids, [&](const TV& x) { return levelset.inside(x); });

    auto dpshn = AttributeName<DruckerPragerStvkHencky<T>>(DruckerPragerStvkHencky<T>::name());
    auto stvkhi = AttributeName<StvkWithHenckyIsotropic<T, dim>>(StvkWithHenckyIsotropic<T, dim>::name());
    T phi = T(0.5) * (1 - std::sqrt((T)5));
    RandomNumber<T> rand;

    for (auto iter = particles.subsetIter({ particle_range }, Particles<T, dim>::X_name(), Particles<T, dim>::mass_name(), dpshn, element_measure_name<T>(), stvkhi); iter; ++iter) {
        TV& x = iter.template get<0>();
        for (const TV& c : grain_centroids) {
            TV v = c - x;
            T alpha = std::min(v.norm() / ((T)1.2 * grain_size), (T)1);
            T scale = phi + 1 / (alpha - phi);
            scale *= rand.randReal(0.75, 1.0);
            x = x + scale * v;
        }
        // scale scale
        T scale_scale = (T)1;
        T scale_scale_min = (T)1;
        if (x(0) > 0.43 && x(1) > 1.05 && (x(2) > -0.2 && x(2) < -0.08) && rand.randReal(0, 1) < 0.5) {
            scale_scale = (T)0.1;
            scale_scale_min = 0.1;
        }
        else if (x(0) > 0.43 && (x(1) > 0.8 && x(1) < 1.05) && (x(2) > 0.1 && x(2) < 0.26) && rand.randReal(0, 1) < 0.5) {
            scale_scale = (T)0.14;
            scale_scale_min = 0.1;
        }
        else if (x(0) > 0.43 && (x(1) > 0.75 && x(1) < 0.905) && (x(2) > -0.26 && x(2) < -0.13) && rand.randReal(0, 1) < 0.7) {
            scale_scale = (T)0.24;
            scale_scale_min = 0.1;
        }

        T min_distance2 = std::numeric_limits<T>::max();
        for (const TV& c : grain_centroids)
            min_distance2 = std::min(min_distance2, (x - c).squaredNorm());
        T& mass = iter.template get<1>();
        DruckerPragerStvkHencky<T>& p = iter.template get<2>();
        T element_measure = iter.template get<3>();
        auto& cm_elastic = iter.template get<4>();
        T scale = 1 - std::min(std::sqrt(min_distance2) * rand.randReal(0.75, 1.25) / (2 * grain_size), (T)1);
        T mass_min = density_min * element_measure;
        mass = scale * (mass - mass_min) + mass_min;
        p.cohesion *= scale * scale_scale;
        cm_elastic.lambda *= std::max(scale * scale_scale, scale_min * scale_scale_min);
        cm_elastic.mu *= std::max(scale * scale_scale, scale_min * scale_scale_min);
    }
}

template <class T, int dim>
template <class TCONST>
void MpmParticleHandleBase<T, dim>::
    addFBasedMpmForceWithMeasure(const TCONST& model, const Range& range, T total_volume)
{
    FBasedMpmForceHelper<TCONST>& helper = mpmforce->getHelper(); // this will let mpmforce create a consitutive model helper
    particles.add(helper.constitutive_model_name(), range, model);
    if (total_volume != 0) {
        particles.add(element_measure_name<T>(), range, total_volume / range.length());
        particles.add(helper.F_name(), range, TM::Identity());
    }
    particles.add(helper.constitutive_model_scratch_name(), range, typename FBasedMpmForceHelper<TCONST>::Scratch());
}

template <class T, int dim>
template <class TCONST>
void MpmParticleHandleBase<T, dim>::
    addFElasticNonequilibratedBasedMpmForce(const TCONST& model, T viscosity_d_input, T viscosity_v_input)
{
    FElasticNonequilibratedBasedMpmForceHelper<TCONST>& helper = mpmforce->getHelper(); // this will let mpmforce create a consitutive model helper
    particles.add(helper.constitutive_model_name(), particle_range, model);
    particles.add(helper.F_name(), particle_range, TM::Identity());
    particles.add(helper.constitutive_model_scratch_name(), particle_range, typename FElasticNonequilibratedBasedMpmForceHelper<TCONST>::Scratch());
    helper.setParameters(viscosity_d_input, viscosity_v_input);
}

template <class T, int dim>
void MpmParticleHandleBase<T, dim>::setMassFromDensity(const T density)
{
    ZIRAN_INFO("Setting mass from densiy: total_volume = ", total_volume, ", particle.count = ", particle_range.length());
    T mp = density * total_volume / particle_range.length();
    for (auto iter = particles.subsetIter(DisjointRanges{ particle_range }, particles.mass_name()); iter; ++iter) {
        iter.template get<0>() = mp;
    }
}
} // namespace ZIRAN
