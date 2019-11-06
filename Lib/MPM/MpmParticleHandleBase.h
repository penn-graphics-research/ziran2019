#ifndef MPM_PARTICLE_HANDLE_BASE_H
#define MPM_PARTICLE_HANDLE_BASE_H

#include <functional>
#include <memory>
#include <string>
#include <Ziran/Math/Geometry/Rotation.h>
#include <Ziran/CS/DataStructure/DisjointRanges.h>
#include <Ziran/CS/DataStructure/Ref.h>
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/Util/AttributeNamesForward.h>
#include <Ziran/CS/Util/AttributeNamesForward.h>
#include <Ziran/Physics/ConstitutiveModel/ConstitutiveModel.h>
#include <MPM/Force/MpmForceBase.h>

namespace ZIRAN {

template <class T, int dim>
class AnalyticLevelSet;
template <class T, int dim>
class Particles;
class PlasticityApplierBase;
template <class T, int dim>
class Scene;
template <class T, int dim>
class MpmSimulation;
template <class T, int dim>
class AnalyticLevelSet;
template <class T, int dim>
class Particles;
template <class T, int dim>
class Scene;

template <class T, int dim>
class MpmParticleHandleBase {
public:
    static constexpr int interpolation_degree = ZIRAN_MPM_DEGREE;
    static const int splat_size = MATH_TOOLS::power(interpolation_degree + 1, dim);

    using TV = Vector<T, dim>;
    typedef Matrix<T, dim, dim> TM;

    Particles<T, dim>& particles;
    Scene<T, dim>& scene;
    MpmForceBase<T, dim>* mpmforce;
    StdVector<std::unique_ptr<PlasticityApplierBase>>& plasticity_appliers;
    StdVector<TV>& scratch_xp;
    T& dt;
    Range particle_range;
    T total_volume;
    int cotangent_manifold_dim;

    MpmParticleHandleBase(Particles<T, dim>& particles, Scene<T, dim>& scene, MpmForceBase<T, dim>* mpmforce,
        StdVector<std::unique_ptr<PlasticityApplierBase>>& plasticity_appliers,
        StdVector<TV>& scratch_xp, T& dt, Range particle_range, T total_volume, int cotangent_manifold_dim = 0);

    // Creates a copy with new particles
    MpmParticleHandleBase copy();

    void transform(const std::function<void(int, Ref<T>, Vector<T, dim>&, Vector<T, dim>&)>& mapping);

    void addVolumeFraction(const T b);

    template <class TCONST>
    void addFBasedMpmForce(const TCONST& model);

    template <class TCONST>
    void addFBasedMpmForceWithPhaseField(const T& percentage, const T& l0, const TCONST& model, bool allow_damage = true, const T random_fluctuation_percentage = 0);

    void addOriginalPositionAsAttribute();

    void randomizePhaseFieldSigmaC(const T random_fluctuation_percentage);

    //void voronoiSigmaC(std::string voronoiPointsFile, T radius, T minPercent, T maxPercent);
    void voronoiSigmaC(std::string voronoiPointsFile, T radius);

    void yDimSigmaCInit(T yMin, T yMax, T maxScale);

    void scaleF(const T scale);

    void scaleF(const TV scale);

    void scaleJ(const T scale);

    void scaleFCurve(int frame, const std::function<T(int)>& growCurve);

    void prescoreWetSand2D(AnalyticLevelSet<T, dim>& levelset, T grain_size, T density_min, T Jp_min, T scale_min = (T)0.5);

    void prescorePhaseFieldSigmaC(AnalyticLevelSet<T, dim>& levelset, T grain_size, T scale_lower_clamp = 0);

    void resetDeformation();

    void addElementMeasure();

    template <class TConst, class TPlastic>
    void addPlasticity(const TConst& cons, const TPlastic& plasticity, std::string strain_name = "F");

    void addDummyPlasticity();

    template <class TCONST>
    void addFBasedMpmForceWithMeasure(const TCONST& model, const Range& range, T total_volume);

    template <class TCONST>
    void addFElasticNonequilibratedBasedMpmForce(const TCONST& model, T viscosity_d_input, T viscosity_v_input);

    void setMassFromDensity(const T density);
};
} // namespace ZIRAN
#endif /* ifndef MPM_PARTICLE_HANDLE */
