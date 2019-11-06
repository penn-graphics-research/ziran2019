#ifndef MPM_INIT_2D_H
#define MPM_INIT_2D_H

#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/CS/Util/AttributeNamesForward.h>
#include <Ziran/Math/Geometry/PartioIO.h>
#include <Ziran/Math/Geometry/MeshConstruction.h>
#include <Ziran/Math/Geometry/AnalyticLevelSet.h>
#include <Ziran/Math/Geometry/CollisionObject.h>
#include <Ziran/Math/MathTools.h>
#include <Ziran/Physics/SoundSpeedCfl.h>

#include <Ziran/Sim/MeshHandle.h>
#include <Ziran/Sim/SceneInitializationCore.h>

#include <MPM/MpmInitializationHelper.h>

#include "FractureSimulation.h"
#include "FractureInit.h"
#include <float.h>

namespace Minchen {
extern double nu, E, rho, v, dx, M, sigmaPercent, bulletFriction;
extern bool usePhaseField;
} // namespace Minchen

namespace ZIRAN {

template <class T>
class FractureInit2D : public FractureInitBase<T, 2> {
public:
    static const int dim = 2;
    using Base = FractureInitBase<T, dim>;
    using TV = Vector<T, dim>;
    using TVI = Vector<int, dim>;

    using Base::init_helper;
    using Base::scene;
    using Base::sim;
    using Base::test_number;
    FractureInit2D(FractureSimulation<T, dim>& sim, const int test_number)
        : Base(sim, test_number)
    {
    }

    void reload() override
    {
    }
};
} // namespace ZIRAN
#endif
