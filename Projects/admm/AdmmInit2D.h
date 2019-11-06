#ifndef ADMM_INIT_2D_H
#define ADMM_INIT_2D_H

#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/Math/MathTools.h>

#include <Ziran/Sim/MeshHandle.h>
#include <Ziran/Sim/SceneInitializationCore.h>

#include "AdmmSimulation.h"
#include "AdmmInit.h"
#include <float.h>

namespace ZIRAN {

template <class T>
class AdmmInit2D : public AdmmInitBase<T, 2> {
public:
    static const int dim = 2;
    using Base = AdmmInitBase<T, dim>;
    using TV = Vector<T, dim>;
    using TVI = Vector<int, dim>;

    using Base::init_helper;
    using Base::scene;
    using Base::sim;
    using Base::test_number;
    AdmmInit2D(AdmmSimulation<T, dim>& sim, const int test_number)
        : Base(sim, test_number)
    {
    }

    void reload() override
    {
    }
};
} // namespace ZIRAN
#endif
