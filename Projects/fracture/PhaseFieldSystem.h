#ifndef PHASE_FIELD_SYSTEM_H
#define PHASE_FIELD_SYSTEM_H

#include <Ziran/Math/Linear/DenseExt.h>
#include <Ziran/Math/Linear/KrylovSolvers.h>
#include <Ziran/Math/Linear/Minres.h>
#include <Ziran/Math/MathTools.h>
#include <Ziran/Sim/SimulationBase.h>
#include <Ziran/Math/Linear/DenseExt.h>
#include <Ziran/Math/Linear/EigenSparseLU.h>

#include <tbb/tbb.h>

namespace ZIRAN {

/**
   Simulation owns Objective and NewtonsMethod.
   This will allow options to do Newton.
   It stores the parameters for them. But it is a little weird.
   This is supposed to be used by both MPM and Elasticity.
 */
template <class Simulation>
class PhaseFieldSystem {
public:
    using T = typename Simulation::Scalar;
    static const int dim = Simulation::dim;
    using Objective = PhaseFieldSystem<Simulation>;
    using TV = Vector<T, dim>;
    using Vec = Vector<T, Eigen::Dynamic>;

    using Scalar = T;

    Simulation& simulation;

    std::function<void(Vec&)> project;
    std::function<void(const Vec&, Vec&)> precondition;
    std::function<void(const Vec&, Vec&)> multiply;

    PhaseFieldSystem(Simulation& simulation)
        : simulation(simulation)
        , project([](Vec&) {})
        , precondition([](const Vec& x, Vec& b) { b = x; })
    {
    }

    PhaseFieldSystem(const PhaseFieldSystem& objective) = delete;

    ~PhaseFieldSystem() {}

    template <class Projection>
    void initialize(const Projection& projection_function)
    {
        ZIRAN_INFO("PhaseFieldSystem only uses matrix free");
        project = projection_function;
    }

    void reinitialize()
    {
    }

    template <class Func>
    void setMultiplier(Func multiply_func)
    {
        multiply = multiply_func;
    }

    template <class Func>
    void setProjection(Func project_func)
    {
        project = project_func;
    }

    template <class Func>
    void setPreconditioner(Func preconditioner_func)
    {
        precondition = preconditioner_func;
    }
};
} // namespace ZIRAN

#endif
