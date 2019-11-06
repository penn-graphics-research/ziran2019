#ifndef MPM_SIMULATION_H
#define MPM_SIMULATION_H

#include <MPM/MpmSimulationBase.h>

namespace ZIRAN {

template <class T, int _dim>
class MpmSimulation : public MpmSimulationBase<T, _dim> {
public:
    using Base = MpmSimulationBase<T, _dim>;

    MpmSimulation()
        : Base()
    {
    }

    const char* name() override { return "mpm"; }
};
} // namespace ZIRAN

#endif
