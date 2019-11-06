#ifndef GLOBAL_SYSTEM_H
#define GLOBAL_SYSTEM_H

#include <Ziran/Math/Linear/DenseExt.h>

namespace ZIRAN {

template <class Simulation, int dim>
class GlobalSystem {
public:
    using T = typename Simulation::Scalar;
    using TVStack = Matrix<T, dim, Eigen::Dynamic>;

    std::function<void(TVStack&)> project;
    std::function<void(const TVStack&, TVStack&)> precondition;
    std::function<void(const TVStack&, TVStack&)> multiply;
    std::function<T(const TVStack&)> residual;

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

    template <class Func>
    void setResidual(Func residual_func)
    {
        residual = residual_func;
    }
};
} // namespace ZIRAN

#endif
