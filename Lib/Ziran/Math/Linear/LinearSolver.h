#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H
#include <Ziran/CS/Util/ErrorContext.h>
#include <Ziran/Math/Linear/KrylovSolvers.h>

namespace ZIRAN {
template <class T, class TM, class TV>
class LinearSolver {
    using TM2 = Matrix<T, 2, 2>;
    using TV2 = Vector<T, 2>;

public:
    T tolerance;
    int max_iterations;
    T relative_tolerance;

    LinearSolver(const int max_it_input)
        : max_iterations(max_it_input)
        , relative_tolerance(1)
    {
        setTolerance();
    }

    virtual ~LinearSolver() {}

    virtual void setTolerance(const T tolerance_input = 16 * std::numeric_limits<T>::epsilon()) { tolerance = tolerance_input; }

    virtual void setRelativeTolerance(const T tolerance_input = 1) { relative_tolerance = tolerance_input; }

    virtual void setMaxIteration(const int max_it_input) { max_iterations = max_it_input; }

    virtual void reinitialize(const TV& b) = 0;

    virtual T dotProduct(const TV& A, const TV& B)
    {
        return (A.array() * B.array()).sum();
    }

    virtual int solve(const TM& A, TV& x, const TV& b, const bool verbose = false) = 0;
};
} // namespace ZIRAN

#endif
