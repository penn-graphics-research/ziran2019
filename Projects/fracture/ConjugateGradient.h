#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include <Ziran/Math/Linear/KrylovSolvers.h>
#include <tbb/tbb.h>

namespace ZIRAN {
template <class T, class TM, class TV>

//
// TODO: This CG seems to be buggy! Fix it using PhysBAM reference.
//
class ConjugateGradient {

    /** All notations adopted from Wikipedia,
         * q denotes A*p in general */
    TV r, p, q, temp;
    TV mr, s;

public:
    T tolerance;
    int max_iterations;

    ConjugateGradient(const int max_it_input)
        : max_iterations(max_it_input)
    {
        setTolerance(std::is_same<T, float>::value ? (T)1e-6 : (T)1e-12);
    }

    ~ConjugateGradient() {}

    void setTolerance(const T tolerance_input = 16 * std::numeric_limits<T>::epsilon()) { tolerance = tolerance_input; }

    void reinitialize(const TV& b)
    {
        r.resizeLike(b);
        p.resizeLike(b);
        q.resizeLike(b);
        temp.resizeLike(b);

        mr.resizeLike(b);
        s.resizeLike(b);
    }

    T dotProduct(const TV& A, const TV& B)
    {
        return (A.array() * B.array()).sum();
    }

    int solve(const TM& A, TV& x, const TV& b, const bool verbose = false)
    {
        reinitialize(x);

#if 1
        // s is search direction
        //system.Set_Boundary_Conditions(x);
        T rho_old = (T)std::numeric_limits<T>::max();
        T convergence_norm = 0;

        r = b;
        A.multiply(x, q); // q=A*x
        r -= q; // r=b-A*x
        A.project(r);
        convergence_norm = std::sqrt(r.squaredNorm());
        T relative_tolerance = convergence_norm * (std::is_same<T, float>::value ? (T)1e-6 : (T)1e-12);
        ZIRAN_INFO("CG starts with convergence_norm = ", convergence_norm);
        int iterations;
        for (iterations = 0;; iterations++) {
            // apply preconditioner
            A.precondition(r, mr);

            // update the search direction s
            T rho = dotProduct(mr, r);
            if (iterations == 0)
                s = mr;
            else
                s = rho / rho_old * s + mr;

            // update x
            A.multiply(s, q);
            A.project(q);
            T s_dot_q = dotProduct(s, q);
            T alpha = s_dot_q ? rho / s_dot_q : (T)std::numeric_limits<T>::max();
            x = alpha * s + x;

            // recompute residual and copy rho_old
            r = -alpha * q + r;
            rho_old = rho;

            // terminate if tolerance is reached
            convergence_norm = std::sqrt(r.squaredNorm());
            ZIRAN_INFO("CG convergence_norm = ", convergence_norm);
            if (convergence_norm < tolerance || convergence_norm < relative_tolerance || iterations == max_iterations) {
                ZIRAN_INFO("CG iterattion ends at ", iterations, " and convergence_norm = ", convergence_norm);
                return iterations;
            }
        }
#else
        ZIRAN_QUIET_TIMER();
        assert(x.size() == b.size());
        reinitialize(b);

        // r = M * (b - A * x) --with assigned dof zeroed out
        A.multiply(x, temp);
        r = b - temp;
        A.project(temp);
        A.precondition(temp, r);

        T r_dot_r = r.squaredNorm(), r_dot_r_new;
        T r_norm = std::sqrt(r_dot_r);
        if (r_norm < tolerance) {
            ZIRAN_VERB_IF(verbose, "Iteration = ", 0);
            ZIRAN_VERB_IF(verbose, "Two norm of the residual = ", r_norm);
            return 0;
        }

        p = r;
        // q = M * A * q;
        A.multiply(p, temp);
        A.precondition(temp, q);

        // alpha = |r|^2 / (p^T * A * p)
        T alpha = r_dot_r / p.dot(q), beta;

        for (int k = 1; k < max_iterations; k++) {
            x += alpha * p;
            r -= alpha * q;

            // zero out the dofs of r
            A.project(r);

            r_dot_r_new = r.squaredNorm();
            r_norm = std::sqrt(r_dot_r_new);

            if (r_norm < tolerance) {
                ZIRAN_VERB("ConjugateGradient iterations ", k);
                return k;

                beta = r_dot_r_new / r_dot_r;
                r_dot_r = r_dot_r_new;
                p = r + beta * p;

                // q = M * A * q;
                A.multiply(p, temp);
                A.precondition(temp, q);

                alpha = r_dot_r / p.dot(q);
            }

            q.setZero();
            r.setZero();
        }
        ZIRAN_VERB_IF(verbose, "ConjugateGradient max iterations reached ", max_iterations);
        return max_iterations;
#endif
    }
};
} // namespace ZIRAN

#endif
