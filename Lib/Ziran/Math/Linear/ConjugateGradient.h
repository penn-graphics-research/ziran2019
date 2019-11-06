#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include <Ziran/Math/Linear/KrylovSolvers.h>
#include <tbb/tbb.h>

#include "LinearSolver.h"

namespace ZIRAN {
template <class T, class TM, class TV>

//
// TODO: This CG seems to be buggy! Fix it using PhysBAM reference.
//
class ConjugateGradient : public LinearSolver<T, TM, TV> {

    using Base = LinearSolver<T, TM, TV>;

    /** All notations adopted from Wikipedia, 
     * q denotes A*p in general */
    TV r, p, q, temp;

public:
    ConjugateGradient(const int max_it_input)
        : Base(max_it_input)
    {
    }

    ~ConjugateGradient() {}

    void reinitialize(const TV& b)
    {
        r.resizeLike(b);
        p.resizeLike(b);
        q.resizeLike(b);
        temp.resizeLike(b);
    }

    int solve_(const TM& A, TV& x, const TV& b, const bool verbose = false)
    {
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
        if (r_norm < Base::tolerance) {
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

        for (int k = 1; k < Base::max_iterations; k++) {
            x += alpha * p;
            r -= alpha * q;

            // zero out the dofs of r
            A.project(r);

            r_dot_r_new = r.squaredNorm();
            r_norm = std::sqrt(r_dot_r_new);

            if (r_norm < Base::tolerance) {
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
        ZIRAN_VERB_IF(verbose, "ConjugateGradient max iterations reached ", Base::max_iterations);
        return Base::max_iterations;
    }

    int solve(const TM& A, TV& x, const TV& b, const bool verbose = false)
    {
        //TODO: adaptive tolerance on unpreconditioned residual norm

        ZIRAN_QUIET_TIMER();
        assert(x.size() == b.size());
        reinitialize(b);
        int cnt = 0;
        T alpha, beta, residual_preconditioned_norm, zTrk, zTrk_last;

        //NOTE: requires that the input x has been projected
        A.multiply(x, temp);
        r = b - temp;
        A.project(r);
        A.precondition(r, q); //NOTE: requires that preconditioning matrix is projected
        p = q;

        zTrk = Base::dotProduct(r, q);
        residual_preconditioned_norm = std::sqrt(zTrk);
        T local_tolerance = std::min(Base::relative_tolerance * residual_preconditioned_norm, Base::tolerance);
        for (cnt = 0; cnt < Base::max_iterations; ++cnt) {
            if (residual_preconditioned_norm < local_tolerance) {
                ZIRAN_VERB_IF(verbose, "\tCG terminates at ", cnt, "; (preconditioned norm) residual = ", residual_preconditioned_norm);
                return cnt;
            }

            if (cnt % 50 == 0) {
                ZIRAN_VERB_IF(verbose, "\tCG iter ", cnt, "; (preconditioned norm) residual = ", residual_preconditioned_norm);
            }

            A.multiply(p, temp);
            A.project(temp);
            alpha = zTrk / Base::dotProduct(temp, p);

            x += p * alpha;
            r -= temp * alpha;
            A.precondition(r, q); //NOTE: requires that preconditioning matrix is projected

            zTrk_last = zTrk;
            zTrk = Base::dotProduct(q, r);
            beta = zTrk / zTrk_last;

            p = q + beta * p;

            residual_preconditioned_norm = std::sqrt(zTrk);
        }
        ZIRAN_VERB_IF(verbose, "ConjugateGradient max iterations reached ", Base::max_iterations);
        return Base::max_iterations;
    }
};
} // namespace ZIRAN

#endif
