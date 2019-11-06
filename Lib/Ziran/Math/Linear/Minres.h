#ifndef MINRES_H
#define MINRES_H
#include <Ziran/CS/Util/ErrorContext.h>
#include <Ziran/Math/Linear/KrylovSolvers.h>

#include <tbb/tbb.h>

#include "LinearSolver.h"

namespace ZIRAN {
template <class T, class TM, class TV>
class Minres : public LinearSolver<T, TM, TV> {

    using Base = LinearSolver<T, TM, TV>;

    using TM2 = Matrix<T, 2, 2>;
    using TV2 = Vector<T, 2>;

    GivensRotation<T> Gk, Gkm1, Gkm2; //The Q in the QR of Hk: Givens rotations

    T gamma, delta, epsilon; //These are the newest entries in R from the QR of Hk
    T beta_kp1, alpha_k, beta_k, sk; //This is the last column in Hk
    TV mk, mkm1, mkm2; //mk is the newest search direction in the memory friendly basis for the Krylov space, you need the other two to generate mk
    TV z, qkp1, qk, qkm1; //These are the Lanczos vectors needed at each iteration, qk is denoted as vk in the notes
    T tk; //This is the step length in the direction of mk
    TV2 last_two_components_of_givens_transformed_least_squares_rhs; //This will track the residual with just a constant number of flops per iteration
    T rhsNorm2;

public:
    Minres(const int max_it_input)
        : Base(max_it_input)
        , Gk(0, 1)
        , Gkm1(0, 1)
        , Gkm2(0, 1)
    {
    }

    ~Minres() {}

    void reinitialize(const TV& b)
    {
        mk.resizeLike(b);
        mkm1.resizeLike(b);
        mkm2.resizeLike(b);
        z.resizeLike(b);
        qkm1.resizeLike(b);
        qk.resizeLike(b);
        qkp1.resizeLike(b);
        mk.setZero();
        mkm1.setZero();
        mkm2.setZero();
        qkm1.setZero();
        qk.setZero();
        qkp1.setZero();
        gamma = 0;
        delta = 0;
        epsilon = 0;
        beta_kp1 = 0;
        alpha_k = 0;
        beta_k = 0;
        tk = 0;
        Gk.setIdentity();
        Gkm1.setIdentity();
        Gkm2.setIdentity();

        rhsNorm2 = b.squaredNorm();
    }

    int solve(const TM& A, TV& x, const TV& b, const bool verbose = false)
    {
        ZIRAN_QUIET_TIMER();
        assert(x.size() == b.size());
        reinitialize(b);

        //qkp1 = b - A * x;
        A.multiply(x, qkp1);
        qkp1 = b - qkp1;

        A.project(qkp1);
        A.precondition(qkp1, z); //z= M_inv*qkp1

        T residual_preconditioned_norm = std::sqrt(Base::dotProduct(z, qkp1));
        beta_kp1 = residual_preconditioned_norm;
        ZIRAN_ASSERT(beta_kp1 == beta_kp1); //check to see that beta_kp1 is not NAN
        T local_tolerance = std::min(Base::relative_tolerance * residual_preconditioned_norm, Base::tolerance);

        //check for convergence on initial guess
        if (residual_preconditioned_norm < local_tolerance) {
            ZIRAN_VERB_IF(verbose, "\tMinres terminates at ", 0, "; (preconditioned norm) residual = ", residual_preconditioned_norm);
            return 0; //Output the number of iterations.
        }

        if (residual_preconditioned_norm > 0) {
            qkp1 /= beta_kp1;
            z /= beta_kp1;
        }
        last_two_components_of_givens_transformed_least_squares_rhs = TV2(residual_preconditioned_norm, 0);

        ZIRAN_VERB_IF(verbose, "\tMinres iter ", 0, "; (preconditioned norm) residual = ", residual_preconditioned_norm);
        for (int k = 0; k < Base::max_iterations; k++) {
            ZIRAN_CONTEXT(residual_preconditioned_norm);
            ZIRAN_CONTEXT("Minres step", k);
            if (residual_preconditioned_norm < local_tolerance) {
                ZIRAN_VERB_IF(verbose, "\tMinres terminates at ", k, "; (preconditioned norm) residual = ", residual_preconditioned_norm);
                return k; //Output the number of iterations.
            }

            //use mk to store zk to save storage
            mkm2.swap(mkm1);
            mkm1.swap(mk);
            mk = z;

            beta_k = beta_kp1;

            // Save the last two Lanczos vectors.
            // Equivalent to qkm1 = qk; qk = qkp1; And qkp1 will be overwritten right after anyway.
            qkm1.swap(qkp1);
            qkm1.swap(qk);

            A.multiply(mk, qkp1); //qpk1 = A * zk; Get the important part of the next Lanczos vector: q_k+1
            A.project(qkp1);
            alpha_k = Base::dotProduct(mk, qkp1);

            qkp1 -= alpha_k * qk;
            qkp1 -= beta_k * qkm1;

            A.precondition(qkp1, z);
            beta_kp1 = std::sqrt(std::max((T)0, Base::dotProduct(z, qkp1)));

            if (beta_kp1 > 0) {
                qkp1 /= beta_kp1;
                z /= beta_kp1;
            }

            residual_preconditioned_norm = applyAllPreviousGivensRotationsAndDetermineNewGivens(); //This determines the newest Givens rotation and applies the previous two where appropriate

            //Three term recurence for the m's, mk stores the old zk which is set at beginning of the iteration
            mk = (mk - delta * mkm1 - epsilon * mkm2) / gamma;
            x += tk * mk;

            if ((k + 1) % 50 == 0) {
                ZIRAN_VERB_IF(verbose, "\tMinres iter ", k + 1, "; (preconditioned norm) residual = ",
                    residual_preconditioned_norm);
            }
        }

        ZIRAN_VERB_IF(verbose, "\tMinres terminates at ", Base::max_iterations, " (max reached); (preconditioned norm) residual = ", residual_preconditioned_norm);
        return Base::max_iterations;
    }

    T applyAllPreviousGivensRotationsAndDetermineNewGivens()
    {
        //QR the LHS: gamma, delta, epsilon
        Gkm2 = Gkm1;
        Gkm1 = Gk;
        TV2 epsilon_k_and_phi_k(0, beta_k);
        Gkm2.rowRotation(epsilon_k_and_phi_k);

        epsilon = epsilon_k_and_phi_k(0);
        TV2 delta_k_and_zsi_k(epsilon_k_and_phi_k(1), alpha_k);
        Gkm1.rowRotation(delta_k_and_zsi_k);
        delta = delta_k_and_zsi_k(0);
        TV2 temp(delta_k_and_zsi_k(1), beta_kp1);

        Gk.compute(temp(0), temp(1));
        Gk.rowRotation(temp);
        gamma = temp(0);

        //Now deal with the RHS: tk and residual (two norm)
        Gk.rowRotation(last_two_components_of_givens_transformed_least_squares_rhs);
        tk = last_two_components_of_givens_transformed_least_squares_rhs(0);
        T residual = last_two_components_of_givens_transformed_least_squares_rhs(1); //This is the two norm of the residual.
        last_two_components_of_givens_transformed_least_squares_rhs = TV2(residual, 0); //Set up for the next iteration
        if (residual < 0)
            return -residual;
        else
            return residual;
    }
};
} // namespace ZIRAN

#endif
