#ifndef GENERALIZED_MINIMAL_RESIDUAL_H
#define GENERALIZED_MINIMAL_RESIDUAL_H
#include <Ziran/CS/Util/ErrorContext.h>
#include <Ziran/Math/Linear/KrylovSolvers.h>

#include <tbb/tbb.h>

namespace ZIRAN {
template <class T, class TM, class TV>
class GeneralizedMinimalResidual {
    using TM2 = Matrix<T, 2, 2>;
    using TV2 = Vector<T, 2>;

public:
    T tolerance;
    int max_iterations;
    T relative_tolerance;
    std::vector<TV> av;

    GeneralizedMinimalResidual(const int max_it_input)
        : max_iterations(max_it_input)
        , relative_tolerance(1)
    {
        setTolerance();
    }

    ~GeneralizedMinimalResidual() {}

    void setTolerance(const T tolerance_input = 16 * std::numeric_limits<T>::epsilon()) { tolerance = tolerance_input; }

    void setRelativeTolerance(const T tolerance_input = 1) { relative_tolerance = tolerance_input; }

    void setMaxIteration(const int max_it_input) { max_iterations = max_it_input; }

    void reinitialize(const TV& b)
    {
        av.clear();
        for (int i = 0; i < max_iterations + 5; ++i)
            av.push_back(b);
    }

    T dotProduct(const TV& A, const TV& B)
    {
        return (A.array() * B.array()).sum();
    }

    TV2 Givens_Transpose_Times(const TV2& u, const TV2& v) const
    {
        return TV2(u[0] * v[0] + u[1] * v[1], u[0] * v[1] - u[1] * v[0]);
    }

    int solve(const TM& A, TV& x, const TV& _b, const bool verbose = false)
    {
        ZIRAN_QUIET_TIMER();
        TV b = _b;
        A.project(b);

        assert(x.size() == b.size());
        reinitialize(b);

        int marker = 1 + 1;
        int next_vector = marker;
        TV& temp_1 = av[0];
        TV& temp_2 = av[1];

        A.precondition(b, temp_2);
        TV prec_b = temp_2;
        T residual = std::sqrt(dotProduct(prec_b, prec_b));
        T local_tolerance = std::min(relative_tolerance * residual, tolerance);

        T phi;
        std::vector<T> t;
        std::vector<std::vector<T>> h, r_2;
        std::vector<TV2> rot;
        int v_cur = 0;

        for (int k = 0, system_size;; k++, system_size++) {
            if (k == 0) {
                system_size = 1;
                next_vector = marker;
                v_cur = next_vector++;
                h.clear();
                r_2.clear();
                rot.clear();
                t.clear();
                {
                    TV temp_3 = x;
                    A.project(temp_3);
                    A.multiply(temp_3, temp_1);
                    A.project(temp_1);
                }
                temp_1 = b - temp_1;
                A.precondition(temp_1, temp_2);
                TV prec = temp_2;
                phi = std::sqrt(dotProduct(prec, prec));
                av[v_cur] = ((T)1 / phi) * prec;
            }

            {
                TV temp_3 = av[v_cur];
                A.project(temp_3);
                A.multiply(temp_3, temp_1);
                A.project(temp_1);
            }
            A.precondition(temp_1, temp_2);
            temp_1 = temp_2;

            r_2.push_back(std::vector<T>());
            std::vector<T>& r_2_cur = r_2[r_2.size() - 1];
            h.push_back(std::vector<T>());
            std::vector<T>& h_cur = h[h.size() - 1];
            assert(r_2_cur.empty());
            assert(h_cur.empty());
            for (int i = 0; i < system_size; ++i) {
                T ip = dotProduct(temp_1, av[marker + i]);
                h_cur.push_back(ip);
                temp_1 += (-ip) * av[marker + i];
            }

            int v_next = next_vector++;
            T ip = std::sqrt(dotProduct(temp_1, temp_1));
            h_cur.push_back(ip);
            av[v_next] = ((T)1 / ip) * temp_1;

            T r = h_cur[0];
            for (int i = 0; i < system_size - 1; ++i) {
                TV2 tmp = Givens_Transpose_Times(rot[i], TV2(r, h_cur[i + 1]));
                r_2_cur.push_back(tmp[0]);
                r = tmp[1];
            }
            TV2 tmp(r, h_cur[system_size]);
            r_2_cur.push_back(tmp.norm());
            rot.push_back(tmp.normalized());

            TV2 ph = phi * rot[system_size - 1];
            t.push_back(ph[0]);
            phi = -ph[1];

            residual = abs(phi);
            if (residual <= local_tolerance || k == max_iterations) {
                for (int i = (int)r_2.size() - 1; i >= 0; --i) {
                    T rhs = t[i];
                    for (int j = i + 1; j < (int)r_2.size(); ++j)
                        rhs -= r_2[j][i] * t[j];
                    t[i] = rhs / r_2[i][i];
                    x += t[i] * av[marker + i];
                }
                A.project(x);
                if (residual <= local_tolerance) {
                    ZIRAN_VERB_IF(true, "GMRES terminates at ", k + 1, "; residual = ", residual);
                    return k; //Output the number of iterations.
                }
                else {
                    ZIRAN_VERB_IF(true, "GMRES terminates at ", max_iterations + 1, " (max reached); residual = ", residual);
                    return max_iterations;
                }
            }
            if ((k + 1) % 50 == 0 || k == 0) {
                ZIRAN_INFO("\tGMRES iter ", k + 1, "; residual = ", residual);
            }

            v_cur = v_next;
        }
    }
};
} // namespace ZIRAN

#endif
