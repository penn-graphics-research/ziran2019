#ifndef DIFF_TEST_H
#define DIFF_TEST_H
#include <Ziran/CS/Util/PrettyPrinting.h>
#include <Ziran/CS/Util/RandomNumber.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <type_traits>
#include <tbb/tbb.h>

namespace ZIRAN {
/**
  Run difftest on a particular simulation
  This is templatized on objective
  You should use double whenever you run this test
  diff_test_perturbation_scale allows you to choose how to scale the perturbation from simple power of 2s
*/
template <class T, int dim, class Objective>
void runDiffTest(const int num_nodes, const Eigen::Matrix<T, dim, Eigen::Dynamic>& dv, Objective& objective, double diff_test_perturbation_scale = 1)
{
    ZIRAN_INFO("Running diff test with perturbation scale: ", diff_test_perturbation_scale);

    bool check_double = std::is_same<double, T>::value;
    ZIRAN_ASSERT(check_double, "Running difftest requires double");
    using TVStack = Eigen::Matrix<T, dim, Eigen::Dynamic>;
    TVStack x, step, f0, df0, f1, df1;
    x.resize(dim, num_nodes);
    step.resize(dim, num_nodes);
    f0.resize(dim, num_nodes);
    df0.resize(dim, num_nodes);
    f1.resize(dim, num_nodes);
    df1.resize(dim, num_nodes);
    x = dv;
    step.setZero();
    RandomNumber<T> rand(123);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_nodes),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(), i_end = range.end(); i < i_end; ++i) {
                for (int j = 0; j < dim; ++j)
                    step.col(i)(j) = rand.randReal(0, 1);
            }
        });
    step.normalize();

    objective.updateState(x);
    if (!objective.matrix_free) {
        objective.buildMatrix();
    }
    T e0 = objective.totalEnergy();
    objective.computeResidual(f0);
    objective.multiply(step, df0);

    const int DIFF_SIZE = 10;
    std::vector<T> energy_difference(DIFF_SIZE);
    std::vector<T> energy_differential(DIFF_SIZE);
    std::vector<T> energy_err(DIFF_SIZE);
    std::vector<T> energy_log_err(DIFF_SIZE);
    std::vector<T> force_difference_norm(DIFF_SIZE);
    std::vector<T> force_differential_norm(DIFF_SIZE);
    std::vector<T> force_err(DIFF_SIZE);
    std::vector<T> force_log_err(DIFF_SIZE);

    std::setprecision(12);
    std::cout << "e0\t=" << std::setw(20) << e0 << std::endl;

    for (int i = 1; i <= 10; ++i) {
        T h = diff_test_perturbation_scale * std::pow((T)(2), -i);
        objective.updateState(x + h * step);
        if (!objective.matrix_free) {
            objective.buildMatrix();
        }
        T e1 = objective.totalEnergy();
        std::cout << "e1\t=" << std::setw(20) << e1 << "\th = " << h << std::endl;
        objective.computeResidual(f1);
        T difference = (e0 - e1) / h;
        TVStack f2 = f0 + f1;
        T differential = tbb::parallel_reduce(tbb::blocked_range<int>(0, f2.cols(), 256), (T)0,
            [&](const tbb::blocked_range<int>& range, T ns) -> T {
                int start = range.begin();
                int length = (range.end() - range.begin());
                const auto& f2_block = f2.middleCols(start, length);
                const auto& step_block = step.middleCols(start, length);
                return ns + (f2_block.array() * step_block.array()).sum();
            },
            [](T a, T b) -> T { return a + b; });
        differential /= 2;

        T err = (difference - differential);
        T log_err = std::log(std::abs(err));
        energy_difference[i - 1] = difference;
        energy_differential[i - 1] = differential;
        energy_err[i - 1] = err;
        energy_log_err[i - 1] = log_err;

        objective.multiply(step, df1);
        TVStack force_difference = (1 / h) * (f0 - f1);
        TVStack force_differential = 0.5 * (df0 + df1);
        T err_force = (force_difference - force_differential).norm();
        T log_err_force = std::log(std::abs(err_force));
        force_difference_norm[i - 1] = force_difference.norm();
        force_differential_norm[i - 1] = force_differential.norm();
        force_err[i - 1] = err_force;
        force_log_err[i - 1] = log_err_force;
    }

    std::cout << std::setprecision(12) << "energy["
              << "i"
              << "] = " << std::setw(20) << "difference"
              << std::setw(20) << "differential"
              << std::setw(20) << "err"
              << std::setw(20) << "log_err"
              << std::endl;
    for (int i = 0; i < DIFF_SIZE; ++i) {
        std::cout << std::setprecision(12) << "energy[" << i << "] = " << std::setw(20) << energy_difference[i]
                  << std::setw(20) << energy_differential[i]
                  << std::setw(20) << energy_err[i]
                  << std::setw(20) << energy_log_err[i]
                  << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::setprecision(12) << "force["
              << "i"
              << "] = " << std::setw(20) << "difference_norm"
              << std::setw(20) << "differential_norm"
              << std::setw(20) << "err"
              << std::setw(20) << "log_err"
              << std::endl;
    for (int i = 0; i < DIFF_SIZE; ++i) {
        std::cout << std::setprecision(12) << "force[" << i << "] = " << std::setw(20) << force_difference_norm[i]
                  << std::setw(20) << force_differential_norm[i]
                  << std::setw(20) << force_err[i]
                  << std::setw(20) << force_log_err[i]
                  << std::endl;
    }
    std::cin.get();
}
} // namespace ZIRAN
#endif
