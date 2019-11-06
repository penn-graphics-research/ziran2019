#include <Ziran/CS/Util/FloatingPointExceptions.h>
#include <Ziran/CS/Util/SignalHandler.h>
#include <Ziran/CS/Util/CommandLineFlags.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/Filesystem.h>
#include <Ziran/CS/Util/PluginManager.h>
#include <tbb/tbb.h>
#include "MpmInit.h"

using namespace ZIRAN;

#define T double
#define dim 3

namespace Minchen {
double p_Jp = 0.95, p_beta = 3, p_nu = 0.39, p_E = 2000, p_vy = -5;
double p_Jp2 = 0.97, p_beta2 = 0.8, p_nu2 = 0.35, p_E2 = 1000, p_rho2 = 2;
bool p_noPlasticity = false;

double v_mu = 0.001;
double v_xxxx = 0.005;

int lin_solver_type;
} // namespace Minchen

int main(int argc, char* argv[])
{
    {
        bool displayHelp = false;
        int test_number = -1; // Non-lua option.
        bool three_d = false; // Non-lua option.
        bool use_double = false; // Non-lua option
        int restart = 0; // Non-lua option
        bool run_diff_test = false; // Non-lua option
        double diff_test_perturbation_scale = 1; // Non-lua option

        // Not checking for nan, because when constitutive model returns that, MpmForceBase is skipping them (treating as zeros)
        // FPE::WatchedScope w(FPE::Mask::Overflow | FPE::Mask::DivZero);
        // Unconmment the following to catch division by 0
        // FPE::WatchedScope w(FPE::Mask::Overflow | FPE::Mask::Invalid | FPE::Mask::DivZero);
        FPE::WatchedScope w(FPE::Mask::Invalid);

        FLAGS::Register helpflag("--help", "Print help (this message) and exit", displayHelp);

        // Non-lua command line options
        FLAGS::Register test_number_flag("-test", "Test number (non-lua test)", test_number);
        FLAGS::Register three_d_flag("--3d", "Dimension is 3(non-lua test)", three_d);
        FLAGS::Register run_diff_test_flag("--run_diff_test", "Run diff test (non-lua test)", run_diff_test);
        FLAGS::Register diff_test_perturbation_scale_flag("-dtps", "diff_test_perturbation_scale (non-lua test)", diff_test_perturbation_scale);
        FLAGS::Register double_flag("--double", "Dimension (non-lua test)", use_double);
        FLAGS::Register restart_flag("-restart", "Restart frame (non-lua test)", restart);
        FLAGS::Register noPlasticity_flag("--noPlasticity", "pure elasticity", Minchen::p_noPlasticity);

        FLAGS::Register p_Jp_flag("-Jp", "Jp", Minchen::p_Jp);
        FLAGS::Register p_beta_flag("-beta", "beta", Minchen::p_beta);
        FLAGS::Register p_nu_flag("-nu", "nu", Minchen::p_nu);
        FLAGS::Register p_E_flag("-E", "E", Minchen::p_E);
        FLAGS::Register p_vy_flag("-vy", "vy", Minchen::p_vy);

        FLAGS::Register p_Jp2_flag("-Jp2", "Jp2", Minchen::p_Jp2);
        FLAGS::Register p_beta2_flag("-beta2", "beta2", Minchen::p_beta2);
        FLAGS::Register p_nu2_flag("-nu2", "nu2", Minchen::p_nu2);
        FLAGS::Register p_E2_flag("-E2", "E2", Minchen::p_E2);
        FLAGS::Register p_rho2_flag("-rho2", "rho2", Minchen::p_rho2);

        FLAGS::Register v_mu_flag("-v_mu", "v_mu", Minchen::v_mu);
        FLAGS::Register v_xxxx_flag("-v_xxxx", "v_xxxx", Minchen::v_xxxx);

        FLAGS::Register lin_solver_type_flag("-linsolvertype", "linear solver type", Minchen::lin_solver_type);

        int num_threads = tbb::task_scheduler_init::automatic;

        FLAGS::Register thread_flag("-t", "Set number of threads", num_threads);

        PluginManager pm;
        try {
            FLAGS::ParseFlags(argc, argv);
            pm.loadAllPlugins();
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            FLAGS::PrintUsage(std::cerr);
            return 1;
        }
        if (displayHelp) {
            std::cout << "Usage:\n";
            FLAGS::PrintUsage(std::cout);
            return 0;
        }
        installSignalHandler();

        tbb::task_scheduler_init init(num_threads);

        MpmSimulation<T, dim> e;
        if (run_diff_test) {
            e.diff_test = true;
            e.diff_test_perturbation_scale = diff_test_perturbation_scale;
        }
        e.logger = LogWorker::initializeLogging();
#if dim == 2
        MpmInit2D<T> h(e, test_number);
#else
        MpmInit3D<T> h(e, test_number);
#endif
        if (!restart)
            h.start();
        else
            h.restart(restart);
    }
    return 0;
}
