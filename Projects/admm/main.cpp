#include <Ziran/CS/Util/FloatingPointExceptions.h>
#include <Ziran/CS/Util/SignalHandler.h>
#include <Ziran/CS/Util/CommandLineFlags.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/Filesystem.h>
#include <Ziran/CS/Util/PluginManager.h>
#include <tbb/tbb.h>
#include "AdmmInit2D.h"
#include "AdmmInit3D.h"

#define T double
#define dim 3

using namespace ZIRAN;

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
        float phase_field_percentage;
        float phase_field_l0_ratio;

        // Not checking for nan, because when constitutive model returns that, MpmForceBase is skipping them (treating as zeros)
        // FPE::WatchedScope w(FPE::Mask::Overflow | FPE::Mask::DivZero);
        // Unconmment the following to catch division by 0
        // FPE::WatchedScope w(FPE::Mask::Overflow | FPE::Mask::Invalid | FPE::Mask::DivZero);
        FPE::WatchedScope w(FPE::Mask::Invalid);

        std::string script_file_name;
        StdVector<std::string> inline_strings;
        FLAGS::Register helpflag("--help", "Print help (this message) and exit", displayHelp);

        // Lua command line options
        FLAGS::Register scriptflag("-script", "Lua script to read for initial data", script_file_name);
        FLAGS::Register iflag("-i", "Append string to script", inline_strings);

        // Non-lua command line options
        FLAGS::Register test_number_flag("-test", "Test number (non-lua test)", test_number);
        FLAGS::Register three_d_flag("--3d", "Dimension is 3(non-lua test)", three_d);
        FLAGS::Register run_diff_test_flag("--run_diff_test", "Run diff test (non-lua test)", run_diff_test);
        FLAGS::Register diff_test_perturbation_scale_flag("-dtps", "diff_test_perturbation_scale (non-lua test)", diff_test_perturbation_scale);
        FLAGS::Register double_flag("--double", "Dimension (non-lua test)", use_double);
        FLAGS::Register restart_flag("-restart", "Restart frame (non-lua test)", restart);

        FLAGS::Register phase_field_percentage_flag("-p", "phase_field_percentage_flag", phase_field_percentage);
        FLAGS::Register phase_field_l0_ratio_flag("-l0", "phase_field_l0_ratio_flag", phase_field_l0_ratio);

        int num_threads = tbb::task_scheduler_init::automatic;

        FLAGS::Register thread_flag("-t", "Set number of threads", num_threads);

        std::stringstream script;
        PluginManager pm;
        try {
            FLAGS::ParseFlags(argc, argv);
            if (!script_file_name.empty()) {
                FILESYSTEM::readFile(script, script_file_name);
            }
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
        for (const auto& s : inline_strings) {
            script << s << "\n";
        }

        tbb::task_scheduler_init init(num_threads);

        if (script_file_name.empty()) {
            ZIRAN_ASSERT(test_number != -1, "No lua script loaded. Either load with --script or set --test");
            AdmmSimulation<T, dim> e;
            if (run_diff_test) {
                e.diff_test = true;
                e.diff_test_perturbation_scale = diff_test_perturbation_scale;
            }
            e.logger = LogWorker::initializeLogging();
#if dim == 2
            AdmmInit2D<T> h(e, test_number);
#else
            AdmmInit3D<T> h(e, test_number);
#endif
            if (!restart)
                h.start();
            else
                h.restart(restart);
        }
    }
    return 0;
}
