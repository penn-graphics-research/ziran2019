#include <MPM/MpmSimulationBase.h>
#include <MPM/Force/FBasedMpmForceHelper.h>
#include <MPM/MpmSimulationDataAnalysis.h>
#include <MPM/MpmParticleHandleBase.h>
#include <Ziran/CS/Util/CommandLineFlags.h>
#include <Ziran/CS/Util/FloatingPointExceptions.h>
#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/CS/Util/SignalHandler.h>
#include <Ziran/CS/Util/AttributeNamesForward.h>
#include <Ziran/Math/Geometry/Grid.h>
#include <Ziran/Math/Geometry/MeshConstruction.h>
#include <Ziran/Math/Geometry/PoissonDisk.h>
#include <Ziran/Math/MathTools.h>
#include <gtest/gtest.h>
#include "MpmSimulation.h"

using namespace ZIRAN;
using namespace ZIRAN::MATH_TOOLS;

int main(int argc, char* argv[])
{
    int r = 0;
    {
        tbb::task_scheduler_init scope;

        openvdb::initialize();
        // FPE::WatchedScope w(FPE::Mask::Overflow | FPE::Mask::Invalid | FPE::Mask::DivZero);
        installSignalHandler();

        bool displayHelp = false;
        FLAGS::Register helpflag("--help", "Print help (this message) and exit", displayHelp);

        testing::InitGoogleTest(&argc, argv);
        // Google Test will remove arguments it recognized
        try {
            FLAGS::ParseFlags(argc, argv);
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            FLAGS::PrintUsage(std::cerr);
            return 1;
        }

        if (displayHelp) {
            std::cout << "\n";
            std::cout << "Custom Flags:\n";
            FLAGS::PrintUsage(std::cout);
        }
        r = RUN_ALL_TESTS();
    }
    sleep(1);
    return r;
}
