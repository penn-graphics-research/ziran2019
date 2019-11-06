#ifndef SIMULATION_BASE_H
#define SIMULATION_BASE_H
#include <Ziran/CS/Util/ErrorContext.h>
#include <Ziran/CS/Util/Filesystem.h>
#include <Ziran/CS/Util/DataDir.h>
#include <Ziran/CS/Util/PrettyPrinting.h>
#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/CS/Util/Timer.h>
#include <Ziran/Sim/TimeStepping.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <vector>
#include <type_traits>
#include <tbb/tbb.h>

namespace ZIRAN {
class RestartException : public std::runtime_error {
public:
    int frame;
    RestartException(int frame)
        : std::runtime_error("Restart Simulation")
        , frame(frame)
    {
    }
};

class SimulationBase {
protected:
    int frame;
    int substep;
    bool restarting;

public:
    TimeStepping step;
    DataDir output_dir;
    bool write_substeps;
    bool write_files;
    bool write_log;
    bool verbose;
    int start_frame, end_frame, restart_frame;
    std::shared_ptr<LogWorker> logger;
    bool diff_test;
    double diff_test_perturbation_scale;
    StdVector<std::function<void(int)>> initialize_callbacks;
    StdVector<std::function<void(int)>> reinitialize_callbacks;
    StdVector<std::function<void(int)>> begin_frame_callbacks;
    StdVector<std::function<void(int)>> end_frame_callbacks;
    StdVector<std::function<void(int, int)>> end_time_step_callbacks;
    StdVector<std::function<void(int, int, double, double)>> begin_time_step_callbacks;
    StdVector<std::function<void(int)>> restart_callbacks;

    StdVector<std::function<void(int, int)>> general_callbacks;

    int binary_ver = 1;

    bool full_implicit;

    enum LinearSolverType {
        MINRES,
        CG
    } linear_solver_type;

    SimulationBase()
        : frame(0)
        , substep(0)
        , restarting(false)
        , output_dir("output")
        , write_substeps(false)
        , write_files(true)
        , write_log(true)
        , verbose(false)
        , start_frame(0)
        , end_frame(0)
        , restart_frame(0)
        , logger(nullptr)
        , diff_test(false)
        , diff_test_perturbation_scale((double)1)
        , full_implicit(false)
        , linear_solver_type(MINRES)
    {
    }

    virtual ~SimulationBase()
    {
        if (logger)
            logger->printTimings();
    }

    /**
      This function will only be called once
    */
    virtual void initialize()
    {
        if (write_files || write_log) {
            output_dir.createPath();
        }
        if (write_log && logger) {
            logger->openLogFile(output_dir.absolutePath("log.txt"), restarting);
        }
        for (auto f : initialize_callbacks) {
            f(frame);
        }
        if (!restarting)
            restart_frame = start_frame;
    }

    /**
      This function can be called more than once on
      the same object, and should be idempotent

      It is used to put the simulation into a consistent
      state given the values of its member variables
    */
    virtual void reinitialize()
    {
        for (auto f : reinitialize_callbacks) {
            f(frame);
        }
    }
    /**
      Callback to advance one time step

      Most of the simulation is done by this callback
      which should advance the state forward by dt
    */
    virtual void advanceOneTimeStep(double dt) = 0;

    /**
      Callback to write the simulation state

      Should be overridden to write the simulation state
    */
    virtual void writeState(std::ostream&) = 0;

    /**
      Callback to read the simulation state

      Should be overridden to read the output of writeState
    */
    virtual void readState(std::istream&) = 0;

    /**
      write to file
      */
    virtual void write(const std::string& filename)
    {
        if (!write_files)
            return;
        std::ofstream file = output_dir.openBinaryOutput(filename);
        ZIRAN_INFO("Writing ", filename);
        writeState(file);
        std::ofstream last_written = output_dir.openTextOutput("last_written.txt");
        last_written << frame;

        std::ofstream binary_ver_file = output_dir.openTextOutput("binary_ver.txt");
        binary_ver_file << binary_ver;
    }

    /**
      read from file
    */
    virtual void read(const std::string& filename)
    {
        int binary_ver_data;
        try {
            std::ifstream binary_ver_file = output_dir.openTextInput("binary_ver.txt");
            binary_ver_file >> binary_ver_data;
        }
        catch (std::exception& e) {
            binary_ver_data = 0;
        }
        ZIRAN_ASSERT(binary_ver_data <= binary_ver, "The data to read is of newer version. Update code base!");
        ZIRAN_INFO("Binary Ver: ", binary_ver_data);
        binary_ver = binary_ver_data;

        std::ifstream last_written = output_dir.openTextInput("last_written.txt");
        int last_written_frame;
        last_written >> last_written_frame;
        ZIRAN_WARN_IF(last_written_frame < frame, "Frame is after the last written frame ", last_written_frame, ". Could be reading stale data");
        std::ifstream file = output_dir.openBinaryInput(filename);
        ZIRAN_INFO("Reading ", filename);
        readState(file);
    }

    /**
      Constructs the output filename
    */
    virtual std::string outputFileName(const std::string& name = "restart", const std::string& suffix = ".dat")
    {
        std::ostringstream filename;
        filename << name << "_" << frame;
        if (substep > 0)
            filename << "_" << substep;
        filename << suffix;
        return filename.str();
    }

    /**
      Callback for beginning of the frame
    */
    virtual void beginFrame(int frame)
    {
        for (auto f : begin_frame_callbacks) {
            f(frame);
        }
    }
    /**
      Callback for end of the frame

      If you override this function
      make sure that you allow for restart to be called
      It's most likely as simple as
      ```
      virtual void endFrame(int frame) override
      {
            SimulationBase::endFrame(frame);
            if(restarting) return;
            // Custom code here
      }
      ``
      {
    */

    virtual void endFrame(int frame)
    {
        for (auto f : end_frame_callbacks) {
            f(frame);
            if (restarting)
                break;
        }
    }
    /**
      Callback for beginning of the time step
    */
    virtual void beginTimeStep(int frame, int substep, double time, double dt)
    {
        for (auto f : begin_time_step_callbacks) {
            f(frame, substep, time, dt);
        }
    }
    /**
      Callback for end of the time step
    */
    virtual void endTimeStep(int frame, int substep)
    {
        if (write_substeps) {
            write(outputFileName());
        }
        for (auto f : end_time_step_callbacks) {
            f(frame, substep);
        }
    }
    /**
      Callback to calculate dt
    */
    virtual double calculateDt()
    {
        return step.max_dt;
    }
    /**
      Call restart when you want to begin the simulation by reading
      the data from a previous simulation

      Safe to call:
      Outside main loop
      During endFrame() callback
      If you need to restart at a different time throw a RestartException

      new_restart_frame should be the frame that we want to read data from
      */
    virtual void restart(int new_restart_frame)
    {
        // Not sure if the simulate loop is running when restart is called
        // If it is not setting start_frame should take care of it
        restart_frame = new_restart_frame;
        frame = restart_frame;
        // Can only restart from frame boundaries not substeps
        substep = 0;
        step.reset(frame);
        restarting = true;
        ZIRAN_INFO("Resetting time = ", step.time);
    }

    virtual void advanceOneFrame()
    {
        ZIRAN_TIMER();
        beginFrame(frame);
        TimeStepping::Event e;
        ZIRAN_EVENT(Color::BLUE, "Frame ", frame, Color::RESET);
        do {
            ZIRAN_CONTEXT(substep);
            bool overwrite_mode = Event == Level::getConsoleLevel();
            ZIRAN_EVENT(ProgressBar(step.time_since_last_frame / step.frame_dt),
                overwrite_mode ? "" : ("frame " + std::to_string(frame)));
            ZIRAN_EVENT(ProgressBar(step.time / (end_frame * step.frame_dt)));
            ZIRAN_EVENT(Color::MAGENTA, "Substep ", std::setw(8), std::right, substep, Color::RESET);
            if (overwrite_mode)
                ZIRAN_EVENT('\r', moveCursorRelative(-4, 0));

            double dt = step.nextDt(calculateDt());
            ZIRAN_CONTEXT(step.time);
            ZIRAN_CONTEXT(dt);
            beginTimeStep(frame, substep, step.time, dt);
            advanceOneTimeStep(dt);
            e = step.advance(dt);
            substep++;
            endTimeStep(frame, substep);
        } while (e == TimeStepping::Substep);
        ZIRAN_EVENT(ProgressBar(1.0));
        ZIRAN_EVENT(ProgressBar(step.time / (end_frame * step.frame_dt)));
        substep = 0;
        endFrame(frame);
        if (!restarting) {
            write(outputFileName());
            ZIRAN_ETA(end_frame - frame);
        }
    }

    /**
      Simulates from start frame to end_frame
    */
    void simulate()
    {
        ZIRAN_TIMER();
        //Don't overwrite the data we will be reading in
        if (!restarting)
            write(outputFileName());

        for (frame = restart_frame + 1; frame <= end_frame; frame++) {
            ZIRAN_CONTEXT(frame);
            if (restarting) {
                ZIRAN_INFO("Restarting frame ", frame);
                frame--; //Need to read the data from the last frame
                read(outputFileName());
                for (auto f : restart_callbacks) {
                    f(frame);
                }
                restart_callbacks.clear();
                restarting = false;
                frame++;
            }
            try {
                advanceOneFrame();
            }
            catch (RestartException& e) {
                restart(e.frame);
            }
        }
    }

    bool currentlyRestarting() const
    {
        return restarting;
    }

    int getCurrentFrame() const
    {
        return frame;
    }

    virtual bool useDouble() = 0;
    virtual int dimension() = 0;
    virtual const char* name() = 0;
};
} // namespace ZIRAN
#endif
