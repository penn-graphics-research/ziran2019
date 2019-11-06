#ifndef TIME_STEPPING_H
#define TIME_STEPPING_H
#include <Ziran/CS/Util/Debug.h>
namespace ZIRAN {
/**

*/
class TimeStepping {
public:
    enum Event {
        Frame,
        Substep
    };

    double time;
    double time_since_last_frame;

    double frame_dt;
    double min_dt;
    double max_dt;

    double max_dt_original;

    TimeStepping()
        : time(0.0)
        , time_since_last_frame(0.0)
        , frame_dt(1.0 / 24)
        , min_dt(1e-6)
        , max_dt(frame_dt)
        , max_dt_original(frame_dt)
    {
    }

    /**
      \brief Calculates the next time step

      The next time step is chosen to be in between min_dt and max_dt,
      and so that the frames times will be hit exactly.  The min_dt
      and max_dt are used as suggestions, and not strictly enforced
      if necessary to make the desired frame timestep.

      \param dt The desired next time step
      \return The next time step
    */
    double nextDt(double dt)
    {
        ZIRAN_ASSERT(dt > 0);
        max_dt = std::min(max_dt, max_dt_original);
        dt = (dt < min_dt) ? min_dt : (dt > max_dt) ? max_dt : dt;

        if (time_since_last_frame + dt >= frame_dt)
            // Final substep in the frame
            dt = frame_dt - time_since_last_frame;
        else if (time_since_last_frame + 2 * dt > frame_dt)
            // Avoid small timesteps to finish off frame
            dt = (frame_dt - time_since_last_frame) / 2;

        return dt;
    }

    /**
      \brief Advances the time

      \param dt The amount of time to advance
      \return Whether the time step finished off a Frame or a Substep
    */
    Event advance(double dt)
    {
        time += dt;
        time_since_last_frame += dt;
        if (time_since_last_frame >= frame_dt) {
            time_since_last_frame = 0;
            return Frame;
        }
        return Substep;
    }

    /**
      \brief Reset the current time to that of the given frame

      \param frame to reset to
      */
    void reset(int frame)
    {
        time = frame * frame_dt;
        time_since_last_frame = 0;
    }
};
} // namespace ZIRAN
#endif
