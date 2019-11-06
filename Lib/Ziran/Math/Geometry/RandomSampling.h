#ifndef RANDOM_SAMPLING_H
#define RANDOM_SAMPLING_H

#include <Ziran/Math/Geometry/AnalyticLevelSet.h>
#include <Ziran/Math/Geometry/VdbLevelSet.h>
#include <Ziran/CS/Util/RandomNumber.h>
#include <vector>

namespace ZIRAN {

//samples uniformly random points
template <class T, int dim>
class RandomSampling {
public:
    typedef Vector<T, dim> TV;
    RandomNumber<T> rnd;
    const TV min_corner, max_corner;
    const int points_number;

    RandomSampling(const RandomNumber<T>& rnd_in, const TV& min_corner_in, const TV& max_corner_in, const int& points_number_in)
        : rnd(rnd_in)
        , min_corner(min_corner_in)
        , max_corner(max_corner_in)
        , points_number(points_number_in)
    {
        ZIRAN_ASSERT(min_corner != max_corner, "min_corner == max_corner in RandomSampling");
    }

    RandomSampling(const int& seed, const TV& min_corner_in, const TV& max_corner_in, const int& points_number_in)
        : min_corner(min_corner_in)
        , max_corner(max_corner_in)
        , points_number(points_number_in)
    {
        rnd.resetSeed(seed);
        ZIRAN_ASSERT(min_corner != max_corner, "min_corner == max_corner in RandomSampling");
    }

    ~RandomSampling()
    {
    }

    template <class Func>
    void sample(StdVector<TV>& samples, Func&& feasible)
    {
        TV new_point;
        while (samples.size() < (size_t)points_number) {
            new_point = rnd.randInBox(min_corner, max_corner);
            if (!feasible(new_point))
                continue;
            samples.emplace_back(new_point);
        }
    }
};
} // namespace ZIRAN
#endif
