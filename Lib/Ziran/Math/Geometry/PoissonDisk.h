#ifndef POISSON_DISK_H
#define POISSON_DISK_H
#include <Ziran/CS/DataStructure/Box.h>
#include <Ziran/CS/DataStructure/MultiArray.h>
#include <Ziran/CS/Util/BinaryIO.h>
#include <Ziran/CS/Util/DataDir.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/Math/MathTools.h>
#include <algorithm>
#include <vector>
namespace ZIRAN {

template <class T, int dim>
class PoissonDisk {
public:
    typedef Vector<T, dim> EigenTV;
    typedef Vector<int, dim> EigenIV;
    RandomNumber<T> rnd;
    T min_distance;
    int max_attempts;
    T h;
    EigenTV min_corner;
    EigenTV max_corner;
    bool periodic;

    PoissonDisk(const RandomNumber<T>& rnd, T min_distance, const EigenTV& min_corner, const EigenTV& max_corner, int max_attempts = 30, bool periodic = false)
        : rnd(rnd)
        , min_distance(min_distance)
        , max_attempts(max_attempts)
        , h(min_distance / std::sqrt((T)dim))
        , min_corner(min_corner)
        , max_corner(max_corner)
        , periodic(periodic)
    {
        ZIRAN_ASSERT(min_corner != max_corner, "min_corner == max_corner in PoissonDisk");
    }

    PoissonDisk(const int seed, T min_distance, const EigenTV& min_corner, const EigenTV& max_corner, int max_attempts = 30, bool periodic = false)
        : min_distance(min_distance)
        , max_attempts(max_attempts)
        , h(min_distance / std::sqrt((T)dim))
        , min_corner(min_corner)
        , max_corner(max_corner)
        , periodic(periodic)
    {
        rnd.resetSeed(seed);
        ZIRAN_ASSERT(min_corner != max_corner, "min_corner == max_corner in PoissonDisk");
    }

    ~PoissonDisk()
    {
    }

private:
    /**
    Convert position in the world space to its corresponding index space.
     */
    EigenIV worldToIndexSpace(const EigenTV& X) const
    {
        EigenIV ijk;
        for (size_t d = 0; d < dim; ++d)
            // ijk(d) = (int)std::floor((X(d) - min_corner(d)) / h);
            ijk(d) = MATH_TOOLS::int_floor((X(d) - min_corner(d)) / h);

        return ijk;
    }

    /**
    Generate one random point around center with a distance between .5 and 1 from center.
     */
    EigenTV generateRandomPointAroundAnnulus(const EigenTV& center)
    {
        while (true) {
            EigenTV v;
            rnd.fill(v, -1, 1);
            T mag2 = v.dot(v);
            if (mag2 >= 0.25 && mag2 <= 1)
                return v * min_distance * 2 + center;
        }
        return EigenTV::Zero();
    }

    /**
    Return true if the distance between the candidate point and any other points in samples are sufficiently far away (at least min_distance away), and false otherwise.
     */
    bool checkDistance(const EigenTV& point, const MultiArray<int, dim>& background_grid, const StdVector<EigenTV>& samples) const
    {
        EigenIV index = worldToIndexSpace(point);
        // Check if we are outside of the background_grid. If so, return false
        for (int d = 0; d < dim; ++d) {
            if (index(d) < 0 || index(d) >= background_grid.size(d))
                return false;
        }
        // If there is already a particle in that cell, return false
        if (background_grid(index) != -1)
            return false;
        T min_distance_sqr = min_distance * min_distance;
        EigenIV local_min_index = index.array() - 2;
        EigenIV local_max_index = index.array() + 3;
        // If not periodic, clamp local_min_index and local_max_index to the size of the background grid
        if (!periodic) {
            local_min_index = local_min_index.cwiseMax(0);
            local_max_index = local_max_index.cwiseMin(background_grid.size);
        }
        // Create local_box for iterator purposes
        Box<int, dim> local_box(local_min_index, local_max_index);
        if (!periodic)
            for (MaxExclusiveBoxIterator<dim> it(local_box); it.valid(); ++it) {
                if (background_grid(it.index) == -1)
                    continue;
                EigenTV x = point - samples[background_grid(it.index)];
                if (x.dot(x) < min_distance_sqr) {
                    return false;
                }
            }
        else {
            for (MaxExclusiveBoxIterator<dim> it(local_box); it.valid(); ++it) {
                EigenIV local_index = it.index;
                // Need to shift point in MultiArray if one of the indices is negative or greater than background_grid.size
                EigenTV shift = EigenTV::Zero();
                for (int d = 0; d < dim; ++d) {
                    // If it.index < 0 update local index to all the way down to the right/top/back
                    // If there is a point in that MultiArray index, we need to shift that point to the left/bottom/front
                    if (it.index(d) < 0) {
                        local_index(d) = local_index(d) % background_grid.size(d) + background_grid.size(d);
                        shift(d) = min_corner(d) - max_corner(d);
                    }
                    // If it.index(d) >= background_grid(d) update local index to all the way down to the left/bottom/front
                    // If there is a point in that MultiArray index, we need to shift that point to the left right/top/back
                    else if (it.index(d) >= background_grid.size(d)) {
                        local_index(d) = local_index(d) % background_grid.size(d);
                        shift(d) = max_corner(d) - min_corner(d);
                    }
                }
                if (background_grid(local_index) == -1)
                    continue;
                EigenTV x = point - (samples[background_grid(local_index)] + shift);
                if (x.dot(x) < min_distance_sqr) {
                    return false;
                }
            }
        }
        return true;
    }

public:
    /**
    Set min_distance between particles based on the number of particles per cell and grid dx.
    ppc = particles per cell. 
     */
    T setDistanceByParticlesPerCell(T dx, T ppc)
    {
        T v = std::pow(dx, dim) / (T)ppc;
        if (dim == 2) {
            min_distance = std::sqrt(v * ((T)2 / 3));
        }
        else if (dim == 3) {
            min_distance = std::pow(v * ((T)13 / 18), (T)1 / 3);
        }
        else {
            ZIRAN_ASSERT(false, "Poisson disk only supports 2D and 3D");
        }
        h = min_distance / std::sqrt((T)dim);
        return min_distance;
    }

    template <class Func>
    void sampleFromPeriodicData(StdVector<Vector<T, 2>>& samples, Func&& feasible)
    {
        ZIRAN_ASSERT(false, "sampleFromPeriodicData not implemented for 2D!");
    }

    /**
    Faster sampling based on particles-1000k.dat file. The function assumes that this file is saved in Data/MpmParticles/particles-1000k.dat.
    min_distance = mimimum distance between particles.
    A few information about particles-1000k.dat:
      - Samples count: 1001436
      - min_distance = 1	h=0.57735
      - min_point =      -60 -59.9999 -59.9998
      - max_point = 59.9999      60 59.9999
    The file is generated by Projects/pdsampler/pdsampler.cpp.
     */
    template <class Func>
    void sampleFromPeriodicData(StdVector<Vector<T, 3>>& samples, Func&& feasible)
    {
        // Figure out the number of offsets
        EigenTV scaled_ref_box_length(120.0, 120.0, 120.0);
        scaled_ref_box_length *= min_distance;
        EigenTV side_length = max_corner - min_corner;
        EigenIV offset_number;
        for (int d = 0; d < dim; ++d)
            offset_number(d) = std::ceil(side_length(d) / scaled_ref_box_length(d)) + 1;

        // Read std vector
        DataDir dir;
        std::ifstream is = dir.openBinaryInput("MpmParticles/particles-1000k.dat");
        size_t size = readEntry<size_t>(is);
        samples.reserve(size * (side_length.prod() / scaled_ref_box_length.prod()));
        readEntry<size_t>(is);
        EigenTV new_point, offset_center, offset_new_point;
        Box<int, dim> index_box(EigenIV::Zero(), offset_number);
        // Read from file as float
        Vector<float, dim> new_point_read;
        for (size_t i = 0; i < size; i++) {
            new_point_read = readEntry<Vector<float, dim>>(is);
            for (int d = 0; d < dim; ++d)
                new_point(d) = min_distance * new_point_read(d) + min_corner(d);
            for (MaxExclusiveBoxIterator<dim> it(index_box); it.valid(); ++it) {
                for (int d = 0; d < dim; ++d)
                    offset_center(d) = (it(d) /*+ (T).5*/) * scaled_ref_box_length(d);
                offset_new_point = new_point + offset_center;
                bool inside = true;
                for (int d = 0; d < dim; ++d)
                    if (offset_new_point(d) < min_corner(d) || offset_new_point(d) > max_corner(d))
                        inside = false;
                if (feasible(offset_new_point) && inside)
                    samples.emplace_back(offset_new_point);
            }
        }
    }

    /**
    Samples particles.
     */
    template <class Func>
    void sample(StdVector<EigenTV>& samples, Func&& feasible)
    {
        /*
      Set up background grid
      dx should be bounded by min_distance / sqrt(dim)
      background_grid is a MultiArray which keeps tracks the indices of points in samples.
      the value of background_grid is initialized to be -1, meaning that there is no particle
      in that cell
      */
        EigenTV cell_numbers_candidate = max_corner - min_corner;
        EigenIV cell_numbers;
        for (int d = 0; d < dim; ++d)
            cell_numbers(d) = std::ceil(cell_numbers_candidate(d) / h);
        MultiArray<int, dim> background_grid(cell_numbers, -1);
        // Set up active list
        StdVector<int> active_list(samples.size());

        /*
      If X already contains points, then set the values of the background grid
      to keep track of these points
      */
        if (samples.size()) {
            for (size_t i = 0; i < samples.size(); ++i) {
                // Change the value of background_grid in that index to be i
                background_grid(worldToIndexSpace(samples[i])) = i;
                active_list.push_back(i);
            }
        }
        else {
            // Generate a random point within the range and append it to samples and active list
            EigenTV first_point = (T).5 * (max_corner + min_corner);
            while (!feasible(first_point)) {
                for (int d = 0; d < dim; ++d) {
                    T r = static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
                    first_point(d) = min_corner(d) + (max_corner(d) - min_corner(d)) * r;
                }
            }
            samples.push_back(first_point);
            background_grid(worldToIndexSpace(first_point)) = 0;
            active_list.push_back(0);
        }

        /*
      While active_list is non-zero, do step 2 in Bridson's proposed algorithm
      */
        while (active_list.size()) {
            // Get a random index from the active list and find the point corresponding to it
            int random_index = rnd.randInt(0, active_list.size() - 1);
            EigenTV current_point = samples[active_list[random_index]];
            // Swap random index with the last element in the active list so that we can pop_back
            // if found_at_least_one is false at the end of this procedure
            iter_swap(active_list.begin() + random_index, active_list.end() - 1);
            // Generate up to max_attempts points in the annulus of radius h and 2h around current_point
            bool found_at_least_one = false;
            for (int i = 0; i < max_attempts; ++i) {
                EigenTV new_point = generateRandomPointAroundAnnulus(current_point);

                // If periodic and new_point is outside of the min_corner, max_corner shift it to be inside
                if (periodic) {
                    for (int d = 0; d < dim; ++d) {
                        if (new_point(d) < min_corner(d))
                            new_point(d) += max_corner(d) - min_corner(d);
                        else if (new_point(d) > max_corner(d))
                            new_point(d) -= max_corner(d) - min_corner(d);
                    }
                }

                if (!feasible(new_point))
                    continue;

                if (checkDistance(new_point, background_grid, samples)) {
                    found_at_least_one = true;
                    // Add new_point to samples
                    samples.push_back(new_point);
                    int index = samples.size() - 1;
                    // Add new_point to active list
                    active_list.push_back(index);
                    // Update background_grid
                    background_grid(worldToIndexSpace(new_point)) = index;
                }
            }
            // If not found at least one, remove random_index from active list
            if (!found_at_least_one) { //active_list.erase(active_list.begin()+random_index);}
                active_list.pop_back();
            }
        }

        // Remove points that are outside of Box(min_corner,max_corner)
        for (int i = samples.size(); i >= 0; i--) {
            EigenTV point = samples.back();
            bool outside = false;
            for (int d = 0; d < dim; ++d) {
                if (point(d) < min_corner(d) || point(d) > max_corner(d)) {
                    outside = true;
                    break;
                }
            }
            if (outside)
                samples.pop_back();
        }
    }
};
} // namespace ZIRAN
#endif
