#ifndef SPATIAL_HASH_H
#define SPATIAL_HASH_H

#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/DataStructure/Box.h>
#include <Ziran/CS/DataStructure/HashTable.h>

namespace ZIRAN {

template <class T, int dim>
class SpatialHash {
public:
    typedef Vector<T, dim> TV;
    typedef Vector<int, dim> IV;

    T h;
    HashTable<IV, StdVector<int>> hash;

    SpatialHash()
        : h(-1)
    {
    }

    void rebuild(const T h_input, const StdVector<TV>& Xs)
    {
        ZIRAN_ASSERT(h_input > 0);

        h = h_input;
        hash.clear();
        for (size_t i = 0; i < Xs.size(); ++i) {
            const TV& X = Xs[i];
            IV cell = IV::Zero();
            for (int d = 0; d < dim; d++)
                cell(d) = (int)(std::floor(X(d) / h));
            hash[cell].push_back(i);
        }
    }

    void oneLayerNeighbors(const TV& X, StdVector<int>& neighbors)
    {
        ZIRAN_ASSERT(h > 0);

        neighbors.clear();
        IV cell = IV::Zero();
        for (int d = 0; d < dim; d++)
            cell(d) = (int)(std::floor(X(d) / h));
        IV local_min_index = cell.array() - 1;
        IV local_max_index = cell.array() + 2;

        Box<int, dim> local_box(local_min_index, local_max_index);

        for (MaxExclusiveBoxIterator<dim> it(local_box); it.valid(); ++it) {
            StdVector<int>* v = hash.get(it.index);
            if (v != nullptr)
                neighbors.insert(neighbors.end(), v->begin(), v->end());
        }
    }

    void oneLayerNeighborsWithinRadius(const TV& X, const StdVector<TV>& points, const T radius, StdVector<int>& neighbors)
    {
        ZIRAN_ASSERT(h > 0);
        const T radius_squared = radius * radius;

        neighbors.clear();
        IV cell = IV::Zero();
        for (int d = 0; d < dim; d++)
            cell(d) = (int)(std::floor(X(d) / h));
        IV local_min_index = cell.array() - 1;
        IV local_max_index = cell.array() + 2;

        Box<int, dim> local_box(local_min_index, local_max_index);

        for (MaxExclusiveBoxIterator<dim> it(local_box); it.valid(); ++it) {
            StdVector<int>* v = hash.get(it.index);
            if (v != nullptr) {
                for (auto i : *v) {
                    if ((points[i] - X).squaredNorm() < radius_squared)
                        neighbors.push_back(i);
                }
            }
        }
    }
};
} // namespace ZIRAN

#endif
