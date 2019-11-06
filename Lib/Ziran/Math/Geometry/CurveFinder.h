#ifndef CURVE_FINDER_H
#define CURVE_FINDER_H

#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/DataStructure/DisjointSets.h>

namespace ZIRAN {

// Helper functions for creating curves
class CurveFinder {
    using CurveID = int;
    using EdgeID = int;
    using VertexID = int;
    using Edge = Eigen::Vector2i;
    const StdVector<Edge>& segments;
    std::vector<std::pair<EdgeID, EdgeID>> adjacent_curve_edges;
    HashTable<VertexID, EdgeID> curve_endpoints;
    DisjointSets<EdgeID> curves;
    std::vector<EdgeID> edge_on_cycle;

    int addCurve(StdVector<int>& pids,
        VertexID& vertex, EdgeID starting_edge, bool is_cycle);

public:
    CurveFinder(const StdVector<Edge>& segments);

    void constructPaths(
        StdVector<int>& pids,
        std::vector<int>& vertex_counts, bool consistent_ordering);

    void constructCycles(StdVector<int>& pids,
        std::vector<int>& vertex_counts);

    void getCurveEndpoints(HashTable<VertexID, EdgeID>& curve_endpoints_in) const;
};
} // namespace ZIRAN

#endif
