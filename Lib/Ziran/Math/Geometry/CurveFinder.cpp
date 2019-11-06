#include "CurveFinder.h"

namespace ZIRAN {

// private
int CurveFinder::addCurve(StdVector<int>& pids,
    VertexID& vertex, EdgeID starting_edge, bool is_cycle)
{
    pids.emplace_back(vertex);
    VertexID& previous_vertex_id = vertex;
    int vertex_count = 1;

    EdgeID previous_edge_id = -1;
    if (is_cycle) {
        std::pair<EdgeID, EdgeID> p = adjacent_curve_edges[starting_edge];
        const Edge& e1 = segments[p.first];
        if (e1[0] == vertex || e1[1] == vertex)
            previous_edge_id = p.first;
        else
            previous_edge_id = p.second;
        assert(segments[previous_edge_id][0] == vertex || segments[previous_edge_id][1] == vertex);
    }
    EdgeID edge_id = starting_edge;
    EdgeID ending_edge = is_cycle ? starting_edge : -1;
    do {
        const Edge& edge = segments[edge_id];
        // Add the vertex that we didn't add last time
        if (edge[0] == previous_vertex_id) {
            previous_vertex_id = edge[1];
        }
        else {
            previous_vertex_id = edge[0];
        }
        pids.emplace_back(previous_vertex_id);
        vertex_count++;

        // Figure out which edge we haven't followed already
        std::pair<EdgeID, EdgeID> p = adjacent_curve_edges[edge_id];
        //ZIRAN_DBUG(p.first, p.second, previous_edge_id);
        assert(is_cycle || (p.first == previous_edge_id) || (p.second == previous_edge_id));
        if (previous_edge_id == p.first) {
            previous_edge_id = edge_id;
            edge_id = p.second;
        }
        else {
            previous_edge_id = edge_id;
            edge_id = p.first;
        }
    } while (edge_id != ending_edge);
    return vertex_count;
}

// public

CurveFinder::CurveFinder(const StdVector<Edge>& segments)
    : segments(segments)
    , adjacent_curve_edges(segments.size(), std::make_pair(-1, -1))
    , curves(segments.size())
{
    for (EdgeID e = 0; e < (EdgeID)segments.size(); e++) {
        Edge edge = segments[e];
        for (int v = 0; v < 2; v++) {
            KeyValuePair<VertexID, EdgeID> ends;
            if (curve_endpoints.remove(edge[v], ends)) {
                auto& q = adjacent_curve_edges[e];
                if (q.first == -1)
                    q.first = ends.value;
                else {
                    assert(q.second == -1);
                    q.second = ends.value;
                }
                auto& p = adjacent_curve_edges[ends.value];
                if (p.first == -1)
                    p.first = e;
                else {
                    assert(p.second == -1);
                    p.second = e;
                }
                EdgeID c0 = curves.getSetId(e);
                EdgeID c1 = curves.getSetId(ends.value);
                if (c0 == c1) {
                    // Cycle created
                    edge_on_cycle.emplace_back(e);
                }
                else {
                    curves.merge(c0, c1);
                }
            }
            else {
                curve_endpoints[edge[v]] = e;
            }
        }
    }
}

void CurveFinder::constructPaths(
    StdVector<int>& pids,
    std::vector<int>& vertex_counts, bool consistent_ordering)
{
    for (auto end : curve_endpoints) {
        if (end.value == -1)
            continue;
        VertexID vertex = end.key;
        EdgeID edge = end.value;
        if (consistent_ordering && segments[edge][1] == vertex)
            continue;
        int num_vertices = addCurve(pids, vertex, end.value, false);
        vertex_counts.emplace_back(num_vertices);
        int* other_end = curve_endpoints.get(vertex);
        //assert(other_end != nullptr);
        if (other_end != nullptr)
            *other_end = -1;
    }
}

void CurveFinder::constructCycles(StdVector<int>& pids,
    std::vector<int>& vertex_counts)
{
    for (EdgeID e : edge_on_cycle) {
        VertexID vertex = segments[e](0);
        int num_vertices = addCurve(pids, vertex, e, true);
        vertex_counts.emplace_back(num_vertices);
    }
}

void CurveFinder::getCurveEndpoints(HashTable<VertexID, EdgeID>& curve_endpoints_in) const
{
    for (auto iter = curve_endpoints.begin(); iter != curve_endpoints.end(); ++iter)
        curve_endpoints_in.insert(iter->key, iter->value);
}
} // namespace ZIRAN
