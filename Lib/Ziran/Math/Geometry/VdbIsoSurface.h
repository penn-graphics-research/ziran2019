#ifdef ZIRAN_WITH_VDB
#ifndef VDB_ISO_SURFACE_H
#define VDB_ISO_SURFACE_H

#include <Ziran/Math/Geometry/SimplexMesh.h>

#include <openvdb/openvdb.h>
#include <openvdb/math/Operators.h>
#include <openvdb/tools/VolumeToMesh.h>
#undef B2

namespace ZIRAN {

/**
   Build the triangle isosurface of a scalar vdb grid 
 */
template <class T, class GridT>
void scalarGridMarchingCube(GridT& grid, const T isovalue, StdVector<Vector<T, 3>>& X, StdVector<Vector<int, 3>>& triangles)
{
    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    openvdb::tools::volumeToMesh(grid, points, quads, isovalue);

    X.resize(points.size());
    triangles.resize(quads.size() * 2);
    for (size_t p = 0; p < points.size(); ++p)
        for (int d = 0; d < 3; d++)
            X[p](d) = points[p](d);
    size_t triangle_count = 0;
    for (size_t q = 0; q < quads.size(); ++q) {
        int a = quads[q](0);
        int b = quads[q](1);
        int c = quads[q](2);
        int d = quads[q](3);
        Vector<int, 3> first_triangle, second_triangle;
        first_triangle << a, b, c;
        second_triangle << a, c, d;
        triangles[triangle_count++] = first_triangle;
        triangles[triangle_count++] = second_triangle;
    }
    ZIRAN_ASSERT(triangle_count == (size_t)2 * (quads.size()));
}
} // namespace ZIRAN

#endif
#endif
