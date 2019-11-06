#ifndef OUTPUT_POLY_MESH_H
#define OUTPUT_POLY_MESH_H
#include <cstdint>
#include <stddef.h>
namespace ZIRAN {
// The OutputPolyMesh class
// is designed to provide a uniform interface
// for polygon meshes for the purpose of
// outputting to a file
// It is not templatized instead it assumes 3d and float
struct OutputPolyMesh {
    const float* positions; //size 3*num_points
    size_t num_points;
    const int32_t* vertices; //size num_vertices
    size_t num_vertices;
    const int32_t* face_vertex_counts; // size num_faces
    size_t num_faces;
};
} // namespace ZIRAN
#endif
