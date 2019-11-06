#ifndef OBJ_IO_H
#define OBJ_IO_H

#include <vector>
#include <fstream>
#include <string>
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {
template <class T, int dim>
class Particles;
/**
    Read points data from obj file (3d) to std::vector (either 2d or 3d).
    The Z value will be ignored if 2D.
*/
template <class T, int dim>
void readPositionObj(std::istream& in, StdVector<Vector<T, dim>>& X);

/**
    Read points data from obj file (3d) to std::vector (either 2d or 3d).
    The Z value will be ignored if 2D.
*/
template <class T, int dim>
void readPositionObj(const std::string& position_file, StdVector<Vector<T, dim>>& X);

/**
    Write points data from std::vector (either 2d or 3d) to obj file (3d).
    The Z value will be 0 if 2D.
*/
template <class T, int dim>
void writePositionObj(std::ostream& os, const StdVector<Vector<T, dim>>& X);

/**
    Write points data from std::vector (either 2d or 3d) to obj file (3d).
    The Z value will be 0 if 2D.
*/
template <class T, int dim>
void writePositionObj(const std::string& filename, const StdVector<Vector<T, dim>>& X);

/**
    Read points data from obj file (3d) to ZIRAN:Particles (either 2d or 3d).
    The Z value will be ignored if ZIRAN:Particles is 2D.
*/
template <class T, int dim>
void readParticleObj(std::istream& in, Particles<T, dim>& particles);

/**
    Read points data from obj file (3d) to ZIRAN:Particles (either 2d or 3d).
    The Z value will be ignored if ZIRAN:Particles is 2D.
*/
template <class T, int dim>
void readParticleObj(const std::string& particle_file, Particles<T, dim>& particles);

/**
    Write points data from ZIRAN:Particles (either 2d or 3d) to obj file (3d).
    The Z value will be 0 if ZIRAN:Particles is 2D.
*/
template <class T, int dim>
void writeParticleObj(std::ostream& os, const Particles<T, dim>& particles);

/**
    Write points data from ZIRAN:Particles (either 2d or 3d) to obj file (3d).
    The Z value will be 0 if ZIRAN:Particles is 2D.
*/
template <class T, int dim>
void writeParticleObj(const std::string& filename, const Particles<T, dim>& particles);

/**
    Write just the faces out in obj format.

    \param[in] triangles the list of triangles.
*/
template <int vertex_per_face>
void writeFacesObj(std::ostream& os, const StdVector<Vector<int, vertex_per_face>>& faces);

/**
    Read a triangle mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] triangles the list of triangles.
*/
template <class T, int dim>
void readTrimeshObj(std::istream& is, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 3>>& triangles);

/**
    Read a triangle mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] triangles the list of triangles.
*/
template <class T, int dim>
void readTrimeshObj(const std::string& particle_file, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 3>>& triangles);

/**
    Write to a 3D triangle mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] triangles the list of triangles.

    The Z value will be 0 if the points are 2D.
*/
template <class T, int dim>
void writeTrimeshObj(std::ostream& os, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 3>>& triangles);

/**
    Write to a 3D triangle mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] triangles the list of triangles.

    The Z value will be 0 if the points are 2D.
*/
template <class T, int dim>
void writeTrimeshObj(const std::string& filename, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 3>>& triangles);

/**
    Read a quadrilateral mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] quads the list of quadrilaterals.
*/
template <class T, int dim>
void readQuadmeshObj(std::istream& is, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 4>>& quads);

/**
    Read a quadrilateral mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] quads the list of quadrilaterals.
*/
template <class T, int dim>
void readQuadmeshObj(const std::string& particle_file, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 4>>& quads);

/**
    Write to a 3D quadrilateral mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] quads the list of quadrilaterals.

    The Z value will be 0 if the points are 2D.
*/
template <class T, int dim>
void writeQuadmeshObj(std::ostream& os, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 4>>& quads);

/**
    Write to a 3D quadrilateral mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] quads the list of quadrilaterals.

    The Z value will be 0 if the points are 2D.
*/
template <class T, int dim>
void writeQuadmeshObj(const std::string& filename, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 4>>& quads);
}; // namespace ZIRAN
#endif
