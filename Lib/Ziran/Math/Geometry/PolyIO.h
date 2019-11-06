#ifndef POLY_IO_H
#define POLY_IO_H

#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/Forward.h>
#include <fstream>

namespace ZIRAN {

/**
    Write points data from std::vector (either 2d or 3d) to poly file (3d).
    The Z value will be 0 if 2D.
*/
template <class T, int dim>
void writePositionPoly(std::ostream& os, const StdVector<Vector<T, dim>>& X)
{
    os << "POINTS\n";
    int count = 0;
    for (auto x : X) {
        os << ++count << ":";
        for (int i = 0; i < dim; i++)
            os << " " << x(i);
        if (dim == 2)
            os << " 0";
        os << "\n";
    }
}

/**
    Write points data from std::vector (either 2d or 3d) to poly file (3d).
    The Z value will be 0 if 2D.
*/
template <class T, int dim>
void writePositionPoly(const std::string& filename, const StdVector<Vector<T, dim>>& X)
{
    std::ofstream fs;
    fs.open(filename);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    writePositionPoly(fs, X);
    fs.close();
}

/**
    Write to a 3D segment mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] segments the list of segments.

    The Z value will be 0 if the points are 2D.
*/
template <class T, int dim>
void writeSegmeshPoly(std::ostream& os, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 2>>& segments)
{
    writePositionPoly(os, X);
    os << "POLYS\n";
    int count = 0;
    for (const Vector<int, 2>& seg : segments)
        os << ++count << ": " << seg(0) + 1 << " " << seg(1) + 1 << "\n"; // poly segment mesh is 1-indexed
    os << "END\n";
}

template <class T, int dim>
void writeSegmeshPoly(const std::string& filename, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 2>>& segments)
{
    std::ofstream fs;
    fs.open(filename);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    writeSegmeshPoly(fs, X, segments);
    fs.close();
}

template <class T, int dim>
void readSegmeshPoly(const std::string& filename, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments);
} // namespace ZIRAN

#endif
