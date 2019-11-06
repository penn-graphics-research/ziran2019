#include <Ziran/Math/Geometry/ObjIO.h>
#include <Ziran/Math/Geometry/Particles.h>

namespace ZIRAN {
/**
    Read points data from obj file (3d) to std::vector (either 2d or 3d).
    The Z value will be ignored if 2D.
*/
template <class T, int dim>
void readPositionObj(std::istream& in, StdVector<Vector<T, dim>>& X)
{
    std::string line;
    Vector<T, dim> position;

    while (std::getline(in, line)) {
        std::stringstream ss(line);
        if (line[0] == 'v' && line[1] == ' ') {
            ss.ignore();
            for (size_t i = 0; i < dim; i++)
                ss >> position(i);
            X.emplace_back(position);
        }
    }
}

/**
    Read points data from obj file (3d) to std::vector (either 2d or 3d).
    The Z value will be ignored if 2D.
*/
template <class T, int dim>
void readPositionObj(const std::string& position_file, StdVector<Vector<T, dim>>& X)
{
    std::ifstream fs;
    fs.open(position_file);
    ZIRAN_ASSERT(fs, "could not open ", position_file);
    readPositionObj(fs, X);
    fs.close();
}

/**
    Write points data from std::vector (either 2d or 3d) to obj file (3d).
    The Z value will be 0 if 2D.
*/
template <class T, int dim>
void writePositionObj(std::ostream& os, const StdVector<Vector<T, dim>>& X)
{
    for (auto x : X) {
        os << "v";
        for (int i = 0; i < dim; i++)
            os << " " << x(i);
        if (dim == 2)
            os << " 0";
        os << "\n";
    }
}

/**
    Write points data from std::vector (either 2d or 3d) to obj file (3d).
    The Z value will be 0 if 2D.
*/
template <class T, int dim>
void writePositionObj(const std::string& filename, const StdVector<Vector<T, dim>>& X)
{
    std::ofstream fs;
    fs.open(filename);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    writePositionObj(fs, X);
    fs.close();
}

/**
    Read points data from obj file (3d) to ZIRAN:Particles (either 2d or 3d).
    The Z value will be ignored if ZIRAN:Particles is 2D.
*/
template <class T, int dim>
void readParticleObj(std::istream& in, Particles<T, dim>& particles)
{
    std::string line;
    Vector<T, dim> position;
    Vector<T, dim> velocity = Vector<T, dim>::Zero();
    T mass = NAN;
    auto ap = particles.appender();

    while (std::getline(in, line)) {
        std::stringstream ss(line);
        if (line[0] == 'v') {
            ss.ignore();
            for (size_t i = 0; i < dim; i++)
                ss >> position(i);
            ap.append(mass, position, velocity);
        }
    }
}

/**
    Read points data from obj file (3d) to ZIRAN:Particles (either 2d or 3d).
    The Z value will be ignored if ZIRAN:Particles is 2D.
*/
template <class T, int dim>
void readParticleObj(const std::string& particle_file, Particles<T, dim>& particles)
{
    std::ifstream fs;
    fs.open(particle_file);
    ZIRAN_ASSERT(fs, "could not open ", particle_file);
    readParticleObj(fs, particles);
    fs.close();
}

/**
    Write points data from ZIRAN:Particles (either 2d or 3d) to obj file (3d).
    The Z value will be 0 if ZIRAN:Particles is 2D.
*/
template <class T, int dim>
void writeParticleObj(std::ostream& os, const Particles<T, dim>& particles)
{
    auto iter = particles.iter(Particles<T, dim>::mass_name(), Particles<T, dim>::X_name(), Particles<T, dim>::V_name());
    for (; iter; ++iter) {
        os << "v";
        Vector<T, dim> X = iter.template get<1>();
        for (int i = 0; i < dim; i++)
            os << " " << X(i);
        if (dim == 2)
            os << " 0";
        os << "\n";
    }
}

/**
    Write points data from ZIRAN:Particles (either 2d or 3d) to obj file (3d).
    The Z value will be 0 if ZIRAN:Particles is 2D.
*/
template <class T, int dim>
void writeParticleObj(const std::string& filename, const Particles<T, dim>& particles)
{
    std::ofstream fs;
    fs.open(filename);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    writeParticleObj(fs, particles);
    fs.close();
}

/**
    Read a triangle mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] triangles the list of triangles.
*/
template <class T, int dim>
void readTrimeshObj(std::istream& is, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 3>>& triangles)
{
    triangles.clear();
    std::string line;
    Vector<T, dim> position;
    Vector<int, 3> tri;
    while (std::getline(is, line)) {
        std::stringstream ss(line);
        if (line[0] == 'v') {
            ss.ignore();
            for (size_t i = 0; i < dim; i++)
                ss >> position(i);
            X.emplace_back(position);
        }
        else if (line[0] == 'f') {
            ss.ignore();
            ss >> tri(0) >> tri(1) >> tri(2);
            tri -= Vector<int, 3>::Ones();
            triangles.emplace_back(tri);
        }
    }
}

/**
    Read a triangle mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] triangles the list of triangles.
*/
template <class T, int dim>
void readTrimeshObj(const std::string& particle_file, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 3>>& triangles)
{
    std::ifstream fs;
    fs.open(particle_file);
    ZIRAN_ASSERT(fs, "could not open ", particle_file);
    readTrimeshObj(fs, X, triangles);
    fs.close();
}

/**
    Write to a 3D triangle mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] triangles the list of triangles.

    The Z value will be 0 if the points are 2D.
*/
template <class T, int dim>
void writeTrimeshObj(std::ostream& os, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 3>>& triangles)
{
    writePositionObj(os, X);
    writeFacesObj(os, triangles);
}

/**
    Write to a 3D triangle mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] triangles the list of triangles.

    The Z value will be 0 if the points are 2D.
*/
template <class T, int dim>
void writeTrimeshObj(const std::string& filename, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 3>>& triangles)
{
    std::ofstream fs;
    fs.open(filename);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    writeTrimeshObj(fs, X, triangles);
    fs.close();
}

/**
    Read a quadrilateral mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] quads the list of quadrilaterals.
*/
template <class T, int dim>
void readQuadmeshObj(std::istream& is, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 4>>& quads)
{
    quads.clear();
    std::string line;
    Vector<T, dim> position;
    Vector<int, 4> quad;
    while (std::getline(is, line)) {
        std::stringstream ss(line);
        if (line[0] == 'v') {
            ss.ignore();
            for (size_t i = 0; i < dim; i++)
                ss >> position(i);
            X.emplace_back(position);
        }
        else if (line[0] == 'f') {
            ss.ignore();
            ss >> quad(0) >> quad(1) >> quad(2) >> quad(3);
            quad -= Vector<int, 4>::Ones();
            quads.emplace_back(quad);
        }
    }
}

/**
    Read a quadrilateral mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] quads the list of quadrilaterals.
*/
template <class T, int dim>
void readQuadmeshObj(const std::string& particle_file, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 4>>& quads)
{
    std::ifstream fs;
    fs.open(particle_file);
    ZIRAN_ASSERT(fs, "could not open ", particle_file);
    readQuadmeshObj(fs, X, quads);
    fs.close();
}
/**
    Write to a 3D quadrilateral mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] quads the list of quadrilaterals.

    The Z value will be 0 if the points are 2D.
*/

template <int vertex_per_face>
void writeFacesObj(std::ostream& os, const StdVector<Vector<int, vertex_per_face>>& faces)
{
    for (const Vector<int, vertex_per_face>& face : faces) {
        os << "f";
        for (int i = 0; i < vertex_per_face; i++)
            os << ' ' << face(i) + 1;
        os << '\n';
    }
}

template <class T, int dim>
void writeQuadmeshObj(std::ostream& os, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 4>>& quads)
{
    writePositionObj(os, X);
    writeFacesObj(os, quads);
}

/**
    Write to a 3D quadrilateral mesh.

    \param[in] particles the list of points (2D or 3D)
    \param[in] quads the list of quadrilaterals.

    The Z value will be 0 if the points are 2D.
*/
template <class T, int dim>
void writeQuadmeshObj(const std::string& filename, const StdVector<Vector<T, dim>>& X, const StdVector<Vector<int, 4>>& quads)
{
    std::ofstream fs;
    fs.open(filename);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    writeQuadmeshObj(fs, X, quads);
    fs.close();
}

template void readParticleObj<float, 2>(std::istream&, Particles<float, 2>&);
template void readParticleObj<float, 3>(std::istream&, Particles<float, 3>&);
template void readParticleObj<float, 2>(std::string const&, Particles<float, 2>&);
template void readParticleObj<float, 3>(std::string const&, Particles<float, 3>&);
template void readPositionObj<float, 1>(std::istream&, StdVector<Vector<float, 1>>&);
template void readPositionObj<float, 2>(std::istream&, StdVector<Vector<float, 2>>&);
template void readPositionObj<float, 3>(std::istream&, StdVector<Vector<float, 3>>&);
template void readPositionObj<float, 1>(std::string const&, StdVector<Vector<float, 1>>&);
template void readPositionObj<float, 2>(std::string const&, StdVector<Vector<float, 2>>&);
template void readPositionObj<float, 3>(std::string const&, StdVector<Vector<float, 3>>&);
template void readTrimeshObj<float, 2>(std::istream&, StdVector<Vector<float, 2>>&, StdVector<Vector<int, 3>>&);
template void readTrimeshObj<float, 3>(std::istream&, StdVector<Vector<float, 3>>&, StdVector<Vector<int, 3>>&);
template void readTrimeshObj<float, 2>(std::string const&, StdVector<Vector<float, 2>>&, StdVector<Vector<int, 3>>&);
template void readTrimeshObj<float, 3>(std::string const&, StdVector<Vector<float, 3>>&, StdVector<Vector<int, 3>>&);
template void writeParticleObj<float, 2>(std::ostream&, Particles<float, 2> const&);
template void writeParticleObj<float, 3>(std::ostream&, Particles<float, 3> const&);
template void writeParticleObj<float, 2>(std::string const&, Particles<float, 2> const&);
template void writeParticleObj<float, 3>(std::string const&, Particles<float, 3> const&);
template void writeTrimeshObj<float, 2>(std::ostream&, StdVector<Vector<float, 2>> const&, StdVector<Vector<int, 3>> const&);
template void writeTrimeshObj<float, 3>(std::ostream&, StdVector<Vector<float, 3>> const&, StdVector<Vector<int, 3>> const&);
template void writeTrimeshObj<float, 2>(std::string const&, StdVector<Vector<float, 2>> const&, StdVector<Vector<int, 3>> const&);
template void writeTrimeshObj<float, 3>(std::string const&, StdVector<Vector<float, 3>> const&, StdVector<Vector<int, 3>> const&);
template void writeTrimeshObj<double, 2>(std::string const&, StdVector<Vector<double, 2>> const&, StdVector<Vector<int, 3>> const&);
template void writeTrimeshObj<double, 3>(std::string const&, StdVector<Vector<double, 3>> const&, StdVector<Vector<int, 3>> const&);
template void readParticleObj<double, 2>(std::istream&, Particles<double, 2>&);
template void readParticleObj<double, 3>(std::istream&, Particles<double, 3>&);
template void readParticleObj<double, 2>(std::string const&, Particles<double, 2>&);
template void readParticleObj<double, 3>(std::string const&, Particles<double, 3>&);
template void readPositionObj<double, 1>(std::istream&, StdVector<Vector<double, 1>>&);
template void readPositionObj<double, 2>(std::istream&, StdVector<Vector<double, 2>>&);
template void readPositionObj<double, 3>(std::istream&, StdVector<Vector<double, 3>>&);
template void readPositionObj<double, 1>(std::string const&, StdVector<Vector<double, 1>>&);
template void readPositionObj<double, 2>(std::string const&, StdVector<Vector<double, 2>>&);
template void readPositionObj<double, 3>(std::string const&, StdVector<Vector<double, 3>>&);
template void readTrimeshObj<double, 2>(std::istream&, StdVector<Vector<double, 2>>&, StdVector<Vector<int, 3>>&);
template void readTrimeshObj<double, 3>(std::istream&, StdVector<Vector<double, 3>>&, StdVector<Vector<int, 3>>&);
template void readTrimeshObj<double, 2>(std::string const&, StdVector<Vector<double, 2>>&, StdVector<Vector<int, 3>>&);
template void readTrimeshObj<double, 3>(std::string const&, StdVector<Vector<double, 3>>&, StdVector<Vector<int, 3>>&);
template void writeParticleObj<double, 2>(std::ostream&, Particles<double, 2> const&);
template void writeParticleObj<double, 3>(std::ostream&, Particles<double, 3> const&);
template void writeTrimeshObj<double, 2>(std::ostream&, StdVector<Vector<double, 2>> const&, StdVector<Vector<int, 3>> const&);
template void writeTrimeshObj<double, 3>(std::ostream&, StdVector<Vector<double, 3>> const&, StdVector<Vector<int, 3>> const&);
template void readQuadmeshObj<float, 2>(std::istream&, StdVector<Vector<float, 2>>&, StdVector<Vector<int, 4>>&);
template void readQuadmeshObj<float, 3>(std::istream&, StdVector<Vector<float, 3>>&, StdVector<Vector<int, 4>>&);
template void readQuadmeshObj<float, 2>(std::string const&, StdVector<Vector<float, 2>>&, StdVector<Vector<int, 4>>&);
template void readQuadmeshObj<float, 3>(std::string const&, StdVector<Vector<float, 3>>&, StdVector<Vector<int, 4>>&);
template void readQuadmeshObj<double, 2>(std::istream&, StdVector<Vector<double, 2>>&, StdVector<Vector<int, 4>>&);
template void readQuadmeshObj<double, 3>(std::istream&, StdVector<Vector<double, 3>>&, StdVector<Vector<int, 4>>&);
template void readQuadmeshObj<double, 2>(std::string const&, StdVector<Vector<double, 2>>&, StdVector<Vector<int, 4>>&);
template void readQuadmeshObj<double, 3>(std::string const&, StdVector<Vector<double, 3>>&, StdVector<Vector<int, 4>>&);
template void writeQuadmeshObj<float, 2>(std::ostream&, StdVector<Vector<float, 2>> const&, StdVector<Vector<int, 4>> const&);
template void writeQuadmeshObj<float, 3>(std::ostream&, StdVector<Vector<float, 3>> const&, StdVector<Vector<int, 4>> const&);
template void writeQuadmeshObj<float, 2>(std::string const&, StdVector<Vector<float, 2>> const&, StdVector<Vector<int, 4>> const&);
template void writeQuadmeshObj<float, 3>(std::string const&, StdVector<Vector<float, 3>> const&, StdVector<Vector<int, 4>> const&);
template void writeQuadmeshObj<double, 2>(std::string const&, StdVector<Vector<double, 2>> const&, StdVector<Vector<int, 4>> const&);
template void writeQuadmeshObj<double, 3>(std::string const&, StdVector<Vector<double, 3>> const&, StdVector<Vector<int, 4>> const&);
template void writeQuadmeshObj<double, 2>(std::ostream&, StdVector<Vector<double, 2>> const&, StdVector<Vector<int, 4>> const&);
template void writeQuadmeshObj<double, 3>(std::ostream&, StdVector<Vector<double, 3>> const&, StdVector<Vector<int, 4>> const&);
template void writeFacesObj<3>(std::ostream& os, const StdVector<Vector<int, 3>>& faces);
template void writeFacesObj<4>(std::ostream& os, const StdVector<Vector<int, 4>>& faces);
template void writePositionObj<float, 2>(std::ostream& os, const StdVector<Vector<float, 2>>& faces);
template void writePositionObj<float, 3>(std::ostream& os, const StdVector<Vector<float, 3>>& faces);
template void writePositionObj<double, 2>(std::ostream& os, const StdVector<Vector<double, 2>>& faces);
template void writePositionObj<double, 3>(std::ostream& os, const StdVector<Vector<double, 3>>& faces);
} // namespace ZIRAN
