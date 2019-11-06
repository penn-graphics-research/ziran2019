#include "PolyIO.h"

namespace ZIRAN {

template <class T, int dim>
void readSegmeshPoly(const std::string& filename, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments)
{
    std::ifstream fs;
    fs.open(filename);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    segments.clear();
    std::string line;
    Vector<T, dim> position;
    Vector<int, 2> seg;
    bool reading_positions = false;
    bool reading_segments = false;
    while (std::getline(fs, line)) {
        std::stringstream ss(line);
        if (line.length() == 0)
            continue; // skip empty line

        if (line.length() >= (size_t)6 && line.substr(0, 6) == "POINTS") {
            reading_positions = true;
            reading_segments = false;
        }
        else if (line.length() >= (size_t)5 && line.substr(0, 5) == "POLYS") {
            reading_segments = true;
            reading_positions = false;
        }
        else if (reading_positions) {
            int id;
            char dummy;
            ss >> id >> dummy;
            for (size_t i = 0; i < dim; i++)
                ss >> position(i);
            X.emplace_back(position);
        }
        else if (reading_segments) {
            int id;
            char dummy;
            ss >> id >> dummy;
            int previous_id = -1;
            int new_point;
            while (ss >> new_point) {
                if (previous_id != -1) {
                    seg(0) = previous_id - 1;
                    seg(1) = new_point - 1;
                    segments.emplace_back(seg);
                }
                previous_id = new_point;
            }
        }
    }
    fs.close();
}

template void readSegmeshPoly<double, 2>(std::string const&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&, std::vector<Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<int, 2, 1, 0, 2, 1>>>&);
template void readSegmeshPoly<double, 3>(std::string const&, std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1, 0, 3, 1>>>&, std::vector<Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<int, 2, 1, 0, 2, 1>>>&);
template void readSegmeshPoly<float, 2>(std::string const&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&, std::vector<Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<int, 2, 1, 0, 2, 1>>>&);
template void readSegmeshPoly<float, 3>(std::string const&, std::vector<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 3, 1, 0, 3, 1>>>&, std::vector<Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<int, 2, 1, 0, 2, 1>>>&);
} // namespace ZIRAN
