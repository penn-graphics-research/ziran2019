#include "UscHairDataIO.h"
#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/CS/Util/BinaryIO.h>

namespace ZIRAN {
template <class T, int dim>
void readUscHairData(std::istream& in, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments, const T percentage_to_keep)
{
    ZIRAN_ASSERT(dim == 3, "Usc Hair data only supports 3d");
    auto num_strands = readEntry<int32_t>(in);
    int particle_index = 0;
    RandomNumber<T> rand;
    for (int i = 0; i < num_strands; i++) {
        bool live = (rand.randReal(0, 1) <= percentage_to_keep);
        auto num_vertices = readEntry<int32_t>(in);
        if (live) {
            X.reserve(X.size() + num_vertices);
            segments.reserve(segments.size() + num_vertices - 1);
            for (int j = 0; j < num_vertices; j++) {
                X.emplace_back(readEntry<Vector<float, 3>>(in).template cast<T>().template head<dim>());
                if (j > 0)
                    segments.emplace_back(particle_index - 1, particle_index);
                particle_index++;
            }
        }
        else {
            for (int j = 0; j < num_vertices; j++) {
                readEntry<Vector<float, 3>>(in);
            }
        }
    }
}

template <class T, int dim>
void readUscHairData(const std::string& filename, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments, const T percentage_to_keep)
{
    std::ifstream fs;
    std::ios_base::iostate exceptionMask = fs.exceptions() | std::ios::failbit;
    fs.exceptions(exceptionMask);
    fs.open(filename, std::ios::in | std::ios::binary);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    readUscHairData(fs, X, segments, percentage_to_keep);
    fs.close();
}
template void readUscHairData<double, 2>(std::string const&, std::vector<Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 1, 0, 2, 1>>>&, std::vector<Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<int, 2, 1, 0, 2, 1>>>&, double);
template void readUscHairData<double, 3>(std::string const&, std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1, 0, 3, 1>>>&, std::vector<Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<int, 2, 1, 0, 2, 1>>>&, double);
template void readUscHairData<float, 2>(std::string const&, std::vector<Eigen::Matrix<float, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 2, 1, 0, 2, 1>>>&, std::vector<Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<int, 2, 1, 0, 2, 1>>>&, float);
template void readUscHairData<float, 3>(std::string const&, std::vector<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<float, 3, 1, 0, 3, 1>>>&, std::vector<Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::aligned_allocator<Eigen::Matrix<int, 2, 1, 0, 2, 1>>>&, float);
} // namespace ZIRAN
