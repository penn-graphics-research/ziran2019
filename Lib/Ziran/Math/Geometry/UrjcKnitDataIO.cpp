#include "UrjcKnitDataIO.h"
#include <Ziran/CS/Util/BinaryIO.h>

namespace ZIRAN {

template <class T, int dim>
void readUrjcKnitData(const std::string& filename, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments)
{
    T offset = 0.001;

    ZIRAN_ASSERT(dim == 3, "Urjc knit data only supports 3d");
    std::ifstream fs;
    fs.open(filename);
    ZIRAN_ASSERT(fs, "could not open ", filename);
    segments.clear();
    std::string line;
    Vector<T, dim> position_A;
    Vector<T, dim> position_B;
    Vector<T, dim> normal;
    ;

    Vector<int, 2> seg;
    bool beginning_of_yarn = true;

    int point_count = 0;
    while (std::getline(fs, line)) {
        std::stringstream ss(line);
        if (line.length() == 0)
            continue; // skip empty line

        if (line.length() >= (size_t)6 && line.substr(0, 6) == "o yarn") {
            beginning_of_yarn = true;
        }
        else if (line.length() >= (size_t)2 && line.substr(0, 2) == "v ") {
            if (beginning_of_yarn) {
                char dummy;
                ss >> dummy;
                for (size_t i = 0; i < dim; i++)
                    ss >> position_A(i);
                X.emplace_back(position_A);
                point_count++;
                beginning_of_yarn = false;
            }
            else {
                char dummy;
                ss >> dummy;
                for (size_t i = 0; i < dim; i++)
                    ss >> position_B(i);
                X.emplace_back(position_B);
                Vector<int, 2> segment;
                segment << point_count - 1, point_count;
                segments.push_back(segment);
                point_count++;
            }
        }
        else if (line.length() >= (size_t)2 && line.substr(0, 2) == "vn") {
            char dummy1, dummy2;
            ss >> dummy1 >> dummy2;
            for (size_t i = 0; i < dim; i++)
                ss >> normal(i);
            if (normal.squaredNorm() > 0)
                X[point_count - 1] += normal.normalized() * offset;
        }
        else
            ZIRAN_ASSERT(false);
    }
    fs.close();
}
} // namespace ZIRAN
