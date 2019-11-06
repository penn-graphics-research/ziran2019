#include <Ziran/Math/Geometry/TetGenIO.h>

namespace ZIRAN {
template <class T, int dim>
void writeNode(std::ostream& out, const StdVector<Vector<T, dim>>& X)
{
    out << "# Node count, " << dim << " dim, no attribute, no boundary marker\n";
    out << X.size() << ' ' << dim << " 0 0\n";
    out << "# Node index, node coordinates\n";
    for (size_t i = 0; i < X.size(); i++) {
        out << i;
        for (int d = 0; d < dim; d++)
            out << ' ' << X[i](d);
        out << '\n';
    }
}
template void writeNode<float, 2>(std::ostream& os, const StdVector<Vector<float, 2>>& faces);
template void writeNode<float, 3>(std::ostream& os, const StdVector<Vector<float, 3>>& faces);
template void writeNode<double, 2>(std::ostream& os, const StdVector<Vector<double, 2>>& faces);
template void writeNode<double, 3>(std::ostream& os, const StdVector<Vector<double, 3>>& faces);
} // namespace ZIRAN
