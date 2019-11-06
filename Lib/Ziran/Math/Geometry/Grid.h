#ifndef _GRID_H_
#define _GRID_H_
#include <Ziran/CS/Util/Forward.h>
#include <Eigen/Core>
#include <Ziran/CS/Util/Debug.h>

namespace ZIRAN {
template <class T, int dim>
class Grid {
public:
    using TV = Vector<T, dim>;
    using TVI = Vector<int, dim>;
    int number_cells_per_dimension[dim];
    int number_nodes_per_dimension[dim];
    T mins[dim];
    T maxs[dim];
    T dx;

    Grid() {}

    Grid(const TVI& number_cells, const T dx_in, const TV& lowerCorner);

    void set(const TVI& number_cells, const T dx_in, const TV& lowerCorner);

    Grid(const int m_input, const T dx_input, const T xmin_input, const T ymin_input);

    Grid(const int m_input, const int n_input, const T dx_input, const T xmin_input, const T ymin_input);

    Grid(const int m_input, const int n_input, const int p_input, const T dx_input, const T xmin_input, const T ymin_input, const T zmin_input);

    void resize(const int m_input, const int n_input, const T dx_input, const T xmin_input, const T ymin_input);

    T Dx() const;

    int numberNodesPerDimension(const int d) const;

    int numberCellsPerDimension(const int d) const;

    int numberNodes() const;

    int numberCells() const;

    T gridVolume() const;

    TV Center_Of_Grid() const;

    void getAllPositions(StdVector<Vector<T, 2>>& x) const;

    void getAllPositions(StdVector<Vector<T, 3>>& x) const;

    inline int nodeIndex(const Vector<int, 2>& index) const
    {
        assert(dim == 2);
        assert(index(0) * number_nodes_per_dimension[1] + index(1) >= 0 && index(0) * number_nodes_per_dimension[1] + index(1) < number_nodes_per_dimension[0] * number_nodes_per_dimension[1]);
        return index(0) * number_nodes_per_dimension[1] + index(1);
    }

    inline int nodeIndex(Vector<int, 3>& index) const
    {
        ZIRAN_ASSERT(dim == 3);
        return index(0) * number_nodes_per_dimension[1] * number_nodes_per_dimension[2] + index(1) * number_nodes_per_dimension[2] + index(2);
    }

    TV node(const TVI& index) const;

    /**
     In the case of 2D, this returns the cell in a grid that covers all of the plane. I.e. the lower left cell in the grid is (0,0), but you could get
     (i,j) with i and j being any integers, negative, larger than number_cells_per_dimension etc. The user can decide how to handle cells not on the grid
     The case is similar for arbitrary dimension.
     */

    void cellContainingPoint(const TV& position, TV& lambda, TVI& index) const;
};
} // namespace ZIRAN
#endif
