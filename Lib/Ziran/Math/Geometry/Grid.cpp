#include "Grid.h"
#include <Ziran/CS/Util/Debug.h>

namespace ZIRAN {

template <class T, int dim>
Grid<T, dim>::Grid(const TVI& number_cells, const T dx_in, const TV& lowerCorner)
{
    set(number_cells, dx_in, lowerCorner);
}

template <class T, int dim>
void Grid<T, dim>::set(const TVI& number_cells, const T dx_in, const TV& lowerCorner)
{
    dx = dx_in;
    for (int d = 0; d < dim; d++) {
        number_cells_per_dimension[d] = number_cells(d);
        number_nodes_per_dimension[d] = number_cells(d) + 1;
        mins[d] = lowerCorner(d);
        maxs[d] = mins[d] + (T)number_cells_per_dimension[d] * dx;
    }
}

template <class T, int dim>
Grid<T, dim>::Grid(const int m_input, const T dx_input, const T xmin_input, const T ymin_input)
    : dx(dx_input)
{
    assert(dim == 2);
    number_cells_per_dimension[0] = m_input;
    number_cells_per_dimension[1] = m_input;
    for (int d = 0; d < dim; d++)
        number_nodes_per_dimension[d] = number_cells_per_dimension[d] + 1;
    mins[0] = xmin_input;
    mins[1] = ymin_input;
    for (int d = 0; d < dim; d++)
        maxs[d] = mins[d] + (T)number_cells_per_dimension[d] * dx;
}

template <class T, int dim>
Grid<T, dim>::Grid(const int m_input, const int n_input, const T dx_input, const T xmin_input, const T ymin_input)
    : dx(dx_input)
{
    assert(dim == 2);
    number_cells_per_dimension[0] = m_input;
    number_cells_per_dimension[1] = n_input;
    for (int d = 0; d < dim; d++)
        number_nodes_per_dimension[d] = number_cells_per_dimension[d] + 1;
    mins[0] = xmin_input;
    mins[1] = ymin_input;
    for (int d = 0; d < dim; d++)
        maxs[d] = mins[d] + (T)number_cells_per_dimension[d] * dx;
}

template <class T, int dim>
Grid<T, dim>::Grid(const int m_input, const int n_input, const int p_input, const T dx_input, const T xmin_input, const T ymin_input, const T zmin_input)
    : dx(dx_input)
{
    ZIRAN_ASSERT(dim == 3);
    number_cells_per_dimension[0] = m_input;
    number_cells_per_dimension[1] = n_input;
    number_cells_per_dimension[2] = p_input;
    for (int d = 0; d < dim; d++)
        number_nodes_per_dimension[d] = number_cells_per_dimension[d] + 1;
    mins[0] = xmin_input;
    mins[1] = ymin_input;
    mins[2] = zmin_input;
    for (int d = 0; d < dim; d++)
        maxs[d] = mins[d] + (T)number_cells_per_dimension[d] * dx;
}

template <class T, int dim>
void Grid<T, dim>::resize(const int m_input, const int n_input, const T dx_input, const T xmin_input, const T ymin_input)
{
    assert(dim == 2);
    number_cells_per_dimension[0] = m_input;
    number_cells_per_dimension[1] = n_input;
    for (int d = 0; d < dim; d++)
        number_nodes_per_dimension[d] = number_cells_per_dimension[d] + 1;
    dx = dx_input;
    mins[0] = xmin_input;
    mins[1] = ymin_input;
    for (int d = 0; d < dim; d++)
        maxs[d] = mins[d] + (T)number_cells_per_dimension[d] * dx;
}

template <class T, int dim>
T Grid<T, dim>::Dx() const
{
    return dx;
}

template <class T, int dim>
int Grid<T, dim>::numberNodesPerDimension(const int d) const
{
    assert(d < dim);
    return number_nodes_per_dimension[d];
}

template <class T, int dim>
int Grid<T, dim>::numberCellsPerDimension(const int d) const
{
    assert(d < dim);
    return number_cells_per_dimension[d];
}

template <class T, int dim>
int Grid<T, dim>::numberNodes() const
{
    int result = 1;
    for (int d = 0; d < dim; d++)
        result *= number_nodes_per_dimension[d];
    return result;
}

template <class T, int dim>
int Grid<T, dim>::numberCells() const
{
    int result = 1;
    for (int d = 0; d < dim; d++)
        result *= number_cells_per_dimension[d];
    return result;
}

template <class T, int dim>
T Grid<T, dim>::gridVolume() const
{
    return numberCells() * std::pow(dx, dim);
}

template <class T, int dim>
typename Grid<T, dim>::TV Grid<T, dim>::Center_Of_Grid() const
{
    TV result;
    for (int d = 0; d < dim; d++)
        result(d) = (T).5 * (mins[d] + maxs[d]);
    return result;
}

template <class T, int dim>
void Grid<T, dim>::getAllPositions(StdVector<Vector<T, 2>>& x) const
{
    x.clear();
    for (int i = 0; i < number_nodes_per_dimension[0]; i++) {
        for (int j = 0; j < number_nodes_per_dimension[1]; j++) {
            Vector<T, 2> xi;
            xi(0) = i * dx + mins[0];
            xi(1) = j * dx + mins[1];
            x.push_back(xi);
        }
    }
}

template <class T, int dim>
void Grid<T, dim>::getAllPositions(StdVector<Vector<T, 3>>& x) const
{
    ZIRAN_ASSERT(dim == 3);
    x.clear();
    for (int i = 0; i < number_nodes_per_dimension[0]; i++) {
        for (int j = 0; j < number_nodes_per_dimension[1]; j++) {
            for (int k = 0; k < number_nodes_per_dimension[2]; k++) {
                Vector<T, 3> xi;
                xi(0) = i * dx + mins[0];
                xi(1) = j * dx + mins[1];
                xi(2) = k * dx + mins[2];
                x.push_back(xi);
            }
        }
    }
}

template <class T, int dim>
typename Grid<T, dim>::TV Grid<T, dim>::node(const TVI& index) const
{
    TV result;
    for (int d = 0; d < dim; d++)
        result(d) = (T)index(d) * dx + mins[d];
    return result;
}

/**
     In the case of 2D, this returns the cell in a grid that covers all of the plane. I.e. the lower left cell in the grid is (0,0), but you could get
     (i,j) with i and j being any integers, negative, larger than number_cells_per_dimension etc. The user can decide how to handle cells not on the grid
     The case is similar for arbitrary dimension.
     */

template <class T, int dim>
void Grid<T, dim>::cellContainingPoint(const TV& position, TV& lambda, TVI& index) const
{

    for (int d = 0; d < dim; d++) {
        T normalized = (position(d) - mins[d]) / dx;
        if (normalized < 0) {
            index(d) = (int)normalized - 1;
        }
        else {
            index(d) = (int)normalized;
        }
    }

    for (int d = 0; d < dim; d++)
        lambda(d) = (position(d) - (dx * (T)index(d) + mins[d])) / dx;
}

template class Grid<double, 2>;
template class Grid<double, 3>;
template class Grid<float, 2>;
template class Grid<float, 3>;
} // namespace ZIRAN
