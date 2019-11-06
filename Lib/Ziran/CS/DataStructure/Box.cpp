#include "Box.h"

namespace ZIRAN {

template <class T, int dim>
Box<T, dim>::Box()
{
    min_corner = EigenTV::Zero();
    max_corner = EigenTV::Zero();
    side_length = EigenTV::Zero();
    volume = 0;
}

template <class T, int dim>
Box<T, dim>::Box(const EigenTV& min_corner, const EigenTV& max_corner)
    : min_corner(min_corner)
    , max_corner(max_corner)
{
    side_length = max_corner - min_corner;
    volume = side_length.prod();
}

template <class T, int dim>
Box<T, dim>::Box(EigenTV&& min_corner, EigenTV&& max_corner)
    : min_corner(min_corner)
    , max_corner(max_corner)
{
    side_length = max_corner - min_corner;
    volume = side_length.prod();
}

template <class T, int dim>
Box<T, dim>::Box(const Box<T, dim>& box)
    : min_corner(box.min_corner)
    , max_corner(box.max_corner)
    , side_length(box.side_length)
    , volume(box.volume)
{
}

template <class T, int dim>
void Box<T, dim>::updateCorners(const EigenTV& new_min_corner, const EigenTV& new_max_corner)
{
    min_corner = new_min_corner;
    max_corner = new_max_corner;
    side_length = max_corner - min_corner;
    volume = side_length.prod();
}

template <int dim>
MaxExclusiveBoxIterator<dim>::MaxExclusiveBoxIterator(const Box<int, dim>& box)
    : box(box)
{
    index = box.min_corner;
}

template <int dim>
MaxExclusiveBoxIterator<dim>::MaxExclusiveBoxIterator(const MaxExclusiveBoxIterator<dim>& rhs)
    : box(rhs.box)
    , index(rhs.index)
{
}

template <int dim>
MaxExclusiveBoxIterator<dim>& MaxExclusiveBoxIterator<dim>::operator++()
{
    if (index == box.max_corner)
        return *this;
    else {
        int d = dim - 1;
        do {
            index(d)++;
            if (valid())
                return *this;
            index(d) = box.min_corner(d);
            d--;
        } while (d > -1);
        index = box.max_corner;
        return *this;
    }
}

template <int dim>
MaxExclusiveBoxIterator<dim> MaxExclusiveBoxIterator<dim>::operator++(int)
{
    MaxExclusiveBoxIterator ret(*this);
    ++(*this);
    return ret;
}

template <int dim>
bool MaxExclusiveBoxIterator<dim>::valid()
{
    bool ret = true;
    for (int d = 0; d < dim; ++d) {
        ret &= (index(d) >= box.min_corner(d) && index(d) < box.max_corner(d));
    }
    return ret;
}

template <int dim>
int& MaxExclusiveBoxIterator<dim>::operator()(int d)
{
    return index(d);
}

template Box<int, 2>::Box(Eigen::Matrix<int, 2, 1, 0, 2, 1> const&, Eigen::Matrix<int, 2, 1, 0, 2, 1> const&);
template Box<int, 3>::Box(Eigen::Matrix<int, 3, 1, 0, 3, 1> const&, Eigen::Matrix<int, 3, 1, 0, 3, 1> const&);
template MaxExclusiveBoxIterator<2>::MaxExclusiveBoxIterator(Box<int, 2> const&);
template MaxExclusiveBoxIterator<2>& MaxExclusiveBoxIterator<2>::operator++();
template bool MaxExclusiveBoxIterator<2>::valid();
template MaxExclusiveBoxIterator<3>::MaxExclusiveBoxIterator(Box<int, 3> const&);
template int& MaxExclusiveBoxIterator<3>::operator()(int);
template MaxExclusiveBoxIterator<3>& MaxExclusiveBoxIterator<3>::operator++();
template bool MaxExclusiveBoxIterator<3>::valid();
} // namespace ZIRAN
