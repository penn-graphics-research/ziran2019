#ifndef BOX_H
#define BOX_H
#include <Ziran/CS/Util/Forward.h>
#include <iostream>
namespace ZIRAN {
/**
Box
*/
template <class T, int dim>
class Box {
public:
    typedef Vector<T, dim> EigenTV;
    // minimum and maximum corner of range
    EigenTV min_corner;
    EigenTV max_corner;
    EigenTV side_length;
    T volume;
    Box();

    Box(const EigenTV& min_corner, const EigenTV& max_corner);

    Box(EigenTV&& min_corner, EigenTV&& max_corner);

    Box(const Box<T, dim>& box);

    void updateCorners(const EigenTV& new_min_corner, const EigenTV& new_max_corner);

    ~Box() {}
};

// to print min corner and max corner of box
template <class T, int dim>
inline std::ostream& operator<<(std::ostream& output, const Box<T, dim>& range)
{
    output << "min_corner = (" << range.min_corner.transpose() << ")\t";
    output << "max_corner = (" << range.max_corner.transpose() << ")\n";
    return output;
}

template <class T, int dim, class Func>
void maxExclusiveBoxOp(Func& func,
    const Box<T, dim>& box)
{
    Vector<T, dim> indices(box.min_corner);
    int d = 0;

    func(indices);
    while (d < dim) {
        if (indices(d) < box.max_corner(d) - 1) {
            indices(d)++;
            func(indices);
            d = 0;
        }
        else {
            indices(d) = box.min_corner(d);
            d++;
        }
    }
}
/**
maxExclusiveBoxOpWithEarlyExit
*/
template <class T, int dim, class Func>
void maxExclusiveBoxOpWithEarlyExit(Func& func,
    const Box<T, dim>& box, Func& exit_func)
{
    Vector<T, dim> indices(box.min_corner);
    int d = 0;

    func(indices);
    if (exit_func(indices))
        return;
    while (d < dim) {
        if (indices(d) < box.max_corner(d) - 1) {
            indices(d)++;
            func(indices);
            if (exit_func(indices))
                return;
            d = 0;
        }
        else {
            indices(d) = box.min_corner(d);
            d++;
        }
    }
}
/**
MaxExclusiveBoxIterator
*/
template <int dim>
class MaxExclusiveBoxIterator {
public:
    typedef Vector<int, dim> EigenIV;
    Box<int, dim> box;
    EigenIV index;
    MaxExclusiveBoxIterator(const Box<int, dim>& box);

    MaxExclusiveBoxIterator(const MaxExclusiveBoxIterator<dim>& rhs);

    // Prefix increment
    MaxExclusiveBoxIterator& operator++();

    // Postfix increment
    MaxExclusiveBoxIterator operator++(int);

    bool valid();

    int& operator()(int d);
};
} // namespace ZIRAN
#endif
