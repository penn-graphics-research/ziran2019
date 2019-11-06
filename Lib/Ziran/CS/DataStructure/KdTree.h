#ifndef KD_TREE_H
#define KD_TREE_H

#include <Ziran/CS/Util/Forward.h>
#include <kdtree++/kdtree.hpp>
#include <iostream>

namespace ZIRAN {
/**
each KdTreeHelper object is a point
*/
template <int dim>
class KdTreeHelper {
    typedef Vector<double, dim> TV;

public:
    int id;
    TV pos;
    KdTreeHelper(int id, const TV& pos)
        : id(id)
        , pos(pos)
    {
    }

    KdTreeHelper(const KdTreeHelper<dim>& x)
    {
        id = x.id;
        pos = x.pos;
    }

    ~KdTreeHelper()
    {
    }

    double distance_to(const KdTreeHelper<dim>& x) const
    {
        return (pos - x.pos).norm();
    }

    inline bool operator==(const KdTreeHelper<dim>& A)
    {
        return this->pos == A.pos && id == A.id;
    }

    std::ostream& operator<<(std::ostream& out)
    {
        return out << " " << this->pos.transpose() << " ";
    }
};

/**
Wrapper class for libkdtree++
*/
template <int dim>
class KdTree {

    typedef Vector<double, dim> TV;

public:
    struct tac {
        typedef double result_type;
        double operator()(KdTreeHelper<dim> const& t, size_t k) const { return t.pos(k); }
    };

    typedef KDTree::KDTree<dim, KdTreeHelper<dim>, tac> TreeType;

    TreeType tree;
    KdTree()
    {
    }

    /**
    add an eigenvec point to the tree
    */
    template <class TV_IN>
    void addPoint(const int i, const TV_IN& new_p);

    void optimize();

    /**
    find the closest point to p in the tree
    */
    template <class TV_IN, class T_IN>
    void findNearest(const TV_IN& new_p, int& id, TV_IN& pos, T_IN& distance);
    /**
    erase point erase_p from the tree
    */
    template <class TV_IN>
    void erasePoint(const int i, const TV_IN& erase_p);
};
} // namespace ZIRAN

#endif
