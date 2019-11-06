#include "KdTree.h"

namespace ZIRAN {
/**
  add an eigenvec point to the tree
  */
template <int dim>
template <class TV_IN>
void KdTree<dim>::addPoint(const int i, const TV_IN& new_p)
{
    TV p;
    for (int i = 0; i < dim; i++)
        p(i) = (double)(new_p(i));
    tree.insert(KdTreeHelper<dim>(i, p));
}

template <int dim>
void KdTree<dim>::optimize()
{
    tree.optimise();
}

/**
  find the closest point to p in the tree
  */
template <int dim>
template <class TV_IN, class T_IN>
void KdTree<dim>::findNearest(const TV_IN& new_p, int& id, TV_IN& pos, T_IN& distance)
{
    TV p;
    for (int i = 0; i < dim; i++)
        p(i) = (double)(new_p(i));
    auto found = tree.find_nearest(KdTreeHelper<dim>(-1, p));
    id = found.first->id;
    for (int i = 0; i < dim; i++)
        pos(i) = (T_IN)((found.first->pos)(i));
    distance = (T_IN)found.second;
}

/**
  erase point erase_p from the tree
  */
template <int dim>
template <class TV_IN>
void KdTree<dim>::erasePoint(const int i, const TV_IN& erase_p)
{
    TV p;
    for (int i = 0; i < dim; i++)
        p(i) = (double)(erase_p(i));
    tree.erase(KdTreeHelper<dim>(i, p));
}
template class KdTree<2>;
template class KdTree<3>;

template void KdTree<2>::addPoint<Eigen::Matrix<float, 2, 1, 0, 2, 1>>(int, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&);
template void KdTree<2>::addPoint<Eigen::Matrix<double, 2, 1, 0, 2, 1>>(int, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&);
template void KdTree<2>::findNearest<Eigen::Matrix<float, 2, 1, 0, 2, 1>, float>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, int&, Eigen::Matrix<float, 2, 1, 0, 2, 1>&, float&);
template void KdTree<2>::findNearest<Eigen::Matrix<double, 2, 1, 0, 2, 1>, double>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, int&, Eigen::Matrix<double, 2, 1, 0, 2, 1>&, double&);
template void KdTree<3>::addPoint<Eigen::Matrix<double, 3, 1, 0, 3, 1>>(int, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&);
template void KdTree<3>::addPoint<Eigen::Matrix<float, 3, 1, 0, 3, 1>>(int, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&);
template void KdTree<3>::erasePoint<Eigen::Matrix<double, 3, 1, 0, 3, 1>>(int, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&);
template void KdTree<3>::findNearest<Eigen::Matrix<double, 3, 1, 0, 3, 1>, double>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, int&, Eigen::Matrix<double, 3, 1, 0, 3, 1>&, double&);
template void KdTree<3>::findNearest<Eigen::Matrix<float, 3, 1, 0, 3, 1>, float>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, int&, Eigen::Matrix<float, 3, 1, 0, 3, 1>&, float&);
} // namespace ZIRAN
