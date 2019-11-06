#ifndef FORWARD_H
#define FORWARD_H
#define EIGEN_MALLOC_ALREADY_ALIGNED 0
#include <Eigen/StdVector>
#include <vector>

namespace ZIRAN {

template <typename T, int dim>
using Vector = Eigen::Matrix<T, dim, 1, 0, dim, 1>;

template <typename T, int n, int m>
using Matrix = Eigen::Matrix<T, n, m, 0, n, m>;

template <typename Type>
using StdVector = std::vector<Type, Eigen::aligned_allocator<Type>>;

} // namespace ZIRAN
#endif
