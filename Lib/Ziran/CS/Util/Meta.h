/**
   Template metaprogramming utilities
*/
#ifndef META_H
#define META_H
#include <Ziran/CS/Util/Forward.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

namespace ZIRAN {

template <typename T>
struct IsEigenMatrix : std::false_type {
};

template <class T, int m, int n>
struct IsEigenMatrix<Matrix<T, m, n>> : std::true_type {
};

template <typename T>
struct IsEigenSparseMatrix : std::false_type {
};

template <class T>
struct IsEigenSparseMatrix<Eigen::SparseMatrix<T>> : std::true_type {
};

template <typename T>
struct IsEigenTriplet : std::false_type {
};

template <class T, class Index>
struct IsEigenTriplet<Eigen::Triplet<T, Index>> : std::true_type {
};

template <typename T>
struct IsUniquePtr : std::false_type {
};

template <typename T>
struct IsUniquePtr<std::unique_ptr<T>> : std::true_type {
};

template <typename T>
struct IsSharedPtr : std::false_type {
};

template <typename T>
struct IsSharedPtr<std::shared_ptr<T>> : std::true_type {
};

template <typename T>
struct IsStdVector : std::false_type {
};

template <typename T, typename Alloc>
struct IsStdVector<std::vector<T, Alloc>> : std::true_type {
};

namespace INTERNAL {
using namespace std;
template <class T, class Enable = void>
struct ScalarTypeHelper {
    using type = typename T::Scalar;
};
template <class T>
struct ScalarTypeHelper<T, enable_if_t<is_arithmetic<T>::value>> {
    using type = T;
};
template <class T, class Enable = void>
struct ValueTypeHelper {
    using type = typename T::value_type;
};
template <class T>
struct ValueTypeHelper<T, enable_if_t<!IsStdVector<T>::value>> {
    using type = T;
};
} // namespace INTERNAL

template <class T>
using ScalarType = typename INTERNAL::ScalarTypeHelper<T>::type;

template <class T>
using ValueType = typename INTERNAL::ValueTypeHelper<T>::type;

template <class MatrixType>
constexpr bool isSize(int m, int n)
{
    return MatrixType::RowsAtCompileTime == m && MatrixType::ColsAtCompileTime == n;
}

/**
  Helper struct for converting a std::vector to a 
  Eigen::Map<Vector<T, Eigen::Dynamic>>
  with the same backing memory.
  */
template <class T>
struct ToVectorMap {
    T* data;
    size_t size;

public:
    ToVectorMap(StdVector<T>& array)
        : data(array.data())
        , size(array.size())
    {
    }

    template <int d>
    ToVectorMap(StdVector<Eigen::Matrix<T, d, 1, 0, d, 1>>& array)
        : data(array.data()->data())
        , size(d * array.size())
    {
    }

    operator Eigen::Map<Vector<T, Eigen::Dynamic>>()
    {
        return Eigen::Map<Vector<T, Eigen::Dynamic>>(data, size, 1);
    }
};

/**
  Helper struct for converting a std::vector to a eigen matrix with the entries as columns
  with the same backing memory.

  converts to Eigen::Map<Matrix<T, Eigen::Dynamic, d>>
  */
template <class T, int d>
struct ToMatrixColMap {
    T* data;
    size_t cols;

public:
    ToMatrixColMap(StdVector<Vector<T, d>>& array)
        : data(array.data()->data())
        , cols(array.size())
    {
    }

    operator Eigen::Map<Matrix<T, d, Eigen::Dynamic>>()
    {
        return Eigen::Map<Matrix<T, d, Eigen::Dynamic>>(data, d, cols);
    }
};

/**
  Helper struct for converting a std::vector to a eigen matrix with the entries as rows
  with the same backing memory.

  converts to Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, d, RowMajor>>
  */
template <class T, int d>
struct ToMatrixRowMap {
    T* data;
    size_t rows;

public:
    ToMatrixRowMap(StdVector<Vector<T, d>>& array)
        : data(array.data()->data())
        , rows(array.size())
    {
    }

    operator Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, d, Eigen::RowMajor>>()
    {
        return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, d, Eigen::RowMajor>>(data, rows, d);
    }
};

/**
  The Call struct is useful for variadic templates
  it's designed for calling a function on all
  values in an argument pack
  `Call{foo(args)...};`
  or (if foo returns void)
  `Call{(foo(args),0)...};`
*/
struct Call {
    template <typename... T>
    /**
        The constructor itself doesn't do anything
        it only exists to evaluate its arguments
    */
    Call(T...)
    {
    }
};

/* 
   Tests whether T has a typedef Type
*/
template <class T>
struct HasType {
    template <class U>
    static char (&test(typename U::Type const*))[1];
    template <class U>
    static char (&test(...))[2];
    static const bool value = (sizeof(test<T>(0)) == 1);
};

template <typename T, typename U>
constexpr inline bool NonSelf()
{
    using DecayedT = typename std::decay<T>::type;
    return !std::is_same<DecayedT, U>::value
        && !std::is_base_of<U, DecayedT>::value;
}
} // namespace ZIRAN
#endif
