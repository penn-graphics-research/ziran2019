#ifndef MULTI_ARRAY_H
#define MULTI_ARRAY_H

#include <Ziran/CS/Util/Forward.h>
#include <Eigen/Core>
#include <type_traits>

namespace ZIRAN {

template <class T, int dim, class enable = void>
class MultiArray;

template <class T, int dim>
class MultiArray<T, dim, std::enable_if_t<dim == 2>> {
public:
    typedef Vector<int, 2> IV;

    Matrix<T, Eigen::Dynamic, Eigen::Dynamic> data;
    Vector<int, 2> size;

    MultiArray()
    {
    }

    MultiArray(const Vector<int, 2>& size_in)
    {
        size = size_in;
        data.resize(size(0), size(1));
    }

    MultiArray(const Vector<int, 2>& size_in, T val)
    {
        size = size_in;
        data.resize(size(0), size(1));
        data.setConstant(val);
    }

    ~MultiArray() {}

    void resize(const Vector<int, 2>& size_in, T val)
    {
        size = size_in;
        data.resize(size(0), size(1));
        data.setConstant(val);
    }

    template <class IV>
    T& operator()(const IV& index)
    {
        return data(index(0), index(1));
    }

    template <class IV>
    const T& operator()(const IV& index) const
    {
        return data(index(0), index(1));
    }

    int count(const T& val)
    {
        int c = 0;
        for (int i = 0; i < size(0); ++i)
            for (int j = 0; j < size(0); ++j)
                c += (data(i, j) == val);
        return c;
    }

    template <class Func>
    void loop(Func&& func)
    {
        for (int i = 0; i < size(0); ++i) {
            for (int j = 0; j < size(1); ++j) {
                IV ij;
                ij << i, j;
                T val = data(i, j);
                func(ij, val);
            }
        }
    }
};

template <class T, int dim>
class MultiArray<T, dim, std::enable_if_t<dim == 3>> {
public:
    typedef Vector<int, 3> IV;
    T* data;
    IV size;

    MultiArray()
        : data(NULL)
    {
    }

    MultiArray(const IV& size_in)
    {
        size = size_in;
        //ZIRAN_ASSERT(size.prod() > 0, "ERROR in MultiArray<T, dim>::MultiArray(const IV& size_in), size.prod() = 0!");
        data = new T[size.prod()];
    }

    MultiArray(const IV& size_in, T val)
    {
        size = size_in;
        //ZIRAN_ASSERT(size.prod() > 0, "ERROR in MultiArray<T, dim>::MultiArray(const IV& size_in, T val), size.prod() = 0!");
        data = new T[size.prod()];
        for (int i = 0; i < size.prod(); ++i)
            data[i] = val;
    }

    MultiArray(const MultiArray<T, dim>& multi)
    {
        size = multi.size;
        //ZIRAN_ASSERT(size.prod() > 0, "ERROR in MultiArray<T, dim>::MultiArray(const MultiArray<T, dim>& multi), size.prod() = 0!");
        data = new T[size.prod()];
        for (int i = 0; i < size.prod(); ++i)
            data[i] = multi.data[i];
    }

    ~MultiArray()
    {
        if (data)
            delete[] data;
    }

    void resize(const IV& size_in, T val)
    {
        size = size_in;
        //ZIRAN_ASSERT(size.prod() > 0, "ERROR in MultiArray<T, dim>::resize(const IV& size_in, T val), size.prod() = 0!");
        data = new T[size.prod()];
        for (int i = 0; i < size.prod(); ++i)
            data[i] = val;
    }

    template <class IV>
    inline int id(const IV& index) const
    {
        return index(0) * size(1) * size(2) + index(1) * size(2) + index(2);
    }

    template <class IV>
    T& operator()(const IV& index)
    {
        return data[id(index)];
    }

    template <class IV>
    const T& operator()(const IV& index) const
    {
        return data[id(index)];
    }

    // number of entries that equal to val
    int count(const T& val)
    {
        int c = 0;
        for (int i = 0; i < size.prod(); ++i)
            c += (data[i] == val);
        return c;
    }

    template <class Func>
    void loop(Func&& func)
    {
        for (int i = 0; i < size(0); ++i) {
            for (int j = 0; j < size(1); ++j) {
                for (int k = 0; k < size(2); ++k) {
                    IV ijk;
                    ijk << i, j, k;
                    T val = data[id(ijk)];
                    func(ijk, val);
                }
            }
        }
    }
};

template <class T, int dim, int res, class enable = void>
class CubeMultiArray;

template <class T, int dim, int res>
class CubeMultiArray<T, dim, res, std::enable_if_t<dim == 2>> {
public:
    T data[res * res];

    template <class IV>
    inline int id(const IV& index) const
    {
        return index[0] * res + index[1];
    }

    template <class IV>
    T& operator()(const IV& index)
    {
        assert((index[0] >= 0) && (index[1] >= 0) && (index[0] < res) && (index[1] < res));
        return data[id(index)];
    }

    template <class IV>
    const T& operator()(const IV& index) const
    {
        assert((index[0] >= 0) && (index[1] >= 0) && (index[0] < res) && (index[1] < res));
        return data[id(index)];
    }

    void fill(const T& z)
    {
        for (int i = 0; i < res * res; i++)
            data[i] = z;
    }
};

template <class T, int dim, int res>
class CubeMultiArray<T, dim, res, std::enable_if_t<dim == 3>> {
public:
    T data[res * res * res];

    template <class IV>
    inline int id(const IV& index) const
    {
        return index[0] * res * res + index[1] * res + index[2];
    }

    template <class IV>
    T& operator()(const IV& index)
    {
        assert((index[0] >= 0) && (index[1] >= 0) && (index[2] >= 0) && (index[0] < res) && (index[1] < res) && (index[2] < res));
        return data[id(index)];
    }

    template <class IV>
    const T& operator()(const IV& index) const
    {
        assert((index[0] >= 0) && (index[1] >= 0) && (index[2] >= 0) && (index[0] < res) && (index[1] < res) && (index[2] < res));
        return data[id(index)];
    }

    void fill(const T& z)
    {
        for (int i = 0; i < res * res * res; i++)
            data[i] = z;
    }
};
} // namespace ZIRAN
#endif
