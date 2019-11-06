#ifndef BINARY_IO_H
#define BINARY_IO_H
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/Meta.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <fstream>

namespace ZIRAN {

template <class Type, class enable = void>
struct RW;

template <class Type>
struct TriviallyCopyableTag {
};
template <class T, int m, int n>
struct EigenMatrixTag {
};
template <class T>
struct EigenSparseMatrixTag {
};
template <class T>
struct UniquePtrTag {
};
template <class T>
struct BasicStringTag {
};
template <class T>
struct CustomTypeTag {
};
template <class T>
struct NoWriteTag {
};

template <class T, size_t size>
struct StdArrayTag {
};

template <class T>
struct RW<T, std::enable_if_t<std::is_trivially_copyable<T>::value>> {
    using Tag = TriviallyCopyableTag<T>;
};

template <class T, int m, int n, int flags, int mm, int nn>
struct RW<Eigen::Matrix<T, m, n, flags, mm, nn>> {
    using Tag = EigenMatrixTag<T, m, n>;
};

template <class T>
struct RW<Eigen::SparseMatrix<T>> {
    using Tag = EigenSparseMatrixTag<T>;
};

template <class T>
struct RW<std::unique_ptr<T>> {
    using Tag = UniquePtrTag<T>;
};

template <class T>
struct RW<std::basic_string<T>> {
    using Tag = BasicStringTag<T>;
};

// template <class T, size_t size>
// struct RW<std::array<T, size>, std::enable_if_t<std::is_trivially_copyable<T>::value>> {
//     using Tag = StdArrayTag<T, size>;
// };

template <size_t size, typename Type>
struct RW<std::array<Type, size>, std::enable_if_t<!std::is_trivially_copyable<Type>::value>> {
    using Tag = StdArrayTag<Type, size>;
};

template <class Type>
void writeEntry(std::ostream& out, const Type& x);

template <class Type>
inline Type readEntry(std::istream& in);

template <class T>
void writeSTDVector(std::ostream& out, const StdVector<T>& a)
{
    writeEntry(out, (uint64_t)a.size());
    writeEntry(out, (uint64_t)sizeof(T));
    for (const auto& x : a)
        writeEntry(out, x);
}

template <class T>
void readSTDVector(std::istream& in, StdVector<T>& a)
{
    uint64_t size = readEntry<uint64_t>(in);
    uint64_t bytes = readEntry<uint64_t>(in);
    ZIRAN_ASSERT(bytes == sizeof(T), "Read error: The size of the types don't match");
    a.clear();
    a.reserve(size);
    for (size_t i = 0; i < size; i++)
        a.emplace_back(readEntry<T>(in));
}

template <class T>
inline void
writeHelper(std::ostream& out, const T a, TriviallyCopyableTag<T>)
{
    out.write(reinterpret_cast<const char*>(&a), sizeof a);
}

template <class T>
inline T
readHelper(std::istream& in, TriviallyCopyableTag<T>)
{
    T a;
    in.read(reinterpret_cast<char*>(&a), sizeof a);
    return a;
}

template <class T, int m, int n, int flags, int mm, int nn>
void writeHelper(std::ostream& out, const Eigen::Matrix<T, m, n, flags, mm, nn>& a, EigenMatrixTag<T, m, n>)
{
    if (m == Eigen::Dynamic)
        writeEntry(out, (uint64_t)a.rows());
    if (n == Eigen::Dynamic)
        writeEntry(out, (uint64_t)a.cols());

    out.write(reinterpret_cast<const char*>(a.data()), sizeof(T) * a.rows() * a.cols());
}

template <class T, int m, int n>
inline Matrix<T, m, n> readHelper(std::istream& in, EigenMatrixTag<T, m, n>)
{
    auto rows = m;
    auto cols = n;
    if (rows == Eigen::Dynamic)
        rows = readEntry<uint64_t>(in);
    if (cols == Eigen::Dynamic)
        cols = readEntry<uint64_t>(in);
    Matrix<T, m, n> a;
    a.resize(rows, cols);
    in.read(reinterpret_cast<char*>(a.data()), sizeof(T) * a.rows() * a.cols());
    return a;
}

template <class CustomType>
void writeHelper(std::ostream& out, const std::unique_ptr<CustomType>& model, UniquePtrTag<CustomType>)
{
    if (model == nullptr)
        writeEntry(out, '\0');
    else {
        writeEntry(out, 'u');
        model->write(out);
    }
}

template <class CustomType>
inline std::unique_ptr<CustomType>
readHelper(std::istream& in, UniquePtrTag<CustomType>)
{
    char c = readEntry<char>(in);
    if (c == '\0')
        return nullptr;
    else {
        return CustomType::readUnique(in);
    }
}

template <class T>
void writeHelper(std::ostream& out, const std::basic_string<T>& a, BasicStringTag<T>)
{
    writeEntry(out, (uint64_t)a.size());
    out.write(a.data(), a.size() * sizeof(T));
}

template <class T>
std::basic_string<T> readHelper(std::istream& in, BasicStringTag<T>)
{
    uint64_t size = readEntry<uint64_t>(in);
    std::basic_string<T> a;
    a.resize(size);
    in.read(&a[0], size * sizeof(T));
    return a;
}

template <class T>
void writeHelper(std::ostream& out, const Eigen::SparseMatrix<T>& a, EigenSparseMatrixTag<T>)
{
    writeEntry(out, a.rows());
    writeEntry(out, a.cols());
    writeEntry(out, (uint64_t)a.nonZeros());
    writeEntry(out, sizeof(Eigen::Triplet<T, typename Eigen::SparseMatrix<T>::Index>));
    for (int k = 0; k < a.outerSize(); ++k)
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(a, k); it; ++it) {
            writeEntry(out, it.row());
            writeEntry(out, it.col());
            writeEntry(out, it.value());
        }
}

template <class T>
inline Eigen::SparseMatrix<T>
readHelper(std::istream& in, EigenSparseMatrixTag<T>)
{
    using TMat = Eigen::SparseMatrix<T>;
    typename TMat::Index rows = readEntry<typename TMat::Index>(in);
    typename TMat::Index cols = readEntry<typename TMat::Index>(in);
    StdVector<Eigen::Triplet<typename TMat::Scalar, typename TMat::Index>> a;
    readSTDVector(in, a);
    TMat mat(rows, cols);
    mat.setFromTriplets(a.begin(), a.end());
    mat.markAsRValue();
    return mat;
}

/**
  Implementation to use for custom types which should make a member function
  called write and a static function called read
  and implement the struct RW
  with member type Tag = CustomTypeTag<MyClass> like so

  struct RW<MyClass> {
    using Tag = CustomTypeTag<MyClass>;
  };
  */
template <class CustomType>
void writeHelper(std::ostream& out, const CustomType& model, CustomTypeTag<CustomType>)
{
    model.write(out);
}

template <class CustomType>
inline CustomType readHelper(std::istream& in, CustomTypeTag<CustomType>)
{
    return CustomType::read(in);
}

template <class T, size_t size>
void writeHelper(std::ostream& out, const std::array<T, size>& a, StdArrayTag<T, size>)
{
    for (const auto& s : a)
        writeEntry(out, s);
}

template <class T, size_t size>
inline std::array<T, size> readHelper(std::istream& in, StdArrayTag<T, size>)
{
    std::array<T, size> a;
    for (auto& s : a)
        s = readEntry<T>(in);
    return a;
}

/**
  Implementation to use for types which should not be written, the read in
  instance will always be default constructed.
  To use this implement the struct RW
  with member type Tag = NoWriteTag<MyClass> like so

  struct RW<MyClass> {
    using Tag = NoWriteTag<MyClass>;
  };
  */
template <class NoWrite>
void writeHelper(std::ostream& out, const NoWrite& model, NoWriteTag<NoWrite>)
{
}

template <class NoWrite>
inline NoWrite readHelper(std::istream& in, NoWriteTag<NoWrite>)
{
    return NoWrite();
}

template <class Type>
void writeEntry(std::ostream& out, const Type& x)
{
    writeHelper(out, x, typename RW<Type>::Tag{});
}

template <class Type>
inline Type readEntry(std::istream& in)
{
    return readHelper(in, typename RW<Type>::Tag{});
}
} // namespace ZIRAN

#endif
