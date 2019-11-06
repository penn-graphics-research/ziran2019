#ifndef DATA_ARRAY_H
#define DATA_ARRAY_H
#include <Ziran/CS/DataStructure/DataArrayBase.h>
#include <Ziran/CS/Util/BinaryIO.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/Meta.h>
#include <string>
#include <vector>

namespace ZIRAN {

namespace DATA {

/**
    DataArray represents a continuous memory block of entry data.
    It has two members.
    1. ranges.
        This is an array of sorted ranges of global entry IDs,
        inclusive on the left, exclusive on the right.
        For example, [0,3), [6,8) means this batch stores data for entries
        0,1,2,6,7
    2. array.
        This is an array of data with any type. E.g., for positions.
        The array stores contiguous data for the entities in the ranges
        Continuing the example if the array is [1, 2, 3, 4, 5]
        then the entries and correspoding values are as follows
        |Entry |   ValueId |   Value |
        |  0   |     0     |     1   |
        |  1   |     1     |     2   |
        |  2   |     2     |     3   |
        |  3   |    NA     |    NA   |
        |  4   |    NA     |    NA   |
        |  5   |    NA     |    NA   |
        |  6   |     3     |     4   |
        |  7   |     4     |     5   |

*/

template <class Type, bool is_const>
class DataArrayIterator;

template <class Type>
class DataArray : public DataArrayBase {
public:
    using DataArrayBase::name;
    using DataArrayBase::ranges;
    /**
      The value array is public for conveniance, it's size should not be changed directly
      unless you can uphold the necessary invariants.
      1. The ranges and entries of the array must be kept consistent
      2. The total number of all possible entries should match the count in the data manager
      */
    StdVector<Type> array;

    DataArray(const std::string& name, const DisjointRanges& ranges, StdVector<Type>&& array)
        : DataArrayBase(name, ranges)
        , array(std::move(array))
    {
        assert(ranges.length() == (int)this->array.size());
    }
    virtual ~DataArray() override {}

    virtual void clear() override
    {
        ranges.ranges.clear();
        array.clear();
    }

    virtual std::unique_ptr<DataArrayBase> makeEmptyCopy() override
    {
        return std::make_unique<DataArray<Type>>(name, DisjointRanges(), StdVector<Type>());
    }

    virtual void reorder(const StdVector<int>& new_to_old) override
    {
        StdVector<Type> array_new;
        DisjointRanges ranges_new;
        bool needs_inverse_map = reordering_callbacks.size() > 0;
        StdVector<int> old_to_new_value_id(array.size());
        array_new.reserve(new_to_old.size());
        int new_value_id = 0;
        for (size_t i = 0; i < new_to_old.size(); ++i) {
            int old_index = new_to_old[i];
            int id_in_array = valueId(old_index);
            if (id_in_array != -1) {
                array_new.push_back(array[id_in_array]);
                ranges_new.append({ (int)i, (int)i + 1 });
                if (needs_inverse_map)
                    old_to_new_value_id[id_in_array] = new_value_id++;
            }
        }
        array_new.swap(array);
        ranges = ranges_new;
        for (const auto& callback : reordering_callbacks)
            callback([&](int old_value_id) {
                return old_to_new_value_id[old_value_id];
            });
    }

    void writeData(std::ostream& out) const override
    {
        writeEntry(out, ranges.lg2_grain_size);
        writeSTDVector(out, ranges.ranges);
        writeSTDVector(out, array);
    }

    void readData(std::istream& in) override
    {
        clear();

        ranges.lg2_grain_size = readEntry<int>(in);
        readSTDVector(in, ranges.ranges);
        readSTDVector(in, array);
        ZIRAN_ASSERT(ranges.length() == (int)array.size());
        ZIRAN_ASSERT(ranges.valid());
    }

    /**
      Resizes the underlying array to have the same size as the disjoint ranges
      Does not preserve the id to value relation
    */
    void lazyResize(const DisjointRanges& dr)
    {
        ranges = dr;
        array.resize(dr.length());
    }

    /**
      Resizes the underlying array to have the same size as the disjoint ranges
      Fill the vector with val
    */
    void lazyResize(const DisjointRanges& dr, const Type& val)
    {
        ranges = dr;
        array.resize(dr.length(), val);
    }

    /**
      Number of elements in the value array
    */
    size_t size() const override
    {
        return array.size();
    }

    /**
      Reserve space for n entries to be appended
    */
    void reserveMore(const size_t n)
    {
        array.reserve(array.size() + n);
    }

    /**
      Indexes into the value array.
      (not according to the entry id)
    */
    Type& operator[](size_t i)
    {
        assert(i < array.size());
        return array[i];
    }

    const Type* getValue(int entry_id) const
    {
        int value_id = valueId(entry_id);
        if (value_id < 0)
            return nullptr;
        else
            return &(array[value_id]);
    }

    Type* getValue(int entry_id)
    {
        int value_id = valueId(entry_id);
        if (value_id < 0)
            return nullptr;
        else
            return &(array[value_id]);
    }

    /**
      Indexes into the value array.
      (not according to the entry id)
    */
    const Type& operator[](size_t i) const
    {
        assert(i < array.size());
        return array[i];
    }

    /**
      Indexes into the value array.
      (not according to the entry id)
    */
    Type& operator()(size_t i)
    {
        assert(i < array.size());
        return array[i];
    }

    /**
      Indexes into the value array.
      (not according to the entry id)
    */
    const Type& operator()(size_t i) const
    {
        assert(i < array.size());
        return array[i];
    }

    bool operator==(const DataArray<Type>& other) const
    {
        return (array == other.array) && (ranges == other.ranges);
    }

    DataArrayIterator<Type, false> iter()
    {
        return DataArrayIterator<Type, false>(*this, ranges);
    }

    DataArrayIterator<Type, true> iter() const
    {
        return DataArrayIterator<Type, true>(*this, ranges);
    }

    DataArrayIterator<Type, false> subsetIter(const DisjointRanges& dr)
    {
        return DataArrayIterator<Type, false>(*this, dr);
    }

    DataArrayIterator<Type, true> subsetIter(const DisjointRanges& dr) const
    {
        return DataArrayIterator<Type, true>(*this, dr);
    }
};

/**
  Allows using the value array as an Eigen Vector
  This doesn't create a copy, it creates an Eigen::Map
  into the underlying array
*/
template <class Type>
Eigen::Map<Vector<ScalarType<Type>, Eigen::Dynamic>> toVectorMap(DataArray<Type>& data)
{
    return ToVectorMap<ScalarType<Type>>(data.array);
}
/**
  Allows using the value array as an Eigen Matrix (entries as cols)
  This doesn't create a copy, it creates an Eigen::Map
  into the underlying array
*/
template <class T, int d>
Eigen::Map<Eigen::Matrix<T, d, Eigen::Dynamic, Eigen::ColMajor>>
toMatrixColMap(DataArray<Vector<T, d>>& data)
{
    return ToMatrixColMap<T, d>(data.array);
}
/**
  Allows using the value array as an Eigen Matrix (entries as rows)
  This doesn't create a copy, it creates an Eigen::Map
  into the underlying array
*/
template <class T, int d>
Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, d, Eigen::RowMajor>>
toMatrixRowMap(DataArray<Vector<T, d>>& data)
{
    return ToMatrixRowMap<T, d>(data.array);
}

/**
  DataArrayIterator is an iterator over the data entries
  It optionally takes a subset which defines which
  rows of entries are to be iterated over
  The subset ranges must be a subset of batch.ranges
  */
template <class Type, bool is_const>
class DataArrayIterator {
    typedef std::conditional_t<is_const,
        typename StdVector<Type>::const_iterator,
        typename StdVector<Type>::iterator>
        IterT;

    //iterator into the ranges
    typename DisjointRanges::ConstIterator r;
    //iterator into the subset
    typename DisjointRanges::ConstIterator s;
    //iterator into the data array
    IterT a;
#ifndef NDEBUG
    IterT aend;
#endif
    // end of subset array
    const typename DisjointRanges::ConstIterator end;
    //current entry id
    int p;

public:
    typedef typename std::conditional<is_const,
        const DataArray<Type>&,
        DataArray<Type>&>::type BatchType;

    typedef typename std::conditional<is_const,
        const Type&,
        Type&>::type ReferenceType;

    DataArrayIterator(BatchType batch, const DisjointRanges& subset)
        : r(batch.ranges.cbegin())
        , s(subset.cbegin())
        , a(batch.array.begin())
#ifndef NDEBUG
        , aend(batch.array.end())
#endif
        , end(subset.cend())
        , p((r != batch.ranges.cend()) ? (r->lower) : -1)
    {
        if (s != end) {
            advanceToEntryId(s->lower);
        }
    }

    DataArrayIterator(BatchType batch)
        : DataArrayIterator(batch, batch.ranges)
    {
    }

    /// Checks if this entry is valid
    explicit operator bool() const
    {
        return (s != end);
    }

    bool operator==(const DataArrayIterator& o) const
    {
        return p == o.p;
    }

    bool operator!=(const DataArrayIterator& o) const
    {
        return p != o.p;
    }

    /// The current entry id
    int entryId() const
    {
        return p;
    }

    /// Advance to the next entry
    DataArrayIterator<Type, is_const>& operator++()
    {
        ++a;
        ++p;
        if (p == s->upper) {
            // we've exhausted the current subset
            ++s;
            if (s != end)
                advanceToEntryId(s->lower);
        }

        return *this;
    }

    /// Advance by n entries
    DataArrayIterator<Type, is_const>& operator+=(int n)
    {
        while (n > 0) {
            int bound = s->upper;
            int diff = std::min(bound - p, n);
            a += diff;
            p += diff;
            n -= diff;
            if (p == bound) {
                // we've exhausted the current subset
                ++s;
                if (s != end)
                    advanceToEntryId(s->lower);
            }
        }

        return *this;
    }

    /// Advance by n entries
    DataArrayIterator<Type, is_const> operator+(int n)
    {
        DataArrayIterator<Type, is_const> iter(*this);
        return iter += n;
    }

    /// Advance the iterator to the given entry id
    void advanceToEntryId(int eid)
    {
        while (r->upper <= eid) {
            // skip all that remains of current r range
            a += r->upper - p;
            ++r;
            p = r->lower;
        }
        if (p < eid) {
            // skip until eid
            a += eid - p;
            p = eid;
        }
    }

    /// Advance to the next entry (postfix)
    DataArrayIterator<Type, is_const> operator++(int)
    {
        DataArrayIterator<Type, is_const> tmp(*this);
        operator++();
        return tmp;
    }

    ReferenceType operator*() const
    {
        assert(a < aend);
        return *a;
    }
};
} // namespace DATA
using DataArrayBase = DATA::DataArrayBase;
template <class Type>
using DataArray = DATA::DataArray<Type>;
} // namespace ZIRAN
#endif
