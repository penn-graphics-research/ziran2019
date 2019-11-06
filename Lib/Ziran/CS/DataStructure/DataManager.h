#ifndef DATA_MANAGER_H
#define DATA_MANAGER_H
#include <Ziran/CS/DataStructure/DataArray.h>
#include <Ziran/CS/Util/BinaryIO.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/Meta.h>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <tbb/tbb.h>
#include <tick/requires.h>

namespace ZIRAN {

namespace DATA {

/**
  An AttributeName stores the name of an attribute, together with its type
  It is used as a convenience to make it so that you only have to state
  the name and the type of the attribute once.
  It also stores the hash value to avoid dynamic memory allocation
  */
template <class Type>
struct AttributeName {
public:
    const std::string name;
    const size_t hash;

    AttributeName(const char* a)
        : name(a)
        , hash(std::hash<std::string>()(name))
    {
    }

    AttributeName(const std::string& a)
        : name(a)
        , hash(std::hash<std::string>()(name))
    {
    }
};

/**
  ValueIdOnly is a marker type used to indicate that the iterator should
  return the corresponding value id of the entry and not a reference to
  the entry itself.
  */
struct ValueIdOnly {
};

template <bool is_const, class... Types>
class DataManagerIterator;

template <class... Types>
class DataAppender;

/**
    DataManager is a database.
    It stores arbitrary number of arrays with arbitrary number of attributes.
    Totally flexible.
    Note that when adding entries, new entries ranges (IDs)
    must be retrived with getNextRange from the database. This is to ensure that
    the entries are stored in the database in a sorted fashion.
    Once created, this class is the only interface of using or modifying entry data.
*/
class DataManager {
public:
    int count; // the total number of entries.
    std::unordered_map<size_t, std::unique_ptr<DataArrayBase>> batches;

    DataManager()
        : count(0)
    {
    }

    // no copying the DataManager
    DataManager(const DataManager& other) = delete;

    /**
      Construct empty arrays with names and types given by the input `AttributeName<T>`
      Note that since `AttributeName<T>` can be freely constructed from a char *
      You can use the form
      `DataManager<int> data("a");`
      or
      `DataManager data(AttributeName<int>("a"));`
      the later is useful when the `AttributeName<int>("a")` is stored in a variable for frequent usage
      */
    template <typename Type, typename... Types>
    DataManager(const AttributeName<Type>& a, const AttributeName<Types>&... b)
        : DataManager(b...)
    {
        add(a);
    }

    virtual ~DataManager() {}

    /**
      Add an empty attribute a
      */
    template <class Type>
    DataArray<Type>& add(const AttributeName<Type>& a)
    {
        return add(a, Range{ count, count }, StdVector<Type>());
    }

    /**
      Add copies of the same value in the specified range
      */
    template <class Type>
    std::enable_if_t<!IsStdVector<Type>::value, DataArray<Type>&>
    add(const char* name, const Range& range, const Type& value)
    {
        return add<Type>(AttributeName<Type>(name), range, value);
    }

    template <class Type, class Type2>
    std::enable_if_t<!IsStdVector<Type2>::value, DataArray<Type>&>
    add(const AttributeName<Type>& a, const Range& range, const Type2& value)
    {
        return add(a, range, StdVector<Type>(range.upper - range.lower, value));
    }

    /**
      Add values for the entries in the corresponding range
      */
    template <class Type>
    DataArray<Type>& add(const char* name, Range range, StdVector<Type>&& array)
    {
        return add<Type>(AttributeName<Type>(name), range, std::move(array));
    }

    template <class Type, class Allocator = std::allocator<Type>>
    DataArray<Type>& add(const AttributeName<Type>& name, Range range, const std::vector<Type, Allocator>& array)
    {
        return add<Type>(AttributeName<Type>(name), range, StdVector<Type>(array.begin(), array.end()));
    }

    template <class Type>
    DataArray<Type>& add(const AttributeName<Type>& a, Range range, StdVector<Type>&& array)
    {
        ZIRAN_ASSERT((range.upper - range.lower) == (int)array.size());
        if (range.upper > count) {
            count = range.upper;
        }

        auto search = batches.find(a.hash);

        if (search == batches.end()) {
            // The attribute doesn't exist create it
            auto r = batches.emplace(a.hash,
                std::make_unique<DataArray<Type>>(
                    a.name, DisjointRanges{ range }, std::move(array)));
            return dynamic_cast<DataArray<Type>&>(*r.first->second);
        }
        else {
            ZIRAN_ASSERT(a.name == search->second->name);
            // The attribute exists, append the entries
            DataArray<Type>& old = dynamic_cast<DataArray<Type>&>(*search->second);
            old.ranges.append(range);
            old.array.insert(old.array.end(),
                std::make_move_iterator(array.begin()),
                std::make_move_iterator(array.end()));
            array.clear();
            return old;
        }
    }

    template <class Type>
    DataArray<Type>& get(const AttributeName<Type>& a) const
    {
        auto search = batches.find(a.hash);
        ZIRAN_ASSERT(search != batches.end(), "Could not find ", a.name);
        ZIRAN_ASSERT(a.name == search->second->name);
        return dynamic_cast<DataArray<Type>&>(*search->second);
    }

    DataArrayBase& get(const AttributeName<ValueIdOnly>& a) const
    {
        auto search = batches.find(a.hash);
        ZIRAN_ASSERT(search != batches.end(), "Could not find ", a.name);
        ZIRAN_ASSERT(a.name == search->second->name);
        return *search->second;
    }

    template <class Type>
    bool exist(const AttributeName<Type>& a) const
    {
        auto search = batches.find(a.hash);
        if (search != batches.end()) {
            ZIRAN_ASSERT(a.name == search->second->name);
        }
        return (search != batches.end());
    }

    bool exist(const char* name) const
    {
        size_t hash(std::hash<std::string>()(name));
        auto search = batches.find(hash);
        if (search != batches.end()) {
            ZIRAN_ASSERT(name == search->second->name);
        }
        return (search != batches.end());
    }

    Range getNextRange(int size) const
    {
        return Range{ count, count + size };
    }

    template <typename... Types>
    DisjointRanges commonRanges(const AttributeName<Types>&... names) const
    {
        return DisjointRanges(get(names).ranges...);
    }

    /// Returns an iterator over the entries which have all the named attributes
    template <typename... Types>
    DataManagerIterator<false, Types...> iter(const AttributeName<Types>&... names)
    {
        return DataManagerIterator<false, Types...>(get(names)...);
    }

    template <typename... Types>
    DataManagerIterator<true, Types...> iter(const AttributeName<Types>&... names) const
    {
        return DataManagerIterator<true, Types...>(get(names)...);
    }

    /**
      Returns an iterator over a subset of the entities,
      the named attributes must have data for all the entries in the subset
      */
    template <typename... Types>
    DataManagerIterator<false, Types...> subsetIter(const DisjointRanges& subset, const AttributeName<Types>&... names)
    {
        return DataManagerIterator<false, Types...>(subset, get(names)...);
    }

    /**
      Returns an iterator over a subset of the entities,
      the named attributes must have data for all the entries in the subset
      */
    template <typename... Types>
    DataManagerIterator<true, Types...> subsetIter(const DisjointRanges& subset, const AttributeName<Types>&... names) const
    {
        return DataManagerIterator<true, Types...>(subset, get(names)...);
    }

    /**
      Returns an appender for the named attributes,
      a helper class for appending data to the entries
      */
    template <typename... Types>
    DataAppender<Types...> appender(const AttributeName<Types>&... names)
    {
        return DataAppender<Types...>(count, add(names)...);
    }

    /*
       Writes the stored data in binary form to out
    */
    virtual void writeData(std::ostream& out) const
    {
        writeEntry(out, count);
        writeEntry(out, (uint64_t)batches.size());
        for (const auto& kv : batches) {
            writeEntry(out, kv.second->name);
            ZIRAN_ASSERT(kv.first == std::hash<std::string>()(kv.second->name),
                "Array ", kv.second->name, " doesn't match hash ", kv.first);
            kv.second->writeData(out);
        }
    }

    /*
       Reads in the data written by writeData
       Empty attribute arrays with the correct types must be created before
       this function is called
    */
    virtual void readData(std::istream& in)
    {
        for (auto& kv : batches)
            kv.second->clear();

        count = readEntry<int>(in);
        uint64_t size = readEntry<uint64_t>(in);
        for (size_t i = 0; i < size; i++) {
            auto name = readEntry<std::string>(in);
            auto search = batches.find(std::hash<std::string>()(name));
            ZIRAN_ASSERT(search != batches.end(), "Array ", name, " was not initialized before reading.");
            search->second->readData(in);
        }
    }

    template <class Func, typename... Types>
    void parallel_for(const Func& f, const AttributeName<Types>&... names)
    {
        auto master = iter(names...);
        tbb::blocked_range<int> range(0, master.common_ranges->length(), 1024);
        return tbb::parallel_for(range,
            [&](const tbb::blocked_range<int>& subrange) {
                auto end = master + subrange.end();
                for (auto i = master + subrange.begin(); i != end; ++i) {
                    i.apply(f);
                }
            });
    }

    template <class Func, typename... Types>
    void parallel_for_subset(const DisjointRanges& subset, const Func& f, const AttributeName<Types>&... names)
    {
        auto master = subsetIter(subset, names...);
        tbb::blocked_range<int> range(0, master.common_ranges->length(), 1024);
        return tbb::parallel_for(range,
            [&](const tbb::blocked_range<int>& subrange) {
                auto end = master + subrange.end();
                for (auto i = master + subrange.begin(); i != end; ++i) {
                    i.apply(f);
                }
            });
    }

    template <class Value, class Func1, class Func2, typename... Types>
    Value map_reduce(const Value& identity, const Func1& map, const Func2& reduction, const AttributeName<Types>&... names) const
    {
        auto master = iter(names...);
        tbb::blocked_range<int> range(0, master.common_ranges->length(), 1024);
        return tbb::parallel_reduce(range, identity,
            [&](const tbb::blocked_range<int>& subrange, const Value& init) -> Value {
                auto end = master + subrange.end();
                Value x = init;
                for (auto i = master + subrange.begin(); i != end; ++i) {
                    x = reduction(x, i.apply(map));
                }
                return x;
            },
            reduction);
    }

    void reorder(const StdVector<int>& new_to_old)
    {
        for (const auto& b : batches) {
            b.second->reorder(new_to_old);
        }
    }
};

/**
  DataManagerIteratorHelper is a helper class for
  the DataManagerIterator it handles the advancement
  for a single DataArray. It is not for external use
  */
template <class Type, bool is_const>
class DataManagerIteratorHelper {
    typedef std::conditional_t<is_const,
        typename StdVector<Type>::const_iterator,
        typename StdVector<Type>::iterator>
        IterT;

    //iterator into the ranges
    typename DisjointRanges::ConstIterator r;
    //iterator into the data array
    IterT a;
#ifndef NDEBUG
    IterT aend;
#endif

public:
    typedef typename std::conditional<is_const,
        const DataArray<Type>&,
        DataArray<Type>&>::type BatchType;

    typedef typename std::conditional<is_const,
        const Type&,
        Type&>::type ReferenceType;

    DataManagerIteratorHelper(BatchType batch, int start)
        : r(batch.ranges.cbegin())
        , a(batch.array.begin())
#ifndef NDEBUG
        , aend(batch.array.end())
#endif
    {
        if (start != -1)
            advance(r->lower, start);
    }

    /// Advance the iterator from the starting entry id to the ending entry id
    void advance(int start, int end)
    {
        int p = start;
        while (r->upper <= end) {
            // skip all that remains of current range
            a += r->upper - p;
            ++r;
            p = r->lower;
        }
        if (p < end) {
            // skip until eid
            a += end - p;
        }
    }

    ReferenceType operator*() const
    {
        assert(a < aend);
        return *a;
    }
};

template <class Type>
AttributeName<ValueIdOnly> valueIdOnly(const AttributeName<Type>& attribute)
{
    return AttributeName<ValueIdOnly>(attribute.name);
}

template <bool is_const>
class DataManagerIteratorHelper<ValueIdOnly, is_const> {
    typedef int IterT;

    //iterator into the ranges
    typename DisjointRanges::ConstIterator r;
    //iterator into the data array
    IterT a;
#ifndef NDEBUG
    IterT aend;
#endif

public:
    typedef const DataArrayBase& BatchType;

    typedef const int& ReferenceType;

    DataManagerIteratorHelper(BatchType batch, int start)
        : r(batch.ranges.cbegin())
        , a(0)
#ifndef NDEBUG
        , aend(batch.ranges.length())
#endif
    {
        if (start != -1)
            advance(r->lower, start);
    }

    /// Advance the iterator from the starting entry id to the ending entry id
    void advance(int start, int end)
    {
        int p = start;
        while (r->upper <= end) {
            // skip all that remains of current range
            a += r->upper - p;
            ++r;
            p = r->lower;
        }
        if (p < end) {
            // skip until eid
            a += end - p;
        }
    }

    ReferenceType operator*() const
    {
        assert(a < aend);
        return a;
    }
};

template <bool is_const, class... Types>
class DataManagerIterator {
public:
    const DisjointRanges* common_ranges;

    bool owns_ranges;
    //iterator into the common_ranges
    typename DisjointRanges::ConstIterator s;
    //current entry id
    int p;

    std::tuple<DataManagerIteratorHelper<Types, is_const>...> iterators;

    DataManagerIterator(typename DataManagerIteratorHelper<Types, is_const>::BatchType... batches)
        : common_ranges(new DisjointRanges(batches.ranges...))
        , owns_ranges(true)
        , s(common_ranges->cbegin())
        , p(s != common_ranges->cend() ? s->lower : -1)
        , iterators(DataManagerIteratorHelper<Types, is_const>(batches, p)...)
    {
    }

    DataManagerIterator(const DisjointRanges& common_ranges_input, typename DataManagerIteratorHelper<Types, is_const>::BatchType... batches)
        : common_ranges(new DisjointRanges(common_ranges_input))
        , owns_ranges(true)
        , s(common_ranges->cbegin())
        , p(s != common_ranges->cend() ? s->lower : -1)
        , iterators(DataManagerIteratorHelper<Types, is_const>(batches, p)...)
    {
    }

    DataManagerIterator(const DataManagerIterator& other)
        : common_ranges(other.common_ranges)
        , owns_ranges(false)
        , s(other.s)
        , p(other.p)
        , iterators(other.iterators)
    {
    }

    DataManagerIterator(DataManagerIterator&& other) = default;

    ~DataManagerIterator()
    {
        if (owns_ranges)
            delete common_ranges;
    }

    explicit operator bool() const
    {
        return s != common_ranges->cend();
    }

    bool operator==(const DataManagerIterator& o) const
    {
        return p == o.p;
    }

    bool operator!=(const DataManagerIterator& o) const
    {
        return p != o.p;
    }

    int entryId() const
    {
        return p;
    }

    DataManagerIterator<is_const, Types...>& operator++()
    {
        return (*this += 1);
    }

    DataManagerIterator<is_const, Types...>& operator+=(int n)
    {
        assert(n >= 0); //TODO n < 0
        while (n > 0) {
            int bound = s->upper;
            int desired = p + n;
            if (desired < bound) {
                advanceToEntryId(desired);
                return *this;
            }
            else {
                ++s;
                n -= bound - p;
                if (s != common_ranges->cend())
                    advanceToEntryId(s->lower);
                else
                    p = bound;
            }
        }
        assert((s == common_ranges->cend()) || ((s->lower <= p) && (p < s->upper)));
        return *this;
    }

    DataManagerIterator<is_const, Types...> operator+(int n)
    {
        DataManagerIterator<is_const, Types...> iter(*this);
        return iter += n;
    }

    std::tuple<typename DataArrayIterator<Types, is_const>::ReferenceType...> operator*()
    {
        return dereferenceHelper(std::index_sequence_for<Types...>());
    }

    void advanceToEntryId(int eid)
    {
        advanceToEntryIdHelper(eid, std::index_sequence_for<Types...>());
    }

    template <class Func>
    auto apply(const Func& f)
    {
        return applyHelper(f, std::index_sequence_for<Types...>());
    }

    template <size_t I>
    auto get() -> decltype(*std::get<I>(iterators))
    {
        return *std::get<I>(iterators);
    }

    template <size_t I>
    auto get() const -> decltype(*std::get<I>(iterators))
    {
        return *std::get<I>(iterators);
    }

private:
    template <std::size_t... I>
    std::tuple<typename DataArrayIterator<Types, is_const>::ReferenceType...> dereferenceHelper(std::index_sequence<I...>)
    {
        return std::tie(*std::get<I>(iterators)...);
    }

    template <std::size_t... I>
    void advanceToEntryIdHelper(int eid, std::index_sequence<I...>)
    {
        Call{ (std::get<I>(iterators).advance(p, eid), 0)... };
        p = eid;
    }

    template <class Func, std::size_t... I>
    auto applyHelper(const Func& f, std::index_sequence<I...>)
    {
        return f(*std::get<I>(iterators)...);
    }
};

// Class indicating to skip appending this attribute
struct Skip {
};

/**
  Helper class for appending data to the manager's arrays.
  Not threadsafe.

  The intended use is as follows
  \snippet Lib/ParticlesTest.h DataAppender example

  */
template <class... Types>
class DataAppender {
public:
    std::tuple<DataArray<Types>&...> batches;
    int& entry_id;

    DataAppender(int& count, DataArray<Types>&... batches)
        : batches(batches...)
        , entry_id(count)
    {
    }

    int entryId() const
    {
        return entry_id;
    }

    template <class... Ts>
    void append(Ts&&... values)
    {
        appendHelper<0>(std::forward<Ts>(values)...);
        entry_id++;
    }

    void reserve(const size_t n)
    {
        reserveHelper(std::index_sequence_for<Types...>(), n);
    }

private:
    template <int index>
    void appendHelper()
    {
    }

    template <int index, class... Ts>
    void appendHelper(const Skip& value, Ts&&... values)
    {
        appendHelper<index + 1>(std::forward<Ts>(values)...);
    }

    template <int index, class T, TICK_REQUIRES(!std::is_same<std::decay<T>, Skip>::value), class... Ts>
    void appendHelper(T&& value, Ts&&... values)
    {
        auto& batch = std::get<index>(batches);
        batch.array.emplace_back(std::forward<T>(value));
        batch.ranges.append({ entry_id, entry_id + 1 });
        appendHelper<index + 1>(std::forward<Ts>(values)...);
    }

    template <std::size_t... I>
    void reserveHelper(std::index_sequence<I...>, const size_t n)
    {
        Call{ (std::get<I>(batches).reserve(n), 0)... };
    }
};
} // namespace DATA
using DATA::DataManager;
template <class Type>
using AttributeName = DATA::AttributeName<Type>;
using DATA::Skip;
using DATA::valueIdOnly;
} // namespace ZIRAN
#endif
