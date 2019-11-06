#ifndef DATA_ARRAY_BASE_H
#define DATA_ARRAY_BASE_H

#include <Ziran/CS/DataStructure/DisjointRanges.h>
#include <memory>
#include <string>
#include <vector>

namespace ZIRAN {

namespace DATA {

/**
    The base class of DataArray. Stores the range vector.
*/
class DataArrayBase {
public:
    std::string name;
    DisjointRanges ranges;
    using ValueIDToNewValueID = std::function<int(int)>;
    using ReorderingCallback = std::function<void(const ValueIDToNewValueID&)>;
    StdVector<ReorderingCallback> reordering_callbacks;

    DataArrayBase(const std::string& name, const DisjointRanges& ranges);
    virtual ~DataArrayBase();

    virtual void clear() = 0;

    virtual std::unique_ptr<DataArrayBase> makeEmptyCopy() = 0;
    virtual void reorder(const StdVector<int>& new_to_old) = 0;

    /**
      Writes the data to out in binary form
      */
    virtual void writeData(std::ostream& out) const = 0;
    /**
      Reads the data from in (that was written using writeData)
      Should replace any existing data
      */
    virtual void readData(std::istream& in) = 0;

    virtual size_t size() const = 0;

    static int valueId(int entry_id, const DisjointRanges& ranges);

    /** 
      Register a function to be called in the event of reordering 
      */
    void registerReorderingCallback(ReorderingCallback&& callback);

    /**
      Calculates the value id of a given entry
      Warning slow, should use iterators if possible
      */
    int valueId(int entry_id) const;

    /**
      Calculates the entry id of a given value
      Warning slow, should use iterators if possible
      */
    int entryId(int valueId) const;
};
}
} // namespace ZIRAN::DATA
#endif
