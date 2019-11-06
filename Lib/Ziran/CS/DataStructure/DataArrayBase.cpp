#include <Ziran/CS/DataStructure/DataArrayBase.h>

namespace ZIRAN {

namespace DATA {
DataArrayBase::DataArrayBase(const std::string& name, const DisjointRanges& ranges)
    : name(name)
    , ranges(ranges)
{
}
DataArrayBase::~DataArrayBase() {}

/** 
  Register a function to be called in the event of reordering 
  */
void DataArrayBase::registerReorderingCallback(ReorderingCallback&& callback)
{
    reordering_callbacks.emplace_back(std::move(callback));
}

/**
  Calculates the value id of a given entry
  Warning slow, should use iterators if possible
  */

int DataArrayBase::valueId(int entry_id, const DisjointRanges& ranges)
{
    int valueId = 0;
    for (const auto& r : ranges) {
        if (entry_id < r.upper) {
            if (entry_id < r.lower)
                return -1;

            valueId += entry_id - r.lower;
            return valueId;
        }
        else {
            valueId += r.upper - r.lower;
        }
    }
    return -1;
}

int DataArrayBase::valueId(int entry_id) const
{
    return valueId(entry_id, ranges);
}

/**
  Calculates the entry id of a given value
  Warning slow, should use iterators if possible
  */
int DataArrayBase::entryId(int valueId) const
{
    for (const auto& r : ranges) {
        int length = r.length();
        if (valueId < length) {
            return r.lower + valueId;
        }
        else {
            valueId -= length;
        }
    }
    return -1;
}
}
} // namespace ZIRAN::DATA
