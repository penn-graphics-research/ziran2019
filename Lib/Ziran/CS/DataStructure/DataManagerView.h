#ifndef DATA_MANAGER_VIEW
#define DATA_MANAGER_VIEW
#include <Ziran/CS/DataStructure/DisjointRanges.h>
#include <Ziran/CS/DataStructure/DataManager.h>
namespace ZIRAN {

struct DataManagerView {
    const DataManager& manager;
    DisjointRanges ranges;
    DataManagerView(const DataManager& manager)
        : manager(manager)
        , ranges{ { 0, manager.count } }
    {
    }

    DataManagerView(const DataManager& manager, const DisjointRanges& ranges)
        : manager(manager)
        , ranges(ranges)
    {
    }
};
} // namespace ZIRAN
#endif
