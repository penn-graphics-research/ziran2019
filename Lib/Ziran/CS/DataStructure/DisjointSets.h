#ifndef DISJOINT_SETS_H
#define DISJOINT_SETS_H
#include <Ziran/CS/DataStructure/HashTable.h>
#include <vector>
namespace ZIRAN {
/**
  Represent disjoint sets of values which collectively form a sequence of
  non negative integers from 0 to n. If you need to represent sets
  of other values, it can be used along with an array or hashtable.

  Common uses: Finding disjoint regions in a graph,
  Kruskal's algorithm (for minimum spanning trees)

  The sets are represented implicitly by the
  function getSetId which for any value
  returns the id of the containing set
  (which is the canonical id of the set)

  Assymptotically optimal with amortized runtime O(alpha(n))
  where alpha is the inverse of the Ackermann function A(n,n)
  i.e. less than 5 for anything that will fit in memory
  */
template <class ID = int>
class DisjointSets {
    StdVector<ID> parent;
    StdVector<uint8_t> rank;
    StdVector<ID> path;
    // ID path[32 * sizeof(ID)];

public:
    DisjointSets(ID starting_number_of_sets = 0)
    {
        reinitialize(starting_number_of_sets);
    }

    void reinitialize(ID starting_number_of_sets)
    {
        parent.clear();
        rank.clear();
        path.clear();

        parent.reserve(starting_number_of_sets);
        rank.reserve(starting_number_of_sets);
        for (int i = 0; i < starting_number_of_sets; i++)
            createNextSet();
    }

    /**
      Creates a new one element set and returns
      the id of the newly created set
      */
    ID createNextSet()
    {
        ID next_set_id = (ID)parent.size();
        parent.emplace_back(next_set_id);
        rank.emplace_back(0);
        return next_set_id;
    }

    /**
      Gets the canonical id of the the set that j belongs to

      (frequently called the find algorithm)
      */
    ID getSetId(ID j)
    {
        size_t i = 0;
        while (j != parent[j]) {
            // assert(i < 32 * sizeof(ID));

            path.emplace_back(j);
            i++;
            // path[i++] = j;

            j = parent[j];
        }
        // Since the set j of all the values seen while searching
        // for the canonical index is the same we can compress the path
        // we searched to find it
        while (i > 0) {
            parent[path[--i]] = j;
        }
        path.clear();
        return j;
    }

    /**
      Merges the two sets together if they are not already
      (frequently called the union algorithm)

      Must be passed in the canonical id of the set (i.e. the 
      one with getSetId(set_id) == set_id)

      returns the set id of the merged set
      which will always be one of the passed in set ids
      */
    ID merge(ID set_id0, ID set_id1)
    {
        assert(parent[set_id0] == set_id0 && parent[set_id1] == set_id1);

        if (set_id0 == set_id1)
            return set_id0;

        // Always merge into the set with the lager rank
        // this keeps the trees well balanced
        if (rank[set_id0] < rank[set_id1]) {
            parent[set_id1] = set_id0;
            return set_id1;
        }
        else if (rank[set_id1] < rank[set_id0]) {
            parent[set_id0] = set_id1;
            return set_id1;
        }
        else {
            parent[set_id1] = set_id0;
            rank[set_id0]++;
            return set_id0;
        }
    }
};
} // namespace ZIRAN
#endif
