#ifndef PRIORITY_QUEUE_H_
#define PRIORITY_QUEUE_H_
#include <Ziran/CS/DataStructure/HashTable.h>
#include <vector>

namespace ZIRAN {
/**
  A (min) priority queue using a combination of a heap and a hashtable
  Values with the smallest priority are stored at the top.
  Supports efficiently updating priorities of existing elements.
  Should not be used to store large values as two copies of 
  the values will be made, to work around this restriction
  store indices or pointers and not the values themselves.

  Also values must be unique and hashable.
  Priorities need not be unique must they must be comparable.
  */
template <class Priority, class Value>
class PriorityQueue {
public:
    using Entry = KeyValuePair<Priority, Value>;

    // private:
    StdVector<Entry> data;
    HashTable<Value, size_t> data_index;

    inline size_t parentIndex(size_t index) const
    {
        return (index - 1) >> 5;
        // return (index - 1) >> 2;
    }

    inline size_t startChildIndex(size_t index) const
    {
        return 32 * index + 1;
        // return 4 * index + 1;
    }

    inline size_t endChildIndex(size_t child_start) const
    {
        using std::min;
        return min(child_start + 32, data.size());
        // return min(child_start + 4, data.size());
    }

    void percolateUp(size_t index)
    {
        using std::swap;
        while (index != 0) {
            auto& current_data = data[index];
            size_t parent_index = parentIndex(index);
            auto& parent_data = data[parent_index];
            if (current_data.key >= parent_data.key)
                break;
            swap(data_index[current_data.value], data_index[parent_data.value]);
            swap(current_data, parent_data);
            index = parent_index;
            ZIRAN_ASSERT(current_data.key >= parent_data.key);
        }
    }

    void percolateDown(size_t index)
    {
        using std::swap;
        size_t child_start = startChildIndex(index);
        while (child_start < data.size()) {
            size_t child_end = endChildIndex(child_start);
            size_t min_child_index = child_start;
            for (size_t i = child_start + 1; i < child_end; i++) {
                if (data[i].key < data[min_child_index].key) {
                    min_child_index = i;
                }
            }
            auto& current_data = data[index];
            auto& min_data = data[min_child_index];
            if (current_data.key <= min_data.key) {
                break;
            }
            swap(data_index[current_data.value], data_index[min_data.value]);
            swap(current_data, min_data);
            index = min_child_index;
            child_start = startChildIndex(index);
        }
    }

public:
    /*
       Insert value with priority, or updates priority if value is already in the queue
       O(log(n))
       */
    bool set(const Value& value, const Priority& priority)
    {
        using std::swap;
        auto result = data_index.insert(value, data.size());
        size_t& current_index = result.value;
        if (result.inserted) {
            ZIRAN_ASSERT(current_index == data.size());
            data.emplace_back(Entry{ priority, value });
            percolateUp(current_index);
            return true;
        }
        else {
            ZIRAN_ASSERT(current_index < data.size());
            auto& kv = data[current_index];
            ZIRAN_ASSERT(kv.value == result.key);
            Priority old_priority = kv.key;
            kv.key = priority;
            if (priority < old_priority)
                percolateUp(current_index);
            else if (priority < old_priority) {
                percolateDown(current_index);
            }
            return false;
        }
    }

    /*
       O(1)
       */
    Priority lookupPriority(const Value& value) const
    {
        size_t* index = data_index.get(value);
        ZIRAN_ASSERT(index != nullptr);
        return data[*index].key;
    }

    /*
       O(1)
       */
    bool empty() const
    {
        return data.size() == 0;
    }

    /*
       O(1)
       */
    size_t size() const
    {
        return data.size();
    }

    /*
       O(log(n))
       */
    bool removeTop(Entry& top)
    {
        using std::swap;
        if (data.size() == 0)
            return false;
        swap(top, data.front());
        swap(data.front(), data.back());
        data.pop_back();
        data_index[data.front().value] = 0;
        percolateDown(0);
        data_index.erase(top.value);
        return true;
    }

    /*
       O(1)
       */
    const Entry* peekTop() const
    {
        if (data.size() == 0)
            return nullptr;
        return data.front();
    }

    void print(std::ostream& out)
    {
        size_t line_end = 1;
        for (size_t i = 0; i < data.size(); i++) {
            if (i == line_end) {
                out << '\n';
                line_end = endChildIndex(startChildIndex(i - 1));
            }
            out << '(' << data[i].key << ", " << data[i].value << ") ";
        }
        out << '\n';
    }
};
} // namespace ZIRAN
#endif
