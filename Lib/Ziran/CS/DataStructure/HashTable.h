#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <Ziran/CS/DataStructure/Hash.h>
#include <Ziran/CS/Util/Debug.h>
#include <Ziran/CS/Util/Forward.h>
#include <algorithm>
#include <assert.h>
#include <functional>
#include <iostream>
#include <memory>
#include <stddef.h>
#include <utility>

namespace ZIRAN {

template <class Key, class Value>
struct KeyValuePair {
    Key key;
    Value value;
};

template <class T>
inline static T makePowerOfTwo(T x)
{
    typename std::make_unsigned<T>::type k = x;
    k--;
    for (size_t i = 1; i < (sizeof(T) * CHAR_BIT); i <<= 1)
        k |= k >> i;
    k++;
    return k;
}

/**
  A hash table
  Uses Robin Hood hashing to deal with collisions
  The entries of the hashtable are stored in a contiguous array
  and will be swapped around.  Therefore large values or keys
  should be boxed i.e. put in unique_ptr or shared_ptr
  */
template <class _Key, class _Value, class Hash = Hash<_Key>,
    class KeyEqual = std::equal_to<_Key>,
    class Allocator = std::allocator<KeyValuePair<_Key, _Value>>>
class HashTable {
public:
    using Key = _Key;
    using Value = _Value;
    using HashType = GetHashType<Hash>;
    using Entry = KeyValuePair<Key, Value>;
    template <bool>
    class Iterator;

    struct InsertResult {
        const Key& key;
        Value& value;
        bool inserted;
        InsertResult(Entry& entry, bool inserted)
            : key(entry.key)
            , value(entry.value)
            , inserted(inserted)
        {
        }
    };

    explicit HashTable(HashType capacity = INITIAL_CAPACITY)
        : m_capacity(makePowerOfTwo(std::max(capacity, (HashType)8)))
        , m_mask(m_capacity - 1u)
        , m_size(0)
        , max_load_factor(DEFAULT_MAX_LOAD_FACTOR)
    {
        entries = alloc.allocate(m_capacity);
        hashes = new HashType[m_capacity];
        for (HashType i = 0; i < m_capacity; i++) {
            hashes[i] = 0;
        }
    }

    HashTable(HashTable&& other)
    {
        *this = std::move(other);
    }

    /* TODO: test
    HashTable(const HashTable& other)
    {
        *this = other;
    }
*/

    ~HashTable()
    {
        clear();
        if (hashes)
            delete[] hashes;
        hashes = nullptr;
        if (entries)
            alloc.deallocate(entries, m_capacity);
        entries = nullptr;
    }

    /**
      Inserts value if the key does not already exist in the hashtable
        returns reference to value in hashtable and whether a new value was added
     */
    template <class K, class V>
    InsertResult insert(K&& key, V&& value)
    {
        if (m_size + 1 >= max_load_factor * m_capacity) {
            rehash(2 * m_capacity);
        }
        HashType hash = computeHash(key);
        return insertHelper(hash, Entry{ std::forward<K>(key), std::forward<V>(value) });
    }

    /**
      Removes the entry from the hashtable for key, returning whether it was removed.
      The removed entry, if any will be moved into the passed in entry
      */
    bool remove(const Key& key, Entry& entry)
    {
        HashType position = getPosition(key);
        if (position == NOT_FOUND)
            return false;
        m_size--;

        for (HashType distance = 0; distance < m_capacity; distance++) {
            HashType next_position = probe(position);
            HashType& next_hash = hashes[next_position];
            if (next_hash == 0)
                break;

            HashType entry_distance = probeDistance(next_position, next_hash);
            if (entry_distance == 0)
                break;

            std::swap(entries[position], entries[next_position]);
            std::swap(hashes[position], hashes[next_position]);
            position = next_position;
        }
        hashes[position] = 0;
        entry = std::move(entries[position]);
        return true;
    }

    /**
      Like remove but deletes the entry
      */
    bool erase(const Key& key)
    {
        Entry entry;
        return remove(key, entry);
    }

    template <class K>
    Value& operator[](K&& key)
    {
        if (m_size + 1 >= max_load_factor * m_capacity) {
            rehash(2 * m_capacity);
        }
        HashType hash = computeHash(key);
        InsertResult new_entry = insertHelper(hash, Entry{ std::forward<K>(key), Value() });
        return new_entry.value;
    }

    // get can safely be used while iterating through the hashtable
    //  because it never inserts or removes entries
    Value* get(const Key& key) const
    {
        HashType position = getPosition(key);
        if (position == NOT_FOUND)
            return nullptr;
        else
            return &entries[position].value;
    }

    void clear()
    {
        if (hashes == nullptr) {
            assert(entries == nullptr && m_size == 0);
            return;
        }
        for (HashType i = 0; i < m_capacity; i++) {
            if (hashes[i] != 0) {
                entries[i].~Entry();
                hashes[i] = 0;
            }
        }
        m_size = 0;
    }

    /* TODO: test
    HashTable& operator=(const HashTable& other)
    {
        if (this != &other) {
            this->~HashTable();
            m_capacity = other.m_capacity;
            m_mask = other.m_mask;
            m_size = other.m_size;
            max_load_factor = other.max_load_factor;
            entries = alloc.allocate(m_capacity);
            hashes = new HashType[m_capacity];
            for(HashType i = 0; i < m_capacity; i++){
                hashes[i] = other.hashes[i];
                if(hashes[i] != 0)
                    alloc.construct(entries + i, other.entries[i]);
            }
        }
        return *this;
    }
*/
    HashTable& operator=(HashTable&& other)
    {
        if (this != &other) {
            this->~HashTable();
            entries = other.entries;
            hashes = other.hashes;
            m_capacity = other.m_capacity;
            m_mask = other.m_mask;
            m_size = other.m_size;
            max_load_factor = other.max_load_factor;
            other.entries = nullptr;
            other.hashes = nullptr;
            other.m_capacity = 0;
            other.m_mask = 0;
            other.m_size = 0;
        }
        return *this;
    }

    void print(std::ostream& os) const
    {
        for (HashType i = 0; i < m_capacity; i++) {
            os << i << ": ";
            HashType hash = hashes[i];
            if (hash != 0) {
                os << entries[i].key << ": ";
                os << entries[i].value << ", ";
                os << "distance: " << probeDistance(i, hash);
                os << " from " << desiredPosition(hash);
            }
            else {
                os << " #empty";
            }
            os << std::endl;
        }
    }

    double loadFactor() const
    {
        return (double)m_size / m_capacity;
    }

    HashType size() const
    {
        return m_size;
    }

    HashType capacity() const
    {
        return m_capacity;
    }

    template <bool is_const>
    HashType probeDistance(const Iterator<is_const>& it) const
    {
        return probeDistance(it.index, hashes[it.index]);
    }

    Iterator<false> begin()
    {
        for (HashType i = 0; i < m_capacity; i++)
            if (hashes[i] != 0)
                return Iterator<false>(*this, i);
        return Iterator<false>(*this, m_capacity);
    }

    Iterator<false> end()
    {
        return Iterator<false>(*this, m_capacity);
    }

    Iterator<true> begin() const
    {
        return cbegin();
    }

    Iterator<true> end() const
    {
        return cend();
    }

    Iterator<true> cbegin() const
    {
        for (HashType i = 0; i < m_capacity; i++)
            if (hashes[i] != 0)
                return Iterator<true>(*this, i);
        return Iterator<true>(*this, m_capacity);
    }

    Iterator<true> cend() const
    {
        return Iterator<true>(*this, m_capacity);
    }

    template <bool is_const>
    class Iterator {
        using HT = std::conditional_t<is_const,
            const HashTable&,
            HashTable&>;
        HT& hash_table;
        HashType index;
        friend class HashTable;

    public:
        using Ptr = std::conditional_t<is_const,
            const KeyValuePair<Key, Value>*,
            KeyValuePair<Key, Value>*>;
        using Ref = std::conditional_t<is_const,
            KeyValuePair<std::reference_wrapper<const Key>, std::reference_wrapper<const Value>>,
            KeyValuePair<std::reference_wrapper<const Key>, std::reference_wrapper<Value>>>;

        Iterator(HT& hash_table, HashType index)
            : hash_table(hash_table)
            , index(index)
        {
        }

        Iterator(const Iterator& other)
            : hash_table(other.hash_table)
            , index(other.index)
        {
        }
        explicit Iterator(const Iterator<false>& other, void** = nullptr)
            : hash_table(other.hash_table)
            , index(other.index)
        {
        }

        Iterator& operator++()
        {
            index++;
            for (; index < hash_table.m_capacity; index++)
                if (hash_table.hashes[index] != 0)
                    break;
            return *this;
        }

        Iterator operator++(int)
        {
            Iterator tmp(*this);
            operator++();
            return tmp;
        }

        TICK_MEMBER_REQUIRES(is_const == false)
        Ref operator*()
        {
            auto& entry = hash_table.entries[index];
            return Ref{ std::cref(entry.key), std::ref(entry.value) };
        }

        TICK_MEMBER_REQUIRES(is_const == true)
        Ref operator*()
        {
            auto& entry = hash_table.entries[index];
            return Ref{ std::cref(entry.key), std::cref(entry.value) };
        }

        Ptr operator->()
        {
            return &hash_table.entries[index];
        }

        bool operator==(const Iterator& other)
        {
            return (index == other.index);
        }

        bool operator!=(const Iterator& other)
        {
            return (index != other.index);
        }
    };

private:
    Entry* entries;
    HashType* hashes;
    HashType m_capacity;
    HashType m_mask;
    HashType m_size;
    float max_load_factor;
    Allocator alloc;
    Hash hasher;
    KeyEqual equals;

    static const HashType INITIAL_CAPACITY = 8;
    static constexpr float DEFAULT_MAX_LOAD_FACTOR = 0.9f;
    static const HashType NOT_FOUND = (HashType)-1;

    HashType computeHash(const Key& key) const
    {
        HashType hash = hasher(key) + 1;
        if (hash == 0)
            hash++;
        return hash;
    }

    void rehash(HashType new_m_capacity)
    {
        HashTable rehashed(new_m_capacity);

        for (HashType i = 0; i < m_capacity; i++) {
            if (hashes[i] != 0) {
                rehashed.insertHelper(hashes[i], std::move(entries[i]));
                entries[i].~Entry();
            }
        }
        assert(m_size == rehashed.m_size);

        delete[] hashes;
        hashes = nullptr;
        alloc.deallocate(entries, m_capacity);
        entries = nullptr;
        m_size = 0;

        *this = std::move(rehashed);
    }

    HashType desiredPosition(HashType hash) const
    {
        return hash & m_mask;
    }

    HashType getPosition(const Key& key) const
    {
        HashType my_hash = computeHash(key);
        HashType position = desiredPosition(my_hash);
        for (HashType distance = 0; distance < m_capacity; distance++) {
            HashType hash = hashes[position];
            if (hash == 0)
                return NOT_FOUND;
            if (my_hash == hash) {
                Entry& entry = entries[position];
                if (equals(entry.key, key))
                    return position;
            }

            HashType entry_distance = probeDistance(position, hash);
            if (entry_distance < distance)
                return NOT_FOUND;
            position = probe(position);
        }
        return NOT_FOUND;
    }

    HashType probe(HashType position) const
    {
        return (position + 1u) & m_mask;
    }

    HashType probeDistance(HashType position, HashType hash) const
    {
        return (position - desiredPosition(hash)) & m_mask;
    }

    InsertResult insertHelper(HashType new_hash, Entry&& new_entry)
    {
        HashType position = desiredPosition(new_hash);
        Entry* inserted_entry = nullptr;
        // pointer to location where we inserted the new entry
        for (HashType distance = 0; distance < m_capacity; distance++) {
            HashType& hash = hashes[position];
            if (hash == 0) {
                // Empty spot was found
                m_size++;
                Entry& entry = entries[position];
                alloc.construct(&entry, std::move(new_entry));
                hash = new_hash;
                if (inserted_entry)
                    return InsertResult(*inserted_entry, true);
                return InsertResult(entry, true);
            }
            if (hash == new_hash) {
                Entry& entry = entries[position];
                if (equals(entry.key, new_entry.key)) {
                    // Found existing entry with same key
                    if (inserted_entry)
                        return InsertResult(*inserted_entry, true);
                    return InsertResult(entry, false);
                }
            }
            HashType entry_distance = probeDistance(position, hash);
            if (entry_distance < distance) {
                // swap entry and new entry to preserve distance invariant
                Entry& entry = entries[position];
                if (inserted_entry == nullptr)
                    inserted_entry = &entry;
                std::swap(entry, new_entry);
                std::swap(hash, new_hash);
                distance = entry_distance;
            }
            position = probe(position);
        }
        ZIRAN_ASSERT(false, "HashTable is full");
    }
};

template <class Key, class Value, class Hash, class KeyEqual, class Allocator>
std::ostream& operator<<(std::ostream& os,
    const HashTable<Key, Value, Hash, KeyEqual, Allocator>& hash)
{
    hash.print(os);
    return os;
}
} // namespace ZIRAN
#endif
