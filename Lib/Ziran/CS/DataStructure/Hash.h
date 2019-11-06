#ifndef HASH_H
#define HASH_H

#include <Ziran/CS/Util/Meta.h>
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/Util/PlatformSpecific.h>
#include <Eigen/Core>
#include <memory>
#include <stdint.h>
#include <tick/requires.h>
//-----------------------------------------------------------------------------
// Platform-specific functions and macros

// Microsoft Visual Studio

#if defined(_MSC_VER)

#include <stdlib.h>

#define ZIRAN_ROTL64(x, y) _rotl64(x, y)

// Other compilers

#else // defined(_MSC_VER)

inline uint64_t rotl64(uint64_t x, int8_t r)
{
    return (x << r) | (x >> (64 - r));
}

#define ZIRAN_ROTL64(x, y) rotl64(x, y)

#endif // !defined(_MSC_VER)

namespace ZIRAN {
template <class T, class = void>
struct Hash;

inline void murmurHash3(const void* key, const int len,
    const uint32_t seed, void* out);

namespace INTERNAL {
template <class T, bool has = HasType<T>::value>
struct GetHashType {
    using Type = size_t;
};

template <class T>
struct GetHashType<T, true> {
    using Type = typename T::Type;
};
} // namespace INTERNAL

template <class T>
using GetHashType = typename INTERNAL::GetHashType<T>::Type;

template <class T, int m, int n>
struct Hash<Matrix<T, m, n>> {
    const uint32_t seed;

    Hash()
        : seed(54943)
    {
    }

    size_t operator()(const Matrix<T, m, n>& A) const
    {
        uint64_t hash[2];
        murmurHash3((const void*)A.data(), sizeof(T) * A.size(), seed, (void*)hash);
        return hash[0] + 7 * hash[1];
    }
};

template <class T>
struct Hash<T, TICK_CLASS_REQUIRES(std::is_arithmetic<T>::value && sizeof(T) <= 4)> {
    using Type = uint32_t;

    Hash()
    {
    }

    Type operator()(T A) const
    {
        uint32_t k;
        char* A_bytes = (char*)&A;
        char* k_bytes = (char*)&k;
        memcpy(k_bytes, A_bytes, sizeof(T));

        k ^= k >> 16;
        k *= 0x85ebca6b;
        k ^= k >> 13;
        k *= 0xc2b2ae35;
        k ^= k >> 16;

        return k;
    }
};

template <class T>
struct Hash<T, TICK_CLASS_REQUIRES(std::is_arithmetic<T>::value && sizeof(T) > 4 && sizeof(T) <= 8)> {

    Hash()
    {
    }

    size_t operator()(T A) const
    {
        uint64_t k;
        char* A_bytes = (char*)&A;
        char* k_bytes = (char*)&k;
        memcpy(k_bytes, A_bytes, sizeof(T));

        k ^= k >> 33;
        k *= ZIRAN_BIG_CONSTANT(0xff51afd7ed558ccd);
        k ^= k >> 33;
        k *= ZIRAN_BIG_CONSTANT(0xc4ceb9fe1a85ec53);
        k ^= k >> 33;

        return k;
    }
};

template <class T>
struct Hash<T, TICK_CLASS_REQUIRES(std::is_arithmetic<T>::value && sizeof(T) > 8)> {
    const uint32_t seed;

    Hash()
        : seed(54943)
    {
    }

    size_t operator()(T A) const
    {
        uint64_t hash[2];
        murmurHash3((const void*)&A, sizeof(T), seed, (void*)hash);
        return hash[0] + 7 * hash[1];
    }
};

template <class T>
struct Hash<std::basic_string<T>> {
    const uint32_t seed;

    Hash()
        : seed(54943)
    {
    }

    size_t operator()(const std::basic_string<T>& A) const
    {
        uint64_t hash[2];
        murmurHash3((const void*)A.data(), sizeof(T) * A.size(), seed, (void*)hash);
        return hash[0] + 7 * hash[1];
    }
};

/** 
  Hash based on the memory address value of the pointer
  */
template <class T>
struct Hash<T*> {
    Hash<uintptr_t> internal;

    Hash()
    {
    }

    size_t operator()(const T* ptr) const
    {
        return internal((uintptr_t)ptr);
    }
};

/** 
  Hash based on the memory address value of the pointer
  */
template <class T>
struct Hash<std::unique_ptr<T>> {
    Hash<uintptr_t> internal;

    Hash()
    {
    }

    size_t operator()(const std::unique_ptr<T>& A) const
    {
        return internal((uintptr_t)A.get());
    }
};

/** 
  Hash based on the memory address value of the pointer
  */
template <class T>
struct Hash<std::shared_ptr<T>> {
    Hash<uintptr_t> internal;

    Hash()
    {
    }

    size_t operator()(const std::shared_ptr<T>& A) const
    {
        return internal((uintptr_t)A.get());
    }
};
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

inline ZIRAN_FORCE_INLINE uint64_t
getBlock64(const uint64_t* p, int i)
{
    return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

inline ZIRAN_FORCE_INLINE uint64_t fmix64(uint64_t k)
{
    k ^= k >> 33;
    k *= ZIRAN_BIG_CONSTANT(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k *= ZIRAN_BIG_CONSTANT(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;

    return k;
}
//-----------------------------------------------------------------------------

// MurmurHash3 was written by Austin Appleby
void murmurHash3(const void* key, const int len,
    const uint32_t seed, void* out)
{
    const uint8_t* data = (const uint8_t*)key;
    const int nblocks = len / 16;

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = ZIRAN_BIG_CONSTANT(0x87c37b91114253d5);
    const uint64_t c2 = ZIRAN_BIG_CONSTANT(0x4cf5ad432745937f);

    //----------
    // body

    const uint64_t* blocks = (const uint64_t*)(key);

    for (int i = 0; i < nblocks; i++) {
        uint64_t k1 = getBlock64(blocks, i * 2 + 0);
        uint64_t k2 = getBlock64(blocks, i * 2 + 1);

        k1 *= c1;
        k1 = ZIRAN_ROTL64(k1, 31);
        k1 *= c2;
        h1 ^= k1;

        h1 = ZIRAN_ROTL64(h1, 27);
        h1 += h2;
        h1 = h1 * 5 + 0x52dce729;

        k2 *= c2;
        k2 = ZIRAN_ROTL64(k2, 33);
        k2 *= c1;
        h2 ^= k2;

        h2 = ZIRAN_ROTL64(h2, 31);
        h2 += h1;
        h2 = h2 * 5 + 0x38495ab5;
    }

    //----------
    // tail

    const uint8_t* tail = (const uint8_t*)(data + nblocks * 16);

    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch (len & 15) {
    case 15:
        k2 ^= ((uint64_t)tail[14]) << 48;
        ZIRAN_FALLTHROUGH;
    case 14:
        k2 ^= ((uint64_t)tail[13]) << 40;
        ZIRAN_FALLTHROUGH;
    case 13:
        k2 ^= ((uint64_t)tail[12]) << 32;
        ZIRAN_FALLTHROUGH;
    case 12:
        k2 ^= ((uint64_t)tail[11]) << 24;
        ZIRAN_FALLTHROUGH;
    case 11:
        k2 ^= ((uint64_t)tail[10]) << 16;
        ZIRAN_FALLTHROUGH;
    case 10:
        k2 ^= ((uint64_t)tail[9]) << 8;
        ZIRAN_FALLTHROUGH;
    case 9:
        k2 ^= ((uint64_t)tail[8]) << 0;
        k2 *= c2;
        k2 = ZIRAN_ROTL64(k2, 33);
        k2 *= c1;
        h2 ^= k2;
        ZIRAN_FALLTHROUGH;
    case 8:
        k1 ^= ((uint64_t)tail[7]) << 56;
        ZIRAN_FALLTHROUGH;
    case 7:
        k1 ^= ((uint64_t)tail[6]) << 48;
        ZIRAN_FALLTHROUGH;
    case 6:
        k1 ^= ((uint64_t)tail[5]) << 40;
        ZIRAN_FALLTHROUGH;
    case 5:
        k1 ^= ((uint64_t)tail[4]) << 32;
        ZIRAN_FALLTHROUGH;
    case 4:
        k1 ^= ((uint64_t)tail[3]) << 24;
        ZIRAN_FALLTHROUGH;
    case 3:
        k1 ^= ((uint64_t)tail[2]) << 16;
        ZIRAN_FALLTHROUGH;
    case 2:
        k1 ^= ((uint64_t)tail[1]) << 8;
        ZIRAN_FALLTHROUGH;
    case 1:
        k1 ^= ((uint64_t)tail[0]) << 0;
        k1 *= c1;
        k1 = ZIRAN_ROTL64(k1, 31);
        k1 *= c2;
        h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len;
    h2 ^= len;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
}

//-----------------------------------------------------------------------------
} // namespace ZIRAN
#endif // _MURMURHASH3_H_
