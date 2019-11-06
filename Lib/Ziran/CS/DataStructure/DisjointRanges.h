#ifndef DISJOINT_RANGES_H
#define DISJOINT_RANGES_H

#include <tbb/tbb.h>

#include <iostream>
#include <vector>

#include <Ziran/CS/Util/Debug.h>

namespace ZIRAN {
struct Range {
    int lower;
    int upper;

    Range()
    {
    }

    Range(int lower, int upper)
        : lower(lower)
        , upper(upper)
    {
    }

    Range(const Range& range) = default;

    bool operator==(const Range& other) const
    {
        return (other.lower == lower) && (other.upper == upper);
    }
    int length() const
    {
        return upper - lower;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Range& m)
{
    return os << "[" << m.lower << ", " << m.upper << ")";
}

class DisjointRanges {
public:
    int lg2_grain_size;
    StdVector<Range> ranges;

    using ConstIterator = StdVector<Range>::const_iterator;
    using Iterator = StdVector<Range>::iterator;

    DisjointRanges(int lg2_grain_size = 7)
        : lg2_grain_size(lg2_grain_size)
    {
    }

    DisjointRanges(std::initializer_list<Range> l, int lg2_grain_size = 7)
        : lg2_grain_size(lg2_grain_size)
    {
        for (auto r : l)
            append(r);
    }

    /**
      Intersection constructor
    */
    template <class... Rest>
    DisjointRanges(const DisjointRanges& dr, const Rest&... rest)
        : lg2_grain_size(dr.lg2_grain_size)
        , ranges(dr.ranges)
    {
        intersect(rest...);
    }

    // tbb splitting constructors are only defined if tbb/tbb_stddef.h was included
    DisjointRanges(DisjointRanges& other, tbb::split)
        : lg2_grain_size(other.lg2_grain_size)
    {
        int length = other.length();
        if (length <= 1 << lg2_grain_size)
            return;
        int num_grains = length >> lg2_grain_size;
        int m = (num_grains - num_grains / 2) << lg2_grain_size;

        auto last = other.cbegin();
        auto iter = other.cbegin();

        for (; m > 0; ++iter) {
            last = iter;
            m -= iter->upper - iter->lower;
        }
        if (m == 0) {
            // Split along existing boundary
            ranges.insert(ranges.end(), iter, other.cend());
            other.ranges.erase(iter, other.cend());
        }
        else {
            ranges.insert(ranges.end(), last, other.cend());
            other.ranges.erase(iter, other.cend());

            other.ranges.back().upper += m;
            ranges.begin()->lower = other.ranges.back().upper;
        }
        assert(length == (other.length() + this->length()));
    }

    void append(const Range& range)
    {
        ZIRAN_ASSERT(ranges.empty() || range.lower >= ranges.back().upper, "New ranges must come after old ones");
        if (!ranges.empty() && ranges.back().upper == range.lower)
            ranges.back().upper = range.upper;
        else {
            if (range.upper > range.lower)
                ranges.push_back(range);
        }
    }

    void append(int lower, int upper)
    {
        append(Range{ lower, upper });
    }

    Iterator begin()
    {
        return ranges.begin();
    }

    Iterator end()
    {
        return ranges.end();
    }

    ConstIterator begin() const
    {
        return ranges.cbegin();
    }

    ConstIterator end() const
    {
        return ranges.cend();
    }

    ConstIterator cbegin() const
    {
        return ranges.cbegin();
    }

    ConstIterator cend() const
    {
        return ranges.cend();
    }

    const Range& operator[](size_t i) const
    {
        return ranges[i];
    }

    Range& operator[](size_t i)
    {
        return ranges[i];
    }

    bool empty() const
    {
        return length() == 0;
    }

    // for tbb
    bool is_divisible() const
    {
        return length() > lg2_grain_size;
    }

    // The number of disjoint ranges
    size_t size() const
    {
        return ranges.size();
    }

    // The total length of the ranges
    int length() const
    {
        int l = 0;
        for (auto range : ranges)
            l += (range.upper - range.lower);
        return l;
    }

    bool operator==(const DisjointRanges& other) const
    {
        return ranges == other.ranges;
    }

    bool valid()
    {
        int last_upper = std::numeric_limits<int>::min();
        for (auto range : ranges) {
            if (range.upper <= range.lower)
                return false;
            if (last_upper >= range.lower)
                return false;
            last_upper = range.upper;
        }
        return true;
    }

    friend void swap(DisjointRanges& a, DisjointRanges& b)
    {
        using std::swap;
        swap(a.lg2_grain_size, b.lg2_grain_size);
        swap(a.ranges, b.ranges);
    }

    void intersect()
    {
    }

    template <class... Rest>
    void intersect(
        const DisjointRanges& ranges1,
        const Rest&... rest)
    {
        auto r0 = cbegin();
        if (r0 == cend())
            return;
        {
            DisjointRanges out(lg2_grain_size);

            auto r1 = ranges1.cbegin();
            if (r1 == ranges1.cend()) {
                ranges.clear();
                return;
            }

            Range p;
            if (r1->lower < r0->lower) {
                p = *r1;
                r1++;
            }
            else {
                p = *r0;
                r0++;
            }
            while (r1 != ranges1.cend() || r0 != cend()) {
                Range q;
                if (r1 == ranges1.cend()) {
                    // Only range0 is left
                    q = *r0;
                    r0++;
                }
                else if (r0 == cend() || r1->lower < r0->lower) {
                    // r1 comes before r0
                    q = *r1;
                    r1++;
                }
                else {
                    q = *r0;
                    r0++;
                }
                if (q.lower < p.upper) {
                    // p and q intersect
                    p.lower = q.lower;
                    if (q.upper < p.upper) {
                        // q is a subset of p
                        out.append(q);
                    }
                    else {
                        out.append(p);
                        // need to check intersection of p with later ranges
                        p = q;
                    }
                }
                else {
                    // need to check intersection of p with later ranges
                    p = q;
                }
            }
            swap(*this, out);
        }
        intersect(rest...);
    }

    /**
      returns the numbers that are in r but not in this
      (as disjoint ranges)
      */
    DisjointRanges complement(Range r)
    {
        DisjointRanges out;
        for (auto dr : *this) {
            if (r.upper <= dr.lower) {
                break;
            }
            if (r.lower <= dr.lower) {
                out.append(Range{ r.lower, dr.lower });
                r.lower = std::min(dr.upper, r.upper);
            }
        }
        out.append(r);
        return out;
    }

    void merge()
    {
    }

    template <class... Rest>
    void merge(
        const DisjointRanges& ranges1,
        const Rest&... rest)
    {
        auto r0 = cbegin();
        auto r1 = ranges1.cbegin();
        if (r0 == cend()) {
            *this = ranges1;
        }
        else if (r1 != ranges1.cend()) {
            DisjointRanges out(lg2_grain_size);
            Range p;
            if (r1->lower < r0->lower) {
                p = *r1;
                r1++;
            }
            else {
                p = *r0;
                r0++;
            }
            while (r1 != ranges1.cend() || r0 != cend()) {
                Range q;
                if (r1 == ranges1.cend()) {
                    // Only range0 is left
                    q = *r0;
                    r0++;
                }
                else if (r0 == cend() || r1->lower < r0->lower) {
                    // r1 comes before r0
                    q = *r1;
                    r1++;
                }
                else {
                    q = *r0;
                    r0++;
                }
                if (q.lower < p.upper) {
                    // p and q intersect
                    p.upper = std::max(q.upper, p.upper);
                }
                else {
                    out.append(p);
                    p = q;
                }
            }
            out.append(p);
            swap(*this, out);
        }
        merge(rest...);
    }
};

inline std::ostream& operator<<(std::ostream& os, const DisjointRanges& m)
{
    for (size_t i = 0; i < m.size(); i++)
        os << m.ranges[i] << " ";
    return os;
}
} // namespace ZIRAN
#endif
