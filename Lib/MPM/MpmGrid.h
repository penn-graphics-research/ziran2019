#ifndef MPM_GRID_H
#define MPM_GRID_H
#include <Ziran/Math/Splines/BSplines.h>
#include <Ziran/CS/Util/Forward.h>
#include <MPM/Forward/MpmForward.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Page_Map.h>
#include <type_traits>

namespace ZIRAN {

using namespace SPGrid;

template <class T, int dim>
class GridState {
public:
    // dont change the order here
    Vector<T, dim> v;
    T m;
    Vector<T, dim> new_v;
    typename std::conditional<std::is_same<T, float>::value, int32_t, int64_t>::type idx;

    Vector<T, dim + dim> padding;
    T phase_field;
    T phase_field_multiplier;

    GridState()
    {
        v = Vector<T, dim>::Zero();
        m = (T)0;
        new_v = Vector<T, dim>::Zero();
        idx = -1;
    }
};

inline constexpr bool is_power_of_two(size_t x)
{
    return x > 0 && (x & (x - 1)) == 0;
}
static_assert(is_power_of_two(sizeof(GridState<float, 2>)), "GridState<float, 2> size must be POT");
static_assert(is_power_of_two(sizeof(GridState<float, 3>)), "GridState<float, 3> size must be POT");
static_assert(is_power_of_two(sizeof(GridState<double, 2>)), "GridState<double, 2> size must be POT");
static_assert(is_power_of_two(sizeof(GridState<double, 3>)), "GridState<double, 3> size must be POT");

inline std::array<int, 2> to_std_array(Vector<int, 2> v)
{
    return std::array<int, 2>{ v[0], v[1] };
}

inline std::array<int, 3> to_std_array(Vector<int, 3> v)
{
    return std::array<int, 3>{ v[0], v[1], v[2] };
}

template <class T, int dim>
class BSplineWeights {
    const T dx;

public:
    static constexpr int interpolation_degree = ZIRAN_MPM_DEGREE;
    Vector<T, interpolation_degree + 1> w[dim];
    Vector<T, interpolation_degree + 1> dw[dim];
    const T one_over_dx;
    Vector<int, dim> base_node;

    BSplineWeights(const Vector<T, dim>& X, T dx)
        : dx(dx), one_over_dx(1 / dx), base_node(Vector<int, dim>::Zero())
    {
        compute(X);
    }

    void compute(const Vector<T, dim>& X)
    {
        Vector<T, dim> X_index_space = one_over_dx * X;
        for (int d = 0; d < dim; d++)
            computeBSplineWeights(X_index_space[d], base_node[d], w[d], &dw[d]);
    }
};

template <class T, int dim>
class MpmGrid {
public:
    static constexpr int log2_page = 12;
    static constexpr int spgrid_size = 4096;
    static constexpr int interpolation_degree = ZIRAN_MPM_DEGREE;
    static constexpr int kernel_size = (dim == 2)
        ? (interpolation_degree + 1) * (interpolation_degree + 1)
        : (interpolation_degree + 1) * (interpolation_degree + 1) * (interpolation_degree + 1);
    using SparseGrid = SPGrid_Allocator<GridState<T, dim>, dim, log2_page>;
    using SparseMask = typename SparseGrid::template Array_type<>::MASK;
    using PageMap = SPGrid_Page_Map<log2_page>;
    std::unique_ptr<SparseGrid> grid;
    std::unique_ptr<PageMap> page_map;
    __m128 __attribute__((aligned(64))) __coord_delta[kernel_size];

    MpmGrid()
    {
        if constexpr (dim == 2) {
            grid = std::make_unique<SparseGrid>(spgrid_size, spgrid_size);
        }
        else {
            grid = std::make_unique<SparseGrid>(spgrid_size, spgrid_size, spgrid_size);
            int cnt = 0;
            for (int i = 0; i < interpolation_degree + 1; ++i)
                for (int j = 0; j < interpolation_degree + 1; ++j)
                    for (int k = 0; k < interpolation_degree + 1; ++k) {
                        __coord_delta[cnt++] = _mm_set_ps((float)0, (float)k, (float)j, (float)i);
                    }
        }
        page_map = std::make_unique<PageMap>(*grid);
    }

    GridState<T, dim>& operator[](const Vector<int, dim>& v)
    {
        return grid->GetArray()(to_std_array(v));
    }

    const GridState<T, dim>& operator[](const Vector<int, dim>& v) const
    {
        return grid->GetArray()(to_std_array(v));
    }

    int getNumNodes()
    {
        // TODO: parallelize
        int total_num_node = 0;
        auto blocks = page_map->Get_Blocks();
        auto grid_array = grid->Get_Array();
        for (int b = 0; b < (int)blocks.second; ++b) {
            GridState<T, dim>* g = reinterpret_cast<GridState<T, dim>*>(&grid_array(blocks.first[b]));
            for (int i = 0; i < (int)SparseMask::elements_per_block; ++i)
                if (g[i].m != 0)
                    g[i].idx = total_num_node++;
        }
        return total_num_node;
    }

    // no guarantee to only iterate valid grid
    template <typename OP>
    void iterateTouchedGrid(const OP& target)
    {
        auto blocks = page_map->Get_Blocks();
        auto grid_array = grid->Get_Array();
        if constexpr (dim == 2) {
            tbb::parallel_for(0, (int)blocks.second, [&](int b) {
                auto base_offset = blocks.first[b];
                auto stdarray_base_coord = SparseMask::LinearToCoord(base_offset);
                Vector<int, dim> base_coord{ stdarray_base_coord[0], stdarray_base_coord[1] };
                auto x = 1 << SparseMask::block_xbits;
                auto y = 1 << SparseMask::block_ybits;
                for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j) {
                        auto offset = SparseMask::Packed_Add(base_offset, SparseMask::Linear_Offset(i, j));
                        GridState<T, dim>& g = reinterpret_cast<GridState<T, dim>&>(grid_array(offset));
                        target(base_coord + Vector<int, dim>{ i, j }, g);
                    }
            });
        }
        else {
            tbb::parallel_for(0, (int)blocks.second, [&](int b) {
                auto base_offset = blocks.first[b];
                auto stdarray_base_coord = SparseMask::LinearToCoord(base_offset);
                Vector<int, dim> base_coord{ stdarray_base_coord[0], stdarray_base_coord[1], stdarray_base_coord[2] };
                auto x = 1 << SparseMask::block_xbits;
                auto y = 1 << SparseMask::block_ybits;
                auto z = 1 << SparseMask::block_zbits;
                for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j)
                        for (int k = 0; k < z; ++k) {
                            auto offset = SparseMask::Packed_Add(base_offset, SparseMask::Linear_Offset(i, j, k));
                            GridState<T, dim>& g = reinterpret_cast<GridState<T, dim>&>(grid_array(offset));
                            target(base_coord + Vector<int, dim>{ i, j, k }, g);
                        }
            });
        }
    }

    // guarantee to only iterate valid grid, whose mass > 0
    template <typename OP>
    void iterateGrid(const OP& target)
    {
        auto blocks = page_map->Get_Blocks();
        auto grid_array = grid->Get_Array();
        if constexpr (dim == 2) {
            tbb::parallel_for(0, (int)blocks.second, [&](int b) {
                auto base_offset = blocks.first[b];
                auto stdarray_base_coord = SparseMask::LinearToCoord(base_offset);
                Vector<int, dim> base_coord{ stdarray_base_coord[0], stdarray_base_coord[1] };
                auto x = 1 << SparseMask::block_xbits;
                auto y = 1 << SparseMask::block_ybits;
                for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j) {
                        auto offset = SparseMask::Packed_Add(base_offset, SparseMask::Linear_Offset(i, j));
                        GridState<T, dim>& g = reinterpret_cast<GridState<T, dim>&>(grid_array(offset));
                        if (g.idx >= 0)
                            target(base_coord + Vector<int, dim>{ i, j }, g);
                    }
            });
        }
        else {
            tbb::parallel_for(0, (int)blocks.second, [&](int b) {
                auto base_offset = blocks.first[b];
                auto stdarray_base_coord = SparseMask::LinearToCoord(base_offset);
                Vector<int, dim> base_coord{ stdarray_base_coord[0], stdarray_base_coord[1], stdarray_base_coord[2] };
                auto x = 1 << SparseMask::block_xbits;
                auto y = 1 << SparseMask::block_ybits;
                auto z = 1 << SparseMask::block_zbits;
                for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j)
                        for (int k = 0; k < z; ++k) {
                            auto offset = SparseMask::Packed_Add(base_offset, SparseMask::Linear_Offset(i, j, k));
                            GridState<T, dim>& g = reinterpret_cast<GridState<T, dim>&>(grid_array(offset));
                            if (g.idx >= 0)
                                target(base_coord + Vector<int, dim>{ i, j, k }, g);
                        }
            });
        }
    }

    template <typename OP>
    inline void iterateKernel(const BSplineWeights<T, dim>& spline, uint64_t base_offset, const OP& target)
    {
        auto grid_array = grid->Get_Array();
        T one_over_dx = spline.one_over_dx;
        auto& w = spline.w;
        auto& dw = spline.dw;
        const Vector<int, dim>& base_coord = spline.base_node;
        Vector<int, dim> coord;
        if constexpr (dim == 2) {
            for (int i = 0; i < interpolation_degree + 1; ++i) {
                T wi = w[0](i);
                T dwidxi = one_over_dx * dw[0](i);
                coord[0] = base_coord[0] + i;
                for (int j = 0; j < interpolation_degree + 1; ++j) {
                    T wj = w[1](j);
                    T wij = wi * wj;
                    T dwijdxi = dwidxi * wj;
                    T dwijdxj = wi * one_over_dx * dw[1](j);
                    coord[1] = base_coord[1] + j;
                    auto offset = SparseMask::Packed_Add(base_offset, SparseMask::Linear_Offset(i, j));
                    GridState<T, dim>& g = reinterpret_cast<GridState<T, dim>&>(grid_array(offset));
                    target(coord, wij, Vector<T, dim>{ dwijdxi, dwijdxj }, g);
                }
            }
        }
        else {
            for (int i = 0; i < interpolation_degree + 1; ++i) {
                T wi = w[0](i);
                T dwidxi = one_over_dx * dw[0](i);
                coord[0] = base_coord[0] + i;
                for (int j = 0; j < interpolation_degree + 1; ++j) {
                    T wj = w[1](j);
                    T wij = wi * wj;
                    T dwijdxi = dwidxi * wj;
                    T dwijdxj = wi * one_over_dx * dw[1](j);
                    coord[1] = base_coord[1] + j;
                    for (int k = 0; k < interpolation_degree + 1; ++k) {
                        coord[2] = base_coord[2] + k;
                        T wk = w[2](k);
                        T wijk = wij * wk;
                        T wijkdxi = dwijdxi * wk;
                        T wijkdxj = dwijdxj * wk;
                        T wijkdxk = wij * one_over_dx * dw[2](k);
                        auto offset = SparseMask::Packed_Add(base_offset, SparseMask::Linear_Offset(i, j, k));
                        GridState<T, dim>& g = reinterpret_cast<GridState<T, dim>&>(grid_array(offset));
                        target(coord, wijk, Vector<T, dim>{ wijkdxi, wijkdxj, wijkdxk }, g);
                    }
                }
            }
        }
    }

    template <typename OP>
    inline void superIterateKernel(const BSplineWeights<T, dim>& spline, const OP& target)
    {
        auto grid_array = grid->Get_Array();
        T one_over_dx = spline.one_over_dx;
        auto& w = spline.w;
        auto& dw = spline.dw;
        const Vector<int, dim>& base_coord = spline.base_node;
        __m128 __base_coord = _mm_set_ps((float)0, (float)base_coord[2], (float)base_coord[1], (float)base_coord[0]);
        int cnt = 0;
        for (int i = 0; i < interpolation_degree + 1; ++i) {
            T wi = w[0](i);
            T dwidxi = one_over_dx * dw[0](i);
            for (int j = 0; j < interpolation_degree + 1; ++j) {
                T wj = w[1](j);
                T wij = wi * wj;
                T dwijdxi = dwidxi * wj;
                T dwijdxj = wi * one_over_dx * dw[1](j);
                for (int k = 0; k < interpolation_degree + 1; ++k) {
                    T wk = w[2](k);
                    T wijk = wij * wk;
                    T wijkdxi = dwijdxi * wk;
                    T wijkdxj = dwijdxj * wk;
                    T wijkdxk = wij * one_over_dx * dw[2](k);
                    __m128 __dw_w = _mm_set_ps(wijk, wijkdxk, wijkdxj, wijkdxi);
                    __m128 __coord = _mm_add_ps(__base_coord, __coord_delta[cnt]);
                    target(cnt, __coord, __dw_w);
                    cnt++;
                }
            }
        }
    }
};

template <class T, int dim>
class CachedOffset {
    using SparseGrid = typename MpmGrid<T, dim>::SparseGrid;
    using SparseMask = typename MpmGrid<T, dim>::SparseMask;
    static constexpr int kernel_size = MpmGrid<T, dim>::kernel_size;
    static constexpr int interpolation_degree = MpmGrid<T, dim>::interpolation_degree;
    uint64_t data[1 << SparseMask::block_bits][kernel_size];

public:
    CachedOffset(uint64_t particle_offset)
    {
        uint64_t base_offset = particle_offset & 0xfffffffffffff000ul;
        for (int delta = 0; delta < 1 << SparseMask::block_bits; ++delta) {
            uint64_t offset = (base_offset | (delta << SparseMask::data_bits));
            int cnt = 0;
            for (int i = 0; i < interpolation_degree + 1; ++i)
                for (int j = 0; j < interpolation_degree + 1; ++j)
                    for (int k = 0; k < interpolation_degree + 1; ++k)
                        data[delta][cnt++] = SparseMask::Packed_Add(offset, SparseMask::Linear_Offset(i, j, k));
        }
    }
    inline uint64_t* getCachedLine(uint64_t particle_offset)
    {
        return data[(particle_offset & 0x0000000000000ffful) >> SparseMask::data_bits];
    }
};

} // namespace ZIRAN

#endif
