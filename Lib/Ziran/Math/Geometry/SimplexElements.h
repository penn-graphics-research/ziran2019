#ifndef SIMPLEX_ELEMENTS_H
#define SIMPLEX_ELEMENTS_H
#include <Ziran/Math/Geometry/ElementManager.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Math/Geometry/SimplexMesh.h>
#include <Ziran/Math/Linear/GivensQR.h>
#include <Ziran/CS/DataStructure/HashTable.h>
#include <Ziran/CS/Util/Timer.h>
#include <algorithm>
#include <metis.h>
#include <tbb/tbb.h>
#include <tick/requires.h>
namespace ZIRAN {

template <class T, int manifold_dim, int dim, class enable = void>
class SimplexElements;

template <class T, int _manifold_dim, int _dim>
class SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>> : public ElementManager<T, _dim> {
public:
    using Scalar = T;
    using ElementManager<T, _dim>::batches;
    using Base = ElementManager<T, _dim>;
    using MeshType = SimplexMesh<_manifold_dim>;
    static constexpr int num_vertices = MeshType::num_vertices;
    static const int manifold_dim = _manifold_dim;
    static const int dim = _dim;
    static_assert(manifold_dim <= dim, "A tetrahedron won't fit in 2d");
    using IV = Eigen::Matrix<int, manifold_dim + 1, 1>;
    using TM = Matrix<T, manifold_dim, manifold_dim>;
    using TM2 = Matrix<T, dim, manifold_dim>;
    using TVStack = Matrix<T, dim, Eigen::Dynamic>;
    using Base::add;
    using Base::appender;
    using Base::count;
    using Base::get;
    using Base::iter;
    using Base::reorder;

    template <class DataType>
    using PerQuadrature = DataType;

    using ShapeGradientType = Eigen::Matrix<T, manifold_dim, num_vertices>;

    DataArray<IV>& indices;
    DataArray<TM>& Dm_inv;
    DataArray<TM>& F;
    DataArray<T>& element_measure;

    // gradient of the isoparametric shape function
    static const ShapeGradientType grad_N_hat;

    StdVector<int> partition_offsets; // with size n+1, where n is # partitions. This is a state and will be read/write.

    // Helpers for parallelization:

    // for splatting to pads
    using Base::pads; // StdVector<TVStack> pads. they are just element partition scratch pads of particles for parallization. no write.
    using Base::pads_dx; // StdVector<TVStack> pads. they are just element partition scratch pads of particles for parallization. no write.
    using Base::particle_l2g; // StdVector<StdVector<int>> particle_l2g; particle index local to global. no write.

    StdVector<IV> local_indices; // mapping from element id to pad wise local particle indices. no write.

    // for gathering
    HashTable<int, int> particle_color; // global to partition id. no write.
    struct NeighborHelper {
        StdVector<int> neighbor_ids;
        StdVector<StdVector<int>> local_particle_ids;
    }; // this is a map from me to my neighbors who have the same colored particles with me. for gathering. no write.
    StdVector<NeighborHelper> neighbors_of_my_color; // This is a state and will be read/write.

    SimplexElements();

    virtual ~SimplexElements();

    inline static AttributeName<IV> indices_name()
    {
        return AttributeName<IV>("indices");
    }

    inline static AttributeName<TM> F_name()
    {
        return AttributeName<TM>("F");
    }

    inline static AttributeName<T> element_measure_name()
    {
        return AttributeName<T>("element measure");
    }

    inline static AttributeName<TM> Dm_inv_name()
    {
        return AttributeName<TM>("Dm inverse");
    }

    inline static AttributeName<TM2> Q_name()
    {
        return AttributeName<TM2>("Q");
    }

    //  template <class DataType>
    //  inline static AttributeName<DataType> attributeName(const std::string& name)
    //  {
    //      return AttributeName<DataType>(name);
    //  }

    void registerParticleReorderingCallback(Particles<T, dim>& particles) override;

    bool renumberParticles(const int particle_count, StdVector<int>& old_to_new, StdVector<int>& new_to_old) override;

    Range addUndeformedMesh(const SimplexMesh<manifold_dim>& mesh, const StdVector<Vector<T, dim>>& X0, size_t particle_offset = 0);

    T totalMeasure(Range element_range) const;

    template <class XArray> // XArray is StdVector<Vector<T, dim>> or TVStack
    void updateF(const XArray& X, const AttributeName<TM2>& F_name) // when called with dF name, writes to dF
    {
        auto master = iter(indices_name(), F_name, Dm_inv_name());
        tbb::blocked_range<int> range(0, master.common_ranges->length(), 1024);
        tbb::parallel_for(range, [&](const tbb::blocked_range<int>& subrange) {
            auto end = master + subrange.end();
            for (auto it = master + subrange.begin(); it != end; ++it) {
                const IV& these_indices = it.template get<0>();
                TM2& F = it.template get<1>();
                const TM& Dm_inverse = it.template get<2>();
                F = dS(these_indices, X) * Dm_inverse;
            }
        });
    }

    inline void copyToPad(const StdVector<Vector<T, dim>>& X, const int k)
    {
        TVStack& pad = pads[k];
        for (int t = 0; t < pad.cols(); t++) // t is pad k's local particle index
            pad.col(t) = X[particle_l2g[k][t]];
    }

    inline void copyToPad(const TVStack& X, const int k)
    {
        TVStack& pad = pads[k];
        for (int t = 0; t < pad.cols(); t++) // t is pad k's local particle index
            pad.col(t) = X.col(particle_l2g[k][t]);
    }

    template <class XArray, TICK_REQUIRES(dim != manifold_dim)> // XArray is StdVector<Vector<T, dim>> or TVStack
    void updateF(const XArray& X, const AttributeName<TM>& F_name) // when called with dF name, writes to dF
    {
        auto master = iter(indices_name(), F_name, Q_name(), Dm_inv_name());
        tbb::blocked_range<int> range(0, master.common_ranges->length(), 1024);
        tbb::parallel_for(range, [&](const tbb::blocked_range<int>& subrange) {
            auto end = master + subrange.end();
            for (auto it = master + subrange.begin(); it != end; ++it) {
                const IV& these_indices = it.template get<0>();
                TM& F = it.template get<1>();
                TM2& Q = it.template get<2>();
                const TM& Dm_inverse = it.template get<3>();
                TM2 F0 = dS(these_indices, X) * Dm_inverse;
                thinGivensQR(F0, Q, F);
            }
        });
    }

    void setMassFromDensity(T density, Range element_range, StdVector<T>& mass)
    {
        // add each element's contribution to its node masses
        for (int e = element_range.lower; e < element_range.upper; e++) {
            T fraction = density / (T)(manifold_dim + 1);
            const IV& these_indices = indices[e];
            T measure = element_measure[e];
            for (int i = 0; i < manifold_dim + 1; i++) {
                mass[these_indices[i]] += measure * fraction;
            }
        }
    }

    TM2 dS(const IV& these_indices, const StdVector<Vector<T, dim>>& X) const;

    TM2 dS(const IV& these_indices, const Matrix<T, dim, Eigen::Dynamic>& X) const;

    auto getVertices(const IV& these_indicies, const StdVector<Vector<T, dim>>& X) -> Matrix<T, dim, manifold_dim + 1> const;

    // need repartition if it is never partitioned, or element count changed.
    bool needRepartitioning() override;

    // use neighbors_of_my_color info to gather pads (no zero out)
    void gatherFromPadsToTVStack(const T scale, TVStack& force) override;

    void resizePadsAndBuildHelpers();

    // build neighbors_of_my_color.
    // when restarting, it is not called. the data is read from disk.
    void buildNeighbors();

    // partition_offsets stores the starting index of elements per partition
    void partitionAndReindexElements(const int n_partitions) override;

    void partitionElements(const int n_partitions, StdVector<StdVector<int>>& buckets);

    // partition_offsets stores n+1 indices, where n is the number of buckets
    void reindexElements(const StdVector<StdVector<int>>& buckets);

    virtual void writeData(std::ostream& out) const override;

    virtual void readData(std::istream& in) override;

    template <class Func, typename... Types>
    static void apply(Func& f, Types&... x)
    {
        f(x...);
    }
};
} // namespace ZIRAN
#endif
