#include "SimplexElements.h"
#include <metis.h>
#include <tbb/tbb.h>
#include <Ziran/CS/Util/Timer.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Math/Geometry/SimplexMesh.h>

namespace ZIRAN {

template <class T, int _manifold_dim, int _dim>
SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::SimplexElements()
    : Base(indices_name(), Dm_inv_name(), F_name(), element_measure_name())
    , indices(get(indices_name()))
    , Dm_inv(get(Dm_inv_name()))
    , F(get(F_name()))
    , element_measure(get(element_measure_name()))
{
    if (dim != manifold_dim)
        add(Q_name());
}

template <class T, int _manifold_dim, int _dim>
SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::~SimplexElements()
{
}

template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::registerParticleReorderingCallback(Particles<T, dim>& particles)
{
    particles.X.registerReorderingCallback(
        [&](const DataArrayBase::ValueIDToNewValueID& old_to_new) {
            for (auto& element : indices.array)
                for (int d = 0; d < manifold_dim + 1; d++)
                    element(d) = old_to_new(element(d));

            for (size_t i = 0; i < particle_l2g.size(); i++) {
                StdVector<int>& zz = particle_l2g[i];
                for (auto& p : zz)
                    p = old_to_new(p);
            }

            HashTable<int, int> particle_color_new;
            for (auto iter : particle_color)
                particle_color_new[old_to_new(iter.key)] = iter.value;

            particle_color = std::move(particle_color_new);
        });
}

template <class T, int _manifold_dim, int _dim>
bool SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::renumberParticles(const int particle_count, StdVector<int>& old_to_new, StdVector<int>& new_to_old)
{
    new_to_old.resize(particle_count, -1);
    old_to_new.resize(particle_count, -1);

    int count = 0;
    for (auto element : indices.array) {
        for (int d = 0; d < manifold_dim + 1; d++) {
            int old_p = element(d);
            if (old_to_new[old_p] == -1) {
                old_to_new[old_p] = count;
                new_to_old[count] = old_p;
                count++;
            }
        }
    }

    for (int i = 0; i < particle_count; i++) {
        if (old_to_new[i] == -1) {
            old_to_new[i] = count;
            new_to_old[count] = i;
            count++;
        }
    }

    // sanity checks
    ZIRAN_ASSERT(count == particle_count);
    for (int i = 0; i < particle_count; i++) {
        ZIRAN_ASSERT(new_to_old[old_to_new[i]] == i);
        ZIRAN_ASSERT(new_to_old[i] >= 0 && new_to_old[i] < count);
        ZIRAN_ASSERT(old_to_new[i] >= 0 && old_to_new[i] < count);
    }

    return true;
}

template <class T, int _manifold_dim, int _dim>
Range SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::addUndeformedMesh(const SimplexMesh<manifold_dim>& mesh, const StdVector<Vector<T, dim>>& X0, size_t particle_offset)
{
    //set up computation of element volume/area
    T d_factorial = (T)1;
    for (int i = 1; i < manifold_dim + 1; i++)
        d_factorial *= i;
    T unit_simplex_measure = 1 / d_factorial;
    //the unit simplex has measure 1 / d_factorial and the element measure is the product of the determinant
    //of the mapping from the unit simplex to the current simplex and the volumn of the unit simplex
    auto ap = appender(indices_name(), Dm_inv_name(), element_measure_name(), F_name());

    Range r;
    r.lower = ap.entryId();
    for (size_t e = 0; e < mesh.numberElements(); e++) {
        IV indices = mesh.indices[e].array() + particle_offset;
        TM2 Ds = dS(indices, X0);
        Matrix<T, dim, dim> Q;
        inplaceGivensQR(Ds, Q); // Ds is dim x manifold_dim  ,  Q is dim x dim
        TM Dm = Ds.template topLeftCorner<manifold_dim, manifold_dim>();
        T element_measure = std::abs(Dm.template triangularView<Eigen::Upper>().determinant() * unit_simplex_measure);

        // Dm_inverse is manifold_dim x manifold_dim
        TM Dm_inverse = Dm.template triangularView<Eigen::Upper>().solve(TM::Identity());

        // For full dimensional mesh, Dm_inverse stores the real Dm_inverse
        if (dim == manifold_dim)
            Dm_inverse *= Q.template topLeftCorner<manifold_dim, manifold_dim>().transpose();

        // For co-dimensional mesh, Dm_inverse stores (R(Ds))^{-1} (e.g. isotropic cloth in 3d)
        // The rotational part is discarded.

        TM F = TM::Identity();
        ap.append(indices, Dm_inverse, element_measure, F);
    }

    r.upper = ap.entryId();

    if (dim != manifold_dim)
        add(Q_name(), r, TM2::Identity()); // Q_name represents world space rotation

    return r;
}

template <class T, int manifold_dim, int dim>
T SimplexElements<T, manifold_dim, dim, std::enable_if_t<manifold_dim <= dim>>::totalMeasure(Range element_range) const
{
    T result = 0;
    for (int e = element_range.lower; e < element_range.upper; e++)
        result += element_measure[e];
    return result;
}

template <class T, int _manifold_dim, int _dim>
typename SimplexElements<T, _manifold_dim, _dim>::ShapeGradientType computeShapeGradients()
{
    typename SimplexElements<T, _manifold_dim, _dim>::ShapeGradientType grad_N_hat;
    grad_N_hat << Vector<T, _manifold_dim>::Constant(-1), Matrix<T, _manifold_dim, _manifold_dim>::Identity();
    return grad_N_hat;
}

template <class T, int _manifold_dim, int _dim>
const typename SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::ShapeGradientType
    SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::grad_N_hat
    = computeShapeGradients<T, _manifold_dim, _dim>();

template <class T, int _manifold_dim, int _dim>
typename SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::TM2 SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::dS(const IV& these_indices, const StdVector<Vector<T, dim>>& X) const
{
    TM2 result;
    for (int i = 1; i < manifold_dim + 1; i++)
        result.col(i - 1) = X[these_indices[i]] - X[these_indices[0]];
    return result;
}

template <class T, int _manifold_dim, int _dim>
auto SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::getVertices(const IV& these_indicies, const StdVector<Vector<T, dim>>& X) -> Matrix<T, dim, manifold_dim + 1> const
{
    Matrix<T, dim, manifold_dim + 1> result;
    for (int i = 0; i < manifold_dim + 1; i++)
        result.col(i) = X[these_indicies[i]];
    return result;
}

template <class T, int _manifold_dim, int _dim>
typename SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::TM2 SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::dS(const IV& these_indices, const Matrix<T, dim, Eigen::Dynamic>& X) const
{
    TM2 result;
    for (int i = 1; i < manifold_dim + 1; i++)
        result.col(i - 1) = X.col(these_indices[i]) - X.col(these_indices[0]);
    return result;
}

// need repartition if it is never partitioned, or element count changed.
template <class T, int _manifold_dim, int _dim>
bool SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::needRepartitioning()
{
    return partition_offsets.size() == 0 || partition_offsets.back() != count;
}

// use neighbors_of_my_color info to gather pads (no zero out)
template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::gatherFromPadsToTVStack(const T scale, TVStack& force)
{
    const int N_partitions = neighbors_of_my_color.size();
    tbb::parallel_for(tbb::blocked_range<int>(0, N_partitions),
        [&](const tbb::blocked_range<int>& subrange) {
            for (int i = subrange.begin(); i < subrange.end(); ++i) {
                NeighborHelper& nh = neighbors_of_my_color[i];
                for (size_t k = 0; k < nh.neighbor_ids.size(); ++k) {
                    int neighbor_id = nh.neighbor_ids[k];
                    for (auto& p : nh.local_particle_ids[k]) {
                        int global = particle_l2g[neighbor_id][p];
                        force.col(global) += scale * pads[neighbor_id].col(p);
                    }
                }
            }
        });
}

template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::resizePadsAndBuildHelpers()
{
    ZIRAN_QUIET_TIMER();
    pads.resize(partition_offsets.size() - 1);
    pads_dx.resize(partition_offsets.size() - 1);

    local_indices.clear();
    particle_l2g.resize(partition_offsets.size() - 1);
    for (size_t i = 0; i < partition_offsets.size() - 1; ++i) // loop over all partitions
    {
        particle_l2g[i].clear();
        int start_element = partition_offsets[i]; // [
        int end_element = partition_offsets[i + 1]; // )
        HashTable<int, int> particles; // key is global particle ids in this parititon, value is local index
        int count = 0;
        for (int e = start_element; e < end_element; e++) {
            IV local_element;
            for (int d = 0; d < manifold_dim + 1; d++) {
                int global_p = indices[e](d);
                auto found = particles.insert(global_p, count);
                if (found.inserted) {
                    particle_l2g[i].emplace_back(global_p);
                    count++;
                }
                local_element(d) = found.value;
            }
            local_indices.emplace_back(local_element);
        }

        // rezize pads
        size_t n_particles = particles.size();
        pads_dx[i].resize(dim, n_particles);
        pads[i].resize(dim, n_particles);
        pads[i].setZero();
    }
}

// build neighbors_of_my_color.
// when restarting, it is not called. the data is read from disk.
template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::buildNeighbors()
{
    ZIRAN_QUIET_TIMER();
    neighbors_of_my_color.resize(particle_l2g.size());
    for (size_t i = 0; i < particle_l2g.size(); ++i) // loop over all partitions
    {
        int my_color = i;
        const auto& my_l2g = particle_l2g[i];
        for (size_t local_pid = 0; local_pid < my_l2g.size(); local_pid++) {
            int global_pid = my_l2g[local_pid];
            int p_color = particle_color[global_pid];
            auto& neighbor = neighbors_of_my_color[p_color];

            auto spot = std::find(neighbor.neighbor_ids.begin(), neighbor.neighbor_ids.end(), my_color);
            size_t index = spot - neighbor.neighbor_ids.begin();
            if (spot == neighbor.neighbor_ids.end()) {
                neighbor.neighbor_ids.emplace_back(my_color);
                neighbor.local_particle_ids.emplace_back();
            }
            neighbor.local_particle_ids[index].emplace_back(local_pid);
        }
    }
}

// partition_offsets stores the starting index of elements per partition
template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::partitionAndReindexElements(const int n_partitions)
{
    if (count == 0)
        return;
    ZIRAN_QUIET_TIMER();
    ZIRAN_ASSERT(n_partitions > 0); // otherwise shouldn't call this function
    StdVector<StdVector<int>> buckets;
    partitionElements(n_partitions, buckets);
    reindexElements(buckets);

    // when restarting, only this function is called.
    resizePadsAndBuildHelpers();

    buildNeighbors();
}

template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::partitionElements(int n_partitions, StdVector<StdVector<int>>& buckets)
{
    if (n_partitions > count)
        n_partitions = count;
    ZIRAN_QUIET_TIMER();
    buckets.resize(n_partitions);
    for (auto& b : buckets)
        b.clear();

    // reindex node indices to a [0, n] integer space
    HashTable<int, int> old_to_new;
    StdVector<int> new_indices;
    int new_id = 0;
    for (auto it = iter(indices_name()); it; ++it) {
        const IV& these_indices = it.template get<0>();
        for (int d = 0; d < manifold_dim + 1; d++) {
            int old_id = these_indices(d);
            auto query = old_to_new.insert(old_id, new_id);
            if (query.inserted) {
                new_id++;
            }
            new_indices.push_back(query.value);
        }
    }
    ZIRAN_ASSERT(new_indices.size() == (size_t)((manifold_dim + 1) * count));

    int ne = count;
    int nn = new_id;
    int* eind = new_indices.data();
    StdVector<int> eptr;
    for (int i = 0; i <= ne; i++)
        eptr.push_back(i * (manifold_dim + 1));
    int objval;
    StdVector<int> epart, npart;
    epart.resize(ne);
    npart.resize(nn);
    int ncommon = manifold_dim;
    int nparts = (int)buckets.size();
    METIS_PartMeshDual(&ne, &nn, eptr.data(), eind, NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval, epart.data(), npart.data());

    for (int i = 0; i < ne; i++) {
        int partition = epart[i];
        buckets[partition].push_back(i);
    }

    // sanity check of buckets
    for (size_t i = 0; i < buckets.size(); i++)
        ZIRAN_ASSERT(buckets[i].size() != 0, "partition failed. some partition is empty.");

    // build particle_color
    particle_color.clear();
    for (auto p : old_to_new)
        particle_color[p.key] = npart[p.value];
}

// partition_offsets stores n+1 indices, where n is the number of buckets
template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::reindexElements(const StdVector<StdVector<int>>& buckets)
{
    ZIRAN_QUIET_TIMER();
    partition_offsets.clear();
    StdVector<int> new_to_old;
    partition_offsets.push_back(0);
    for (size_t i = 0; i < buckets.size(); ++i) {
        for (size_t j = 0; j < buckets[i].size(); ++j)
            new_to_old.push_back(buckets[i][j]);
        partition_offsets.push_back(partition_offsets.back() + buckets[i].size());
    }

    reorder(new_to_old);
}

template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::writeData(std::ostream& out) const
{
    Base::writeData(out);
    writeSTDVector(out, partition_offsets);
    if (partition_offsets.size() > 1) {
        writeEntry(out, (uint64_t)neighbors_of_my_color.size());
        for (auto& i : neighbors_of_my_color) {
            writeSTDVector(out, i.neighbor_ids);
            writeEntry(out, (uint64_t)i.local_particle_ids.size());
            for (auto& k : i.local_particle_ids) {
                writeSTDVector(out, k);
            }
        }
    }
}

template <class T, int _manifold_dim, int _dim>
void SimplexElements<T, _manifold_dim, _dim, std::enable_if_t<_manifold_dim <= _dim>>::readData(std::istream& in)
{
    Base::readData(in);
    readSTDVector(in, partition_offsets);
    if (partition_offsets.size() > (size_t)0) {
        ZIRAN_ASSERT(partition_offsets.size() > (size_t)1);
        resizePadsAndBuildHelpers();

        uint64_t s = readEntry<uint64_t>(in);

        neighbors_of_my_color.resize(s);
        for (uint64_t i = 0; i < s; i++) {
            readSTDVector(in, neighbors_of_my_color[i].neighbor_ids);
            uint64_t z = readEntry<uint64_t>(in);
            neighbors_of_my_color[i].local_particle_ids.resize(z);
            for (uint64_t k = 0; k < z; k++) {
                readSTDVector(in, neighbors_of_my_color[i].local_particle_ids[k]);
            }
        }
    }
}

template class SimplexElements<double, 1, 2, void>;
template class SimplexElements<double, 1, 3, void>;
template class SimplexElements<double, 2, 2, void>;
template class SimplexElements<double, 2, 3, void>;
template class SimplexElements<double, 3, 3, void>;
template class SimplexElements<float, 1, 2, void>;
template class SimplexElements<float, 1, 3, void>;
template class SimplexElements<float, 2, 2, void>;
template class SimplexElements<float, 2, 3, void>;
template class SimplexElements<float, 3, 3, void>;
} // namespace ZIRAN
