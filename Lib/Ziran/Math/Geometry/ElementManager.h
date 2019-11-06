#ifndef ELEMENT_MANAGER_H
#define ELEMENT_MANAGER_H
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/CS/DataStructure/DataManager.h>
namespace ZIRAN {

//*****************************************************************************
// class ElementManager : public DataManager
//
// ElementManager is a date manager for finite elements.
// It contains a bunch of empty virtual functions.
// It also contains three scratch data structures for partition-based parallization for scatter operations on a mesh.
//*****************************************************************************

template <class T, int dim>
class ElementManager : public DataManager {
public:
    using TVStack = Matrix<T, dim, Eigen::Dynamic>;

    template <typename... Types>
    ElementManager(const AttributeName<Types>&... b)
        : DataManager(b...)
    {
    }

    //
    // These three data structures are for partition-based parallization for scatter operations on a mesh.
    //
    StdVector<TVStack> pads_dx; // they are just element partition scratch pads of particles for parallization. no write.
    StdVector<TVStack> pads; // they are just element partition scratch pads of particles for parallization. no write.
    StdVector<StdVector<int>> particle_l2g; // particle index local to global. no write.

    virtual void partitionAndReindexElements(const int n_partitions) = 0;

    virtual bool needRepartitioning() = 0;

    virtual void gatherFromPadsToTVStack(const T scale, TVStack& force) = 0;

    virtual bool renumberParticles(const int particle_count, StdVector<int>& old_to_new, StdVector<int>& new_to_old) = 0;

    virtual void registerParticleReorderingCallback(Particles<T, dim>& particles) = 0;
};
} // namespace ZIRAN
#endif
