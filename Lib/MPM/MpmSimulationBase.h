#ifndef MPM_SIMULATION_BASE_H
#define MPM_SIMULATION_BASE_H

#include <Ziran/Math/Nonlinear/NewtonsMethod.h>
#include <Ziran/Sim/BackwardEuler.h>
#include <Ziran/Sim/Scene.h>

#include <MPM/Force/MpmForceBase.h>
#include <MPM/MpmGrid.h>

namespace ZIRAN {
template <class T, int dim>
class AnalyticCollisionObject;
template <class T, int dim>
class SourceCollisionObject;
template <class T, int dim>
struct CollisionNode;
class PlasticityApplierBase;

template <class T, int _dim>
class MpmSimulationBase : public SimulationBase {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef T Scalar;
    static const int dim = _dim;
    static constexpr int interpolation_degree = ZIRAN_MPM_DEGREE;
    using Base = SimulationBase;
    typedef Vector<T, dim> TV;
    typedef Vector<int, dim> IV;
    using TVStack = Matrix<T, dim, Eigen::Dynamic>;
    typedef Vector<T, 4> TV4;
    typedef Matrix<T, dim, dim> TM;
    typedef Matrix<T, 4, 4> TM4;
    typedef tbb::concurrent_vector<CollisionNode<T, dim>> CollisionNodeArray;
    using Simulation = MpmSimulationBase<T, dim>;
    using Objective = BackwardEulerLagrangianForceObjective<Simulation>;
    using Base::diff_test;
    typedef typename MpmGrid<T, dim>::SparseGrid SparseGrid;
    typedef typename MpmGrid<T, dim>::SparseMask SparseMask;

    enum DebugEvents : int {
        BeforeNewtonSolve,
        AfterNewtonSolve
    };

    enum TransferScheme {
        APIC_blend_RPIC,
        FLIP_blend_PIC,
        OTHER
    };
    TransferScheme transfer_scheme;

    T apic_rpic_ratio;
    T flip_pic_ratio;
    T flip_rpic_ratio;

    bool print_stats;
    bool write_partio;
    bool write_meshes;
    bool symplectic;
    bool quasistatic;
    bool autorestart;
    bool ignoreCollisionObject;
    bool useTrialCollision;

    bool mls_mpm;

    TV gravity;
    T cfl;

    T D_inverse;
    T dx;
    T dt;

    int element_partitions;

    StdVector<std::function<void()>> before_euler_callbacks;
    StdVector<std::function<void()>> before_p2g_callbacks;

    // scene
    Scene<T, dim> scene;
    Particles<T, dim>& particles;

    // particle data
    DataArray<T>& element_measure;

    // For Lagrangian Forces
    StdVector<TV> scratch_xp;
    TVStack scratch_vp;
    TVStack scratch_fp;

    StdVector<TM> scratch_stress;
    StdVector<TM> scratch_gradV;

    // collision
    StdVector<std::unique_ptr<AnalyticCollisionObject<T, dim>>> collision_objects;
    CollisionNodeArray collision_nodes;

    // external body force field: fext(x,t)
    StdVector<std::function<TV(const TV&, const T)>> fext;

    // TODO: be more elegant?
    // vector<bool> have concurrent writes issue
    StdVector<int> explicit_collision_nodes;

    // implicit related
    Objective objective;
    NewtonsMethod<Objective> newton;
    Vector<T, Eigen::Dynamic> mass_matrix;
    TVStack dv; // Newton is formed from g(dv)=0
    TVStack vn;

    std::unique_ptr<LagrangianForce<T, dim>> inertia;

    StdVector<std::unique_ptr<MpmForceBase<T, dim>>> forces;
    std::unique_ptr<MpmForceBase<T, dim>>& force;

    //For plasticity
    StdVector<std::unique_ptr<PlasticityApplierBase>> plasticity_appliers;

    std::function<void(T&, const TV&, TV&)> explicit_velocity_field;

    // check spgrid
    StdVector<uint64_t> particle_sorter;
    StdVector<uint64_t> particle_base_offset;
    StdVector<int> particle_order;
    std::vector<std::pair<int, int>> particle_group;
    std::vector<uint64_t> block_offset;
    MpmGrid<T, dim> grid;
    int num_nodes = 0;

public:
    MpmSimulationBase();

    virtual ~MpmSimulationBase() override;

    virtual void setDx(const T new_dx);

    virtual void addCollisionObject(AnalyticCollisionObject<T, dim>&& object);

    virtual int addSourceCollisionObject(SourceCollisionObject<T, dim>&& object);

    virtual void initialize() override;

    virtual void reinitialize() override;

    virtual void trialCollision();

    virtual void beginTimeStep(int frame, int substep, double time, double dt) override;

    virtual void advanceOneTimeStep(double dt) override;

    virtual void buildAndPassGradVnToMpmForceHelpers();

    virtual void buildGradVn(StdVector<TM>& grad_vn);

    template <bool USE_APIC_BLEND_RPIC, bool USE_MPM_DEGREE_ONE, bool USE_MLS_MPM>
    void particlesToGridWithForceHelper();

    template <bool USE_APIC_BLEND_RPIC, bool USE_MPM_DEGREE_ONE>
    void particlesToGridHelper();

    virtual void performRpicVelocityFilteringAtBeginningOfTimeStep();

    virtual void particlesToGrid();

    virtual void startBackwardEuler();

    virtual void backwardEulerStep();

    T st_tau = 0;

    void addSurfaceTension();

    // Only called by symplecticEulerStep
    virtual void forcesUpdateParticleState();

    virtual T totalEnergy(std::string info);

    virtual void moveNodes(const TVStack& dv_in);

    /**
       Callback to write the simulation state

       Should be overridden to write the simulation state
    */
    virtual void writeState(std::ostream& out) override;

    /**
       Callback to read the simulation state

       Should be overridden to read the output of writeState
    */
    virtual void readState(std::istream& in) override;

    /**
       Callback to calculate dt based on CFL
    */
    virtual double calculateDt() override;

    // Build and mass_matrix (for implicit)
    virtual void buildMassMatrix();

    virtual void addScaledForces(const T scale, TVStack& f); // called BE Newton objective (for computing residual)

    virtual void addScaledForceDifferentials(const T scale, const TVStack& x, TVStack& f); // called BE Newton objective (for computing residual)

    virtual void gridVelocityExplicitUpdate(double dt);

    virtual void constructNewVelocityFromNewtonResult();

    template <bool USE_APIC_BLEND_RPIC, bool USE_MPM_DEGREE_ONE, bool USE_MLS_MPM>
    void gridToParticlesHelper(double dt);

    virtual void gridToParticles(double dt);

    virtual void applyPlasticity();

    virtual void sortParticlesAndPolluteGrid();

    virtual void buildInitialDvAndVnForNewton();

    /**
       result(0) is purely particle speed maximum.
       result(1) is APIC enhanced 'speed' maximum.
    */
    Eigen::Array<T, 2, 1> evalMaxParticleSpeed(TV& min_corner, TV& max_corner);

    virtual void writeSimulationInformation();

    // Called by writeSimulationInformation
    virtual void writeParticlePerCellHistogram();

    // Accessors. Mainly used in MpmInitializationBase to pass as arguments for MpmParticleHandleBase.

    inline Particles<T, dim>& getParticles() { return particles; }

    inline Scene<T, dim>& getScene() { return scene; }

    virtual MpmForceBase<T, dim>* getMpmForce();

    inline StdVector<std::unique_ptr<PlasticityApplierBase>>& getPlasticityAppliers() { return plasticity_appliers; }

    inline StdVector<TV>& getScratchXp() { return scratch_xp; }

    inline T& getDt() { return dt; }

    bool useDouble() override { return std::is_same<T, double>::value; }
    int dimension() override { return dim; }
    const char* name() override { return "mpm_base"; }

    template <typename OP>
    inline void parallel_for_updating_grid(const OP& target)
    {
        for (uint64_t color = 0; color < (1 << dim); ++color) {
            tbb::parallel_for(0, (int)particle_group.size(), [&](int group_idx) {
                if ((block_offset[group_idx] & ((1 << dim) - 1)) != color)
                    return;
                for (int idx = particle_group[group_idx].first; idx <= particle_group[group_idx].second; ++idx) {
                    int i = particle_order[idx];
                    target(i);
                }
            });
        }
    }
};
} // namespace ZIRAN
#endif
