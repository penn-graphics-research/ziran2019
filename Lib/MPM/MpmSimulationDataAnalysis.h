#ifndef MPM_SIMULATION_DATA_ANALYSIS_H
#define MPM_SIMULATION_DATA_ANALYSIS_H

#include "MpmSimulationBase.h"

#include <Ziran/Math/Geometry/Rotation.h>
#include <float.h>
#include <mutex>

namespace ZIRAN {

template <class T, int dim>
class MpmSimulationDataAnalysis {

public:
    static constexpr int interpolation_degree = ZIRAN_MPM_DEGREE;

    typedef Vector<T, dim> TV;
    typedef Vector<double, dim> DoubleV;
    typedef Matrix<double, dim, dim> DoubleM;
    typedef Matrix<T, dim, dim> TM;

    MpmSimulationBase<T, dim>& mpm;

public:
    explicit MpmSimulationDataAnalysis(MpmSimulationBase<T, dim>& mpm_in)
        : mpm(mpm_in)
    {
    }

    ~MpmSimulationDataAnalysis() {}

    /** result(0) is purely particle speed maximum.
        result(1) is APIC enhanced 'speed' maximum. */
    // remained in MpmSimulation class because needed by calculateDt()
    Eigen::Array<T, 2, 1> evalMaxParticleSpeed()
    {
        Eigen::Array<T, 2, 1> init((T)0, (T)0);
        if (mpm.transfer_scheme == MpmSimulationBase<T, dim>::APIC_blend_RPIC) {
            T scale = (6 * std::sqrt((T)dim)) / (mpm.dx * mpm.D_inverse);
            return mpm.particles.map_reduce(
                init,
                [scale](const TV& velocity, const TM& C) {
                    T velocity_norm = velocity.norm();
                    return Eigen::Array<T, 2, 1>(velocity_norm,
                        velocity_norm + C.norm() * scale);
                },
                [](const Eigen::Array<T, 2, 1>& a, const Eigen::Array<T, 2, 1>& b) {
                    return a.max(b);
                },
                mpm.particles.V_name(), mpm.C_name());
        }
        else
            return mpm.particles.map_reduce(
                init,
                [](const TV& velocity) {
                    T velocity_norm = velocity.norm();
                    return Vector<T, 2>(velocity_norm, velocity_norm);
                },
                [](const Eigen::Array<T, 2, 1>& a, const Eigen::Array<T, 2, 1>& b) {
                    return a.max(b);
                },
                mpm.particles.V_name());
    }

    /** Evaluate total linear momentum from particles */
    TV evalTotalLinearMomentumParticles() // summing mass*velocity over particles
    {
        DoubleV zero = DoubleV::Zero();

        DoubleV p = mpm.particles.map_reduce(
            zero,
            [](T m, const TV& v) {
                return double(m) * v.template cast<double>();
            },
            [](const DoubleV& a, const DoubleV& b) {
                return a + b;
            },
            mpm.particles.mass_name(), mpm.particles.V_name());
        return p.template cast<T>();
    }

    /** Evaluate total angular momentum from particles */
    AngularVelocity<T, dim> evalTotalAngularMomentumParticles() // summing mass*(velocity `cross_product` position) over particles
    {
        const AngularVelocity<double, dim> init;
        AngularVelocity<double, dim> result;

        if (mpm.transfer_scheme == MpmSimulationBase<T, dim>::APIC_blend_RPIC) {
            double D = (double)1 / double(mpm.D_inverse);
            result = mpm.particles.map_reduce(
                init,
                [D](const T& mass, const TV& position, const TV& velocity, const TM& C) {
                    DoubleM tempDoubleM = (D * C.transpose().template cast<double>()) + (position.template cast<double>() * velocity.transpose().template cast<double>());
                    return LeviCivitaContract(tempDoubleM) * double(mass);
                },
                [](const AngularVelocity<double, dim>& a, const AngularVelocity<double, dim>& b) { return a + b; },
                mpm.particles.mass_name(), mpm.particles.X_name(), mpm.particles.V_name(), C_name<TM>());
        }
        else {
            result = mpm.particles.map_reduce(
                init,
                [](const T& mass, const TV& position, const TV& velocity) {
                    DoubleM tempDoubleM = position.template cast<double>() * velocity.transpose().template cast<double>();
                    return LeviCivitaContract(tempDoubleM) * double(mass);
                },
                [](const AngularVelocity<double, dim>& a, const AngularVelocity<double, dim>& b) { return a + b; },
                mpm.particles.mass_name(), mpm.particles.X_name(), mpm.particles.V_name());
        }

        return result.template cast<T>();
    }

    /** Evaluate total kinetic energy from particles */
    T evalTotalKineticEnergyParticles() // 0.5*(summing mass*velocity.squaredNorm() over particles)
    {
        if (mpm.particles.get(mpm.particles.mass_name()).array.size() == 0)
            return (T)0;
        auto ranges = mpm.particles.get(mpm.particles.mass_name()).ranges;
        return (T)(0.5 * tbb::parallel_reduce(ranges, (double)0, [&](const DisjointRanges& subset, double init) -> double {
            for (auto iter = mpm.particles.subsetIter(subset, mpm.particles.mass_name(), mpm.particles.V_name()); iter; ++iter) {
                auto& mass = iter.template get<0>();
                auto& velocity = iter.template get<1>().template cast<double>();
                init += (double)mass * velocity.squaredNorm();
            }
            return init; }, [](double a, double b) -> double { return a + b; }));
    }

    T evalTotalMassParticles()
    {
        return mpm.particles.map_reduce(
            (T)0,
            [](T m) -> T { return m; },
            [](const T& a, const T& b) -> T {
                return a + b;
            },
            mpm.particles.mass_name());
    }

    T evalTotalMassGrid()
    {
        std::mutex mtx;
        double result = 0;
        mpm.grid.iterateGrid([&](Vector<int, dim> node, GridState<T, dim>& g) {
            mtx.lock();
            result += (double)g.m;
            mtx.unlock();
        });
        return (T)result;
    }

    T evalMaxSpeedGrid()
    {
        std::mutex mtx;
        double result = 0;
        mpm.grid.iterateGrid([&](Vector<int, dim> node, GridState<T, dim>& g) {
            mtx.lock();
            result = std::max(result, g.v.template cast<T>().length());
            mtx.unlock();
        });
        return (double)result;
    }

    Vector<T, dim> evalTotalLinearMomentumGrid(bool new_v = false)
    {
        std::mutex mtx;
        Vector<double, dim> result = Vector<double, dim>::Zero();
        mpm.grid.iterateGrid([&](Vector<int, dim> node, GridState<T, dim>& g) {
            Vector<double, dim> momentum = Vector<double, dim>::Zero();
            if (new_v) {
                momentum = (double)g.m * g.new_v.template cast<double>();
            }
            else {
                momentum = (double)g.m * g.v.template cast<double>();
            }
            mtx.lock();
            result += momentum;
            mtx.unlock();
        });
        return result.template cast<T>();
    }

    AngularVelocity<T, dim> evalTotalAngularMomentumGrid(bool new_v = false)
    {
        std::mutex mtx;
        AngularVelocity<double, dim> result;
        result.setZero();
        mpm.grid.iterateGrid([&](Vector<int, dim> node, GridState<T, dim>& g) {
            Vector<double, dim> x = node.template cast<double>() * (double)mpm.dx;
            double m = (double)g.m;
            Vector<double, dim> v = Vector<double, dim>::Zero();
            if (new_v)
                v = g.new_v.template cast<double>();
            else
                v = g.v.template cast<double>();
            DoubleM tempDoubleM = x * v.transpose();
            AngularVelocity<double, dim> momentum = LeviCivitaContract(tempDoubleM) * m;
            mtx.lock();
            result += momentum;
            mtx.unlock();
        });
        return result.template cast<T>();
    }

    AngularVelocity<T, dim> evalIncompressibleTotalAngularMomentumGrid(bool new_v = false)
    {
        std::mutex mtx;
        AngularVelocity<double, dim> result = AngularVelocity<double, dim>::Zero();
        mpm.grid.iterateGrid([&](Vector<int, dim> node, GridState<T, dim>& g) {
            Vector<double, dim> x = node.template cast<double>() * (double)mpm.dx;
            double m = (double)g.m;
            Vector<double, dim> v = Vector<double, dim>::Zero();
            if (new_v)
                v = g.new_v.template cast<double>();
            else
                v = g.v.template cast<double>();
            AngularVelocity<double, dim> momentum = LeviCivitaContract(x * v.transpose()) * 1; // assume constant mass
            mtx.lock();
            result += momentum;
            mtx.unlock();
        });
        return result.template cast<T>();
    }

    double evalTotalKineticEnergyGrid(bool new_v = false)
    {
        std::mutex mtx;
        double result = 0;
        mpm.grid.iterateGrid([&](Vector<int, dim> node, GridState<T, dim>& g) {
            double m = (double)g.m;
            Vector<double, dim> v = Vector<double, dim>::Zero();
            if (new_v)
                v = g.new_v.template cast<double>();
            else
                v = g.v.template cast<double>();
            double energy = 0.5 * m * v.squaredNorm();
            mtx.lock();
            result += energy;
            mtx.unlock();
        });
        return result;
    }

    T evalMaxMassParticles()
    {
        return mpm.particles.map_reduce(
            (T)0,
            [](T m) -> T { return m; },
            [](const T& a, const T& b) -> T {
                return std::max(a, b);
            },
            mpm.particles.mass_name());
    }

    T evalMinMassParticles()
    {
        return mpm.particles.map_reduce(
            (T)FLT_MAX,
            [](T x) -> T { return x; },
            [](const T& a, const T& b) -> T {
                return std::min(a, b);
            },
            mpm.particles.mass_name());
    }
};
} // namespace ZIRAN
#endif
