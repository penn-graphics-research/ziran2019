#ifndef MPM_INIT_3D_H
#define MPM_INIT_3D_H

#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/Math/Geometry/MeshConstruction.h>
#include <Ziran/Physics/SoundSpeedCfl.h>
#include <Ziran/Sim/MeshHandle.h>
#include <Ziran/Sim/SceneInitializationCore.h>
#include <Ziran/Math/Geometry/AnalyticLevelSet.h>
#include "MpmSimulation.h"

namespace Minchen {
extern double p_Jp, p_beta, p_nu, p_E, p_vy;
extern double p_Jp2, p_beta2, p_nu2, p_E2, p_rho2;
extern bool p_noPlasticity;

extern double v_mu;
extern double v_xxxx;

extern int lin_solver_type;
} // namespace Minchen

namespace ZIRAN {

template <class T, int dim>
class MpmInitBase;

template <class T>
class MpmInit3D : public MpmInitBase<T, 3> {
public:
    static const int dim = 3;
    using Base = MpmInitBase<T, dim>;
    using TV2 = Vector<T, 2>;
    using TVI2 = Vector<int, 2>;
    using TV = Vector<T, dim>;
    using TM = Eigen::Matrix<T, dim, dim>;
    using TVI = Vector<int, dim>;

    using Base::init_helper;
    using Base::scene;
    using Base::sim;
    using Base::test_number;

    MpmInit3D(MpmSimulation<T, dim>& sim, const int test_number)
        : Base(sim, test_number)
    {
    }

    void reload() override
    {
        // Watermelon smash
        // ./mpm -test 1 --3d -Jp 0.99 -Jp2 0.97 -beta 2 -beta2 1 -nu 0.4 -nu2 0.35 -E 2000 -E2 1000 -rho 2 -vy -6.5
        if (test_number == 1) {
            //            std::string Jp_str = "Jp_" + std::to_string(Minchen::p_Jp);
            //            Jp_str.erase(Jp_str.find_last_not_of('0') + 1, std::string::npos);
            //            Jp_str += "_" + std::to_string(Minchen::p_Jp2);
            //            Jp_str.erase(Jp_str.find_last_not_of('0') + 1, std::string::npos);
            //
            //            std::string beta_str = "beta_" + std::to_string(Minchen::p_beta);
            //            beta_str.erase(beta_str.find_last_not_of('0') + 1, std::string::npos);
            //            beta_str += "_" + std::to_string(Minchen::p_beta2);
            //            beta_str.erase(beta_str.find_last_not_of('0') + 1, std::string::npos);
            //
            //            std::string nu_str = "nu_" + std::to_string(Minchen::p_nu);
            //            nu_str.erase(nu_str.find_last_not_of('0') + 1, std::string::npos);
            //            nu_str += "_" + std::to_string(Minchen::p_nu2);
            //            nu_str.erase(nu_str.find_last_not_of('0') + 1, std::string::npos);
            //
            //            std::string E_str = "E_" + std::to_string(Minchen::p_E);
            //            E_str.erase(E_str.find_last_not_of('0') + 1, std::string::npos);
            //            E_str += "_" + std::to_string(Minchen::p_E2);
            //            E_str.erase(E_str.find_last_not_of('0') + 1, std::string::npos);
            //
            //            std::string rho_str = "rho_" + std::to_string(Minchen::p_rho2);
            //            rho_str.erase(rho_str.find_last_not_of('0') + 1, std::string::npos);
            //
            //            std::string vy_str = "vy_" + std::to_string(Minchen::p_vy);
            //            vy_str.erase(vy_str.find_last_not_of('0') + 1, std::string::npos);
            //
            //            std::cout << "Jp = " << Minchen::p_Jp << ", beta = " << Minchen::p_beta << std::endl;
            //            std::cout << "nu = " << Minchen::p_nu << ", E = " << Minchen::p_E << std::endl;
            //            std::cout << "Jp2 = " << Minchen::p_Jp2 << ", beta2 = " << Minchen::p_beta2 << std::endl;
            //            std::cout << "nu2 = " << Minchen::p_nu2 << ", E2 = " << Minchen::p_E2 << std::endl;
            //            std::cout << "rho2 = " << Minchen::p_rho2 << std::endl;
            //            std::cout << "vy = " << Minchen::p_vy << std::endl;
            //
            //            std::string path("output/watermelon_smash_3d_" + Jp_str + "_" + beta_str + "_" + nu_str + "_" + E_str + "_" + rho_str + "_" + vy_str);
            //            if (Minchen::p_noPlasticity) {
            //                std::cout << "pure elasticity" << std::endl;
            //                path = "output/watermelon_smash_3d_pureElasticity";
            //            }
            //
            //            sim.output_dir.path = path.c_str();
            //
            //            if (!Minchen::p_noPlasticity) {
            //                // extra attribute dump (logJp)
            //                sim.end_frame_callbacks.emplace_back(
            //                    [&](int frame) {
            //                        std::string filename = sim.output_dir.absolutePath(sim.outputFileName("logJp", ".bgeo"));
            //                        writePartioAttribute(filename, sim.particles, AttributeName<NonAssociativeCamClay<T>>(NonAssociativeCamClay<T>::name()), (T)-1);
            //                    });
            //            }

            T Jp = 0.99;
            T Jp2 = 0.97;
            T beta = 2;
            T beta2 = 1;
            T nu = 0.4;
            T nu2 = 0.35;
            T E = 2000;
            T E2 = 1000;
            //vy = -6.5;

            sim.output_dir.path = "output/timingWatermelon";

            sim.end_frame = 80;
            sim.dx = 0.007; //2.7 million particles with dx = 0.004 * 1.7
            sim.gravity = -3 * TV::Unit(1);
            sim.symplectic = true;
            // sim.transfer_scheme = MpmSimulationBase<T, dim>::APIC_blend_RPIC;
            // sim.apic_rpic_ratio = 0;
            sim.transfer_scheme = MpmSimulationBase<T, dim>::FLIP_blend_PIC;
            sim.flip_pic_ratio = 0.99;
            sim.objective.matrix_free = true;
            sim.verbose = true;
            sim.cfl = 0.6;
            init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx, AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for SPGrid.

            //Add watermelon shell modeled with pumpkin parameters and add watermelon innards modeled with cam clay snow params
            T rho = 2;
            T suggested_dt = evaluateTimestepLinearElasticityAnalysis(E, nu, rho, sim.dx, sim.cfl);
            if (sim.symplectic) { sim.step.max_dt = suggested_dt * 0.6; }

            int particlesPerCell = 8;
            MpmParticleHandleBase<T, dim> particles_handle_outside = init_helper.sampleFromVdbFile("LevelSets/watermelon_outside.vdb", rho, particlesPerCell);

            auto initial_translation = [](int index, Ref<T> mass, TV& X, TV& V) {
                X = X + TV(4, 1.5, 4);
                V = TV(0, -6.5, 0);
            };
            particles_handle_outside.transform(initial_translation);

            NeoHookeanBorden<T, dim> modelOutside(E, nu);
            NonAssociativeCamClay<T> pOutside(/*logJp*/ std::log(Jp), /*friction_angle*/ 45, /*beta*/ beta, /*xi*/ 3, dim); // test 19 pumpkin params
            particles_handle_outside.addFBasedMpmForce(modelOutside); //.97, 45, .8, 3
            particles_handle_outside.addPlasticity(modelOutside, pOutside, "F");

            ZIRAN_INFO("peel particle #: ", sim.particles.count);

            MpmParticleHandleBase<T, dim> particles_handle_inside = init_helper.sampleFromVdbFile("LevelSets/watermelon_inside.vdb", Minchen::p_rho2, particlesPerCell);
            particles_handle_inside.transform(initial_translation);

            NeoHookeanBorden<T, dim> modelInside(E2, nu2);
            NonAssociativeCamClay<T> pInside(/*logJp*/ std::log(Jp2), /*friction_angle*/ 45, /*beta*/ beta2, /*xi*/ 3, dim); // test 9 cam clay snow params (M = 0.6 --> phi = 15.83 deg)
            particles_handle_inside.addFBasedMpmForce(modelInside);
            particles_handle_inside.addPlasticity(modelInside, pInside, "F");

            //Add ground
            T groundFriction = 0.2;
            TV ground_origin(0, .2, 0);
            TV ground_normal(0, 1, 0);
            HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::SLIP);
            ground_object.setFriction(groundFriction);
            init_helper.addAnalyticCollisionObject(ground_object); // add ground
        }

        // Pumpkin crash --> One pumpkin smashing down on another pumpkin
        // ./mpm -test 2 --3d -Jp 0.96 -beta 2 -nu 0.39 -E 2000
        if (test_number == 2) {

            T Jp = 0.96;
            T beta = 2;
            T nu = 0.39;
            T E = 2000;

            //            std::string Jp_str = std::to_string(Minchen::p_Jp);
            //            Jp_str.erase(Jp_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::string beta_str = std::to_string(Minchen::p_beta);
            //            beta_str.erase(beta_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::string nu_str = std::to_string(Minchen::p_nu);
            //            nu_str.erase(nu_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::string E_str = std::to_string(Minchen::p_E);
            //            E_str.erase(E_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::cout << "Jp = " << Minchen::p_Jp << ", beta = " << Minchen::p_beta << std::endl;
            //            std::cout << "nu = " << Minchen::p_nu << ", E = " << Minchen::p_E << std::endl;
            //            std::string path("output/pumpkin_smash_3d_" + Jp_str + "_" + beta_str + "_" + nu_str + "_" + E_str);
            //            if (Minchen::p_noPlasticity) {
            //                std::cout << "pure elasticity" << std::endl;
            //                path = "output/pumpkin_smash_3d_pureElasticity";
            //            }

            // extra attribute dump (logJp)
            //            if (!Minchen::p_noPlasticity) {
            //                sim.end_frame_callbacks.emplace_back(
            //                    [&](int frame) {
            //                        std::string filename = sim.output_dir.absolutePath(sim.outputFileName("logJp", ".bgeo"));
            //                        writePartioAttribute(filename, sim.particles, AttributeName<NonAssociativeCamClay<T>>(NonAssociativeCamClay<T>::name()), (T)-1);
            //                    });
            //            }

            //sim.output_dir.path = path.c_str();

            sim.output_dir.path = "output/timingPumpkin";
            sim.end_frame = 80;
            sim.dx = 0.007;
            //            sim.dx = 0.006; //2.4 million particles
            sim.gravity = -3 * TV::Unit(1);

            sim.symplectic = true;
            // sim.transfer_scheme = MpmSimulationBase<T, dim>::APIC_blend_RPIC;
            // sim.apic_rpic_ratio = 0;
            sim.transfer_scheme = MpmSimulationBase<T, dim>::FLIP_blend_PIC;
            sim.flip_pic_ratio = 0.98;
            sim.objective.matrix_free = true;
            sim.verbose = false;
            sim.cfl = 0.6;
            init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx, AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for SPGrid.

            T rho = 2;
            T suggested_dt = evaluateTimestepLinearElasticityAnalysis(E, nu, rho, sim.dx, sim.cfl);

            if (sim.symplectic) { sim.step.max_dt = suggested_dt * 0.6; }

            int particlesPerCell = 8;
            MpmParticleHandleBase<T, dim> particles_handle_pumpkin1 = init_helper.sampleFromVdbFile("LevelSets/pumpkinShiftedUp.vdb", rho, particlesPerCell);
            MpmParticleHandleBase<T, dim> particles_handle_pumpkin2 = init_helper.sampleFromVdbFile("LevelSets/pumpkinShiftedUp.vdb", rho, particlesPerCell);

            auto initial_translation_pumpkin1 = [](int index, Ref<T> mass, TV& X, TV& V) {
                X[1] *= 2;
                X = X + TV(4, 0.2, 4);
                V = TV(0, 0, 0);
            };
            particles_handle_pumpkin1.transform(initial_translation_pumpkin1);

            auto initial_translation_pumpkin2 = [](int index, Ref<T> mass, TV& X, TV& V) {
                Rotation<T, dim> rotation(TV(0, 1, 0), TV(-1, 1, 0));
                X = rotation.rotation * X + TV(4.05, 1.8, 4);
                V = TV(-0.1, -5.6, 0);
            };
            particles_handle_pumpkin2.transform(initial_translation_pumpkin2);

            NeoHookeanBorden<T, dim> model(Minchen::p_E, Minchen::p_nu);
            NonAssociativeCamClay<T> p(/*logJp*/ std::log(Jp), /*friction_angle*/ 45, /*beta*/ beta, /*xi*/ 3, dim); // .97, 45, .8, 3
            particles_handle_pumpkin1.addFBasedMpmForce(model);
            if (!Minchen::p_noPlasticity) {
                particles_handle_pumpkin1.addPlasticity(model, p, "F");
            }
            particles_handle_pumpkin2.addFBasedMpmForce(model);
            if (!Minchen::p_noPlasticity) {
                particles_handle_pumpkin2.addPlasticity(model, p, "F");
            }

            //Add ground, and a wall to prevent negative position in each direction (x and z)!
            T wallFriction = 0.1;
            T groundFriction = 0.4;

            TV ground_origin(0.0, 0.2, 0.0);
            TV ground_normal(0, 1, 0);
            HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::SLIP);
            ground_object.setFriction(groundFriction);
            init_helper.addAnalyticCollisionObject(ground_object); // add ground
        }

        // drop round cookie
        // ./mpm -test 3 --3d
        if (test_number == 3) {
            sim.output_dir.path = "output/high_oreo_with_heart";
            sim.end_frame = 240;
            sim.dx = 0.004;
            sim.gravity = -3 * TV::Unit(1);

            sim.symplectic = true;
            // sim.transfer_scheme = MpmSimulationBase<T, dim>::APIC_blend_RPIC;
            // sim.apic_rpic_ratio = 0;
            sim.transfer_scheme = MpmSimulationBase<T, dim>::FLIP_blend_PIC;
            sim.flip_pic_ratio = 0.98;
            sim.objective.matrix_free = true;
            sim.verbose = false;
            sim.cfl = 0.6;
            //sim.step.frame_dt = sim.step.max_dt;
            init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx, AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for SPGrid.

            TV domain_center = 2048 * sim.dx * TV(1, 0, 1);
            TV initial_V = TV(0, -.2, 0);
            T rho = 2;
            T rho_heart = rho * .2;
            T nu = 0.35, youngs = 1000;
            T youngs_crust = youngs;
            T suggested_dt = evaluateTimestepLinearElasticityAnalysis(youngs_crust, nu, rho, sim.dx, sim.cfl);
            sim.step.max_dt = suggested_dt * 0.7;
            if (sim.symplectic) { ZIRAN_ASSERT(sim.step.max_dt <= suggested_dt, suggested_dt); }

            { // crust
                MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleFromVdbFile("LevelSets/oreo_crust.vdb", rho);
                auto initial_translation = [=](int index, Ref<T> mass, TV& X, TV& V) {
                    Rotation<T, dim> rotation(TV(0, 1, 0), TV(0.1, 1, 0.2));
                    X = rotation.rotation * X + TV(.0, .9, .0) + domain_center;
                    V = initial_V;
                };
                particles_handle.transform(initial_translation);
                NeoHookeanBorden<T, dim> model(youngs_crust, nu);
                NonAssociativeCamClay<T> p(/*logJp*/ std::log((T)1.01), /*friction_angle*/ 60, /*beta*/ 3, /*xi*/ 3, dim); // .97, 45, .8, 3
                particles_handle.addFBasedMpmForce(model);
                particles_handle.addPlasticity(model, p, "F");
            }

            ZIRAN_INFO("crust particle #: ", sim.particles.count);

            { // heart
                T radius = 0.2;
                T height = 0.02;
                T theta = M_PI / 2;
                Vector<T, 4> cylinder_rotation(std::cos(theta / 2), 0, 0, std::sin(theta / 2));
                TV translation(0, 0, 0);
                CappedCylinder<T, dim> cylinder(radius, height, cylinder_rotation, translation);
                MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleInAnalyticLevelSet(cylinder, rho_heart, 8);
                auto initial_translation = [=](int index, Ref<T> mass, TV& X, TV& V) {
                    Rotation<T, dim> rotation(TV(0, 1, 0), TV(0.1, 1, 0.2));
                    X = rotation.rotation * X + TV(.0, .9, .0) + domain_center;
                    V = initial_V;
                };
                particles_handle.transform(initial_translation);
                NeoHookeanBorden<T, dim> model(youngs, nu);
                NonAssociativeCamClay<T> p(/*logJp*/ std::log((T).99), /*friction_angle*/ 45, /*beta*/ .7, /*xi*/ 1, dim); // .97, 45, .8, 3
                particles_handle.addFBasedMpmForce(model);
                particles_handle.addPlasticity(model, p, "F");
            }

            TV ground_origin(0, 0.5, 0);
            TV ground_normal(0, 1, 0);
            HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::SLIP);
            ground_object.setFriction(.15);
            init_helper.addAnalyticCollisionObject(ground_object); // add ground
        }
    }
};
} // namespace ZIRAN
#endif
