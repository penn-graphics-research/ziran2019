#ifndef MPM_INIT_3D_H
#define MPM_INIT_3D_H

#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/Math/Geometry/MeshConstruction.h>
#include <Ziran/Math/Geometry/VoronoiNoise.h>
#include <Ziran/Physics/SoundSpeedCfl.h>
#include <Ziran/Sim/MeshHandle.h>
#include <Ziran/Sim/SceneInitializationCore.h>
#include "FractureSimulation.h"
#include "FractureInit.h"

namespace Minchen {
extern double nu, E, rho, v, dx, M, sigmaPercent, bulletFriction;
extern bool usePhaseField;
} // namespace Minchen

namespace ZIRAN {

template <class T>
class FractureInit3D : public FractureInitBase<T, 3> {
public:
    static const int dim = 3;
    using Base = FractureInitBase<T, dim>;
    using TV2 = Vector<T, 2>;
    using TVI2 = Vector<int, 2>;
    using TV = Vector<T, dim>;
    using TM = Eigen::Matrix<T, dim, dim>;
    using TVI = Vector<int, dim>;

    using Base::init_helper;
    using Base::scene;
    using Base::sim;
    using Base::test_number;

    FractureInit3D(FractureSimulation<T, dim>& sim, const int test_number)
        : Base(sim, test_number)
    {
    }

    void reload() override
    {

        // scripted bullet shooting a jello dinosaur
        // ./mpm -test 24 --3d --usePhaseField -nu 0.4 -E 20 -rho 0.1 -v 10 -dx 0.012 -M 15 -sigmaPercent 0.015 -bulletFriction 0.1
        // ./mpm -test 24 --3d --usePhaseField -nu 0.4 -E 20 -rho 0.1 -v 10 -dx 0.012 -M 10 -sigmaPercent 0.02 -bulletFriction 0.1
        // ./mpm -test 24 --3d -nu 0.4 -E 20 -rho 0.1 -v 10 -dx 0.012 -bulletFriction 0.1
        if (test_number == 1) {
            //            std::string name = "Dino";
            //
            //            std::string nu_str = "nu" + std::to_string(Minchen::nu);
            //            nu_str.erase(nu_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::string E_str = "E" + std::to_string(Minchen::E);
            //            E_str.erase(E_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::string rho_str = "rho" + std::to_string(Minchen::rho);
            //            rho_str.erase(rho_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::string v_str = "v" + std::to_string(Minchen::v);
            //            v_str.erase(v_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::string dx_str = "dx" + std::to_string(Minchen::dx);
            //            dx_str.erase(dx_str.find_last_not_of('0') + 1, std::string::npos);
            //            std::string bulletFriction_str = "bulletFric" + std::to_string(Minchen::bulletFriction);
            //            bulletFriction_str.erase(bulletFriction_str.find_last_not_of('0') + 1, std::string::npos);
            //            if (Minchen::usePhaseField) {
            //                std::string M_str = "M" + std::to_string(Minchen::M);
            //                M_str.erase(M_str.find_last_not_of('0') + 1, std::string::npos);
            //                std::string sigmaPercent_str = "sigma" + std::to_string(Minchen::sigmaPercent);
            //                sigmaPercent_str.erase(sigmaPercent_str.find_last_not_of('0') + 1, std::string::npos);
            //
            //                sim.output_dir.path = "output/bullet_jello" + name + "_pf_" + nu_str + "_" + E_str + "_" + rho_str + "_" + v_str + "_" + dx_str + "_" + bulletFriction_str + "_" + M_str + "_" + sigmaPercent_str;
            //            }
            //            else {
            //                sim.output_dir.path = "output/bullet_jello" + name + "_" + nu_str + "_" + E_str + "_" + rho_str + "_" + v_str + "_" + dx_str + "_" + bulletFriction_str;
            //            }

            T jelloE = 20;
            T jelloNu = 0.4;
            T jelloRho = 0.1;
            sim.use_phase_field = 1;
            T percent = 0.015;
            sim.parabolic_M = 15;
            sim.lumping = true;
            T bulletFriction = 0.1;

            sim.output_dir.path = "output/timingDinosaur2";

            sim.end_frame = 100; //100 for final demo
            sim.dx = 0.012;
            sim.gravity = -3 * TV::Unit(1);
            sim.step.max_dt = 1e-4;
            sim.step.frame_dt = (T)1 / 30;
            sim.symplectic = false;
            sim.transfer_scheme = MpmSimulationBase<T, dim>::APIC_blend_RPIC;
            sim.apic_rpic_ratio = 1;
            sim.objective.matrix_free = true;
            sim.verbose = false;
            sim.cfl = 0.8;

            sim.newton.tolerance = 1e-4;
            sim.newton.max_iterations = 200;
            sim.objective.minres.max_iterations = 20;
            sim.objective.minres.tolerance = 1e-5;

            T groundY = sim.dx * 3;

            //Setup the dino
            int particlesPerCellPane = 8;
            NeoHookeanBorden<T, dim> model(jelloE, jelloNu); //Elasticity Handler

            MpmParticleHandleBase<T, dim> particles_handle_top = init_helper.sampleFromVdbFile(
                "LevelSets/dino.vdb",
                jelloRho, particlesPerCellPane);
            auto initial_translation_kb = [&](int index, Ref<T> mass, TV& X, TV& V) {
                X = X + TV(5, groundY, 5.5);
            };
            particles_handle_top.transform(initial_translation_kb);

            particles_handle_top.addFBasedMpmForceWithPhaseField(percent, (T)sim.dx * 0.5, model);

#define COLLISION_OBJECT_BULLET
#ifdef COLLISION_OBJECT_BULLET
            //Setup the moving scripted bullet
            VdbLevelSet<T, dim> bullet_ls("LevelSets/bullet.vdb");
            auto bulletTransform = [](T time, AnalyticCollisionObject<T, dim>& object) {
                TV translation_velocity(10, 0,
                    0); //based on the partial derivatives of the parametric equations for this parabola
                //                TV translation(4 + Minchen::v * time, 0.45, 5.5); //multiply each velocity by dt to get dx!
                TV translation(4 + 10 * time, 0.975, 5.5); //multiply each velocity by dt to get dx!
                object.setTranslation(translation, translation_velocity);

                //                object.setScaling(0.5, 0);

                //                // rotation
                //                Vector<T, 3> rotAxis(0, 0, 1);
                //                T angle = -90;
                //
                //                T angle_rad = angle * M_PI / 180;
                //                T cosine = std::cos(angle_rad / 2), sine = std::sin(angle_rad / 2);
                //                Vector<T, 4> rotation(cosine, sine * rotAxis[0], sine * rotAxis[1], sine * rotAxis[2]);
                //                object.setRotation(rotation);
                Vector<T, 3> omega(0, 0, 0);
                object.setAngularVelocity(omega);
            };
            AnalyticCollisionObject<T, dim> bullet_object(bulletTransform, bullet_ls, AnalyticCollisionObject<T, dim>::STICKY);
            bullet_object.setFriction(bulletFriction);
            init_helper.addAnalyticCollisionObject(bullet_object);
#else
            //Setup the moving elastic bullet
            MpmParticleHandleBase<T, dim> particles_handle_bullet = init_helper.sampleFromVdbFile("/home/joshwolper/Desktop/Fracture Media/Models/bullet.vdb",
                bulletRho, particlesPerCellPane);
            auto initial_translation_bullet = [](int index, Ref<T> mass, TV& X, TV& V) {
                X = X + TV(4 + 10 * time, 0.975, 5.5);
                V = TV(10, 0, 0);
            };
            particles_handle_bullet.transform(initial_translation_bullet);
            NeoHookeanBorden<T, dim> bulletModel(bulletE, bulletNu);
            particles_handle_bullet.addFBasedMpmForceWithPhaseField(/* resistence to break */ 0.0075, (T)sim.dx * 0.5, bulletModel, false);
#endif

            //Set up a ground plane
            TV ground_origin(0, groundY, 0);
            TV ground_normal(0, 1, 0);
            HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::STICKY);
            ground_object.setFriction(1000.0);
            init_helper.addAnalyticCollisionObject(ground_object); // add ground

            //Set two walls
            TV xWall_origin(groundY, 0, 0);
            TV xWall_normal(1, 0, 0);
            HalfSpace<T, dim> xWall_ls(xWall_origin, xWall_normal);
            AnalyticCollisionObject<T, dim> xWall_object(xWall_ls, AnalyticCollisionObject<T, dim>::STICKY);
            xWall_object.setFriction(.2);
            init_helper.addAnalyticCollisionObject(xWall_object); // add ground

            TV zWall_origin(0, 0, groundY);
            TV zWall_normal(0, 0, 1);
            HalfSpace<T, dim> zWall_ls(zWall_origin, zWall_normal);
            AnalyticCollisionObject<T, dim> zWall_object(zWall_ls, AnalyticCollisionObject<T, dim>::STICKY);
            zWall_object.setFriction(.2);
            init_helper.addAnalyticCollisionObject(zWall_object); // add ground
        }

        // Tear bread
        // ./fracture -test 25 --3d
        if (test_number == 2) {
            sim.output_dir.path = "output/timingBread";
            sim.end_frame = 3600; //3600 for final demo
            sim.dx = 0.01 * 0.5 / 1.31 / 1.35;
            sim.gravity = -2 * TV::Unit(1);
            sim.step.max_dt = 0.00005 / 1.31;
            sim.step.frame_dt = (T)1 / (T)480;
            sim.symplectic = true; //explicit dynamics
            sim.transfer_scheme = MpmSimulationBase<T, dim>::FLIP_blend_PIC;
            sim.flip_pic_ratio = 0;
            sim.objective.matrix_free = true;
            sim.verbose = false;
            sim.cfl = 0.4;
            init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx, AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for SPGrid.

            //Solver params
            sim.newton.tolerance = 1e-3;
            sim.newton.max_iterations = 5;
            sim.objective.minres.max_iterations = 10000;
            sim.objective.minres.tolerance = 1e-4;

            sim.begin_time_step_callbacks.push_back(
                [this](int frame, int substep, double time, double dt) {
                    if (frame > 38) {
                        sim.performRpicVelocityFilteringAtBeginningOfTimeStep();
                        sim.performRpicVelocityFilteringAtBeginningOfTimeStep();
                        sim.performRpicVelocityFilteringAtBeginningOfTimeStep();
                    }
                });

            //Phase field params
            sim.use_phase_field = true;
            sim.parabolic_M = 1; // smaller M -> slower propagation
            sim.lumping = true;
            T pfPercent = 0.01;
            T l0 = (1 / (T)128) / (T)2;
            bool allowDamage = true;
            {
                sim.delete_particle_threshold = -1;
                Sphere<T, dim> ballLevelSet(TV(4, 4), sim.dx);
                AnalyticCollisionObject<T, dim> ground_object(ballLevelSet, AnalyticCollisionObject<T, dim>::STICKY);
                init_helper.addAnalyticCollisionObject(ground_object);
            }

            //Neo Hookean Borden Variables
            T rho = 2;
            T nu = 0.4, youngs = 500;

            int ppc = 8;
            MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleFromVdbFile("LevelSets/breadxxx.vdb", rho, ppc);
            TV domain_center = 2048 * sim.dx * TV(0, 0, 0);
            auto initial_translation = [=](int index, Ref<T> mass, TV& X, TV& V) {
                X = X + domain_center;
            };
            particles_handle.transform(initial_translation);

            NeoHookeanBorden<T, dim> model(youngs, nu); //neohookean borden elasticity
            //NonAssociativeCamClay<T> p(/*logJp*/ 0.99, /*friction_angle*/ 45, /*beta*/ 60, /*xi*/ 6, dim);
            //            NonAssociativeVonMises<T> p(/*tauY*/ 1000 * 1.0, /*alpha*/ 0, /*hardeningCoefficient*/ 0, dim); //von mises plasticity

            particles_handle.addFBasedMpmForceWithPhaseField(pfPercent, l0, model, allowDamage);
            //            particles_handle.addPlasticity(model, p, "F");

            for (auto iter = particles_handle.particles.subsetIter(DisjointRanges{ particles_handle.particle_range }, phase_field_name<PhaseField<T, dim>>()); iter; ++iter) {
                auto& pf = iter.template get<0>();
                auto& X = particles_handle.particles.X(iter.entryId());
                pf.residual_phase = 0.001 * 2;
            }

            //Set up collision objects
            // original sphere -> left, right
            T sphereRadius = 0.2;

            TV leftSphereCenter = TV(3.2, 3.8, 3.5);
            TV rightSphereCenter = TV(3.8, 3.8, 3.5);

            auto leftSphereTransform = [](T time, AnalyticCollisionObject<T, dim>& object) {
                TV translation_velocity(-0.125, 0., 0.);
                TV translation(-0.125 * time, 0, 0.);
                object.setTranslation(translation, translation_velocity);
            };
            auto rightSphereTransform = [](T time, AnalyticCollisionObject<T, dim>& object) {
                TV translation_velocity(0.125, 0., 0.);
                TV translation(0.125 * time, 0, 0.);
                object.setTranslation(translation, translation_velocity);
            };

            Sphere<T, dim> leftSphere(leftSphereCenter, sphereRadius);
            Sphere<T, dim> rightSphere(rightSphereCenter, sphereRadius);

            AnalyticCollisionObject<T, dim> leftObject(leftSphereTransform, leftSphere, AnalyticCollisionObject<T, dim>::STICKY);
            AnalyticCollisionObject<T, dim> rightObject(rightSphereTransform, rightSphere, AnalyticCollisionObject<T, dim>::STICKY);

            init_helper.addAnalyticCollisionObject(leftObject);
            init_helper.addAnalyticCollisionObject(rightObject);
        }

        //Armadillo stretch (full elastic but also trying metal after)
        if (test_number == 3) {
            //NOTE: These parameters are FINAL
            sim.output_dir.path = "output/timingArmadilloStretch";
            sim.end_frame = 360; //360 for final demo
            T frameRate = 120;
            sim.step.frame_dt = (T)1 / frameRate;
            sim.dx = 0.005; //use 0.005 for high res
            sim.gravity = -0.2 * TV::Unit(1);
            sim.step.max_dt = /*1e-3*/ 3e-4;
            sim.newton.max_iterations = 5;
            sim.newton.tolerance = 1e-3;
            sim.objective.minres.max_iterations = 10000;
            sim.objective.minres.tolerance = 1e-4;
            sim.quasistatic = false;
            sim.symplectic = true;
            sim.objective.matrix_free = true;
            sim.verbose = false;
            sim.cfl = 0.4;
            sim.transfer_scheme = MpmSimulationBase<T, dim>::APIC_blend_RPIC;
            sim.apic_rpic_ratio = 0; //full RPIC
            init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx, AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for SPGrid.

            //Elasticity Params
            T Youngs = 50;
            T nu = 0.4;
            T rho = 2;

            //Phase Field Params
            sim.use_phase_field = true;
            sim.parabolic_M = 1; //was -1 before to match physbam, changed to tune pure elasticity with PF
            sim.lumping = true;
            {
                sim.delete_particle_threshold = -1; //0.0075 looks decent, but lose detail at the fracture interface
                Sphere<T, dim> ballLevelSet(TV(10, 10), sim.dx);
                AnalyticCollisionObject<T, dim> ground_object(ballLevelSet, AnalyticCollisionObject<T, dim>::STICKY);
                init_helper.addAnalyticCollisionObject(ground_object);
            }
            T percentage = 0.01; //Fanfu's param is 0.01
            T l0 = (1 / (T)128) / (T)2;
            bool allowDamage = true;

            //Setup the armadillo material
            int ppc = 40;
            MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleFromVdbFile("LevelSets/armadillo.vdb", rho, ppc);

            //Elasticity Handler
            NeoHookeanBorden<T, dim> model(Youngs, nu);
            particles_handle.addFBasedMpmForceWithPhaseField(percentage, l0, model, allowDamage);

            //Sigma C Voronoi Initialization Routine
            if (sim.use_phase_field) {
                //std::string filename = "/home/joshwolper/Desktop/Fracture Media/Models/Voronoi Files/armadilloVoro.obj";
                //T radius = 0.1; //0.1 seems good
                //T magnitude = 10; //magnitude = 1 keeps the function continuous
                //T minP = -0.9999999; //should not go below -1
                //T maxP = 3000;
                //T zeta = 75;
                //particles_handle.voronoiSigmaC(filename, radius);
            }

            //Set residual phase
            for (auto iter = particles_handle.particles.subsetIter(DisjointRanges{ particles_handle.particle_range }, phase_field_name<PhaseField<T, dim>>()); iter; ++iter) {
                auto& pf = iter.template get<0>();
                auto& X = particles_handle.particles.X(iter.entryId());
                pf.residual_phase = 0.001;
            }

            //Now setup all the spheres (1 static, 4 moving)
            T radA = 0.08;
            T radB = 0.12;
            T radC = 0.12;
            T radD = 0.11;
            T radE = 0.14;
            TV centerA(2.01, 2.31, 1.94);
            TV centerB(1.72128, 2.39828, 2.15566);
            TV centerC(2.30214, 2.35291, 2.21367);
            TV centerD(1.89017, 1.70982, 1.90085);
            TV centerE(2.1572, 1.71448, 1.94296);
            Sphere<T, dim> sphereA(centerA, radA);
            Sphere<T, dim> sphereB(centerB, radB);
            Sphere<T, dim> sphereC(centerC, radC);
            Sphere<T, dim> sphereD(centerD, radD);
            Sphere<T, dim> sphereE(centerE, radE);

            //SPHERE A
            AnalyticCollisionObject<T, dim> sphereObjectA(sphereA, AnalyticCollisionObject<T, dim>::STICKY); //this sphere does not move
            init_helper.addAnalyticCollisionObject(sphereObjectA);

            //SPHERE B
            auto sphereTransformB = [](T time, AnalyticCollisionObject<T, dim>& object) {
                T sphereSpeed = 0.2;
                TV translation_velocity(-1 * sphereSpeed, sphereSpeed, 0);
                TV translation(-1 * sphereSpeed * time, sphereSpeed * time, 0); //multiply each velocity by dt to get dx!
                object.setTranslation(translation, translation_velocity);
                Vector<T, 3> omega(0, 0, 0);
                object.setAngularVelocity(omega);
            };
            AnalyticCollisionObject<T, dim> sphereObjectB(sphereTransformB, sphereB, AnalyticCollisionObject<T, dim>::STICKY); //this sphere moves
            init_helper.addAnalyticCollisionObject(sphereObjectB);

            //SPHERE C
            auto sphereTransformC = [](T time, AnalyticCollisionObject<T, dim>& object) {
                T sphereSpeed = 0.2;
                TV translation_velocity(sphereSpeed, sphereSpeed, 0);
                TV translation(sphereSpeed * time, sphereSpeed * time, 0); //multiply each velocity by dt to get dx!
                object.setTranslation(translation, translation_velocity);
                Vector<T, 3> omega(0, 0, 0);
                object.setAngularVelocity(omega);
            };
            AnalyticCollisionObject<T, dim> sphereObjectC(sphereTransformC, sphereC, AnalyticCollisionObject<T, dim>::STICKY); //this sphere moves
            init_helper.addAnalyticCollisionObject(sphereObjectC);

            //SPHERE D
            auto sphereTransformD = [](T time, AnalyticCollisionObject<T, dim>& object) {
                T sphereSpeed = 0.2;
                TV translation_velocity(-1 * sphereSpeed, -1 * sphereSpeed, 0);
                TV translation(-1 * sphereSpeed * time, -1 * sphereSpeed * time, 0); //multiply each velocity by dt to get dx!
                object.setTranslation(translation, translation_velocity);
                Vector<T, 3> omega(0, 0, 0);
                object.setAngularVelocity(omega);
            };
            AnalyticCollisionObject<T, dim> sphereObjectD(sphereTransformD, sphereD, AnalyticCollisionObject<T, dim>::STICKY); //this sphere moves
            init_helper.addAnalyticCollisionObject(sphereObjectD);

            //SPHERE E
            auto sphereTransformE = [](T time, AnalyticCollisionObject<T, dim>& object) {
                T sphereSpeed = 0.2;
                TV translation_velocity(sphereSpeed, -1 * sphereSpeed, 0);
                TV translation(sphereSpeed * time, -1 * sphereSpeed * time, 0); //multiply each velocity by dt to get dx!
                object.setTranslation(translation, translation_velocity);
                Vector<T, 3> omega(0, 0, 0);
                object.setAngularVelocity(omega);
            };
            AnalyticCollisionObject<T, dim> sphereObjectE(sphereTransformE, sphereE, AnalyticCollisionObject<T, dim>::STICKY); //this sphere moves
            init_helper.addAnalyticCollisionObject(sphereObjectE);

            //Set up a ground plane
            TV ground_origin(0, 0, 0);
            TV ground_normal(0, .1, 0);
            HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::SLIP);
            ground_object.setFriction(.1575);
            init_helper.addAnalyticCollisionObject(ground_object); // add ground
        }
    }
};
} // namespace ZIRAN
#endif
