#ifndef ADMM_INIT_3D_H
#define ADMM_INIT_3D_H

#include <Ziran/CS/Util/RandomNumber.h>
#include <Ziran/Math/Geometry/MeshConstruction.h>
#include <Ziran/Physics/SoundSpeedCfl.h>
#include <Ziran/Sim/MeshHandle.h>
#include <Ziran/Sim/SceneInitializationCore.h>
#include "AdmmSimulation.h"
#include "AdmmInit.h"

namespace ZIRAN {

template <class T, int dim>
class AdmmInitBase;

template <class T>
class AdmmInit3D : public AdmmInitBase<T, 3> {
public:
    static const int dim = 3;
    using Base = AdmmInitBase<T, dim>;
    using TV2 = Vector<T, 2>;
    using TVI2 = Vector<int, 2>;
    using TV = Vector<T, dim>;
    using TM = Eigen::Matrix<T, dim, dim>;
    using TVI = Vector<int, dim>;

    using Base::init_helper;
    using Base::scene;
    using Base::sim;
    using Base::test_number;

    AdmmInit3D(AdmmSimulation<T, dim>& sim, const int test_number)
        : Base(sim, test_number)
    {
    }

    void reload() override
    {
        if (test_number == 1) {
            sim.output_dir.path = "output/twist";
            sim.end_frame = 150;
            sim.dx = 0.01;
            sim.gravity = 0 * TV::Unit(1);
            sim.step.max_dt = sim.step.frame_dt / 4;
            sim.symplectic = false;
            sim.cfl = 0.6;
            sim.rho_scale = 0.5;
            sim.autorestart = false;

            sim.use_ruiz = false;
            sim.admm_max_iterations = 30;
            sim.local_tolerance = 1e-8;
            sim.global_tolerance = 1e-14;
            sim.global_max_iteration = 15;

            sim.enable_visco = false;

            {
                AxisAlignedAnalyticBox<T, dim> box_out(TV(4.25, 4.9, 4.9), TV(5.75, 5.1, 5.1));
                AxisAlignedAnalyticBox<T, dim> box_in(TV(4.25, 4.925, 4.925), TV(5.75, 5.075, 5.075));
                DifferenceLevelSet<T, dim> materialRegionLevelSet;
                materialRegionLevelSet.add(box_out, box_in);
                StdVector<TV> meshed_points;
                readPositionObj("metal.obj", meshed_points);
                MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleInAnalyticLevelSetWithExistingPoints(meshed_points, materialRegionLevelSet, 2, 8);
                //// Elasticity
                // StvkWithHencky<T, dim> model(100, 0.4);
                // particles_handle.addFBasedMpmForce(model);
                // particles_handle.addDummyPlasticity();
                //// Metal
                StvkWithHencky<T, dim> model(2000, 0.4);
                particles_handle.addFBasedMpmForce(model);
                VonMisesStvkHencky<T, dim> p(1, FLT_MAX, 0);
                particles_handle.addPlasticity(model, p, "F");
            }
            { // floor
                TV ground_origin(5, 4.25, 5);
                TV ground_normal(0, 1, 0);
                HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
                AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::STICKY);
                init_helper.addAnalyticCollisionObject(ground_object);
            }

            auto sphere_transform1 = [](T time, AnalyticCollisionObject<T, dim>& object) {
                object.setTranslation(TV(4.25, 5, 5), TV(0, 0, 0));
            };
            auto sphere_transform2 = [](T time, AnalyticCollisionObject<T, dim>& object) {
                T stretch_time = 3; // 3
                T release_time = 4; // 4
                if (time < stretch_time) {
                    T f = 0.1;
                    T theta = 2 * M_PI * f;
                    TV omega(-theta, 0, 0);
                    object.setRotation(Vector<T, 4>(1, 0, 0, 0));
                    object.setAngularVelocity(omega);
                    object.setTranslation(TV(5.75, 5, 5), TV(0, 0, 0));
                }
                else if (time < release_time) {
                    TV omega(0, 0, 0);
                    object.setRotation(Vector<T, 4>(1, 0, 0, 0));
                    object.setAngularVelocity(omega);
                    object.setTranslation(TV(5.75, 5, 5), TV(0, 0, 0));
                }
                else {
                    object.setTranslation(TV(15, 5, 5), TV(0, 0, 0));
                }
            };

            TV sphere_center1(0, 0, 0);
            Sphere<T, dim> sphere1(sphere_center1, 0.2);
            AnalyticCollisionObject<T, dim> sphere_object1(sphere_transform1, sphere1, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(sphere_object1);
            TV sphere_center2(0, 0, 0);
            Sphere<T, dim> sphere2(sphere_center2, 0.2);
            AnalyticCollisionObject<T, dim> sphere_object2(sphere_transform2, sphere2, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(sphere_object2);
        }

        if (test_number == 2) {
            sim.output_dir.path = "output/crash";
            sim.end_frame = 120;
            sim.dx = 0.02;
            sim.gravity = 0 * TV::Unit(1);
            sim.step.max_dt = sim.step.frame_dt / 12;
            sim.symplectic = false;
            sim.cfl = 0.6;
            sim.rho_scale = 1.0;
            sim.autorestart = false;

            sim.use_ruiz = false;
            sim.admm_max_iterations = 30;
            sim.local_tolerance = 1e-8;
            sim.global_tolerance = 1e-14;
            sim.global_max_iteration = 25;

            sim.enable_visco = false;

            {
                StdVector<TV> meshed_points;
                readPositionObj("car_left.obj", meshed_points);
                MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleFromVdbFileWithExistingPoints(meshed_points, "LevelSets/car_left.vdb", 2, 8);
                StvkWithHencky<T, dim> model(2000, 0.4);
                particles_handle.addFBasedMpmForce(model);
                VonMisesStvkHencky<T, dim> p(2, FLT_MAX, 0);
                particles_handle.addPlasticity(model, p, "F");
                auto initial_translation = [=](int index, Ref<T> mass, TV& X, TV& V) {
                    V = TV(0, 0, 1);
                    X -= TV(0, 0, 2);
                };
                particles_handle.transform(initial_translation);
            }
            {
                StdVector<TV> meshed_points;
                readPositionObj("car_right.obj", meshed_points);
                MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleFromVdbFileWithExistingPoints(meshed_points, "LevelSets/car_right.vdb", 2, 8);
                StvkWithHencky<T, dim> model(2000, 0.4);
                particles_handle.addFBasedMpmForce(model);
                VonMisesStvkHencky<T, dim> p(2, FLT_MAX, 0);
                particles_handle.addPlasticity(model, p, "F");
            }
            { // floor
                TV ground_origin(5, 4.86, 5);
                TV ground_normal(0, 1, 0);
                HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
                AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::SLIP);
                init_helper.addAnalyticCollisionObject(ground_object);
            }
            { // floor
                TV ground_origin(5, 5, 6.4);
                TV ground_normal(0, 0, -1);
                HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
                AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::SLIP);
                init_helper.addAnalyticCollisionObject(ground_object);
            }
        }

        if (test_number == 3) {
            sim.output_dir.path = "output/collapse";
            sim.end_frame = 120;
            sim.dx = 0.004;
            sim.gravity = -0.5 * TV::Unit(1);
            sim.step.max_dt = sim.step.frame_dt / 6;
            sim.symplectic = false;
            sim.cfl = 0.4;
            sim.rho_scale = 0.1;
            sim.autorestart = false;

            sim.use_ruiz = false;
            sim.admm_max_iterations = 15;
            sim.local_tolerance = 1e-8;
            sim.global_tolerance = 1e-14;
            sim.global_max_iteration = 15;

            sim.enable_visco = false;

            T rho = 1400;
            T nu = 0.3, youngs = 1e7;

            T radius = 0.1;
            T height = 0.4;
            T theta = 0;
            Vector<T, 4> rotation(std::cos(theta / 2), std::sin(theta / 2), 0, 0);
            TV translation(5, height / 2 + 5 + sim.dx, 5);
            CappedCylinder<T, dim> cylinder(radius, height, rotation, translation);
            MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleInAnalyticLevelSet(cylinder, rho, 8);
            StvkWithHencky<T, dim> model(youngs, nu);
            particles_handle.addFBasedMpmForce(model);
            DruckerPragerStvkHencky<T> p(30);
            particles_handle.addPlasticity(model, p, "F");

            TV ground_origin(5, 5, 5);
            TV ground_normal(0, 1, 0);
            HalfSpace<T, dim> ground_ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground_object(ground_ls, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(ground_object);
        }

        if (test_number == 4) {
            sim.output_dir.path = "output/lion";
            sim.viscosity_v = 10;
            sim.viscosity_d = 10;
            sim.end_frame = 600;
            sim.dx = 0.0075 * 2;
            sim.gravity = TV(0, -.3, 0);
            sim.step.max_dt = 0.0005 * 2 * 16; //*1 for tet mesh
            sim.step.frame_dt = 1. / 24.;
            sim.quasistatic = false;
            sim.symplectic = false;
            sim.objective.matrix_free = true;
            sim.rho_scale = 0.8;
            sim.verbose = false;
            sim.cfl = 0.6;
            sim.use_elasticity_plasticity = false;
            init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx, AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for SPGrid.

            sim.use_ruiz = false;
            sim.admm_max_iterations = 20;
            sim.local_tolerance = 1e-8;
            sim.global_tolerance = 1e-14;
            sim.global_max_iteration = 20;

            sim.enable_visco = true;

            // ground is at 0
            TV ground_origin(0, 0.6, 0);
            TV ground_normal(0, 1, 0);
            HalfSpace<T, dim> ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground([&](T, AnalyticCollisionObject<T, dim>&) {}, ls, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(ground);

            T rho = 2;
            T ppc = 500;
            T E = 50;

#if 0 // pure elasticity tet mesh


            StvkWithHenckyIsotropic<T, dim> equilibrated_model(E, 0.4);
            equilibrated_model.mu *= 0.1;

            using TMeshTet = SimplexMesh<3>;
            auto reader = MeshReader<T, dim, TMeshTet>("TetMesh/lion_minchen.vtk");
            MeshHandle<T, dim, TMeshTet> mesh = scene.createMesh((std::function<void(TMeshTet&, StdVector<TV>&)>)reader);
            DeformableObjectHandle<T, dim, TMeshTet> deformable = scene.addDeformableObject(mesh);
            deformable.setMassFromDensity(rho);

            deformable.transform(
                    [&](int index, Ref<T> mass, TV& X, TV& V) {
                        X = X + TV(5, 5 , 5);
                    });
            deformable.addFemHyperelasticForce(equilibrated_model);
#else // viscoelastic MPM

            // Embedding a mesh for rendering.
            StdVector<TV> meshed_points;
            readPositionObj("lion.obj", meshed_points);
            MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleFromVdbFileWithExistingPoints(meshed_points, "LevelSets/lion.vdb", rho, ppc);
            auto initial_translation = [=](int index, Ref<T> mass, TV& X, TV& V) {
                X = X + TV(5, 5, 5);
            };
            particles_handle.transform(initial_translation);

            // Not embedding a mesh for rendering.
            // MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleInAnalyticLevelSet(cylinder, rho, ppc);

            StvkWithHenckyIsotropic<T, dim> equilibrated_model(E, 0.4);
            equilibrated_model.mu *= 0.1;
            particles_handle.addFBasedMpmForce(equilibrated_model);
            StvkWithHencky<T, dim> nonequilibrated_model(E, 0.4);
            particles_handle.addFElasticNonequilibratedBasedMpmForce(nonequilibrated_model, (T)10, (T)10);

#endif

            TV translation_freeze = TV(5, 5.345, 5) + TV(0, 24 * 0.0025, 0) * 50. / 24.;
            T rotation_freeze = (T)24 * 5 / 180 * M_PI * 50. / 24.;

            auto cylinder_transform = [translation_freeze, rotation_freeze](T time, AnalyticCollisionObject<T, dim>& object) {
                T frame_dt = 1. / 24.;
                T frame_move = 50;
                T frame_freeze = 100;
                if (time < frame_move * frame_dt) {
                    // rotate
                    T theta = (T)24 * 5 / 180 * M_PI;
                    Vector<T, 4> rotation(std::cos(theta * time / 2), 0, std::sin(theta * time / 2), 0);
                    object.setRotation(rotation);
                    TV omega(0, theta, 0);
                    object.setAngularVelocity(omega);

                    TV translation_fixed(5, 5.345, 5);

                    TV translation_velocity(0, 24 * 0.0025, 0);
                    TV translation = translation_fixed + translation_velocity * time;
                    object.setTranslation(translation, translation_velocity);
                }
                else if (time < frame_freeze * frame_dt) {

                    Vector<T, 4> rotation(std::cos(rotation_freeze / 2), 0, std::sin(rotation_freeze / 2), 0);
                    object.setRotation(rotation);
                    TV omega(0, 0, 0);
                    object.setAngularVelocity(omega);

                    TV translation = translation_freeze;
                    TV translation_velocity(0, 0, 0);
                    object.setTranslation(translation, translation_velocity);
                }
                else {
                    TV translation(0, 0, 0);
                    TV translation_velocity(0, 0, 0);
                    object.setTranslation(translation, translation_velocity);
                }
            };

            T radius = .09;
            T height = .2;
            T theta = 0;
            Vector<T, 4> rotation(std::cos(theta / 2), std::sin(theta / 2), 0, 0);
            TV translation(0, 0, 0);
            CappedCylinder<T, dim> cylinder(radius, height, rotation, translation);
            AnalyticCollisionObject<T, dim> sphere_object1(cylinder_transform, cylinder, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(sphere_object1);

            TV box_center(5, 5, 5);
            TV box_size(.2, .1, 2);
            TV box_min_corner = box_center - box_size * .5;
            TV box_max_corner = box_center + box_size * .5;
            AxisAlignedAnalyticBox<T, dim> box(box_min_corner, box_max_corner);
            AnalyticCollisionObject<T, dim> box_object(box, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(box_object);
        }

        // viscoelastic silicone rubber
        if (test_number == 5) {
            sim.output_dir.path = "output/silicone";
            sim.viscosity_v = 10;
            sim.viscosity_d = 10;

            sim.end_frame = 400;
            sim.dx = 0.0075;
            sim.gravity = -0.3 * TV::Unit(1);
            sim.step.max_dt = 8.5e-3 / 2;
            sim.symplectic = false;
            sim.cfl = 0.6;
            sim.rho_scale = 0.8;
            sim.autorestart = false;

            sim.use_ruiz = false;
            sim.admm_max_iterations = 40;
            sim.local_tolerance = 1e-8;
            sim.global_tolerance = 1e-8;
            sim.global_max_iteration = 100000;

            sim.enable_visco = true;
            sim.use_elasticity_plasticity = false;

            init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx, AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for SPGrid.

            // ground is at 0
            TV ground_origin(0, 0.6, 0);
            TV ground_normal(0, 1, 0);
            HalfSpace<T, dim> ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground([&](T, AnalyticCollisionObject<T, dim>&) {}, ls, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(ground);

            T rho = 2;
            T ppc = 20;
            T E = 50;

            T radius = .3;
            T height = .05;
            T theta = 0;
            Vector<T, 4> rotation(std::cos(theta / 2), std::sin(theta / 2), 0, 0);
            TV translation(1, 1, 1);
            CappedCylinder<T, dim> cylinder(radius, height, rotation, translation);

            // Embedding a mesh for rendering.
            StdVector<TV> meshed_points;
            readPositionObj("silicone_rubber_tri.obj", meshed_points);
            MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleInAnalyticLevelSetWithExistingPoints(meshed_points, cylinder, rho, ppc);

            // Not embedding a mesh for rendering.
            // MpmParticleHandleBase<T, dim> particles_handle = init_helper.sampleInAnalyticLevelSet(cylinder, rho, ppc);

            StvkWithHenckyIsotropic<T, dim> equilibrated_model(E, 0.4);
            equilibrated_model.mu *= 0.1;
            particles_handle.addFBasedMpmForce(equilibrated_model);
            StvkWithHencky<T, dim> nonequilibrated_model(E, 0.4);
            particles_handle.addFElasticNonequilibratedBasedMpmForce(nonequilibrated_model, (T)10, (T)10);

            auto sphere_transform1 = [](T time, AnalyticCollisionObject<T, dim>& object) {
                T v = 0.1;
                T stretch_time = 2;
                T release_time = 4;
                if (time < stretch_time) {
                    TV translation(v * time, 0, 0);
                    TV translation_velocity(v, 0, 0);
                    object.setTranslation(translation, translation_velocity);
                }
                else if (time < release_time) {
                    TV translation(v * stretch_time, 0, 0);
                    TV translation_velocity(0, 0, 0);
                    object.setTranslation(translation, translation_velocity);
                }
                else {
                    TV translation(10, 0, 0);
                    TV translation_velocity(0, 0, 0);
                    object.setTranslation(translation, translation_velocity);
                }
            };
            auto sphere_transform2 = [](T time, AnalyticCollisionObject<T, dim>& object) {
                T v = -0.1;
                T stretch_time = 2;
                T release_time = 4;
                if (time < stretch_time) {
                    TV translation(v * time, 0, 0);
                    TV translation_velocity(v, 0, 0);
                    object.setTranslation(translation, translation_velocity);
                }
                else if (time < release_time) {
                    TV translation(v * stretch_time, 0, 0);
                    TV translation_velocity(0, 0, 0);
                    object.setTranslation(translation, translation_velocity);
                }
                else {
                    TV translation(10, 0, 0);
                    TV translation_velocity(0, 0, 0);
                    object.setTranslation(translation, translation_velocity);
                }
            };

            TV sphere_center1(1.3, 1, 1);
            Sphere<T, dim> sphere1(sphere_center1, 0.1);
            AnalyticCollisionObject<T, dim> sphere_object1(sphere_transform1, sphere1, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(sphere_object1);
            TV sphere_center2(.7, 1, 1);
            Sphere<T, dim> sphere2(sphere_center2, 0.1);
            AnalyticCollisionObject<T, dim> sphere_object2(sphere_transform2, sphere2, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(sphere_object2);
        }
        // honey coil
        // ./mpm -test 2008 -E 1e-4
        if (test_number == 6) {
            T p_E = 1e-4;
            sim.output_dir.path = "output/coil_1e-4";
            sim.viscosity_v = 0.4;
            sim.viscosity_d = 0.4;
            sim.end_frame = 240;
            sim.dx = 0.01;
            sim.gravity = TV(0, -1, 0);
            sim.step.max_dt = 5e-4;
            sim.quasistatic = false;
            sim.symplectic = false;
            sim.objective.matrix_free = true;
            sim.verbose = false;
            sim.cfl = 0.6;
            init_helper.addAllWallsInDomain(4096 * sim.dx, 5 * sim.dx, AnalyticCollisionObject<T, dim>::STICKY); // add safety domain walls for SPGrid.
            sim.rho_scale = 0.8;

            sim.use_ruiz = false;
            sim.admm_max_iterations = 30;
            sim.local_tolerance = 1e-8;
            sim.global_tolerance = 1e-8;
            sim.global_max_iteration = 100000;

            sim.enable_visco = true;
            sim.use_elasticity_plasticity = false;

            // ground is at 0
            TV ground_origin(0, 0.095, 0);
            TV ground_normal(0, 1, 0);
            HalfSpace<T, dim> ls(ground_origin, ground_normal);
            AnalyticCollisionObject<T, dim> ground([&](T, AnalyticCollisionObject<T, dim>&) {}, ls, AnalyticCollisionObject<T, dim>::STICKY);
            init_helper.addAnalyticCollisionObject(ground);

            // create source collision object
            T rho = 2;
            T E = 100;
            int ppc = 8;
            Sphere<T, dim> sphere(TV(1, 1, 1), .03);
            TV material_speed(0.0, -0.8, 0);
            SourceCollisionObject<T, dim> sphere_source(sphere, material_speed);
            int source_id = init_helper.addSourceCollisionObject(sphere_source);
            init_helper.sampleSourceAtTheBeginning(source_id, rho, ppc);

            // initialize empty particle handle for restarting
            MpmParticleHandleBase<T, dim> empty_particle_handle = init_helper.getZeroParticle();
            StvkWithHenckyIsotropic<T, dim> equilibrated_model(E, 0.4);
            equilibrated_model.lambda *= 10;
            empty_particle_handle.addFBasedMpmForce(equilibrated_model);
            VonMisesStvkHencky<T, dim> p(p_E, FLT_MAX, 0);
            empty_particle_handle.addPlasticity(equilibrated_model, p, "F");
            StvkWithHencky<T, dim> nonequilibrated_model(E * .2, 0.4);
            empty_particle_handle.addFElasticNonequilibratedBasedMpmForce(nonequilibrated_model, (T).4, (T).4);

            sim.end_time_step_callbacks.push_back(
                [this, source_id, rho, ppc, E, p_E](int frame, int substep) {
                    if (frame < 2400) {
                        // add more particles from source Collision object
                        int N = init_helper.sourceSampleAndPrune(source_id, rho, ppc);
                        if (N) {
                            MpmParticleHandleBase<T, dim> source_particles_handle = init_helper.getParticlesFromSource(source_id, rho, ppc);
                            StvkWithHenckyIsotropic<T, dim> equilibrated_model(E, 0.4);
                            equilibrated_model.lambda *= 10;
                            source_particles_handle.addFBasedMpmForce(equilibrated_model);
                            VonMisesStvkHencky<T, dim> p(p_E, FLT_MAX, 0);
                            source_particles_handle.addPlasticity(equilibrated_model, p, "F");
                            StvkWithHencky<T, dim> nonequilibrated_model(E * .2, 0.4);
                            source_particles_handle.addFElasticNonequilibratedBasedMpmForce(nonequilibrated_model, (T).4, (T).4);
                        }
                    }
                });
        }
    }
};
} // namespace ZIRAN
#endif
