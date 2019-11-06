#ifndef PARTIO_IO_H
#define PARTIO_IO_H

#include <Eigen/Dense>
#include <Partio.h>
#include <vector>
#include <Ziran/CS/Util/Forward.h>
#include <Ziran/CS/Util/AttributeNamesForward.h>
#include <Ziran/Math/Geometry/Particles.h>
#include <Ziran/Physics/PlasticityApplier.h>
#include <Ziran/Physics/ConstitutiveModel/ConstitutiveModel.h>

namespace ZIRAN {

/**
   read particles positions from .bgeo file
   fill in mass with NAN
   fill in velocity with 0
*/
template <class T, int dim>
void readPartio(const std::string& particleFile, Particles<T, dim>& particles)
{
    Partio::ParticlesData* data = Partio::read(particleFile.c_str());
    ZIRAN_ASSERT(data, "could not open ", particleFile);

    Partio::ParticleAttribute posAttr;
    ZIRAN_ASSERT(data->attributeInfo("position", posAttr) && posAttr.type == Partio::VECTOR && posAttr.count == 3, "could not get position as vector of size 3");

    Partio::ParticleAccessor posAcc(posAttr);
    Partio::ParticlesData::const_iterator it = data->begin();
    it.addAccessor(posAcc);

    auto ap = particles.appender();

    Vector<T, dim> velocity = Vector<T, dim>::Zero();
    T mass = NAN;
    for (; it != data->end(); ++it) {
        float* pos = posAcc.raw<float>(it);
        Vector<T, dim> position = Vector<T, dim>::Zero();
        for (size_t i = 0; i < dim; i++)
            position(i) = pos[i];
        ap.append(mass, position, velocity);
    }
    data->release();
}

/**
   Suppose that you have StdVector<Vector<T, dim>> particle_positions that
   store the position of the particles
   Use this function to write these positions into bgeo.
   Example:
   StdVector<EigenTV> samples;
   std::string filename = "test.bgeo";
   writePositionVectorToPartio(filename, samples);
*/
template <class T, int dim>
void writePositionVectorToPartio(const std::string& particleFile,
    const StdVector<Vector<T, dim>>& particles)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (size_t iter = 0; iter < particles.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        for (int k = 0; k < dim; k++)
            p[k] = particles[iter](k);
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

template <class T, int dim>
void writePositionVectorToPartio(const std::string& particleFile,
    const StdVector<Vector<T, dim>>& particles, const StdVector<int>& inside_out)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH;
    Partio::ParticleAttribute boolH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    boolH = parts->addAttribute("inside", Partio::INT, 1);
    std::cout << "samples.size() = " << particles.size() << std::endl;
    std::cout << "inside out.size() = " << inside_out.size() << std::endl;
    for (size_t iter = 0; iter < particles.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        for (int k = 0; k < dim; k++)
            p[k] = particles[iter](k);

        int* b = parts->dataWrite<int>(boolH, idx);
        int bb = (int)(inside_out[iter]);
        b[0] = bb;
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

template <class T, int dim>
void writeMatrixVectorToPartio(
    const std::string& particleFile,
    const StdVector<Matrix<T, dim, dim>>& matrices)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH;
    Partio::ParticleAttribute col1, col2, col3;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    col1 = parts->addAttribute("firstcolumn", Partio::VECTOR, 3);
    col2 = parts->addAttribute("secondcolumn", Partio::VECTOR, 3);
    col3 = parts->addAttribute("thirdcolumn", Partio::VECTOR, 3);
    for (size_t iter = 0; iter < matrices.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* r1 = parts->dataWrite<float>(col1, idx);
        float* r2 = parts->dataWrite<float>(col2, idx);
        float* r3 = parts->dataWrite<float>(col3, idx);
        for (int k = 0; k < 3; k++) {
            p[k] = 0;
            r1[k] = 0;
            r2[k] = 0;
            r3[k] = 0;
        }
        for (int k = 0; k < dim; k++) {
            r1[k] = matrices[iter](k, 0);
            r2[k] = matrices[iter](k, 1);
            if (dim == 3)
                r3[k] = matrices[iter](k, 2);
        }
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

//#########################################################################
// Function: writePartio
//
// Dump out a bgeo file that contains particles' X, V, m.
//#########################################################################
template <class T, int dim>
void writePartio(const std::string& particleFile,
    const Particles<T, dim>& particles)
{
    // AttributeName<T> mass = Particles<T, dim>::mass_name();
    AttributeName<Vector<T, dim>> x_name = Particles<T, dim>::X_name();
    // AttributeName<Vector<T, dim>> v_name = Particles<T, dim>::V_name();

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    // mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    // vH = parts->addAttribute("v", Partio::VECTOR, 3);
    // auto iter = particles.iter(mass, x_name, v_name);
    auto iter = particles.iter(x_name);
    for (; iter; ++iter) {
        int idx = parts->addParticle();
        // float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        // float* v = parts->dataWrite<float>(vH, idx);

        // T mass = iter.template get<0>();
        // m[0] = mass;

        Vector<T, dim> X = iter.template get<0>();
        // Vector<T, dim> V = iter.template get<2>();

        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        // for (int k = 0; k < 3; k++)
        //    v[k] = (T)0;

        for (int k = 0; k < dim; k++)
            p[k] = X(k);
        // for (int k = 0; k < dim; k++)
        //    v[k] = V(k);
    }

    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

//#########################################################################
// Function: writePartioAttribute （scalar)
//
// Dump out a bgeo file that contains a certain scalar particle attribute.
// For particles without the attribute, put a background value there.
//#########################################################################
template <class T, int dim>
void writePartioAttribute(const std::string& particleFile, const Particles<T, dim>& particles, AttributeName<T> scalar_attribute, const T background_value)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);

    // gather data as a flat array
    StdVector<T> data(particles.count, background_value);
    for (auto iter = particles.iter(scalar_attribute); iter; ++iter)
        data[iter.entryId()] = iter.template get<0>();

    // write to partio structure
    for (int k = 0; k < particles.count; k++) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        p[0] = (float)data[k];
        p[1] = (float)data[k];
        p[2] = (float)data[k];
    }

    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

//#########################################################################
// Function: writePartioAttribute (plasticity paramater in a Vec3)
//
// Dump out a bgeo file that contains a certain attribute from a plasticity model.
// For particles without the attribute, put a background value there.
// The plasticity model needs to provide function fillAttributesToVec3
//#########################################################################
template <class T, class T_PLASTICITY, int dim>
void writePartioAttribute(const std::string& particleFile, const Particles<T, dim>& particles, AttributeName<T_PLASTICITY> plasticity, const T background_value)
{
    using TV3 = Vector<T, 3>;
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);

    // gather data as a flat array
    StdVector<TV3> data(particles.count, TV3::Ones() * background_value);
    for (auto iter = particles.iter(plasticity); iter; ++iter) {
        auto pp = iter.template get<0>();
        pp.fillAttributesToVec3(data[iter.entryId()]);
    }

    // write to partio structure
    for (int k = 0; k < particles.count; k++) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        p[0] = (float)data[k](0);
        p[1] = (float)data[k](1);
        p[2] = (float)data[k](2);
    }

    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

template <class T, int dim>
void writePartioWaterDebugReversible(const std::string& particleFile, const Particles<T, dim>& particles)
{
    AttributeName<T> mass = Particles<T, dim>::mass_name();
    AttributeName<Vector<T, dim>> x_name = Particles<T, dim>::X_name();
    AttributeName<Vector<T, dim>> v_name = Particles<T, dim>::V_name();

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH, water_pressureH, water_grad_porosityH;
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);
    water_pressureH = parts->addAttribute("pressure", Partio::VECTOR, 1);
    water_grad_porosityH = parts->addAttribute("grad_porosity", Partio::VECTOR, 3);
    auto iter = particles.iter(mass, x_name, v_name, water_pressure_name<T>(), water_grad_porosity_name<T, dim>());
    for (; iter; ++iter) {
        int idx = parts->addParticle();
        float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        float* water_pressure = parts->dataWrite<float>(water_pressureH, idx);
        float* water_grad_porosity = parts->dataWrite<float>(water_grad_porosityH, idx);

        T mass = iter.template get<0>();
        Vector<T, dim> X = iter.template get<1>();
        Vector<T, dim> V = iter.template get<2>();
        T pressure = iter.template get<3>();
        Vector<T, dim> grad_porosity = iter.template get<4>();

        m[0] = mass;

        p[2] = (T)0;
        v[2] = (T)0;
        water_grad_porosity[2] = 0;

        for (int k = 0; k < dim; k++)
            p[k] = X(k);
        for (int k = 0; k < dim; k++)
            v[k] = V(k);

        water_pressure[0] = pressure;

        for (int k = 0; k < dim; k++)
            water_grad_porosity[k] = grad_porosity(k);
    }

    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}
} // namespace ZIRAN
#endif
