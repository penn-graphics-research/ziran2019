#ifndef ATTRIBUTE_NAMES_FORWARD_H
#define ATTRIBUTE_NAMES_FORWARD_H
#include <Ziran/CS/DataStructure/DataManager.h>

namespace ZIRAN {

// Previously on Projects/mpm/JBasedMpmForceHelper
template <class T>
inline static AttributeName<T> element_measure_name()
{
    return AttributeName<T>("element measure");
}

template <class T>
inline static AttributeName<T> J_name()
{
    return AttributeName<T>("J");
}

template <class T>
inline static AttributeName<T> ether_drag_name()
{
    return AttributeName<T>("ether");
}

template <class TM>
inline static AttributeName<TM> C_name()
{
    return AttributeName<TM>("C");
}

template <class TM>
inline static AttributeName<TM> phase_field_name()
{
    return AttributeName<TM>("phase field");
}

// For multiphase mpm, i.e. beach
template <class T>
inline static AttributeName<T> volume_fraction_name()
{
    return AttributeName<T>("volume fraction");
}

template <class T>
inline static AttributeName<T> surface_tension_coefficient_name()
{
    return AttributeName<T>("surface tension");
}

template <class T>
inline static AttributeName<T> density_summation_bulk_name()
{
    return AttributeName<T>("density summation bulk");
}

template <class T>
inline static AttributeName<T> initial_density_name()
{
    return AttributeName<T>("initial density");
}

template <class T>
inline static AttributeName<T> density_name()
{
    return AttributeName<T>("density");
}

template <class T>
inline static AttributeName<T> interpolated_mass_name()
{
    return AttributeName<T>("interpolated mass");
}

template <class T, int dim>
inline static AttributeName<Vector<T, dim>> normal_name()
{
    return AttributeName<Vector<T, dim>>("normal");
}

template <class T, int dim>
inline static AttributeName<Vector<T, dim>> X0_name()
{
    return AttributeName<Vector<T, dim>>("X0");
}

template <class T, int manifold_dim>
inline static AttributeName<Vector<T, manifold_dim + 1>> bary_weights_name()
{
    return AttributeName<Vector<T, manifold_dim + 1>>("bary_weights");
}

inline static AttributeName<int> parent_element_name()
{
    return AttributeName<int>("parent_element");
}

template <class T>
inline static AttributeName<T> saturation_name()
{
    return AttributeName<T>("saturation");
}

template <class T>
inline static AttributeName<T> cohesion_base_name()
{
    return AttributeName<T>("cohesion_base");
}

template <class T, int dim>
inline static AttributeName<Eigen::Matrix<T, dim, dim>> F_name()
{
    return AttributeName<Eigen::Matrix<T, dim, dim>>("F");
}

template <class T, int dim>
inline static AttributeName<Eigen::Matrix<T, dim, dim>> FEN_name()
{
    return AttributeName<Eigen::Matrix<T, dim, dim>>("Fe_Nonequilibrated");
}

template <class T>
inline static AttributeName<T> water_pressure_name()
{
    return AttributeName<T>("water_pressure");
}

template <class T, int dim>
inline static AttributeName<Eigen::Matrix<T, dim, 1>> water_grad_porosity_name()
{
    return AttributeName<Eigen::Matrix<T, dim, 1>>("water_grad_porosity");
}

template <class T>
inline static AttributeName<T> conductivity_name()
{
    return AttributeName<T>("conductivity");
}

template <class T>
inline static AttributeName<T> capacity_name()
{
    return AttributeName<T>("capacity");
}

template <class T>
inline static AttributeName<T> temperature_name()
{
    return AttributeName<T>("temperature");
}

template <class TV>
inline static AttributeName<TV> temperature_gradient_name()
{
    return AttributeName<TV>("temperature gradient");
}

template <class TM>
inline static AttributeName<TM> vtau_name()
{
    return AttributeName<TM>("vtau");
}

template <class T>
inline static AttributeName<T> yield_stress_name()
{
    return AttributeName<T>("yield_stress");
}

// template <class T>
// inline static AttributeName<T> p_name()
// {
//     return AttributeName<T>("p");
// }

// template <class T>
// inline static AttributeName<T> q_name()
// {
//     return AttributeName<T>("q");
// }

template <class TV>
inline static AttributeName<TV> fp_name()
{
    return AttributeName<TV>("fp");
}

static inline AttributeName<int> species_name()
{
    return AttributeName<int>("species");
}

template <class T>
inline static AttributeName<T> F_Dilational_name()
{
    return AttributeName<T>("F_Dilational");
}

template <class TM>
inline static AttributeName<TM> F_Distortional_name()
{
    return AttributeName<TM>("F_Distortional");
}
} // namespace ZIRAN
#endif
