/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_FREEZING_FREEZINGPROCESSDATA_H
#define PROCESSLIB_FREEZING_FREEZINGPROCESSDATA_H

namespace MeshLib
{
    class Element;
}


namespace ProcessLib
{

template <typename T>
struct Parameter;

namespace Freezing
{

// copy that and change names first still only thermal conductivity, later on more parameters will be added.
struct FreezingProcessData
{
    FreezingProcessData(
            ProcessLib::Parameter<double> const& porosity_,
            ProcessLib::Parameter<double> const& intrinsic_permeability_,
        //    ProcessLib::Parameter<double> const& storage_,
            ProcessLib::Parameter<double> const& viscosity_,
            ProcessLib::Parameter<double> const& density_solid_,
            ProcessLib::Parameter<double> const& density_fluid_,
            ProcessLib::Parameter<double> const& density_ice_,
            ProcessLib::Parameter<double> const& specific_heat_capacity_solid_,
            ProcessLib::Parameter<double> const& specific_heat_capacity_fluid_,
            ProcessLib::Parameter<double> const& specific_heat_capacity_ice_,
            ProcessLib::Parameter<double> const& thermal_conductivity_solid_,
            ProcessLib::Parameter<double> const& thermal_conductivity_fluid_,
            ProcessLib::Parameter<double> const& thermal_conductivity_ice_,
            ProcessLib::Parameter<double> const& thermal_expansion_coefficient_,
            ProcessLib::Parameter<double> const& sigmoid_coefficient_,
            ProcessLib::Parameter<double> const& latent_heat_,
            ProcessLib::Parameter<double> const& compressibility_,
            ProcessLib::Parameter<double> const& reference_temperature_,
            ProcessLib::Parameter<double> const& specific_body_force_,
            bool const has_gravity_)
            : porosity(porosity_),
              intrinsic_permeability(intrinsic_permeability_),
        //      storage(storage_),
              viscosity(viscosity_),
              density_solid(density_solid_),
              density_fluid(density_fluid_),
              density_ice(density_ice_),
              specific_heat_capacity_solid(specific_heat_capacity_solid_),
              specific_heat_capacity_fluid(specific_heat_capacity_fluid_),
              specific_heat_capacity_ice(specific_heat_capacity_ice_),
              thermal_conductivity_solid(thermal_conductivity_solid_),
              thermal_conductivity_fluid(thermal_conductivity_fluid_),
              thermal_conductivity_ice(thermal_conductivity_ice_),
              thermal_expansion_coefficient(thermal_expansion_coefficient_),
              sigmoid_coefficient(sigmoid_coefficient_),
              latent_heat(latent_heat_),
              compressibility(compressibility_),
              reference_temperature(reference_temperature_),
              specific_body_force(specific_body_force_),
              has_gravity(has_gravity_)
    {
    }

    FreezingProcessData(FreezingProcessData&& other)
        : porosity(other.porosity),
          intrinsic_permeability(other.intrinsic_permeability),
       //   storage(other.storage),
          viscosity(other.viscosity),
          density_solid(other.density_solid),
          density_fluid(other.density_fluid),
          density_ice(other.density_ice),
          specific_heat_capacity_solid(other.specific_heat_capacity_solid),
          specific_heat_capacity_fluid(other.specific_heat_capacity_fluid),
          specific_heat_capacity_ice(other.specific_heat_capacity_ice),
          thermal_conductivity_solid(other.thermal_conductivity_solid),
          thermal_conductivity_fluid(other.thermal_conductivity_fluid),
          thermal_conductivity_ice(other.thermal_conductivity_ice),
          thermal_expansion_coefficient(other.thermal_expansion_coefficient),
          sigmoid_coefficient(other.sigmoid_coefficient),
          latent_heat(other.latent_heat),
          compressibility(other.compressibility),
          reference_temperature(other.reference_temperature),
          specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity)
    {
    }

    //! Copies are forbidden.
    FreezingProcessData(FreezingProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(FreezingProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(FreezingProcessData&&) = delete;


    Parameter<double> const& porosity;
    Parameter<double> const& intrinsic_permeability;
//    Parameter<double> const& storage;
    Parameter<double> const& viscosity;
    Parameter<double> const& density_solid;
    Parameter<double> const& density_fluid;
    Parameter<double> const& density_ice;
    Parameter<double> const& specific_heat_capacity_solid;
    Parameter<double> const& specific_heat_capacity_fluid;
    Parameter<double> const& specific_heat_capacity_ice;
    Parameter<double> const& thermal_conductivity_solid;
    Parameter<double> const& thermal_conductivity_fluid;
    Parameter<double> const& thermal_conductivity_ice;
    Parameter<double> const& thermal_expansion_coefficient;
    Parameter<double> const& sigmoid_coefficient;
    Parameter<double> const& latent_heat;
    Parameter<double> const& compressibility;
    Parameter<double> const& reference_temperature;
    Parameter<double> const& specific_body_force;
    bool const has_gravity;
};

} // namespace Freezing
} // namespace ProcessLib

#endif // PROCESSLIB_GROUNDWATERFLOW_GROUNDWATERFLOWPROCESSDATA_H
