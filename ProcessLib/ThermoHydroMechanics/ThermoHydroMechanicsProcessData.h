/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <int DisplacementDim>
struct ThermoHydroMechanicsProcessData
{
    ThermoHydroMechanicsProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& storage_coefficient_,
        Parameter<double> const& fluid_viscosity_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& porosity_,
        Parameter<double> const& solid_density_,
        Parameter<double> const& fluid_density_,
        Parameter<double> const& solid_linear_thermal_expansion_coefficient_,
        Parameter<double> const& fluid_volumetric_thermal_expansion_coefficient_,
        Parameter<double> const& fluid_specific_heat_capacity_,
        Parameter<double> const& solid_specific_heat_capacity_,
        Parameter<double> const& solid_thermal_conductivity_,
        Parameter<double> const& fluid_thermal_conductivity_,
        Parameter<double> const& reference_temperature_,
        Eigen::Matrix<double, DisplacementDim, 1>
            specific_body_force_)
        : material{std::move(material_)},
          intrinsic_permeability(intrinsic_permeability_),
          storage_coefficient(storage_coefficient_),
          fluid_viscosity(fluid_viscosity_),
          biot_coefficient(biot_coefficient_),
          porosity(porosity_),
          solid_density(solid_density_),
          fluid_density(fluid_density_),
          solid_linear_thermal_expansion_coefficient(solid_linear_thermal_expansion_coefficient_),
          fluid_volumetric_thermal_expansion_coefficient(fluid_volumetric_thermal_expansion_coefficient_),
          fluid_specific_heat_capacity(fluid_specific_heat_capacity_),
          solid_specific_heat_capacity(solid_specific_heat_capacity_),
          solid_thermal_conductivity(solid_thermal_conductivity_),
          fluid_thermal_conductivity(fluid_thermal_conductivity_),
          reference_temperature(reference_temperature_),
          specific_body_force(std::move(specific_body_force_))
    {
    }

    ThermoHydroMechanicsProcessData(ThermoHydroMechanicsProcessData&& other)
        : material{std::move(other.material)},
          intrinsic_permeability(other.intrinsic_permeability),
          storage_coefficient(other.storage_coefficient),
          fluid_viscosity(other.fluid_viscosity),
          biot_coefficient(other.biot_coefficient),
          porosity(other.porosity),
          solid_density(other.solid_density),
          fluid_density(other.fluid_density),
          solid_linear_thermal_expansion_coefficient(other.solid_linear_thermal_expansion_coefficient),
          fluid_volumetric_thermal_expansion_coefficient(other.fluid_volumetric_thermal_expansion_coefficient),
          fluid_specific_heat_capacity(other.fluid_specific_heat_capacity),
          solid_specific_heat_capacity(other.solid_specific_heat_capacity),
          solid_thermal_conductivity(other.solid_thermal_conductivity),
          fluid_thermal_conductivity(other.fluid_thermal_conductivity),
          reference_temperature(other.reference_temperature),
          specific_body_force(other.specific_body_force),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    ThermoHydroMechanicsProcessData(ThermoHydroMechanicsProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(ThermoHydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermoHydroMechanicsProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& storage_coefficient;
    Parameter<double> const& fluid_viscosity;
    Parameter<double> const& biot_coefficient;
    Parameter<double> const& porosity;
    Parameter<double> const& solid_density;
    Parameter<double> const& fluid_density;
    Parameter<double> const& solid_linear_thermal_expansion_coefficient;
    Parameter<double> const& fluid_volumetric_thermal_expansion_coefficient;
    Parameter<double> const& fluid_specific_heat_capacity;
    Parameter<double> const& solid_specific_heat_capacity;
    Parameter<double> const& solid_thermal_conductivity;
    Parameter<double> const& fluid_thermal_conductivity;
    Parameter<double> const& reference_temperature;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0.0;
    double t = 0.0;
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
