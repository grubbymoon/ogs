/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_
#define PROCESSLIB_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
struct HydroMechanicsProcessData
{
    HydroMechanicsProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& storage_coefficient_,
        Parameter<double> const& fluid_viscosity_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& porosity_,
        Parameter<double> const& solid_density_,
        Parameter<double> const& fluid_density_,
        Parameter<double> const& specific_body_force_)
        : material{std::move(material_)},
          intrinsic_permeability(intrinsic_permeability_),
          storage_coefficient(storage_coefficient_),
          fluid_viscosity(fluid_viscosity_),
          biot_coefficient(biot_coefficient_),
          porosity(porosity_),
          solid_density(solid_density_),
          fluid_density(fluid_density_),
          specific_body_force(specific_body_force_)
    {
    }

    HydroMechanicsProcessData(HydroMechanicsProcessData&& other)
        : material{std::move(other.material)},
          intrinsic_permeability(other.intrinsic_permeability),
          storage_coefficient(other.storage_coefficient),
          fluid_viscosity(other.fluid_viscosity),
          biot_coefficient(other.biot_coefficient),
          porosity(other.porosity),
          solid_density(other.solid_density),
          fluid_density(other.fluid_density),
          specific_body_force(other.specific_body_force),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    HydroMechanicsProcessData(HydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicsProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& storage_coefficient;
    Parameter<double> const& fluid_viscosity;
    Parameter<double> const& biot_coefficient;
    Parameter<double> const& porosity;
    Parameter<double> const& solid_density;
    Parameter<double> const& fluid_density;
    Parameter<double> const& specific_body_force;
    double dt;
    double t;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib

#endif  // PROCESSLIB_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_
