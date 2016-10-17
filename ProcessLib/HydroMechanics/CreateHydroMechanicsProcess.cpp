/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHydroMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"

#include "HydroMechanicsProcess.h"
#include "HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
class HydroMechanicsProcess;

extern template class HydroMechanicsProcess<2>;
extern template class HydroMechanicsProcess<3>;

template <int DisplacementDim>
std::unique_ptr<Process> createHydroMechanicsProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "HYDRO_MECHANICS");
    DBUG("Create HydroMechanicsProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//! \ogs_file_param_special{process__HYDRO_MECHANICS_process_variables__process_variable}
          "pressure", "displacement"});

    DBUG("Associate displacement with process variable \'%s\'.",
         process_variables[1].get().getName().c_str());

    if (process_variables[1].get().getNumberOfComponents() !=
        DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            process_variables[1].get().getName().c_str(),
            process_variables[1].get().getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate pressure with process variable \'%s\'.",
         process_variables[0].get().getName().c_str());
    if (process_variables[0].get().getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            process_variables[0].get().getName().c_str(),
            process_variables[0].get().getNumberOfComponents(),
            DisplacementDim);
    }


    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{process__HYDRO_MECHANICS_constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropic")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
                parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Intrinsic permeability
    auto& intrinsic_permeability = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HYDRO_MECHANICS_intrinsic_permeability}
        "intrinsic_permeability",
        parameters, 1);

    DBUG("Use \'%s\' as hydraulic conductivity parameter.",
         intrinsic_permeability.name.c_str());

    // Storage coefficient
    auto& storage_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HYDRO_MECHANICS_storage_coefficient}
        "storage_coefficient",
        parameters, 1);

    DBUG("Use \'%s\' as storage coefficient parameter.",
         storage_coefficient.name.c_str());


    // Fluid viscosity
    auto& fluid_viscosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HYDRO_MECHANICS_fluid_viscosity}
        "fluid_viscosity",
        parameters, 1);
    DBUG("Use \'%s\' as fluid viscosity parameter.",
         fluid_viscosity.name.c_str());

    // Biot coefficient
    auto& biot_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HYDRO_MECHANICS_biot_coefficient}
        "biot_coefficient",
        parameters, 1);
    DBUG("Use \'%s\' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Porosity
    auto& porosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HYDRO_MECHANICS_porosity}
        "porosity",
        parameters, 1);
    DBUG("Use \'%s\' as porosity parameter.",
         porosity.name.c_str());

    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HYDRO_MECHANICS_solid_density}
        "solid_density",
        parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.",
         solid_density.name.c_str());

    // Fluid density
    auto& fluid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HYDRO_MECHANICS_fluid_density}
        "fluid_density",
        parameters, 1);
    DBUG("Use \'%s\' as fluid density parameter.",
         fluid_density.name.c_str());

    // Specific body force
    auto& specific_body_force = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HYDRO_MECHANICS_specific_body_force}
        "specific_body_force",
        parameters, DisplacementDim);
    DBUG("Use \'%s\' as specific body force parameter.",
         specific_body_force.name.c_str());

    HydroMechanicsProcessData<DisplacementDim> process_data{
        std::move(material), intrinsic_permeability, storage_coefficient,
        fluid_viscosity,     biot_coefficient,       porosity,
        solid_density,       fluid_density,          specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HydroMechanics_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<HydroMechanicsProcess<DisplacementDim>>{
        new HydroMechanicsProcess<DisplacementDim>{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data),
            std::move(secondary_variables), std::move(named_function_caller)}};
}

template std::unique_ptr<Process> createHydroMechanicsProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createHydroMechanicsProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace HydroMechanics
}  // namespace ProcessLib
