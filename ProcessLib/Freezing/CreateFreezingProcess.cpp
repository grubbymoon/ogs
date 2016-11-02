/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateFreezingProcess.h"

#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "FreezingProcess.h"
#include "FreezingProcessData.h"

namespace ProcessLib
{
namespace Freezing
{
std::unique_ptr<Process> createFreezingProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "FREEZING");

    DBUG("Create FreezingProcess.");

    // Process variable.
    auto process_variables =
        findProcessVariables(variables, config, {
            "temperature_variable" , "pressure_variable"});  // configure two Pcs

        // Porosity parameter.
        auto& porosity = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__porosity}
            "porosity", parameters, 1);
        DBUG("Use \'%s\' as porosity parameter.", porosity.name.c_str());

        // Parameter for the intrinsic permeability.
        auto& intrinsic_permeability = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__intrinsic_permeability}
            "intrinsic_permeability", parameters, 1);
        DBUG("Use \'%s\' as intrinsic_permeability parameter.",
             intrinsic_permeability.name.c_str());

 /*       // Parameter for the storage.
        auto& storage = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__storage}
            "storage", parameters, 1);
        DBUG("Use \'%s\' as storage parameter.", storage.name.c_str());  */

        // Parameter for the reference_temperature.
        auto& reference_temperature = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__reference_temperature}
            "reference_temperature", parameters, 1);
        DBUG("Use \'%s\' as reference_temperature parameter.",
             reference_temperature.name.c_str());

        // Parameter for the viscosity.
        auto& viscosity = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__viscosity}
            "viscosity", parameters, 1);
        DBUG("Use \'%s\' as viscosity parameter.", viscosity.name.c_str());

        // Parameter for the density of the solid.
        auto& density_solid = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__density_solid}
            "density_solid", parameters, 1);
        DBUG("Use \'%s\' as density_solid parameter.",
             density_solid.name.c_str());

        // Parameter for the density of the fluid.
        auto& density_fluid = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__density_fluid}
            "density_fluid", parameters, 1);
        DBUG("Use \'%s\' as density_fluid parameter.",
             density_fluid.name.c_str());

        // Parameter for the density of the ice.
        auto& density_ice = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__density_ice}
            "density_ice", parameters, 1);
        DBUG("Use \'%s\' as density_ice parameter.",
             density_ice.name.c_str());

        // Parameter for the specific heat capacity of the solid.
        auto& specific_heat_capacity_solid = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__specific_heat_capacity_solid}
            "specific_heat_capacity_solid", parameters, 1);
        DBUG("Use \'%s\' as specific_heat_capacity_solid parameter.",
             specific_heat_capacity_solid.name.c_str());

        // Parameter for the specific heat capacity of the fluid.
        auto& specific_heat_capacity_fluid = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__specific_heat_capacity_fluid}
            "specific_heat_capacity_fluid", parameters, 1);
        DBUG("Use \'%s\' as specific_heat_capacity_fluid parameter.",
             specific_heat_capacity_fluid.name.c_str());


        // Parameter for the specific heat capacity of the ice.
        auto& specific_heat_capacity_ice = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__specific_heat_capacity_ice}
            "specific_heat_capacity_ice", parameters, 1);
        DBUG("Use \'%s\' as specific_heat_capacity_ice parameter.",
             specific_heat_capacity_ice.name.c_str());

        // Parameter for the thermal conductivity of the solid.
        auto& thermal_conductivity_solid = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__thermal_conductivity_solid}
            "thermal_conductivity_solid", parameters, 1);
        DBUG("Use \'%s\' as thermal_conductivity_solid parameter.",
             thermal_conductivity_solid.name.c_str());

        // Parameter for the thermal conductivity of the fluid.
        auto& thermal_conductivity_fluid = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__thermal_conductivity_fluid}
            "thermal_conductivity_fluid", parameters, 1);
        DBUG("Use \'%s\' as thermal_conductivity_fluid parameter.",
             thermal_conductivity_fluid.name.c_str());

        // Parameter for the thermal conductivity of the fluid.
        auto& thermal_conductivity_ice = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__thermal_conductivity_ice}
            "thermal_conductivity_ice", parameters, 1);
        DBUG("Use \'%s\' as thermal_conductivity_ice parameter.",
             thermal_conductivity_ice.name.c_str());

        // Parameter for the thermal expansion coefficient.
        auto& thermal_expansion_coefficient = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__thermal_expansion_coefficient}
            "thermal_expansion_coefficient", parameters, 1);
        DBUG("Use \'%s\' as thermal_expansion_coefficient parameter.",
             thermal_expansion_coefficient.name.c_str());

        // Parameter for the sigmoid coefficient.
        auto& sigmoid_coefficient = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__sigmoid_coefficient}
            "sigmoid_coefficient", parameters, 1);
        DBUG("Use \'%s\' as sigmoid_coefficient parameter.",
             sigmoid_coefficient.name.c_str());

        // Parameter for the latent heat.
        auto& latent_heat = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__latent_heat}
            "latent_heat", parameters, 1);
        DBUG("Use \'%s\' as latent heat parameter.",
             latent_heat.name.c_str());

        // Parameter for the compressibility.
        auto& compressibility = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__compressibility}
            "compressibility", parameters, 1);
        DBUG("Use \'%s\' as compressibility parameter.",
             compressibility.name.c_str());

        // Specific body force parameter.
        auto& specific_body_force = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__Freezing__specific_body_force}
            "specific_body_force", parameters, mesh.getDimension());
        DBUG("Use \'%s\' as specific body force parameter.",
             specific_body_force.name.c_str());

        // Assume constant parameter, then check the norm at arbitrary
        // SpatialPosition and time.
        assert(dynamic_cast<ConstantParameter<double>*>(&specific_body_force));
        bool const has_gravity =
    MathLib::toVector(specific_body_force(0, SpatialPosition{})).norm() > 0;


    FreezingProcessData process_data{porosity,
                                     intrinsic_permeability,
                                //     storage,
                                     viscosity,
                                     density_solid,
                                     density_fluid,
                                     density_ice,
                                     specific_heat_capacity_solid,
                                     specific_heat_capacity_fluid,
                                     specific_heat_capacity_ice,
                                     thermal_conductivity_solid,
                                     thermal_conductivity_fluid,
                                     thermal_conductivity_ice,
                                     thermal_expansion_coefficient,
                                     sigmoid_coefficient,
                                     latent_heat,
                                     compressibility,
                                     reference_temperature,
                                     specific_body_force,
                                     has_gravity};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
            {"Freezing_temperature_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<Process>{new FreezingProcess{
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller)}};
}

}  // namespace Freezing
}  // namespace ProcessLib
