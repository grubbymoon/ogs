/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateFreezingProcess.h"

#include "FreezingProcess.h"
#include "FreezingProcessData.h"

namespace ProcessLib
{
namespace Freezing
{
std::unique_ptr<Process> createFreezingProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "FREEZING");

    DBUG("Create FreezingProcess.");

    // Process variable.
    auto process_variables =
        findProcessVariables(variables, config, {
            "temperature_variable" , "pressure_variable"});  // configure two Pcs

    // Thermal conductivity parameter.
    auto& thermal_conductivity = findParameter<double,
                                                 MeshLib::Element const&>(
        config,
        //! \ogs_file_param_special{process__Freezing__thermal_conductivity}
        "thermal_conductivity",
        parameters);

    DBUG("Use \'%s\' as thermal conductivity parameter.",
         thermal_conductivity.name.c_str());

    FreezingProcessData process_data{thermal_conductivity};

    SecondaryVariableCollection secondary_variables{
        //! \ogs_file_param{process__secondary_variables}
        config.getConfigSubtreeOptional("secondary_variables"),
        {
            /*//! \ogs_file_param_special{process__Freezing__secondary_variables__heat_flux_x}
         "heat_flux_x",
         //! \ogs_file_param_special{process__Freezing__secondary_variables__heat_flux_y}
         "heat_flux_y",
         //! \ogs_file_param_special{process__Freezing__secondary_variables__heat_flux_z}
         "heat_flux_z" */
         //! \ogs_file_param_special{process__Freezing__secondary_variables__darcy_velocity_x}
         "darcy_velocity_x",
         //! \ogs_file_param_special{process__Freezing__secondary_variables__darcy_velocity_y}
         "darcy_velocity_y",
         //! \ogs_file_param_special{process__Freezing__secondary_variables__darcy_velocity_z}
         "darcy_velocity_z"}};
    ProcessOutput
        //! \ogs_file_param{process__output}
        process_output{config.getConfigSubtree("output"), process_variables,
                       secondary_variables};

    return std::unique_ptr<Process>{new FreezingProcess{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(process_output)}};
}

}  // namespace Freezing
}  // namespace ProcessLib
