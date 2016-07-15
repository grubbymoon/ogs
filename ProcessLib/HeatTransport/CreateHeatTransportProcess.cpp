/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHeatTransportProcess.h"

#include "HeatTransportProcess.h"
#include "HeatTransportProcessData.h"

namespace ProcessLib
{
namespace HeatTransport
{
std::unique_ptr<Process> createHeatTransportProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "GROUNDWATER_FLOW");

    DBUG("Create HeatTransportProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//! \ogs_file_param_special{process__GROUNDWATER_FLOW__process_variables__process_variable}
         "process_variable"});

    // thermal conductivity parameter.
    auto& thermal_conductivity = findParameter<double,
                                                 MeshLib::Element const&>(
        config,
        //! \ogs_file_param_special{process__GROUNDWATER_FLOW__thermal_conductivity}
        "thermal_conductivity",
        parameters);

    DBUG("Use \'%s\' as thermal conductivity parameter.",
         thermal_conductivity.name.c_str());

    HeatTransportProcessData process_data{thermal_conductivity};

    SecondaryVariableCollection secondary_variables{
        //! \ogs_file_param{process__secondary_variables}
        config.getConfigSubtreeOptional("secondary_variables"),
        {//! \ogs_file_param_special{process__HeatTransport__secondary_variables__heat_flux_x}
         "heat_flux_x",
         //! \ogs_file_param_special{process__HeatTransport__secondary_variables__heat_flux_y}
         "heat_flux_y",
         //! \ogs_file_param_special{process__HeatTransport__secondary_variables__heat_flux_z}
         "heat_flux_z"}};

    ProcessOutput
        //! \ogs_file_param{process__output}
        process_output{config.getConfigSubtree("output"), process_variables,
                       secondary_variables};

    return std::unique_ptr<Process>{new HeatTransportProcess{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(process_output)}};
}

}  // namespace HeatTransport
}  // namespace ProcessLib
