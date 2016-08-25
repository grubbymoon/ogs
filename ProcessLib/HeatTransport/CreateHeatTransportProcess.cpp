/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHeatTransportProcess.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
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
    std::unique_ptr<NumLib::ConvergenceCriterion>&& convergence_criterion,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "HEAT_TRANSPORT");

    DBUG("Create HeatTransportProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//!
         //\ogs_file_param_special{process__HEAT_TRANSPORT__process_variables__process_variable}
         "process_variable"});

    // thermal conductivity parameter.
    auto& thermal_conductivity = findParameter<double, MeshLib::Element const&>(
        config,
        //!
        //\ogs_file_param_special{process__HEAT_TRANSPORT__thermal_conductivity}
        "thermal_conductivity",
        parameters);

    DBUG("Use \'%s\' as thermal conductivity parameter.",
         thermal_conductivity.name.c_str());

    // heat capacity parameter.
    auto& heat_capacity = findParameter<double, MeshLib::Element const&>(
        config,
        //! \ogs_file_param_special{process__HEAT_TRANSPORT__heat_capacity}
        "heat_capacity",
        parameters);

    DBUG("Use \'%s\' as heat capacity parameter.", heat_capacity.name.c_str());

    // density parameter.
    auto& density = findParameter<double, MeshLib::Element const&>(
        config,
        //! \ogs_file_param_special{process__HEAT_TRANSPORT__density}
        "density",
        parameters);

    DBUG("Use \'%s\' as density parameter.", density.name.c_str());

    HeatTransportProcessData process_data{thermal_conductivity, heat_capacity,
                                          density};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HeatConduction_temperature"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    //! \ogs_file_param{process__output}
    ProcessOutput process_output{config.getConfigSubtree("output")};

    return std::unique_ptr<Process>{new HeatTransportProcess{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(convergence_criterion), std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(process_output), std::move(named_function_caller)}};
}

}  // namespace HeatTransport
}  // namespace ProcessLib
