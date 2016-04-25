/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_FREEZINGPROCESS_H_
#define PROCESS_LIB_FREEZINGPROCESS_H_

#include <cassert>

#include "AssemblerLib/VectorMatrixAssembler.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/ProcessUtil.h"

#include "FreezingFEM.h"
#include "FreezingProcessData.h"

namespace MeshLib
{
    class Element;
}

namespace ProcessLib
{

namespace Freezing
{

template<typename GlobalSetup>
class FreezingProcess final
        : public Process<GlobalSetup>
{

    using Base = Process<GlobalSetup>;
    using GlobalMatrix = typename GlobalSetup::MatrixType;
    using GlobalVector = typename GlobalSetup::VectorType;

public:
    FreezingProcess(
        MeshLib::Mesh& mesh,
        typename Base::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<typename Base::TimeDiscretization>&& time_discretization,
        std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
        FreezingProcessData&& process_data)
        : Process<GlobalSetup>(mesh, nonlinear_solver, std::move(time_discretization),
                      std::move(process_variables))
        , _process_data(std::move(process_data))
  {
        if (dynamic_cast<NumLib::ForwardEuler<GlobalVector>*>(
                    &Base::getTimeDiscretization()) != nullptr)
        {
            ERR("FreezingProcess can not be solved with the ForwardEuler"
                " time discretization scheme. Aborting");
            // Because the M matrix is not assembled. Thus, the linearized system
            // would be singular. The same applies to CrankNicolson with theta = 0.0,
            // but this case is not checked here.
            // Anyway, the FreezingProcess shall be transferred to a simpler
            // ODESystemTag in the future.
            std::abort();
        }
    }

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override
    {
        return false;
    }

    //! @}

private:
    using LocalAssemblerInterface =
        ProcessLib::LocalAssemblerInterface<GlobalMatrix, GlobalVector>;

    using GlobalAssembler = AssemblerLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector, LocalAssemblerInterface,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    void createAssemblers(AssemblerLib::LocalToGlobalIndexMap const& dof_table,
                          MeshLib::Mesh const& mesh,
                          unsigned const integration_order) override
    {
        DBUG("Create global assembler.");
        _global_assembler.reset(new GlobalAssembler(dof_table));

        ProcessLib::createLocalAssemblers<GlobalSetup, LocalAssemblerData>(
                    mesh.getDimension(), mesh.getElements(),
                    dof_table, integration_order, _local_assemblers,
                    _process_data);
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override
    {
        DBUG("Assemble FreezingProcess.");

        // Call global assembler for each local assembly item.
        GlobalSetup::executeMemberDereferenced(
            *_global_assembler, &GlobalAssembler::assemble,
            _local_assemblers, t, x, M, K, b);
    }


    FreezingProcessData _process_data;

    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;
};

template <typename GlobalSetup>
std::unique_ptr<FreezingProcess<GlobalSetup>>
createFreezingProcess(
    MeshLib::Mesh& mesh,
    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    config.checkConfParam("type", "FREEZING");

    DBUG("Create FreezingProcess.");

    // Process variable.
    auto process_variables =
        findProcessVariables(variables, config, {
            "temperature_variable" , "pressure_variable"});  // configure two Pcs

    // Thermal conductivity parameter.
    auto& thermal_conductivity =
            findParameter<double, MeshLib::Element const&>(
                config, "thermal_conductivity", parameters);

        DBUG("Use \'%s\' as thermal conductivity parameter.",
             thermal_conductivity.name.c_str());

        FreezingProcessData process_data {
            thermal_conductivity
        };

    return std::unique_ptr<FreezingProcess<GlobalSetup>>{
        new FreezingProcess<GlobalSetup>{
            mesh, nonlinear_solver,std::move(time_discretization),
            std::move(process_variables),
            std::move(process_data)
    }

  };
}

}   // namesapce Freezing
}   // namespace ProcessLib

#endif  // PROCESS_LIB_FREEZINGPROCESS_H_
