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

#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"
#include "FreezingFEM.h"
#include "FreezingProcessData.h"


namespace ProcessLib
{
namespace Freezing
{
class FreezingProcess final : public Process
{

public:
    FreezingProcess(
        MeshLib::Mesh& mesh,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        FreezingProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    FreezingProcessData _process_data;

    std::vector<std::unique_ptr<FreezingLocalAssemblerInterface>>
_local_assemblers;
};

}   // namespace Freezing
}   // namespace ProcessLib

#endif // PROCESS_LIB_FreezingPROCESS_H_
