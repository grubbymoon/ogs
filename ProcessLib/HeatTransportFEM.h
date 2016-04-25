/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_HEATTRANSPORT_FEM_H_
#define PROCESS_LIB_HEATTRANSPORT_FEM_H_

#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Parameter.h"
#include "ProcessLib/ProcessUtil.h"
#include "HeatTransportProcessData.h"

namespace ProcessLib
{

namespace HeatTransport
{

template <typename ShapeFunction,
         typename IntegrationMethod,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData : public ProcessLib::LocalAssemblerInterface<GlobalMatrix, GlobalVector>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

public:
    /// The thermal_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       HeatTransportProcessData const& process_data)
        : _element(element)
        , _shape_matrices(
              initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod, GlobalDim>(
                  element, integration_order))
        , _process_data(process_data)
        , _localK(local_matrix_size, local_matrix_size) // TODO narrowing conversion
        , _localM(local_matrix_size, local_matrix_size)
        , _integration_order(integration_order)
    {}

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/) override
    {
        _localK.setZero();
        _localM.setZero();
        const double density = 1000 ;
        const double heat_capacity = 3000 ;

        IntegrationMethod integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            auto const& sm = _shape_matrices[ip];

            auto const& wp = integration_method.getWeightedPoint(ip);
            auto const k = _process_data.thermal_conductivity(_element);

            _localK.noalias() += sm.dNdx.transpose() *
                                  k * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _localM.noalias() += sm.N *
                                  density*heat_capacity*sm.N.transpose() *
                                  sm.detJ * wp.getWeight();
        }
    }

    void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& /*b*/)
        const override
    {
        K.add(indices, _localK);
        M.add(indices, _localM);
    }

private:
    MeshLib::Element const& _element;
    std::vector<ShapeMatrices> _shape_matrices;
    HeatTransportProcessData const& _process_data;

    NodalMatrixType _localK;
    NodalMatrixType _localM;

    unsigned const _integration_order;
};


}   // namespace HeatTransport
}   // namespace ProcessLib

#endif  // PROCESS_LIB_HEATTRANSPORT_FEM_H_
