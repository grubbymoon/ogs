/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_FREEZINGPROCESS_FEM_H_
#define PROCESS_LIB_FREEZINGPROCESS_FEM_H_

#include <memory>
#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"

#include "Parameter.h"
#include "ProcessUtil.h"


namespace ProcessLib
{

namespace Freezing
{

template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void assemble(double const t, std::vector<double> const& local_x) = 0;

    virtual void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    /// The thermal_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(
            MeshLib::Element const& e,
            std::size_t const local_matrix_size,
            unsigned const integration_order,
            Parameter<double, MeshLib::Element const&> const&
                thermal_conductivity
            /* ... */)
    {
        _integration_order = integration_order;

        _shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod_, GlobalDim>(
                e, integration_order);

        _thermal_conductivity = [&thermal_conductivity, &e]()
        {
            return thermal_conductivity(e);
        };

        _localK.reset(new NodalMatrixType(local_matrix_size, local_matrix_size));
        _localM.reset(new NodalMatrixType(local_matrix_size, local_matrix_size));
    }

    void assemble(double const /*t*/, std::vector<double> const& local_x) override
    {
        _localK->setZero();
        _localM->setZero();
        const double density = 1000 ;
        const double heat_capacity = 3000 ;

        IntegrationMethod_ integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();

        double T_int_pt = 0.0;
        std::array<double*, 1> const int_pt_array = { &T_int_pt };

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            auto const& sm = _shape_matrices[ip];

            NumLib::shapeFunctionInterpolate(local_x, sm.N, int_pt_array);

            // use T_int_pt here ...

            auto const& wp = integration_method.getWeightedPoint(ip);
            _localK->noalias() += sm.dNdx.transpose() *
                                  _thermal_conductivity() * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _localM->noalias() += sm.N *
                                  density*heat_capacity*sm.N.transpose() *
                                  sm.detJ * wp.getWeight();
        }
    }

    void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& /*b*/)
        const override
    {
        K.add(indices, *_localK);
        M.add(indices, *_localM);
    }

private:
    std::vector<ShapeMatrices> _shape_matrices;
    std::function<double(void)> _thermal_conductivity;

    std::unique_ptr<NodalMatrixType> _localK;
    std::unique_ptr<NodalMatrixType> _localM;

    unsigned _integration_order = 2;
};


}   // namespace Freezing
}   // namespace ProcessLib

#endif  // PROCESS_LIB_FREEZING_FEM_H_
