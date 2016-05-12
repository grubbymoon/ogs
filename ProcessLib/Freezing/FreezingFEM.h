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

#include <Eigen/Dense>
#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Parameter.h"
#include "ProcessLib/ProcessUtil.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "FreezingProcessData.h"
#include "FreezingMaterialModel.h"
using namespace Eigen;

namespace ProcessLib
{

namespace Freezing
{

template <typename ShapeFunction,
         typename IntegrationMethod,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData : public ProcessLib::LocalAssemblerInterface<GlobalMatrix, GlobalVector>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LAT = LocalAssemblerTraits<ShapeMatricesType, ShapeFunction::NPOINTS,
        2 /* number of pcs vars */, GlobalDim>;
    // Two Pcs variable T and P
    using NodalMatrixType = typename LAT::LocalMatrix;
    using NodalVectorType = typename LAT::LocalVector;

public:
    /// The thermal_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       FreezingProcessData const& process_data)
        : _element(element)
        , _shape_matrices(
              initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod, GlobalDim>(
                  element, integration_order))
        , _process_data(process_data)
        , _localK(local_matrix_size, local_matrix_size)
        , _localM(local_matrix_size, local_matrix_size)
        , _integration_order(integration_order)
   {}

    void assemble(double const /*t*/, std::vector<double> const& local_x) override
    {
        _localK.setZero();
        _localM.setZero();

        const double density_water = 1000.0 ;
        const double density_ice = 1000.0 ; // for the mass balance
        const double density_soil = 2000.0 ;
        const double specific_heat_capacity_soil = 1000.0 ;
        const double specific_heat_capacity_ice = 2000.0 ;
        const double specific_heat_capacity_water = 4200.0 ;
        const double thermal_conductivity_ice = 2 ;
        const double thermal_conductivity_soil = 3.2 ;
        const double thermal_conductivity_water = 0.6 ;
        const double specific_storage = 2e-4 ;       // m-1
        const double hydraulic_conductivity = 8e-7 ; // m/s
        double porosity = 0.5 ;
        double phi_i = 0.0 ;
        double sigmoid_coeff = 5.0 ;
        double latent_heat = 334000 ;
        double sigmoid_derive = 0.0 ;
        double thermal_conductivity = 0.0 ;
        double heat_capacity = 0.0 ;
        double Real_hydraulic_conductivity = 0.0 ;

        IntegrationMethod integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();

        double T_int_pt = 0.0;
        double p_int_pt = 0.0;

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            auto const num_nodes = ShapeFunction::NPOINTS;
            typedef Matrix<double, num_nodes, num_nodes> MatrixNN;
            MatrixNN _Ktt;
            MatrixNN _Mtt;
            MatrixNN _Kpp;
            MatrixNN _Mpp;
            MatrixNN _Ktp;
            MatrixNN _Mtp;
            MatrixNN _Kpt;
            MatrixNN _Mpt;
            auto const& sm = _shape_matrices[ip];
            int n = n_integration_points;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T_int_pt, p_int_pt);

            // use T_int_pt here ...

            phi_i = CalcIceVolFrac(T_int_pt, sigmoid_coeff, porosity);

          Real_hydraulic_conductivity = KozenyKarman(hydraulic_conductivity, porosity, phi_i);
          //    Real_hydraulic_conductivity = hydraulic_conductivity ;
            sigmoid_derive = Calcsigmoidderive(phi_i, sigmoid_coeff, porosity);

            thermal_conductivity = TotalThermalConductivity(porosity, phi_i, thermal_conductivity_ice,
thermal_conductivity_soil, thermal_conductivity_water);

            heat_capacity = EquaHeatCapacity(phi_i, density_water, density_soil, density_ice,
 specific_heat_capacity_soil, specific_heat_capacity_ice,
 specific_heat_capacity_water, porosity, sigmoid_derive, latent_heat);

            auto const p_nodal_values =
                    Eigen::Map<const Eigen::VectorXd>(&local_x[num_nodes], num_nodes);

            auto const& wp = integration_method.getWeightedPoint(ip);
            _Ktt.noalias() += sm.dNdx.transpose() *
                                  thermal_conductivity * sm.dNdx *
                                  sm.detJ * wp.getWeight() + sm.N*density_water
                                  *specific_heat_capacity_water*Real_hydraulic_conductivity*
                                  (sm.dNdx*p_nodal_values).transpose()*sm.dNdx*sm.detJ*wp.getWeight();
            _Kpp.noalias() += sm.dNdx.transpose() *
                                  Real_hydraulic_conductivity * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Kpt.noalias() += sm.dNdx.transpose() *
                                  0 * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Ktp.noalias() += sm.dNdx.transpose() *
                                  0 * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Mtt.noalias() += sm.N *heat_capacity*sm.N.transpose() *
                                  sm.detJ * wp.getWeight();
            _Mpp.noalias() += sm.N *specific_storage*sm.N.transpose() *
                                  sm.detJ * wp.getWeight();
            _Mpt.noalias() += sm.N *0*sm.N.transpose() *
                                  sm.detJ * wp.getWeight();
            _Mtp.noalias() += sm.N *0*sm.N.transpose() *
                                  sm.detJ * wp.getWeight();


            _localK.block<num_nodes,num_nodes>(0,0).noalias() += _Ktt;
            _localM.block<num_nodes,num_nodes>(0,0).noalias() += _Mtt;
            _localK.block<num_nodes,num_nodes>(num_nodes,num_nodes).noalias() += _Kpp;
            _localM.block<num_nodes,num_nodes>(num_nodes,num_nodes).noalias() += _Mpp;
            _localK.block<num_nodes,num_nodes>(num_nodes,0).noalias() += _Ktp;
            _localM.block<num_nodes,num_nodes>(num_nodes,0).noalias() += _Mtp;
            _localK.block<num_nodes,num_nodes>(0,num_nodes).noalias() += _Kpt;
            _localM.block<num_nodes,num_nodes>(0,num_nodes).noalias() += _Mpt;
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
    FreezingProcessData const& _process_data;

    NodalMatrixType _localK;
    NodalMatrixType _localM;

    unsigned const _integration_order;
};


}   // namespace Freezing
}   // namespace ProcessLib

#endif  // PROCESS_LIB_FREEZING_FEM_H_
