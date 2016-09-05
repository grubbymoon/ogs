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

#include <iostream>
#include <Eigen/Dense>
#include <vector>

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "FreezingProcessData.h"
#include "FreezingMaterialModel.h"
using namespace Eigen;

namespace ProcessLib
{

namespace Freezing
{

const unsigned NUM_NODAL_DOF = 2;

class FreezingLocalAssemblerInterface
        : public ProcessLib::LocalAssemblerInterface
        , public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtHeatFluxX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxZ(
std::vector<double>& /*cache*/) const = 0;
};


template <typename ShapeFunction,
         typename IntegrationMethod,
         unsigned GlobalDim>
class LocalAssemblerData
        : public FreezingLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
    ShapeMatricesType, ShapeFunction::NPOINTS,
        NUM_NODAL_DOF, /* number of pcs vars */ GlobalDim>;
    // Two Pcs variable T and P
    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;

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
        , _localRhs(local_matrix_size)
        , _integration_order(integration_order)
   {
        // This assertion is valid only if all nodal d.o.f. use the same shape matrices.
       assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
   }

    /* commit t? */
    void assembleConcrete(double const /*t*/, std::vector<double> const& local_x,
                  NumLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
                  GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override
    {
        _localK.setZero();
        _localM.setZero();
        _localRhs.setZero();

        const double density_water = 1000.0 ;
        const double temperature0 = 273.15 ;
        const double beta = 2.14e-4 ;  // thermal expansion coefficient
        const double density_ice = 1000.0 ; // for the mass balance
        const double density_soil = 2000.0 ;
        const double specific_heat_capacity_soil = 250.0 ;
        const double specific_heat_capacity_ice = 2000.0 ;
        const double specific_heat_capacity_water = 1100.0 ;
        const double thermal_conductivity_ice = 2.0 ;
        const double thermal_conductivity_soil = 5.0 ;
        const double thermal_conductivity_water = 1.0 ;
        const double specific_storage = 0 ;       // m-1 if 0 then constant velocity
        const double hydraulic_conductivity = 5.787e-12 ; // m/s permeability/viscosity
        const double g = 9.81 ;
        double porosity = 0.1 ;
        double phi_i = 0.0 ;
        double sigmoid_coeff = 5.0 ;
        double latent_heat = 334000 ;
        double sigmoid_derive = 0.0 ;
        double thermal_conductivity = 0.0 ;
        double heat_capacity = 0.0 ;
        double Real_hydraulic_conductivity = 0.0 ;
        //_gravity_v(element.getDimension()) = 0.;
        //_gravity_v[(element.getDimension() - 1)] = 1.;


        IntegrationMethod integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNumberOfPoints();

        double T_int_pt = 0.0;
        double p_int_pt = 0.0;

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            auto const num_nodes = ShapeFunction::NPOINTS; /* difference with n_integration_points?*/
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
           // auto const k = _process_data.thermal_conductivity(_element);
            typedef Matrix<double, num_nodes, num_nodes> MatrixNN;
            MatrixNN _Ktt;
            MatrixNN _Mtt;
            MatrixNN _Kpp;
            MatrixNN _Mpp;
            MatrixNN _Ktp;
            MatrixNN _Mtp;
            MatrixNN _Kpt;
            MatrixNN _Mpt;
            typedef Matrix<double, num_nodes, 1> VecterNN;
            VecterNN _Bpp;

            int n = n_integration_points;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T_int_pt, p_int_pt);

            // use T_int_pt here ...

            double density_water_T = DensityWater_T(density_water, T_int_pt, temperature0, beta);
            std::cout << density_water_T << std::endl;

            phi_i = CalcIceVolFrac(T_int_pt, sigmoid_coeff, porosity);

          Real_hydraulic_conductivity = KozenyKarman(hydraulic_conductivity, porosity, phi_i);
          //    Real_hydraulic_conductivity = hydraulic_conductivity ;
            sigmoid_derive = Calcsigmoidderive(phi_i, sigmoid_coeff, porosity);

          /*  thermal_conductivity = TotalThermalConductivity(porosity, phi_i, thermal_conductivity_ice,
thermal_conductivity_soil, thermal_conductivity_water); */
            thermal_conductivity = thermal_conductivity_soil*(1 - porosity)
                    + thermal_conductivity_water*porosity; // mixed thermal conductivity

           /* heat_capacity = EquaHeatCapacity(phi_i, density_water, density_soil, density_ice,
 specific_heat_capacity_soil, specific_heat_capacity_ice,
 specific_heat_capacity_water, porosity, sigmoid_derive, latent_heat); */
            heat_capacity = density_soil*specific_heat_capacity_soil*(1 - porosity)
                    + density_water_T*specific_heat_capacity_water*porosity;  // mixed heat capacity
            auto const detJ_w_NT = (sm.detJ * wp.getWeight() * sm.N.transpose()).eval();
            auto const p_nodal_values =
                    Eigen::Map<const Eigen::VectorXd>(&local_x[num_nodes], num_nodes);
            auto const velocity =  (-hydraulic_conductivity*
                    sm.dNdx*p_nodal_values).eval();
            //std::cout << velocity << std::endl;
            auto const detJ_w_NT_vT_dNdx =
                (detJ_w_NT * velocity.transpose() * sm.dNdx).eval();
            // matrix assembly
            _Ktt.noalias() += sm.dNdx.transpose() *
                                  thermal_conductivity * sm.dNdx *
                                  sm.detJ * wp.getWeight() + detJ_w_NT_vT_dNdx*density_water_T*specific_heat_capacity_water;
 /*sm.N*density_water*specific_heat_capacity_water*Real_hydraulic_conductivity*
            (sm.dNdx*p_nodal_values).transpose()*sm.dNdx*sm.detJ*wp.getWeight(); */
            _Kpp.noalias() += sm.dNdx.transpose() *
                                 /* Real_hydraulic_conductivity */ hydraulic_conductivity* sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Kpt.noalias() += sm.dNdx.transpose() *
                                  0 * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Ktp.noalias() += sm.dNdx.transpose() *
                                  0 * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Mtt.noalias() += sm.N.transpose() *heat_capacity*sm.N *
                                  sm.detJ * wp.getWeight();
            _Mpp.noalias() += sm.N.transpose() *specific_storage*sm.N *
                                  sm.detJ * wp.getWeight();
            _Mpt.noalias() += sm.N.transpose() *0*sm.N *
                                  sm.detJ * wp.getWeight();
            _Mtp.noalias() += sm.N.transpose() *(-beta)*sm.N *
                                  sm.detJ * wp.getWeight();
            _Bpp.noalias() += hydraulic_conductivity* sm.detJ * wp.getWeight()*sm.dNdx.transpose().col(_element.getDimension()-1)*g*density_water_T;


            _localK.block<num_nodes,num_nodes>(0,0).noalias() += _Ktt;
            _localM.block<num_nodes,num_nodes>(0,0).noalias() += _Mtt;
            _localK.block<num_nodes,num_nodes>(num_nodes,num_nodes).noalias() += _Kpp;
            _localM.block<num_nodes,num_nodes>(num_nodes,num_nodes).noalias() += _Mpp;
            _localK.block<num_nodes,num_nodes>(num_nodes,0).noalias() += _Ktp;
            _localM.block<num_nodes,num_nodes>(num_nodes,0).noalias() += _Mtp;
            _localK.block<num_nodes,num_nodes>(0,num_nodes).noalias() += _Kpt;
            _localM.block<num_nodes,num_nodes>(0,num_nodes).noalias() += _Mpt;
            _localRhs.block<num_nodes,1>(num_nodes,0).noalias() -= _Bpp  ;
            // heat flux only computed for output.
            auto const heat_flux = (-thermal_conductivity * sm.dNdx *
                Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes)
                ).eval();

            for (unsigned d=0; d<GlobalDim; ++d) {
                _heat_fluxes[d][ip] = heat_flux[d];
        }
    }


        K.add(indices, _localK);
        M.add(indices, _localM);
        b.add(indices.rows, _localRhs);
  }

    Eigen::Map<const Eigen::RowVectorXd>
    getShapeMatrix(const unsigned integration_point) const override
    {
       auto const& N = _shape_matrices[integration_point].N;

       // assumes N is stored contiguously in memory
       return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

   std::vector<double> const&
   getIntPtHeatFluxX(std::vector<double>& /*cache*/) const override
   {
       assert(_heat_fluxes.size() > 0);
       return _heat_fluxes[0];
   }

   std::vector<double> const&
   getIntPtHeatFluxY(std::vector<double>& /*cache*/) const override
   {
       assert(_heat_fluxes.size() > 1);
       return _heat_fluxes[1];
   }

   std::vector<double> const&
   getIntPtHeatFluxZ(std::vector<double>& /*cache*/) const override
   {
       assert(_heat_fluxes.size() > 2);
       return _heat_fluxes[2];
}

private:
    MeshLib::Element const& _element;
    std::vector<ShapeMatrices> _shape_matrices;
    FreezingProcessData const& _process_data;

    NodalMatrixType _localK;
    NodalMatrixType _localM;
    NodalVectorType _localRhs;
    NodalVectorType _gravity_v;

    unsigned const _integration_order;
    std::vector<std::vector<double>> _heat_fluxes
            = std::vector<std::vector<double>>(
    GlobalDim, std::vector<double>(ShapeFunction::NPOINTS));
};


}   // namespace Freezing
}   // namespace ProcessLib

#endif  // PROCESS_LIB_FREEZING_FEM_H_
