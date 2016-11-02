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
#include "ProcessLib/Parameter/Parameter.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "FreezingProcessData.h"
#include "FreezingMaterialModel.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
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
    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
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
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    /// The thermal_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       FreezingProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _darcy_velocities(
              GlobalDim,
              std::vector<double>(_integration_method.getNumberOfPoints()))
   {
        // This assertion is valid only if all nodal d.o.f. use the same shape matrices.
       assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
       (void)local_matrix_size;
   }

    /* commit t? */
    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
               // This assertion is valid only if all nodal d.o.f. use the same shape
               // matrices.
               assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

               auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
                   local_M_data, local_matrix_size, local_matrix_size);
               auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
                   local_K_data, local_matrix_size, local_matrix_size);
               auto local_b = MathLib::createZeroedVector<NodalVectorType>(
                   local_b_data, local_matrix_size);

   /*     const double density_water = 1000.0 ;
        const double density_water_T = 1000.0 ;  // thermal expansion is ignored
        //  const double temperature0 = 20.0 ;   // reference temperature for thermal expansion
        const double beta = 4.3e-4 ;  // thermal expansion coefficient
        const double density_ice = 1000.0 ; // for the mass balance
        const double density_soil = 2650.0;
        const double specific_heat_capacity_soil = 835;
        const double specific_heat_capacity_ice = 2060.0 ;
        const double specific_heat_capacity_water = 4182.0 ;
        const double thermal_conductivity_ice = 2.14 ;
        const double thermal_conductivity_soil = 9.0 ;
        const double thermal_conductivity_water = 0.6 ;
        double specific_storage = 0.0 ;       // m-1 if 0 then no stoarge term
        const double permeability = 1.3e-10 ;
        const double viscosity = 1.793e-3 ;
        // const double temperature_con = 75 ;    // reference temperature for viscosity
        const double g = 0 ;  // if 0 then do not consider gravity
        double porosity = 0.37 ; // 0.001 is easier for thermal convection converge
        double phi_i = 0.0 ;
        double sigmoid_coeff = 5.0 ;
        double latent_heat = 334000 ;
        double sigmoid_derive = 0.0 ;
        double thermal_conductivity = 0.0 ;
        double heat_capacity = 0.0 ;
        double Real_hydraulic_conductivity = 0.0 ;
        double compressibility = 1e-8 ;
    //    typedef Matrix<double, 2, 1> Vecterg;
     //   Vecterg  _gravity_v;
     //   _gravity_v[0] = 0.;
     //   _gravity_v[1] = 1.;   */


      //  IntegrationMethod integration_method(_integration_order);
        unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const temperature0 =
                    _process_data.reference_temperature(t, pos)[0];  // 20.0;
        auto const beta =
                    _process_data.thermal_expansion_coefficient(t, pos)[0];

        auto const sigmoid_coefficient =
                    _process_data.sigmoid_coefficient(t, pos)[0];

        auto const latent_heat =
                    _process_data.latent_heat(t, pos)[0];

        auto const body_force = _process_data.specific_body_force(t, pos);
        assert(body_force.size() == GlobalDim);
        auto const b = MathLib::toVector<GlobalDimVectorType>(body_force, GlobalDim);

        double T_int_pt = 0.0;
        double p_int_pt = 0.0;

        for (unsigned ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const density_fluid = _process_data.density_fluid(t, pos)[0];
            auto const density_solid = _process_data.density_solid(t, pos)[0];
            auto const density_ice = _process_data.density_ice(t, pos)[0];
            auto const intrinsic_permeability =
                _process_data.intrinsic_permeability(t, pos)[0];
            auto const viscosity0 = _process_data.viscosity(t, pos)[0];  // 1e-3

            auto const porosity = _process_data.porosity(t, pos)[0];
            auto const compressibility = _process_data.compressibility(t, pos)[0];

            auto const specific_heat_capacity_solid =
                _process_data.specific_heat_capacity_solid(t, pos)[0];
            auto const specific_heat_capacity_fluid =
                _process_data.specific_heat_capacity_fluid(t, pos)[0];
            auto const specific_heat_capacity_ice =
                _process_data.specific_heat_capacity_ice(t, pos)[0];
            auto const thermal_conductivity_solid =
                _process_data.thermal_conductivity_solid(t, pos)[0];
            auto const thermal_conductivity_fluid =
                _process_data.thermal_conductivity_fluid(t, pos)[0];
            auto const thermal_conductivity_ice =
                _process_data.thermal_conductivity_ice(t, pos)[0];


            auto const num_nodes = ShapeFunction::NPOINTS; /* difference with n_integration_points?*/
            auto const& sm = _shape_matrices[ip];
            auto const& wp = _integration_method.getWeightedPoint(ip);
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
            VecterNN _Bp;

            // int n = n_integration_points;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, sm.N, T_int_pt, p_int_pt);

            // use T_int_pt here ...

          //  double density_water_T =  DensityWater_T(density_water, T_int_pt, temperature0, beta);
          //  double density_water_T = density_water*(1 - beta*(T_int_pt - temperature0));
          //  std::cout << density_water_T << std::endl;
          //  double viscosity = Viscosity(viscosity0, T_int_pt, temperature_con, temperature0); // calculate the dynamic viscosity
            double hydraulic_conductivity = intrinsic_permeability/viscosity0 ; // m/s permeability/viscosity
            double phi_i = CalcIceVolFrac(T_int_pt, sigmoid_coefficient, porosity);

          //  Real_hydraulic_conductivity = KozenyKarman(hydraulic_conductivity, porosity, phi_i);
          //    Real_hydraulic_conductivity = hydraulic_conductivity ;
            double Kr = Relative_permeability(phi_i, porosity);
            double Real_hydraulic_conductivity = Hydraulic_conductivity(intrinsic_permeability, Kr, viscosity0);
            double sigmoid_derive = Calcsigmoidderive(sigmoid_coefficient, porosity, T_int_pt);
            double specific_storage = (porosity -phi_i)*compressibility ;

            double thermal_conductivity = TotalThermalConductivity(porosity, phi_i, thermal_conductivity_ice,
                   thermal_conductivity_solid, thermal_conductivity_fluid);
          //  thermal_conductivity = thermal_conductivity_soil*(1 - porosity)
           //         + thermal_conductivity_water*porosity; // mixed thermal conductivity

            double heat_capacity = EquaHeatCapacity(phi_i, density_fluid, density_solid, density_ice,
                   specific_heat_capacity_solid, specific_heat_capacity_ice,
                   specific_heat_capacity_fluid, porosity, sigmoid_derive, latent_heat);

           /* heat_capacity = density_soil*specific_heat_capacity_soil*(1 - porosity)
                    + density_water*specific_heat_capacity_water*porosity;*/  // mixed heat capacity
            auto const detJ_w_NT = (sm.detJ * wp.getWeight() * sm.N.transpose()).eval();
            auto p_nodal_values =
                    Eigen::Map<const Eigen::VectorXd>(&local_x[num_nodes], num_nodes);

            auto velocity =  (-Real_hydraulic_conductivity*
                    sm.dNdx*p_nodal_values /*- g*_gravity_v*density_water_T*hydraulic_conductivity*/).eval();

            if (_process_data.has_gravity)
            velocity += density_fluid*Real_hydraulic_conductivity*b ;
            //std::cout << velocity << std::endl;
            auto detJ_w_NT_vT_dNdx =
                (detJ_w_NT * velocity.transpose() * sm.dNdx).eval();
            // matrix assembly
            _Ktt.noalias() += sm.dNdx.transpose() *
                                  thermal_conductivity * sm.dNdx *
                                  sm.detJ * wp.getWeight() + detJ_w_NT_vT_dNdx*density_fluid*specific_heat_capacity_fluid;
 /*sm.N*density_water*specific_heat_capacity_water*Real_hydraulic_conductivity*
            (sm.dNdx*p_nodal_values).transpose()*sm.dNdx*sm.detJ*wp.getWeight(); */
            _Kpp.noalias() += sm.dNdx.transpose() *
                                 Real_hydraulic_conductivity  /*hydraulic_conductivity*/* sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Kpt.noalias() += sm.dNdx.transpose() *
                                  0.0 * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Ktp.noalias() += sm.dNdx.transpose() *
                                  0.0 * sm.dNdx *
                                  sm.detJ * wp.getWeight();
            _Mtt.noalias() += sm.N.transpose() *heat_capacity*sm.N *
                                  sm.detJ * wp.getWeight();
            _Mpp.noalias() += sm.N.transpose() *specific_storage*sm.N *
                                  sm.detJ * wp.getWeight();
            _Mpt.noalias() += sm.N.transpose() *0.0*sm.N *
                                  sm.detJ * wp.getWeight();
            _Mtp.noalias() += sm.N.transpose() *(-beta*porosity*0)*sm.N *    // beta or -beta or phi*beta or phi*(-beta)
                                  sm.detJ * wp.getWeight();
        //    _Bpp.noalias() += hydraulic_conductivity* sm.detJ * wp.getWeight()*sm.dNdx.transpose().col(_element.getDimension()-1)*g*density_water_T;
            if (_process_data.has_gravity)
            _Bp.noalias() += hydraulic_conductivity* sm.detJ * wp.getWeight()*sm.dNdx.transpose().col(1)*b*density_fluid;
            /* with Oberbeck-Boussing assumption density difference only exists in buoyancy effects */

            local_K.template block<num_nodes, num_nodes>(0, 0).noalias() += _Ktt;
            local_M.template block<num_nodes, num_nodes>(0, 0).noalias() += _Mtt;
            local_K.template block<num_nodes, num_nodes>(num_nodes, num_nodes)
                .noalias() += _Kpp;
            local_M.template block<num_nodes, num_nodes>(num_nodes, num_nodes)
                .noalias() += _Mpp;
            local_K.template block<num_nodes, num_nodes>(num_nodes, 0).noalias() += _Ktp;
            local_M.template block<num_nodes, num_nodes>(num_nodes, 0).noalias() += _Mtp;
            local_K.template block<num_nodes, num_nodes>(0, num_nodes).noalias() += _Kpt;
            local_M.template block<num_nodes, num_nodes>(0, num_nodes).noalias() += _Mpt;
            local_b.template block<num_nodes, 1>(num_nodes, 0).noalias() -= _Bp;
            // heat flux only computed for output.
         /*   auto const heat_flux = (-thermal_conductivity * sm.dNdx *
                Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes)
                ).eval();

            for (unsigned d=0; d<GlobalDim; ++d) {
                _heat_fluxes[d][ip] = heat_flux[d];
        } */
           // velocity computed for output.
            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = velocity[d];
            }
    }
  }

    Eigen::Map<const Eigen::RowVectorXd>
    getShapeMatrix(const unsigned integration_point) const override
    {
       auto const& N = _shape_matrices[integration_point].N;

       // assumes N is stored contiguously in memory
       return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

   std::vector<double> const&
   getIntPtDarcyVelocityX(std::vector<double>& /*cache*/) const override
   {
       assert(_darcy_velocities.size() > 0);
       return _darcy_velocities[0];
   }

   std::vector<double> const&
   getIntPtDarcyVelocityY(std::vector<double>& /*cache*/) const override
   {
       assert(_darcy_velocities.size() > 1);
       return _darcy_velocities[1];
   }

   std::vector<double> const&
   getIntPtDarcyVelocityZ(std::vector<double>& /*cache*/) const override
   {
       assert(_darcy_velocities.size() > 2);
       return _darcy_velocities[2];
   }

private:
    MeshLib::Element const& _element;
    FreezingProcessData const& _process_data;

/*    NodalMatrixType _localK;
    NodalMatrixType _localM;
    NodalVectorType _localRhs;
    NodalVectorType _gravity_v; */

    //unsigned const _integration_order;
    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices> _shape_matrices;
 //   std::vector<std::vector<double>> _heat_fluxes
  //          = std::vector<std::vector<double>>(
 //   GlobalDim, std::vector<double>(ShapeFunction::NPOINTS));
    std::vector<std::vector<double>> _darcy_velocities ;
};


}   // namespace Freezing
}   // namespace ProcessLib

#endif  // PROCESS_LIB_FREEZING_FEM_H_
