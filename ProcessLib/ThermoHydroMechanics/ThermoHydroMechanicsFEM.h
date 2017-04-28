/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// currently lambda, beta, K are considered to be constant

#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "FreezingMaterialModel.h"

#include "ThermoHydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointData(IntegrationPointData&& other)
        : b_matrices(std::move(other.b_matrices)),
          sigma_eff(std::move(other.sigma_eff)),
          sigma_eff_prev(std::move(other.sigma_eff_prev)),
          eps(std::move(other.eps)),
          eps_prev(std::move(other.eps_prev)),
          solid_material(other.solid_material),
          material_state_variables(std::move(other.material_state_variables)),
          //C_solid(std::move.C_solid)
          //C_ice(std::move.C_ice)
          C(std::move(other.C)),
          integration_weight(std::move(other.integration_weight))
    {
    }
#endif  // _MSC_VER

    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        N_u;
    typename BMatricesType::BMatrixType b_matrices;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    typename ShapeMatricesTypePressure::NodalRowVectorType _N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType _dNdx_p;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C;
   // typename BMatricesType::KelvinMatrixType C_solid;
   // typename BMatricesType::KelvinMatrixType C_ice;
    double integration_weight;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(double const t,
                                    SpatialPosition const& x_position,
                                    double const dt,
                                    DisplacementVectorType const& u
                                    /*double phi_i*/)
    {
        eps.noalias() = b_matrices * u;
        solid_material.computeConstitutiveRelation(
            t, x_position, dt, eps_prev, eps, sigma_eff_prev, sigma_eff, /*C_solid, C_ice,*/ C,/* phi_i,*/
            *material_state_variables);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

struct ThermoHydroMechanicsLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
    virtual std::vector<double> const& getIntPtSigmaXX(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtSigmaYY(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtSigmaZZ(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtSigmaXY(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtSigmaXZ(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtSigmaYZ(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtEpsilonXX(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtEpsilonYY(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtEpsilonZZ(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtEpsilonXY(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtEpsilonXZ(
            std::vector<double>& cache) const = 0;

        virtual std::vector<double> const& getIntPtEpsilonYZ(
            std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
class ThermoHydroMechanicsLocalAssembler
    : public ThermoHydroMechanicsLocalAssemblerInterface
{
public:
/*    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for displacement.
    // (Higher order elements = ShapeFunctionDisplacement).
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename BMatricesType::NodalForceVectorType;

    // Types for pressure.
    // (Lower order elements = Order(ShapeFunctionDisplacement) - 1).
    using ShapeFunctionPressure =
        typename LowerDim<ShapeFunctionDisplacement>::type;
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>; */

    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;
    using ShapeMatrices = typename ShapeMatricesTypeDisplacement::ShapeMatrices;


    ThermoHydroMechanicsLocalAssembler(ThermoHydroMechanicsLocalAssembler const&) = delete;
    ThermoHydroMechanicsLocalAssembler(ThermoHydroMechanicsLocalAssembler&&) = delete;

    ThermoHydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices_u =
            initShapeMatrices<ShapeFunctionDisplacement, ShapeMatricesTypeDisplacement,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        auto const shape_matrices_p =
            initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices_u[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() * sm.integralMeasure *
                shape_matrices_u[ip].detJ;
            ip_data.b_matrices.resize(
                kelvin_vector_size,
                ShapeFunctionDisplacement::NPOINTS * DisplacementDim);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesTypeDisplacement>(
                    e, shape_matrices_u[ip].N);
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS>(
                shape_matrices_u[ip].dNdx, ip_data.b_matrices,
                is_axially_symmetric, shape_matrices_u[ip].N, x_coord);

            ip_data.sigma_eff.resize(kelvin_vector_size);
            ip_data.sigma_eff_prev.resize(kelvin_vector_size);
            ip_data.eps.resize(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.C.resize(kelvin_vector_size, kelvin_vector_size);
            //ip_data.C_ice.resize(kelvin_vector_size, kelvin_vector_size);
            //ip_data.C_solid.resize(kelvin_vector_size, kelvin_vector_size);

            //ip_data._N_u = shape_matrices_u[ip].N;
            ip_data.N_u = ShapeMatricesTypeDisplacement::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);
            for (int i = 0; i < DisplacementDim; ++i)
                ip_data.N_u
                    .template block<1, displacement_size / DisplacementDim>(
                        i, i * displacement_size / DisplacementDim)
                    .noalias() = shape_matrices_u[ip].N;

            ip_data._N_p = shape_matrices_p[ip].N;
            ip_data._dNdx_p = shape_matrices_p[ip].dNdx;
            _secondary_data.N[ip] = shape_matrices_u[ip].N;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "ThermoHydroMechanicsLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_matrix_size = local_x.size();
        assert(local_x.size() == temperature_size + pressure_size + displacement_size);

        auto T =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_x.data() + temperature_index,
                                      temperature_size);

        auto p =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                pressure_size> const>(local_x.data() + pressure_index,
                                      pressure_size);

        auto u = Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

        auto T_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                temperature_size> const>(local_xdot.data() + temperature_index,
                                      temperature_size);

        auto p_dot =
            Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                pressure_size> const>(local_xdot.data() + pressure_index,
                                      pressure_size);

        auto u_dot = Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
                                      displacement_size);

     /*   auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_rhs =
            MathLib::createZeroedVector<NodalDisplacementVectorType>(
                local_rhs_data, local_matrix_size); */

        // Jacobian Matrix
        auto local_Jac = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                displacement_size + pressure_size + temperature_size,
                displacement_size + pressure_size + temperature_size>>(
            local_Jac_data, displacement_size + pressure_size + temperature_size,
            displacement_size + pressure_size + temperature_size);

        // Jacobian Matrix numerical
        Eigen::MatrixXd local_Jac_numerical = Eigen::MatrixXd::Zero(
            local_matrix_size, local_matrix_size);
      /*  auto local_Jac_numerical = MathLib::createZeroedMatrix<
            typename ShapeMatricesTypeDisplacement::template MatrixType<
                displacement_size + pressure_size + temperature_size,
                displacement_size + pressure_size + temperature_size>>(
            local_Jac_data, displacement_size + pressure_size + temperature_size,
            displacement_size + pressure_size + temperature_size);*/

        // Residual vector for numerical jacobian
        Eigen::MatrixXd local_b_p = Eigen::VectorXd::Zero(local_matrix_size);

        Eigen::MatrixXd local_b_m = Eigen::VectorXd::Zero(local_matrix_size);

        // Residual
        auto local_rhs = MathLib::createZeroedVector<
            typename ShapeMatricesTypeDisplacement::template VectorType<
                displacement_size + pressure_size + temperature_size>>(
            local_rhs_data, displacement_size + pressure_size + temperature_size);
      /*  auto local_rhs_num = MathLib::createZeroedVector<
            typename ShapeMatricesTypeDisplacement::template VectorType<
                displacement_size + pressure_size + temperature_size>>(
            local_rhs_data, displacement_size + pressure_size + temperature_size); */

        // coeff matrix is used for the residual
        typename ShapeMatricesTypePressure::NodalMatrixType MTT;
        MTT.setZero(temperature_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType MTT_coeff;
        MTT_coeff.setZero(temperature_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType KTT;
        KTT.setZero(temperature_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType KTT_coeff;
        KTT_coeff.setZero(temperature_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType KTp;
        KTp.setZero(temperature_size, pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType Kpp;
        Kpp.setZero(pressure_size, pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType Mpp;
        Mpp.setZero(pressure_size, pressure_size);

        typename ShapeMatricesTypePressure::NodalMatrixType MpT;
        MpT.setZero(pressure_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType MpT_coeff;
        MpT_coeff.setZero(pressure_size, temperature_size);

        typename ShapeMatricesTypePressure::NodalMatrixType Mpu;
        Mpu.setZero(pressure_size, displacement_size);

        typename ShapeMatricesTypeDisplacement::template MatrixType<displacement_size,
                                                        displacement_size>
            Kuu;
        Kuu.setZero(displacement_size,  displacement_size);

        typename ShapeMatricesTypeDisplacement::template MatrixType<displacement_size,
                                                        displacement_size>
            Kuu_coeff;
        Kuu_coeff.setZero(displacement_size,  displacement_size);

   //     typename ShapeMatricesTypeDisplacement::template MatrixType<displacement_size,
   //                                                     pressure_size>
   //         Kup;
   //     Kup.setZero(displacement_size, pressure_size);

        typename ShapeMatricesTypeDisplacement::template MatrixType<displacement_size,
                                                        pressure_size>
            Kup_coeff;
        Kup_coeff.setZero(displacement_size, pressure_size);

        typename ShapeMatricesTypeDisplacement::template MatrixType<displacement_size,
                                                        temperature_size>
            KuT;
        KuT.setZero(displacement_size, temperature_size);

        typename ShapeMatricesTypeDisplacement::template MatrixType<displacement_size,
                                                        temperature_size>
            KuT_coeff;
        KuT_coeff.setZero(displacement_size, temperature_size);

        double const& dt = _process_data.dt;

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        auto const newton_factor = _process_data.newton_factor(t, x_position)[0];

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N_p = _ip_data[ip]._N_p;
            auto const& N_u = _ip_data[ip].N_u;
            auto const& dNdx_p = _ip_data[ip]._dNdx_p;

            auto const& N_T = N_p ;
            auto const& dNdx_T = dNdx_p ;

            auto const& B = _ip_data[ip].b_matrices;
            auto const& sigma_eff = _ip_data[ip].sigma_eff;

            auto const& C = _ip_data[ip].C;
            //auto const& C_ice = _ip_data[ip].C_ice;
            //auto const& C_solid = _ip_data[ip].C_solid;

            double S =
                _process_data.storage_coefficient(t, x_position)[0];
            //double const K_over_mu =
              //  _process_data.intrinsic_permeability(t, x_position)[0] /
              //  _process_data.fluid_viscosity(t, x_position)[0];
            double const intrinsic_permeability = _process_data.intrinsic_permeability(t, x_position)[0] ;
            double const viscosity0 = _process_data.fluid_viscosity(t, x_position)[0];
            double const beta_s = _process_data.beta_solid(t, x_position)[0];
            double const beta_f = _process_data.beta_fluid(t, x_position)[0];
            double const lambda_f = _process_data.lambda_f(t, x_position)[0];
            double const lambda_s = _process_data.lambda_s(t, x_position)[0];
            double const C_f = _process_data.fluid_heat_capacity(t, x_position)[0];
            double const C_s = _process_data.solid_heat_capacity(t, x_position)[0];
            double const T0 = _process_data.reference_temperature(t, x_position)[0];
            auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
            auto const rho_sr = _process_data.solid_density(t, x_position)[0];
            auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
            auto const porosity = _process_data.porosity(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;


            auto const rho_ir = 920 ;
            auto rho_i = 0.0;
            auto rho_f = 0.0;
            auto rho_s = 0.0;
            auto const C_i = 2090 ;
            auto const lambda_i = 2.2;
            double phi_i = 0.0;
            auto const sigmoid_coeff = 2;
            auto const sigmoid_coeff1 = 2 ;
            auto const latent_heat = 334000;
            //auto const beta_i = 2.1e-5 ;
            auto beta_i = 0 ;
            double const beta_IF = 0.09 ; // Freezing strain

            auto sigmoid_derive = 0.0 ;
            //auto sigmoid_derive_second = 0.0 ;
            //double lambda = 0.0 ;
            double heat_capacity = 0.0 ;
         //   assert(body_force.size() == DisplacementDim);
         //   auto const b =
         //       Eigen::Map<typename ShapeMatricesType::template VectorType<
         //           DisplacementDim> const>(body_force.data(), DisplacementDim);

            // get the element values for T and p
         //   auto T_int_pt = 0.0 ;
        //    double p_int_pt = 0.0 ;
        //    auto const num_nodes = 4;  // only for quad, need to change for Shapefunction::NPOINTS
            // order matters First T then P
        //    double const * const ptr = local_x[0] ;
        //    Eigen::Map<Eigen::VectorXd> local_Tp(ptr,num_nodes*2) ;
        //    std::vector<double>::const_iterator first = local_x.begin() ;
        //    std::vector<double>::const_iterator last = local_x.begin() + num_nodes*2;
        //    std::vector<double> local_Tp(first, last);
        //    NumLib::shapeFunctionInterpolate(local_Tp, N_p, T_int_pt, p_int_pt);
            auto T_int_pt = N_T * T ;
            double delta_T(T_int_pt - T0);

            rho_f = rho_fr*(1 - beta_f * delta_T);
            rho_s = rho_sr*(1 - beta_s * delta_T);
            rho_i = rho_ir*(1 - beta_i * delta_T);

           // Freezing process
            phi_i = CalcIceVolFrac(T_int_pt, sigmoid_coeff, porosity);
            //std::cout << "phi_i: " << phi_i <<std::endl;
            double multipl1 = beta_IF/3*phi_i ;
            //std::cout << "1: " << multipl1 <<std::endl;

            sigmoid_derive = Calcsigmoidderive(sigmoid_coeff1, porosity, T_int_pt); // dphi_i/dT
            double Kr = Relative_permeability(phi_i, porosity);
            double Real_hydraulic_conductivity = Hydraulic_conductivity(intrinsic_permeability, Kr, viscosity0);
            double K_over_mu = Real_hydraulic_conductivity ;
            //std::cout << K_over_mu << std::endl;

            //sigmoid_derive_second = Calcsigmoidsecondderive(sigmoid_coeff, porosity, T_int_pt);
            // lambda to be constant for now
            //lambda = (porosity - phi_i) * lambda_f + (1 - porosity) * lambda_s + phi_i * lambda_i;
            heat_capacity =
                    EquaHeatCapacity(phi_i, rho_fr, rho_sr, rho_ir,
            C_s, C_i, C_f, porosity, sigmoid_derive, latent_heat);
            //std::cout << "C_ice: "<< C_ice << std::endl;
            //std::cout << "C_solid: "<< C_solid << std::endl;
           // std::cout << "C: "<< C << std::endl;
            // dC(T)/dT
    //        auto const C_derive =  sigmoid_derive * C_ice;

      /*      auto p_nodal_values = Eigen::Map<const Eigen::VectorXd>(
            &local_x[num_nodes], num_nodes);   */

      /*      auto T_nodal_values = Eigen::Map<const Eigen::VectorXd>(
            &local_x[0], num_nodes); */

            auto velocity =
                 (-K_over_mu * dNdx_p * p).eval();
            velocity += K_over_mu * rho_fr * b;  // do not consider the density change here

            // mapping matrix (m)
            auto const& identity2 = MaterialLib::SolidModels::Invariants<
                kelvin_vector_size>::identity2;

            //
            // displacement equation, displacement part (K_uu) coeff matrix and Jacobian matrix are same here
            //
    /*        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u, phi_i);    */

            Kuu.noalias() += B.transpose() * C * B * w;
            Kuu_coeff.noalias() += B.transpose() * C * B * w;

           /* typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);
            for (int i = 0; i < DisplacementDim; ++i)
                N_u.template block<1, displacement_size / DisplacementDim>(
                       i, i * displacement_size / DisplacementDim)
                    .noalias() = _ip_data[ip]._N_u; */

            // todo rho should also be a function of phi_i
            double rho= (porosity - phi_i) * rho_f + (1 - porosity) * rho_s + phi_i * rho_i;
            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() -=
                (B.transpose() * sigma_eff - N_u.transpose() * rho * b) * w;
            local_rhs  // cancelling out the reference temperature and add the freezing expansion term (check whether N_T is necessary))
                .template segment<displacement_size>(displacement_index)
                .noalias() -=
            //    B.transpose()* (identity2 * beta_s/3 * T0) ;
                 B.transpose() * (C * identity2 * (beta_s/3) * T0 - multipl1* C * identity2) * w ;

            //
            // displacement equation, temperature part (K_uT) Jacobian Matrix TODO beta is considered to be constant here, in fact is function of T
            //
            double beta = beta_s * (1 - porosity) + porosity * beta_f;
          /*  KuT.noalias() += B.transpose() * (C_derive * B * u - C * identity2 * beta / 3 * N_T -
                                              C_derive*identity2*N_T*beta/3*(N_T*T) - C_derive*identity2*N_T*beta_IF)*w; */

            // displacement equation, temperature part (K_uT) Coeff Matrix
            KuT_coeff.noalias() += B.transpose() * C * identity2 * (beta_s/3) * N_T * w ;

            //
            // displacement equation, pressure part (K_up) Coeff Matrix and Jacobian Matrix are the same
            //
            // used for numerical Jacobian
        //    Kup.noalias() = B.transpose() * alpha * identity2 * N_p * w ;
            // add cryo suction effect
            Kup_coeff.noalias() += B.transpose() * alpha * identity2 * N_p * w ;


            //
            // pressure equation, pressure part (K_pp and M_pp). coeff matrix and Jacobian matrix are the same (todo permeability is a function of T)
            //
            Kpp.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w ;

            //S =  0.2/2e9;
            Mpp.noalias() += N_p.transpose() * S * N_p * w;
            //
            // fp
            //
            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() += dNdx_p.transpose() * rho_fr * K_over_mu * b * w ;
            //
            // pressure equation, temperature part (M_pT)
            //
         /*   MpT.noalias() += N_p.transpose() * (rho_ir*sigmoid_derive*(1/rho_fr - 1/rho_ir) - beta) * N_p * w / dt + N_p.transpose()*(rho_ir*sigmoid_derive*(1/rho_fr - 1/rho_ir))* N_p * w *T_dot ; */
            MpT_coeff.noalias() += N_p.transpose() * (rho_ir*sigmoid_derive*(1/rho_fr - 1/rho_ir) - beta) * N_p  * w ;


            //
            // pressure equation, displacement part.  (M_pu)
            //
            // Reusing Kup.transpose().
            //Mpu.noalias() += Kup_coeff.transpose() ;

            //
            // temperature equation, temperature part.
            //
            double lambda = (porosity - phi_i) * lambda_f + (1 - porosity) * lambda_s + phi_i * lambda_i;
            //KTT.noalias() += (dNdx_T.transpose() * lambda * dNdx_T + dNdx_T.transpose() * velocity * N_p * rho_fr * C_f * 0) * w ;
            // TODO lambda is constant here, in fact it should be function of T
            KTT_coeff.noalias() += (dNdx_T.transpose() * lambda * dNdx_T + N_T.transpose() * velocity.transpose() * dNdx_T * rho_fr * C_f *0)* w ;
            //double heat_capacity = porosity * C_f * rho_fr + (1 - porosity) * C_s * rho_sr;
            MTT_coeff.noalias() += N_T.transpose() * heat_capacity * N_T * w;

            // Jacobian matrix
       /*     KTT.noalias() += MTT_coeff/dt + KTT_coeff + sigmoid_derive*(rho_ir*C_i +-rho_fr*C_f - C_i*latent_heat)*T_dot; */

            //
            // temperature equation, pressure part !!!!!!.positive or negative no coeff matrix
            //
           // KTp.noalias() +=  K_over_mu * rho_fr * C_f * dNdx_p.transpose() * (dNdx_T * T) * N_T * w ;
      /*      KTp.noalias() +=  K_over_mu * rho_fr * C_f * N_T.transpose() * (dNdx_T * T).transpose() * dNdx_T * w ;  */


          //  KTp.noalias() +=  N_T.transpose() * heat_capacity * N_T * w;

            // construction of Residual without Perturbation
   /*         double const T_int_pt = N_T*T;
            // update the constitutive relation for displacement becasue of the increment
            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u, phi_i);
            // thermal expansion
            double delta_T(T_int_pt - T0);
            rho_fr = rho_fr*(1 - beta_f * delta_T);
            rho_sr = rho_sr*(1 - beta_s * delta_T);
            rho_ir = rho_ir*(1 - beta_i * delta_T);
            // Freezing process
            phi_i = CalcIceVolFrac(T_int_pt, sigmoid_coeff, porosity);
            // Permeability change due to freezing is not considered currenty in order to compare with analytical solution
            sigmoid_derive = Calcsigmoidderive(sigmoid_coeff, porosity, T_int_pt); // dphi_i/dT
            //sigmoid_derive_second = Calcsigmoidsecondderive(sigmoid_coeff, porosity, T_int_pt_m);
            // lambda to be constant for now
            lambda = (porosity - phi_i) * lambda_f + (1 - porosity) * lambda_s + phi_i * lambda_i;
            heat_capacity = EquaHeatCapacity(phi_i, rho_fr, rho_sr, rho_ir,
            C_s, C_i, C_f, porosity, sigmoid_derive, latent_heat);
            local_rhs
                .template segment<temperature_size>(temperature_index)
                .noalias() =
                -N_T.transpose() * heat_capacity * N_T * w * T_dot -
                    (dNdx_T.transpose() * lambda * dNdx_T - N_T.transpose() * velocity.transpose() * dNdx_T * rho_fr * C_f * 0)* T * w ;
            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() = -(N_p.transpose() * S * N_p * w * p_dot + dNdx_p.transpose() * K_over_mu * dNdx_p * w * p +
                    N_p.transpose()*(rho_ir*sigmoid_derive*(1/rho_fr-1/rho_ir)-beta)*N_p* w * T_dot + Kup_coeff.transpose()* u_dot)
                    + dNdx_p.transpose() * rho_fr * K_over_mu * b * w;
            local_rhs
                .template segment<displacement_size>(displacement_index)
                .noalias() =
                -(B.transpose() * sigma_eff - N_u.transpose() * rho * b) * w -
                    B.transpose() * (C * identity2 * (beta_s/3) * T0 + phi_i * C_ice * identity2 * beta_IF) * w +
                 B.transpose() * C * identity2 * (beta_s/3) * N_T * T* w + B.transpose() * alpha * identity2 * N_p *p * w;   */


            // calculate numerical Jac
                        double num_p = 1e-7;
                        //double _Theta = 0.5;
                        Eigen::VectorXd num_vec = Eigen::VectorXd::Zero(local_matrix_size);
                        std::vector<double> local_perturbed = local_x;
                        for (Eigen::MatrixXd::Index i = 0; i < local_matrix_size; i++)
                        {
                            num_vec[i] = num_p*(1+std::abs(local_x[i]));
                            local_perturbed[i] += num_vec[i];
                            auto T_p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                                temperature_size> const>(local_perturbed.data() + temperature_index,
                                                        temperature_size);
                            auto p_p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                                pressure_size> const>(local_perturbed.data() + pressure_index,
                                                        pressure_size);
                            auto u_p = Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
                                displacement_size> const>(local_perturbed.data() + displacement_index,
                                                          displacement_size);
                            typename ShapeMatricesTypeDisplacement::template MatrixType<displacement_size,
                                                                            pressure_size>
                                Kup;
                            Kup.setZero(displacement_size, pressure_size);
                            Kup = B.transpose() * alpha * identity2 * N_p * w ;
                            double T_int_pt_p = N_T*T_p;
                            phi_i = CalcIceVolFrac(T_int_pt_p, sigmoid_coeff, porosity);
                            Kr = Relative_permeability(phi_i, porosity);
                            Real_hydraulic_conductivity = Hydraulic_conductivity(intrinsic_permeability, Kr, viscosity0);
                            K_over_mu = Real_hydraulic_conductivity ;
                            //S = 0.2/2e9 ;
                            double multipl2 = phi_i*beta_IF/3;
                            //std::cout << "2: " << multipl2 <<std::endl;
                            // update the constitutive relation for displacement becasue of the increment
                            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u_p/*, phi_i*/);
                            // thermal expansion
                            double delta_T_p(T_int_pt_p - T0);
                            rho_f = rho_fr*(1 - beta_f * delta_T_p);
                            rho_s = rho_sr*(1 - beta_s * delta_T_p);
                            rho_i = rho_ir*(1 - beta_i * delta_T_p);
                            rho = rho_s * (1 - porosity) + (porosity - phi_i) * rho_f + phi_i * rho_i;
                            // Freezing process
                            // Permeability change due to freezing is not considered currenty in order to compare with analytical solution
                            sigmoid_derive = Calcsigmoidderive(sigmoid_coeff1, porosity, T_int_pt_p); // dphi_i/dT
                            //sigmoid_derive_second = Calcsigmoidsecondderive(sigmoid_coeff, porosity, T_int_pt);
                            // lambda to be constant for now
                            lambda = (porosity - phi_i) * lambda_f + (1 - porosity) * lambda_s + phi_i * lambda_i;
                            heat_capacity =
                                    EquaHeatCapacity(phi_i, rho_fr, rho_sr, rho_ir,
                            C_s, C_i, C_f, porosity, sigmoid_derive, latent_heat);
                            // dC(T)/dT
                            // auto const C_derive =  sigmoid_derive * C_ice;
                            local_b_p
                                .template block<temperature_size,1>(temperature_index,0)
                                .noalias() =
                                N_T.transpose() * heat_capacity * N_T * w * T_p/dt +
                                   (dNdx_T.transpose() * lambda * dNdx_T + N_T.transpose() * velocity.transpose() * dNdx_T * rho_f * C_f *0)* T_p * w ;
                            local_b_p.template block<pressure_size,1>(pressure_index,0)
                                .noalias() = N_p.transpose() * S * N_p * w * p_p/dt + dNdx_p.transpose() * K_over_mu * dNdx_p * w * p_p
                                 + N_p.transpose()*(rho_ir*sigmoid_derive*(1/rho_fr-1/rho_ir)-beta)*N_p* w * T_p/dt   + Kup.transpose() * u_p/dt
                                    - dNdx_p.transpose() * rho_f * K_over_mu * b * w;
                            local_b_p
                                .template block<displacement_size,1>(displacement_index,0)
                                .noalias() =
                                (B.transpose() * sigma_eff - N_u.transpose() * rho * b) * w +
                                    B.transpose() * (C * identity2 * (beta_s/3) * T0 + multipl2*C * identity2 )* w -
                                 B.transpose() * C * identity2 * (beta_s/3) * N_T * T_p * w - B.transpose() * alpha * identity2 * N_p * p_p * w;
                            //std::cout << "C: " << C << std::endl ;
                            local_perturbed[i] = local_x[i] - num_vec[i];
                            auto T_m = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                                temperature_size> const>(local_perturbed.data() + temperature_index,
                                                        temperature_size);
                            auto p_m = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                                pressure_size> const>(local_perturbed.data() + pressure_index,
                                                        pressure_size);
                            auto u_m = Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
                                displacement_size> const>(local_perturbed.data() + displacement_index,
                                                          displacement_size);
                            double T_int_pt_m = N_T*T_m;
                            phi_i = CalcIceVolFrac(T_int_pt_m, sigmoid_coeff, porosity);
                            Kr = Relative_permeability(phi_i, porosity);
                            Real_hydraulic_conductivity = Hydraulic_conductivity(intrinsic_permeability, Kr, viscosity0);
                            K_over_mu = Real_hydraulic_conductivity ;
                            //S = 0.2/2e9 ;
                            double multipl3 = phi_i*beta_IF/3 ;
                            //std::cout << "3: " << multipl3 <<std::endl;
                            // update the constitutive relation for displacement becasue of the increment
                            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u_m/*, phi_i*/);
                            // thermal expansion
                            double delta_T_m(T_int_pt_m - T0);
                            rho_f = rho_fr*(1 - beta_f * delta_T_m);
                            rho_s = rho_sr*(1 - beta_s * delta_T_m);
                            rho_i = rho_ir*(1 - beta_i * delta_T_m);
                            rho = rho_s * (1 - porosity) + (porosity - phi_i) * rho_f + phi_i * rho_i;
                            // Freezing process                            
                            // Permeability change due to freezing is not considered currenty in order to compare with analytical solution
                            sigmoid_derive = Calcsigmoidderive(sigmoid_coeff1, porosity, T_int_pt_m); // dphi_i/dT
                            //sigmoid_derive_second = Calcsigmoidsecondderive(sigmoid_coeff, porosity, T_int_pt_m);
                            // lambda to be constant for now
                            lambda = (porosity - phi_i) * lambda_f + (1 - porosity) * lambda_s + phi_i * lambda_i;
                            heat_capacity =
                                    EquaHeatCapacity(phi_i, rho_fr, rho_sr, rho_ir,
                            C_s, C_i, C_f, porosity, sigmoid_derive, latent_heat);
                            local_b_m
                                .template block<temperature_size,1>(temperature_index,0)
                                .noalias() =
                                N_T.transpose() * heat_capacity * N_T * w * T_m/dt +
                                    (dNdx_T.transpose() * lambda * dNdx_T + N_T.transpose() * velocity.transpose() * dNdx_T * rho_fr * C_f*0 )* T_m * w ;
                            local_b_m.template block<pressure_size,1>(pressure_index,0)
                                .noalias() = N_p.transpose() * S * N_p * w * p_m/dt + dNdx_p.transpose() * K_over_mu * dNdx_p * w * p_m  +
                                    N_p.transpose()*(rho_ir*sigmoid_derive*(1/rho_fr-1/rho_ir)-beta)*N_p* w * T_m/dt  + Kup.transpose()* u_m/dt
                                    - dNdx_p.transpose() * rho_f * K_over_mu * b * w;
                            local_b_m
                                .template block<displacement_size,1>(displacement_index,0)
                                .noalias() =
                                (B.transpose() * sigma_eff - N_u.transpose() * rho * b) * w +
                                    B.transpose() * (C * identity2 * (beta_s/3) * T0 +  C * identity2 *multipl3) * w -
                                 B.transpose() * C * identity2 * (beta_s/3) * N_T * T_m* w - B.transpose() * alpha * identity2 * N_p *(p_m) * w;
                            local_perturbed[i] = local_x[i];
                            local_Jac_numerical.col(i).noalias() += (local_b_p - local_b_m) / (2.0 * num_vec[i]);



                        }
        }
 /*       // temperature equation, temperature part
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += KTT + MTT / dt ;

        // temperature equation, pressure part
        local_Jac
            .template block<temperature_size, pressure_size>(
                temperature_index, pressure_index)
            .noalias() -= KTp ;

        // displacement equation, displacment part
        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += Kuu ;

        // displacement equation, temperature part
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= KuT ;

        // displacement equation, pressure part
        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() -= Kup;

        // pressure equation, temperature part.
        local_Jac
            .template block<pressure_size, temperature_size>(pressure_index,
                                                          temperature_index)
            .noalias() -= MpT ;

        // pressure equation, pressure part.
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += Kpp + Mpp / dt;

        // pressure equation, displacement part.
        local_Jac
            .template block<pressure_size, displacement_size>(
                pressure_index, displacement_index)
            .noalias() += Mpu / dt; */

        // pressure equation (f_p)
        local_rhs.template segment<pressure_size>(pressure_index)
          .noalias() -=
            Kpp * p + Mpp * p_dot + MpT_coeff * T_dot + Kup_coeff.transpose() * u_dot;


        // displacement equation (f_u) ref = 0 case  (freezing cyro suction pressure)
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() += Kup_coeff* p + KuT_coeff * T  ;

        // temperature equation (f_T)
        local_rhs.template segment<temperature_size>(temperature_index)
            .noalias() -= KTT_coeff * T + MTT_coeff * T_dot;


        local_Jac = local_Jac_numerical*newton_factor  ;
       // local_rhs = local_rhs ;


    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }


    void postTimestepConcrete(std::vector<double> const& local_x) override
         {
             double const& t = _process_data.t;

             auto p =
                 Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                     pressure_size> const>(local_x.data() + pressure_index,
                                           pressure_size);
             auto T =
                 Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
                     temperature_size> const>(local_x.data() + temperature_index,
                                           temperature_size);

             using GlobalDimVectorType =
                 typename ShapeMatricesTypePressure::GlobalDimVectorType;

             unsigned const n_integration_points =
                 _integration_method.getNumberOfPoints();

             SpatialPosition x_position;
             x_position.setElementID(_element.getID());
             for (unsigned ip = 0; ip < n_integration_points; ip++)
             {
                 x_position.setIntegrationPoint(ip);
               //  double const K_over_mu =
               //      _process_data.intrinsic_permeability(t, x_position)[0] /
               //      _process_data.fluid_viscosity(t, x_position)[0];
                 double const intrinsic_permeability = _process_data.intrinsic_permeability(t, x_position)[0] ;
                 double const viscosity0 = _process_data.fluid_viscosity(t, x_position)[0];

                 auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
                 auto const& b = _process_data.specific_body_force;
                 auto const porosity = _process_data.porosity(t, x_position)[0];

                 // Compute the velocity
                 auto const& N_T = _ip_data[ip]._N_p;
                 auto const sigmoid_coeff = 2;
                 auto T_int_pt = N_T * T ;
                 double phi_i = CalcIceVolFrac(T_int_pt, sigmoid_coeff, porosity);
                 double Kr = Relative_permeability(phi_i, porosity);
                 double Real_hydraulic_conductivity = Hydraulic_conductivity(intrinsic_permeability, Kr, viscosity0);
                 auto K_over_mu = Real_hydraulic_conductivity ;

                 auto const& dNdx_p = _ip_data[ip]._dNdx_p;
                 GlobalDimVectorType const darcy_velocity = -K_over_mu * dNdx_p * p
                 -K_over_mu* rho_fr* b;
                 for (unsigned d = 0; d < DisplacementDim; ++d)
                 {
                     _darcy_velocities[d][ip] = darcy_velocity[d];
                 }
             }
         }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSigmaXX(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 0);
    }

    std::vector<double> const& getIntPtSigmaYY(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 1);
    }

    std::vector<double> const& getIntPtSigmaZZ(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 2);
    }

    std::vector<double> const& getIntPtSigmaXY(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 3);
    }

    std::vector<double> const& getIntPtSigmaXZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 4);
    }

    std::vector<double> const& getIntPtSigmaYZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 5);
    }

    std::vector<double> const& getIntPtEpsilonXX(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 0);
    }

    std::vector<double> const& getIntPtEpsilonYY(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 1);
    }

    std::vector<double> const& getIntPtEpsilonZZ(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 2);
    }

    std::vector<double> const& getIntPtEpsilonXY(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 3);
    }

    std::vector<double> const& getIntPtEpsilonXZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 4);
    }

    std::vector<double> const& getIntPtEpsilonYZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 5);
    }

    std::vector<double> const& getIntPtDarcyVelocityX(
             std::vector<double>& /*cache*/) const override
    {
         assert(_darcy_velocities.size() > 0);
         return _darcy_velocities[0];
    }

    std::vector<double> const& getIntPtDarcyVelocityY(
            std::vector<double>& /*cache*/) const override
    {
         assert(_darcy_velocities.size() > 1);
         return _darcy_velocities[1];
    }

    std::vector<double> const& getIntPtDarcyVelocityZ(
    std::vector<double>& /*cache*/) const override
    {
         assert(_darcy_velocities.size() > 2);
         return _darcy_velocities[2];
    }


private:

    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                                 std::size_t const component) const
        {
            cache.clear();
            cache.reserve(_ip_data.size());

            for (auto const& ip_data : _ip_data) {
                if (component < 3)  // xx, yy, zz components
                    cache.push_back(ip_data.sigma_eff[component]);
                else    // mixed xy, yz, xz components
                    cache.push_back(ip_data.sigma_eff[component] / std::sqrt(2));
            }

            return cache;
        }
        std::vector<double> const& getIntPtEpsilon(
            std::vector<double>& cache, std::size_t const component) const
        {
            cache.clear();
            cache.reserve(_ip_data.size());

            for (auto const& ip_data : _ip_data) {
                cache.push_back(ip_data.eps[component]);
            }

            return cache;
        }

    ThermoHydroMechanicsProcessData<DisplacementDim>& _process_data;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    std::vector<std::vector<double>> _darcy_velocities =
    std::vector<std::vector<double>>(
              DisplacementDim,
    std::vector<double>(_integration_method.getNumberOfPoints()));

    static const int temperature_index = 0;
    static const int temperature_size = ShapeFunctionPressure::NPOINTS;
    static const int pressure_index = ShapeFunctionPressure::NPOINTS;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS*2;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerData final
    : public ThermoHydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                                ShapeFunctionPressure ,
                                                IntegrationMethod, DisplacementDim>
{
public:
    LocalAssemblerData(LocalAssemblerData const&) = delete;
    LocalAssemblerData(LocalAssemblerData&&) = delete;

    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       ThermoHydroMechanicsProcessData<DisplacementDim>& process_data)
        : ThermoHydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                             ShapeFunctionPressure, IntegrationMethod,
                                             DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
