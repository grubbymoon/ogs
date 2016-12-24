/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_SOLIDMODELS_FREEZINGLINEARELASTICISOTROPIC_H_
#define MATERIALLIB_SOLIDMODELS_FREEZINGLINEARELASTICISOTROPIC_H_

#include "MechanicsFreezingBase.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
class FreezingLinearElasticIsotropic final : public MechanicsFreezingBase<DisplacementDim>
{
public:
    /// Variables specific to the material model
    class MaterialProperties
    {
        using P = ProcessLib::Parameter<double>;
        using X = ProcessLib::SpatialPosition;

    public:
        MaterialProperties(P const& youngs_modulus_solid, P const& poissons_ratio_solid, P const& youngs_modulus_ice, P const& poissons_ratio_ice)
            : _youngs_modulus_solid(youngs_modulus_solid), _poissons_ratio_solid(poissons_ratio_solid),
              _youngs_modulus_ice(youngs_modulus_ice), _poissons_ratio_ice(poissons_ratio_ice)
        {
        }

        /// Lamé's first parameter for solid.
        double lambda_solid(double const t, X const& x) const
        {
            return _youngs_modulus_solid(t, x)[0] * _poissons_ratio_solid(t, x)[0] /
                   (1 + _poissons_ratio_solid(t, x)[0]) /
                   (1 - 2 * _poissons_ratio_solid(t, x)[0]);
        }

        /// Lamé's second parameter, the shear modulus for solid.
        double mu_solid(double const t, X const& x) const
        {
            return _youngs_modulus_solid(t, x)[0] /
                   (2 * (1 + _poissons_ratio_solid(t, x)[0]));
        }

        /// Lamé's first parameter for ice.
        double lambda_ice(double const t, X const& x) const
        {
            return _youngs_modulus_ice(t, x)[0] * _poissons_ratio_ice(t, x)[0] /
                   (1 + _poissons_ratio_ice(t, x)[0]) /
                   (1 - 2 * _poissons_ratio_ice(t, x)[0]);
        }

        /// Lamé's second parameter, the shear modulus for ice.
        double mu_ice(double const t, X const& x) const
        {
            return _youngs_modulus_ice(t, x)[0] /
                   (2 * (1 + _poissons_ratio_ice(t, x)[0]));
        }


    private:
        P const& _youngs_modulus_solid;
        P const& _poissons_ratio_solid;
        P const& _youngs_modulus_ice;
        P const& _poissons_ratio_ice;
    };

    struct MaterialStateVariables
        : public MechanicsFreezingBase<DisplacementDim>::MaterialStateVariables
    {
        void pushBackState() {}
    };

    std::unique_ptr<
        typename MechanicsFreezingBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::unique_ptr<
            typename MechanicsFreezingBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables};
    }

public:
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    explicit FreezingLinearElasticIsotropic(
        MaterialProperties const& material_properties)
        : _mp(material_properties)
    {
    }

    bool computeFreezingConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const /*dt*/,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        KelvinVector& sigma,
        KelvinMatrix& C_solid,
        KelvinMatrix& C_ice,
        KelvinMatrix& C,
        double phi_i,
        typename MechanicsFreezingBase<DisplacementDim>::MaterialStateVariables&
        /*material_state_variables*/) override
    {
        //phi_i = 0.5;
        C_solid.setZero();

        C_solid.template topLeftCorner<3, 3>().setConstant(_mp.lambda_solid(t, x));
        C_solid.noalias() += 2 * _mp.mu_solid(t, x) * KelvinMatrix::Identity();

        C_ice.setZero();

        C_ice.template topLeftCorner<3, 3>().setConstant(_mp.lambda_ice(t, x));
        C_ice.noalias() += 2 * _mp.mu_ice(t, x) * KelvinMatrix::Identity();

        C = phi_i*C_ice + (1-phi_i)*C_solid;
        sigma.noalias() = sigma_prev + C * (eps - eps_prev);
        return true;
    }

private:
    MaterialProperties _mp;
};

extern template class FreezingLinearElasticIsotropic<2>;
extern template class FreezingLinearElasticIsotropic<3>;

}  // namespace Solids
}  // namespace MaterialLib

#endif  // MATERIALLIB_SOLIDMODELS_LINEARELASTICISOTROPIC_H_
