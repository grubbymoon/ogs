/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_SOLIDMODELS_CREATEFREEZINGLINEARELASTICISOTROPIC_H_
#define MATERIALLIB_SOLIDMODELS_CREATEFREEZINGLINEARELASTICISOTROPIC_H_

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter
#include "FreezingLinearElasticIsotropic.h"

namespace MaterialLib
{
namespace Solids
{

template <int DisplacementDim>
std::unique_ptr<FreezingLinearElasticIsotropic<DisplacementDim>>
createFreezingLinearElasticIsotropic(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__SMALL_DEFORMATION__constitutive_relation__type}
    config.checkConfigParameter("type", "FreezingLinearElasticIsotropic");
    DBUG("Create FreezingLinearElasticIsotropic material");

    // Youngs modulus for solid
    auto& youngs_modulus_solid = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{process__SMALL_DEFORMATION__constitutive_relation__LinearElasticIsotropic__youngs_modulus}
        config, "youngs_modulus_solid", parameters, 1);

    DBUG("Use '%s' as youngs_modulus_solid parameter.", youngs_modulus_solid.name.c_str());

    // Poissons ratio for solid
    auto& poissons_ratio_solid = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{process__SMALL_DEFORMATION__constitutive_relation__LinearElasticIsotropic__poissons_ratio}
        config, "poissons_ratio_solid", parameters, 1);

    DBUG("Use '%s' as poissons_ratio_ice parameter.", poissons_ratio_solid.name.c_str());

    // Youngs modulus for ice
    auto& youngs_modulus_ice = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{process__SMALL_DEFORMATION__constitutive_relation__LinearElasticIsotropic__youngs_modulus}
        config, "youngs_modulus_ice", parameters, 1);

    DBUG("Use '%s' as youngs_modulus_ice parameter.", youngs_modulus_ice.name.c_str());

    // Poissons ratio for ice
    auto& poissons_ratio_ice = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{process__SMALL_DEFORMATION__constitutive_relation__LinearElasticIsotropic__poissons_ratio}
        config, "poissons_ratio_ice", parameters, 1);

    DBUG("Use '%s' as poissons_ratio_ice parameter.", poissons_ratio_ice.name.c_str());

    typename FreezingLinearElasticIsotropic<DisplacementDim>::MaterialProperties mp{
        youngs_modulus_solid, poissons_ratio_solid, youngs_modulus_ice, poissons_ratio_ice};

    return std::unique_ptr<FreezingLinearElasticIsotropic<DisplacementDim>>{
        new FreezingLinearElasticIsotropic<DisplacementDim>{mp}};
}

}  // namespace Solids
}  // namespace MaterialLib

#endif  // MATERIALLIB_SOLIDMODELS_CREATELINEARELASTICISOTROPIC_H_
