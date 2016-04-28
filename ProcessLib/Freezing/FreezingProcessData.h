/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_FREEZING_FREEZINGPROCESSDATA_H
#define PROCESSLIB_FREEZING_FREEZINGPROCESSDATA_H

namespace MeshLib
{
    class Element;
}


namespace ProcessLib
{

template <typename ReturnType, typename... Args>
struct Parameter;

namespace Freezing
{

// copy that and change names
struct FreezingProcessData
{
    FreezingProcessData(
            ProcessLib::Parameter<double, MeshLib::Element const&> const&
            thermal_conductivity_
            )
        : thermal_conductivity(thermal_conductivity_)
    {}

    FreezingProcessData(FreezingProcessData&& other)
        : thermal_conductivity(other.thermal_conductivity)
    {}

    //! Copies are forbidden.
    FreezingProcessData(FreezingProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(FreezingProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(FreezingProcessData&&) = delete;

    Parameter<double, MeshLib::Element const&> const& thermal_conductivity;
};

} // namespace Freezing
} // namespace ProcessLib

#endif // PROCESSLIB_GROUNDWATERFLOW_GROUNDWATERFLOWPROCESSDATA_H
