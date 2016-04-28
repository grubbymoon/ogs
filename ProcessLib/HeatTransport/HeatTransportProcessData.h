/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_HEATTRANSPORT_HEATTRANSPORTPROCESSDATA_H
#define PROCESSLIB_HEATTRANSPORT_HEATTRANSPORTPROCESSDATA_H

namespace MeshLib
{
    class Element;
}


namespace ProcessLib
{

template <typename ReturnType, typename... Args>
struct Parameter;

namespace HeatTransport
{

// copy that and change names
struct HeatTransportProcessData
{
    HeatTransportProcessData(
            ProcessLib::Parameter<double, MeshLib::Element const&> const&
            thermal_conductivity_
            )
        : thermal_conductivity(thermal_conductivity_)
    {}

    HeatTransportProcessData(HeatTransportProcessData&& other)
        : thermal_conductivity(other.thermal_conductivity)
    {}

    //! Copies are forbidden.
    HeatTransportProcessData(HeatTransportProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HeatTransportProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HeatTransportProcessData&&) = delete;

    Parameter<double, MeshLib::Element const&> const& thermal_conductivity;
};

} // namespace HeatTransport
} // namespace ProcessLib

#endif // PROCESSLIB_GROUNDWATERFLOW_GROUNDWATERFLOWPROCESSDATA_H
