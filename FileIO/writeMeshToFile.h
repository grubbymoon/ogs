/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef WRITEMESHTOFILE_H_
#define WRITEMESHTOFILE_H_

#include <string>

namespace MeshLib
{
class Mesh;
}

namespace FileIO
{
void writeMeshToFile(const MeshLib::Mesh &mesh, const std::string &file_name);
}

#endif // WRITEMESHTOFILE_H_
