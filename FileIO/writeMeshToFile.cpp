/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "writeMeshToFile.h"

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "MeshLib/Mesh.h"

#include "Legacy/MeshIO.h"
#include "VtkIO/VtuInterface.h"

namespace FileIO
{
void writeMeshToFile(const MeshLib::Mesh &mesh, const std::string &file_name)
{
	if (BaseLib::hasFileExtension("msh", file_name))
	{
		Legacy::MeshIO meshIO;
		meshIO.setMesh(&mesh);
		meshIO.writeToFile(file_name);
	} else if (BaseLib::hasFileExtension("vtu", file_name)) {
		FileIO::VtuInterface writer(&mesh);
		writer.writeToFile(file_name);
	} else {
		ERR("writeMeshToFile(): Unknown mesh file format in file %s.", file_name.c_str());
	}
}

} // end namespace FileIO
