/**
 * @brief Extracts the surface from the given mesh.
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "tclap/CmdLine.h"

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/StringTools.h"
#include "BaseLib/FileTools.h"

#include "FileIO/readMeshFromFile.h"
#include "FileIO/writeMeshToFile.h"

#include "MathLib/Vector3.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSurfaceExtraction.h"

int main (int argc, char* argv[])
{
	ApplicationsLib::LogogSetup logog_setup;

	TCLAP::CmdLine cmd("Tool extracts the surface of the given mesh.", ' ',
	                   "0.1");
	TCLAP::ValueArg<std::string> mesh_in(
	    "i", "mesh-input-file",
	    "the name of the file containing the input mesh", true, "",
	    "file name of input mesh");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out(
	    "o", "mesh-output-file",
	    "the name of the file the surface mesh should be written to", false, "",
	    "file name of output mesh");
	cmd.add(mesh_out);
	TCLAP::ValueArg<double> x("x", "x-component", "x component of the normal",
	                          false, 0, "floating point value");
	cmd.add(x);
	TCLAP::ValueArg<double> y("y", "y-component", "y component of the normal",
	                          false, 0, "floating point value");
	cmd.add(y);
	TCLAP::ValueArg<double> z("z", "z-component", "z component of the normal",
	                          false, -1.0, "floating point value");
	cmd.add(z);
	TCLAP::ValueArg<double> angle_arg(
	    "a", "angle", "angle between given normal and element normal", false,
	    90, "floating point value");
	cmd.add(angle_arg);

	cmd.parse(argc, argv);

	std::unique_ptr<MeshLib::Mesh const> mesh(
	    FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %u nodes, %u elements.", mesh->getNNodes(), mesh->getNElements());

	// extract surface
	MathLib::Vector3 const dir(x.getValue(), y.getValue(), z.getValue());
	double const angle(angle_arg.getValue());
	std::unique_ptr<MeshLib::Mesh> surface_mesh(
	    MeshLib::MeshSurfaceExtraction::getMeshSurface(
	        *mesh, dir, angle, "OriginalSubsurfaceNodeIDs"));

	std::string out_fname(mesh_out.getValue());
	if (out_fname.empty())
		out_fname = BaseLib::dropFileExtension(mesh_in.getValue()) + "_sfc.vtu";
	FileIO::writeMeshToFile(*surface_mesh, out_fname);

	return EXIT_SUCCESS;
}
