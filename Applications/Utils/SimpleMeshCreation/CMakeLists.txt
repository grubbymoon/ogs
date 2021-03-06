
include_directories(
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/Legacy
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/Gui/VtkVis
)

if(OGS_BUILD_GUI)
	add_executable(createMeshElemPropertiesFromASCRaster createMeshElemPropertiesFromASCRaster.cpp)
	set_target_properties(createMeshElemPropertiesFromASCRaster PROPERTIES FOLDER Utilities)
	target_link_libraries(createMeshElemPropertiesFromASCRaster FileIO)
	ADD_VTK_DEPENDENCY(createMeshElemPropertiesFromASCRaster)
endif() # OGS_BUILD_GUI (VtkVis-target is existing)


add_executable(generateStructuredMesh generateStructuredMesh.cpp)
target_link_libraries(generateStructuredMesh FileIO)
set_target_properties(generateStructuredMesh PROPERTIES FOLDER Utilities)

####################
### Installation ###
####################
install(TARGETS
	generateStructuredMesh
	RUNTIME DESTINATION bin
	COMPONENT Utilities
)
if(TARGET VtkVis)
	install(TARGETS
		createMeshElemPropertiesFromASCRaster
		RUNTIME DESTINATION bin
		COMPONENT Utilities
	)
endif()
