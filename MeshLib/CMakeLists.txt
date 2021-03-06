# Source files
GET_SOURCE_FILES(SOURCES_MESHLIB)

# It should be removed too for other MPI based DDC approach in future.
if(NOT OGS_USE_PETSC)
	list(REMOVE_ITEM SOURCES_MESHLIB NodePartitionedMesh.h)
endif()

GET_SOURCE_FILES(SOURCES_ELEMENTS Elements)
GET_SOURCE_FILES(SOURCES_EDITING MeshEditing)
GET_SOURCE_FILES(SOURCES_GENERATORS MeshGenerators)
GET_SOURCE_FILES(SOURCES_SEARCH MeshSearch)

GET_SOURCE_FILES(SOURCES_QUALITY MeshQuality)
set(SOURCES ${SOURCES_MESHLIB} ${SOURCES_ELEMENTS} ${SOURCES_EDITING} ${SOURCES_GENERATORS} ${SOURCES_QUALITY} ${SOURCES_SEARCH})

# Create the library
add_library(MeshLib STATIC ${SOURCES})

target_link_libraries(MeshLib INTERFACE
	BaseLib
	FileIO
	GeoLib
	MathLib
	${VTK_LIBRARIES}
	logog
)

ADD_VTK_DEPENDENCY(MeshLib)

if(TARGET Eigen)
    add_dependencies(MeshLib Eigen)
endif()

if(TARGET Boost)
    add_dependencies(MeshLib Boost)
endif()
