set(SOURCES
	VtkPickCallback.cpp
	VtkCustomInteractorStyle.cpp
)

set(HEADERS
	VtkPickCallback.h
	VtkCustomInteractorStyle.h
)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/../../BaseLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../GeoLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../FileIO
	${CMAKE_CURRENT_SOURCE_DIR}/../../FileIO/GmshIO
	${CMAKE_CURRENT_SOURCE_DIR}/../../MeshLib
	${CMAKE_CURRENT_SOURCE_DIR}/../../MathLib
	${CMAKE_CURRENT_SOURCE_DIR}/../VtkVis
)

add_library( VtkAct STATIC
	${SOURCES}
	${HEADERS}
)

ADD_VTK_DEPENDENCY(VtkAct)

target_link_libraries( VtkAct PUBLIC Qt4::QtCore )
target_link_libraries( VtkAct INTERFACE ${VTK_LIBRARIES} )

set_property(TARGET VtkAct PROPERTY FOLDER "DataExplorer")
