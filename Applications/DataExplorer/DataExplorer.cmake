# Source files
set(SOURCES
	mainwindow.cpp
	${CMAKE_SOURCE_DIR}/Applications/Utils/OGSFileConverter/OGSFileConverter.cpp
	${CMAKE_SOURCE_DIR}/Applications/Utils/OGSFileConverter/FileListDialog.cpp
)

# Moc Header files
set(MOC_HEADERS
	mainwindow.h
	${CMAKE_SOURCE_DIR}/Applications/Utils/OGSFileConverter/OGSFileConverter.h
	${CMAKE_SOURCE_DIR}/Applications/Utils/OGSFileConverter/FileListDialog.h
)

# UI files
set(UIS
	mainwindow.ui
	${CMAKE_SOURCE_DIR}/Applications/Utils/OGSFileConverter/OGSFileConverter.ui
	${CMAKE_SOURCE_DIR}/Applications/Utils/OGSFileConverter/FileList.ui
)


# Run Qts user interface compiler uic on .ui files
qt4_wrap_ui(UI_HEADERS ${UIS} )

qt4_add_resources(QTRESOURCES ./Img/icons.qrc )

# Run Qts meta object compiler moc on header files
qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS} )

# Include the headers which are generated by uic and moc
# and include additional header
set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)
include_directories(
	${SOURCE_DIR_REL}/BaseLib
	${SOURCE_DIR_REL}/MathLib
	${SOURCE_DIR_REL}/GeoLib
	${SOURCE_DIR_REL}/FileIO
	${SOURCE_DIR_REL}/MeshLib
	${SOURCE_DIR_REL}/MeshLibGEOTOOLS
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_BINARY_DIR}/Base
	${CMAKE_CURRENT_BINARY_DIR}/DataView
	${CMAKE_CURRENT_BINARY_DIR}/DataView/StratView
	${CMAKE_CURRENT_BINARY_DIR}/DataView/DiagramView
	${CMAKE_CURRENT_BINARY_DIR}/VtkVis
	${CMAKE_CURRENT_BINARY_DIR}/VtkAct
	${CMAKE_CURRENT_BINARY_DIR}/Applications/Utils/OGSFileConverter
	${CMAKE_CURRENT_SOURCE_DIR}/Base
	${CMAKE_CURRENT_SOURCE_DIR}/DataView
	${CMAKE_CURRENT_SOURCE_DIR}/DataView/StratView
	${CMAKE_CURRENT_SOURCE_DIR}/DataView/DiagramView
	${CMAKE_CURRENT_SOURCE_DIR}/VtkVis
	${CMAKE_CURRENT_SOURCE_DIR}/VtkAct
)

# Put moc files in a project folder
source_group("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")

# Application icon
set(APP_ICON ${CMAKE_SOURCE_DIR}/scripts/packaging/ogs-de-icon.icns)

# Create the executable
add_executable(DataExplorer MACOSX_BUNDLE
	main.cpp
	${SOURCES}
	${MOC_HEADERS}
	${MOC_SOURCES}
	${UIS}
	${QTRESOURCES}
	${APP_ICON}
	exe-icon.rc
)

target_link_libraries(DataExplorer
	${QT_LIBRARIES}
	ApplicationsLib
	BaseLib
	GeoLib
	FileIO
	InSituLib
	MeshLib
	#MSHGEOTOOLS
	QtBase
	QtDataView
	QtStratView
	VtkVis
	VtkAct
	${Boost_LIBRARIES}
	${CATALYST_LIBRARIES}
	zlib
	shp
)

if(VTK_NETCDF_FOUND)
	target_link_libraries(DataExplorer vtkNetCDF vtkNetCDF_cxx )
else()
	target_link_libraries(DataExplorer ${Shapelib_LIBRARIES} )
endif () # Shapelib_FOUND

if (GEOTIFF_FOUND)
	target_link_libraries(DataExplorer ${GEOTIFF_LIBRARIES} )
endif () # GEOTIFF_FOUND

add_dependencies (DataExplorer VtkVis)

if(MSVC)
	# Set linker flags
	set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:MSVCRT /IGNORE:4099")
	target_link_libraries(DataExplorer winmm)
endif()

### OpenSG support ###
if (VTKOSGCONVERTER_FOUND)
	USE_OPENSG(DataExplorer)
	include_directories(${VTKOSGCONVERTER_INCLUDE_DIRS})
	target_link_libraries(DataExplorer ${VTKOSGCONVERTER_LIBRARIES})
endif ()

if(VTKFBXCONVERTER_FOUND)
	target_link_libraries(DataExplorer ${VTKFBXCONVERTER_LIBRARIES})
endif()

include(AddCatalystDependency)
ADD_CATALYST_DEPENDENCY(DataExplorer)

set_property(TARGET DataExplorer PROPERTY FOLDER "DataExplorer")


####################
### Installation ###
####################
if(APPLE)
	include(packaging/PackagingMacros)
	ConfigureMacOSXBundle(DataExplorer ${APP_ICON})

	install(TARGETS DataExplorer DESTINATION .)
	set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION .)
	include(InstallRequiredSystemLibraries)
	include(DeployQt4)
	INSTALL_QT4_EXECUTABLE(DataExplorer.app "" "" "" "" "" ogs_gui)
else()
	install(TARGETS DataExplorer RUNTIME DESTINATION bin COMPONENT ogs_gui)
endif()

cpack_add_component(ogs_gui
	DISPLAY_NAME "OGS Data Explorer"
	DESCRIPTION "The graphical user interface for OpenGeoSys."
	GROUP Applications
)
set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "DataExplorer" "OGS Data Explorer" PARENT_SCOPE)
set(CPACK_NSIS_MENU_LINKS ${CPACK_NSIS_MENU_LINKS} "bin/DataExplorer.exe" "Data Explorer" PARENT_SCOPE)
if(APPLE)
	return()
endif()

if(MSVC)
	set(OGS_GUI_EXE ${EXECUTABLE_OUTPUT_PATH}/Release/DataExplorer.exe)
else()
	set(OGS_GUI_EXE ${EXECUTABLE_OUTPUT_PATH}/DataExplorer)
endif()

include(GetPrerequisites)
if(EXISTS ${OGS_GUI_EXE})
	if(MSVC)
		get_prerequisites(${OGS_GUI_EXE} OGS_GUI_DEPENDENCIES 1 1 "" "")
	else()
		get_prerequisites(${OGS_GUI_EXE} OGS_GUI_DEPENDENCIES 0 1 "/usr/local/lib;/;${VTK_DIR};/usr/lib64;" "")
	endif()
	message(STATUS "DataExplorer depends on:")
	foreach(DEPENDENCY ${OGS_GUI_DEPENDENCIES})
		if(NOT ${DEPENDENCY} STREQUAL "not") # Some bug on Linux?
			gp_resolve_item("/" "${DEPENDENCY}" ${OGS_GUI_EXE} "/usr/local/lib;/;${VTK_DIR};/usr/lib64;" DEPENDENCY_PATH)
			get_filename_component(RESOLVED_DEPENDENCY_PATH "${DEPENDENCY_PATH}" REALPATH)
			string(TOLOWER ${DEPENDENCY} DEPENDENCY_LOWER)
			if("${DEPENDENCY_LOWER}" MATCHES "tiff|blas|lapack|proj|jpeg|qt|gfortran|vtk|boost|png")
				set(DEPENDENCY_PATHS ${DEPENDENCY_PATHS} ${RESOLVED_DEPENDENCY_PATH} ${DEPENDENCY_PATH})
				message("    ${DEPENDENCY}")
			endif()
		endif()
	endforeach()
	install(FILES ${DEPENDENCY_PATHS} DESTINATION bin COMPONENT ogs_gui)
	if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
		install(PROGRAMS data-explorer.sh DESTINATION . COMPONENT ogs_gui)
	endif()
	add_custom_command(TARGET DataExplorer POST_BUILD COMMAND ;) # For caching: excetuting empty command
else()
	# Run CMake after DataExplorer was built to run GetPrerequisites on executable
	add_custom_command(TARGET DataExplorer POST_BUILD COMMAND ${CMAKE_COMMAND}
		ARGS ${CMAKE_SOURCE_DIR} WORKING_DIRECTORY ${CMAKE_BINARY_DIR} VERBATIM)
endif()