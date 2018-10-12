APPEND_SOURCE_FILES(SOURCES)

add_library(ThermoHydroMechanics ${SOURCES})
target_link_libraries(ThermoHydroMechanics PUBLIC ProcessLib)

include(Tests.cmake)
