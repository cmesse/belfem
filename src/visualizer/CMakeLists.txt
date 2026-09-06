set( LIBNAME visualize )

set( SOURCES
        cl_VTK_MeshView.cpp
        cl_VTK_Curve.cpp
)

# VTK 9.6 headers emit trailing ';' after macros, which trips the project's
# global -pedantic-errors. Silence -Wpedantic for this VTK-only module so the
# (otherwise harmless) VTK header noise does not fail the build.
add_compile_options( -Wno-pedantic )

include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/io )
include_directories( ${BELFEM_SOURCE_DIR}/containers )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        integration
        mesh
        visualize
)

set( EXECNAME visualize )
set( MAIN     visualize.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

