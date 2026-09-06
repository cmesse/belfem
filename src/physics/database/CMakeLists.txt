set( LIBNAME database )

set( SOURCES
        cl_Database.cpp
        cl_DatabaseProjector.cpp
)

include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        mesh
        integration
        interpolation
        kernel )

# db2exo reads HDF5 and writes exodus. Without either it would still build and
# link -- both backends compile to stubs -- and then abort at runtime, so the
# target is gated instead. A name in BELFEM_INSTALL_EXECUTABLES with no target
# is ignored ( config/globals.cmake ), so the install entry needs no guard
if( USE_HDF5 AND USE_EXODUS )
    set( EXECNAME db2exo )
    set( MAIN     db2exo.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()
