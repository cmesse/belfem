set( LIBNAME edgeinterpolation )

set( SOURCES
        cl_EF_TRI3.cpp
        cl_EF_TRI6.cpp
        cl_EF_TET4.cpp
        cl_EF_TET10.cpp
        cl_EdgeFunctionFactory.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
set ( LIBLIST
        mesh
        edgeinterpolation
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )