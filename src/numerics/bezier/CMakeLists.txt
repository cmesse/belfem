set( LIBNAME bezier )

set( SOURCES
        cl_Bezier.cpp )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/bezier )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )