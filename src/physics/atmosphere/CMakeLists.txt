set( LIBNAME atmosphere )

set( SOURCES
        cl_AtmosphereModel_ISA1976.cpp
        cl_Atmosphere.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        sparse
        spline
        gastables
        gasmodels
        atmosphere )

if( USE_EXAMPLES )
    set( EXECNAME atmos )
    set( MAIN     atmos.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()
