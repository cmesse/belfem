set( LIBNAME iwg )

set( SOURCES
        cl_IWG.cpp
        cl_IWG_Timestep.cpp
        cl_IWG_StaticHeatConduction.cpp
        cl_IwgFactory.cpp
        )

if ( USE_GASMODELS )
    set(SOURCES
            ${SOURCES}
            cl_IWG_TATCAD.cpp)
endif()

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )

if ( USE_GASMODELS )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )

    set ( LIBLIST
            sparse
            mesh
            integration
            interpolation
            materials
            gasmodels
            gastables
            iwg
            kernel
            )
else()
    set ( LIBLIST
        sparse
        mesh
        integration
        interpolation
        materials
        iwg
        kernel
        )
endif()

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )