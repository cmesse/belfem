set( LIBNAME iwg )

set( SOURCES
        cl_IwgFactory.cpp
        cl_IWG.cpp
        cl_IWG_Poisson.cpp
        cl_IWG_StaticHeatConduction.cpp
        cl_IWG_TransientHeatConduction.cpp
        cl_IWG_Timestep.cpp
        cl_TimestepMatrices.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )

set ( LIBLIST
        sparse
        mesh
        integration
        interpolation
        materials
        iwg
        kernel
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )