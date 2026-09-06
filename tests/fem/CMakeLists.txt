set( TESTNAME fem )
set( TESTLABELS fast )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )

set( SOURCES
        test_FacetIntegrationPoints.cpp
        test_LagrangeInterpolation.cpp
        test_EdgeFunctions.cpp
        test_Integration.cpp
        test_IntegrationExactness.cpp
        test_AndersonMixing.cpp
        test_BdfTimestepMethod.cpp
        test_DofSeeding.cpp
        test_GhostSwitch.cpp
        )

set( LIBLIST
        mesh
        integration
        interpolation
        kernel
        iwg
        postproc
        )

# the interface-orientation battery exercises the Maxwell IWGs; the maxwell
# library only exists when the module is enabled
if( USE_MAXWELL )
    include_directories( ${BELFEM_SOURCE_DIR}/fem/maxwell )
    include_directories( ${BELFEM_SOURCE_DIR}/fem/maxwell/matrices )
    list( APPEND SOURCES test_InterfaceOrientation.cpp )
    list( APPEND LIBLIST maxwell )
endif()

list( APPEND LIBLIST materials )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
