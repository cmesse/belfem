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
include_directories( ${BELFEM_SOURCE_DIR}/fem/thermal )

set( SOURCES
        test_FacetIntegrationPoints.cpp
        test_LagrangeInterpolation.cpp
        test_EdgeFunctions.cpp
        test_Integration.cpp
        test_IntegrationExactness.cpp
        test_AndersonMixing.cpp
        test_BdfTimestepMethod.cpp
        test_DofSeeding.cpp
        test_DirichletPin.cpp
        test_Pipette.cpp
        test_IwgNodeCoords.cpp
        test_GhostSwitch.cpp
        test_BlockRebuild.cpp
        test_Tmatrix.cpp
        test_OwnershipContracts.cpp
        test_SolverParametersDeck.cpp
        test_RhoWindow.cpp
        test_IwgDefaults.cpp
        test_ElementMapper.cpp
        test_DofSelfContainment.cpp
        test_QoldTable.cpp
        test_IntegrationAuto.cpp
        test_DeckIntegers.cpp
        test_PenaltyUnits.cpp
        test_Hex8FlipEdges.cpp
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
    # cl_MaxwellFactory.hpp includes cl_Topology.hpp
    include_directories( ${BELFEM_SOURCE_DIR}/homology )
    list( APPEND SOURCES test_InterfaceOrientation.cpp )
    list( APPEND SOURCES test_SideSetRebuild.cpp )
    list( APPEND SOURCES test_ThinShellSelection.cpp )
    list( APPEND SOURCES test_PostprocessorRecovery.cpp )
    list( APPEND SOURCES test_AdmissionRefusals.cpp )
    list( APPEND SOURCES test_IwgPenaltyDefaults.cpp )
    list( APPEND SOURCES test_FreeCutCurrents.cpp )
    list( APPEND SOURCES test_Hex8TsVolume.cpp )
    list( APPEND SOURCES test_BackgroundField.cpp )
    list( APPEND SOURCES test_ThinShellTerminalCurves.cpp )
    list( APPEND LIBLIST maxwell )
endif()

list( APPEND LIBLIST materials )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
