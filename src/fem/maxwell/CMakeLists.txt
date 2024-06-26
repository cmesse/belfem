set( LIBNAME maxwell )

set( SOURCES
        cl_BhExtrapolationOde.cpp
        cl_Maxwell_FieldList.cpp
        cl_IWG_Maxwell.cpp
        #cl_IWG_Maxwell_2D.cpp
        cl_IWG_Maxwell_HPhi_Tri3.cpp
        cl_IWG_Maxwell_HPhi_Tri6.cpp
        cl_IWG_Maxwell_L2_Current.cpp
        cl_IWG_Maxwell_L2_Magfield.cpp
        cl_IWG_Maxwell_Thermal2D.cpp
        cl_MaxwellMaterial.cpp
        cl_MaxwellJob.cpp
        cl_FEM_MaxwellBoundaryConditionCurrent.cpp
        cl_FEM_MaxwellBoundaryConditionMagfield.cpp
        cl_MaxwellMeshSynch.cpp
        cl_Maxwell_ThermalMeshExtractor.cpp
        cl_MaxwellFactory.cpp
        fn_FEM_compute_brandt.cpp
        fn_FEM_compute_normb.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/edgeinterpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/maxwell )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/ode )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )

set ( LIBLIST
        sparse
        ode
        mesh
        integration
        interpolation
        edgeinterpolation
        materials
        iwg
        maxwell
        postproc
        kernel
        postproc
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

if( USE_EXAMPLES )
    set( EXECNAME test_geometrytools )
    set( MAIN test_geometrytools.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()