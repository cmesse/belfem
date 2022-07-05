set( LIBNAME maxwell )

set( SOURCES
        cl_Maxwell_FieldList.cpp
        cl_IWG_Maxwell.cpp
        cl_IWG_Maxwell_HA_Tri3.cpp
        cl_IWG_Maxwell_HA_Tri6.cpp
        cl_IWG_Maxwell_HPhi_Tri3.cpp
        cl_IWG_Maxwell_HPhi_Tri6.cpp
        cl_IWG_Maxwell_HPhi_Tet4.cpp
        cl_IWG_Maxwell_HPhi_Tet10.cpp
        cl_IWG_Maxwell_L2_Current.cpp
        cl_IWG_Maxwell_L2_Magfield.cpp
        cl_MaxwellMaterial.cpp
        cl_MaxwellJob.cpp
        cl_FEM_MaxwellBoundaryConditionCurrent.cpp
        cl_FEM_MaxwellBoundaryConditionMagfield.cpp
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
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )

set ( LIBLIST
        sparse
        mesh
        integration
        interpolation
        materials
        iwg
        maxwell
        kernel
        postproc
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )