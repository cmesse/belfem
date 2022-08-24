set( LIBNAME postproc )

set( SOURCES
        fn_FEM_mises_planestress.cpp
        fn_Mesh_compute_surface_normals.cpp
        fn_Mesh_compute_edge_lengths.cpp
        fn_Mesh_compute_volume.cpp
        fn_Mesh_compute_surface.cpp
        fn_Mesh_integrate_scalar_over_sidesets.cpp
        cl_OneDMapper.cpp
        cl_Pipette.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )

set ( LIBLIST
        sparse
        mesh
        integration
        interpolation
        materials
        iwg
        postproc
        kernel
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set( EXECNAME pipette )
set( MAIN pipette.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )


if( USE_EXAMPLES )
set( EXECNAME normaltest )
set( MAIN normaltest.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

set( EXECNAME pentatest )
set( MAIN pentatest.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

set( EXECNAME mappertest )
set( MAIN mappertest.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

endif()