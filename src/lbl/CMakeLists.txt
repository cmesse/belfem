include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
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

if ( USE_GASMODELS )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )

    set ( LIBLIST
            sparse
            bezier
            mesh
            integration
            ode
            materials
            interpolation
            integration
            edgeinterpolation
            gastables
            gasmodels
            iwg
            kernel
            postproc
            maxwell
             )
else()
    set ( LIBLIST
            sparse
            bezier
            mesh
            integration
            ode
            materials
            interpolation
            integration
            edgeinterpolation
            iwg
            kernel
            postproc
            maxwell )
endif()


#set( EXECNAME test_dofmgr )
#set( MAIN test_dofmgr.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )




#set( EXECNAME petsctst )
#set( MAIN petsctst.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_cut )
#set( MAIN test_cut.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

set( EXECNAME hphirun )
set( MAIN hphirun.cpp)
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )