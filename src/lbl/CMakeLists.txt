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
            interpolation
            edgeinterpolation
            gastables
            gasmodels
            iwg
            maxwell
            materials
            kernel
            postproc
            integration )
else()
    set ( LIBLIST
            sparse
            bezier
            mesh
            integration
            ode
            interpolation
            edgeinterpolation
            iwg
            maxwell
            materials
            kernel
            postproc
            integration )
endif()


#set( EXECNAME test_dofmgr )
#set( MAIN test_dofmgr.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_n2ha )
#set( MAIN test_n2ha.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_e2ha )
#set( MAIN test_e2ha.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

set( EXECNAME petsctst )
set( MAIN petsctst.cpp)
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_cut )
#set( MAIN test_cut.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

set( EXECNAME test_input )
set( MAIN test_input.cpp)
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_curved )
#set( MAIN test_curved.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_tape )
#set( MAIN test_tape.cpp )
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_l2 )
#set( MAIN test_l2.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_phiback )
#set( MAIN test_phiback.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME test_tet10vis )
#set( MAIN test_tet10vis.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME ndlc_l2test )
#set( MAIN ndlc_l2test.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME ndlc_htest )
#set( MAIN ndlc_htest.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )



#set( EXECNAME surfacetest )
#set( MAIN surfacetest.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )