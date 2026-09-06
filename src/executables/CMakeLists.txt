include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/physics/database )
include_directories( ${BELFEM_SOURCE_DIR}/circuit )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/ode )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/maxwell )
include_directories( ${BELFEM_SOURCE_DIR}/fem/thermal )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/homology )

set ( LIBLIST
        sparse
        bezier
        mesh
        sources
        integration
        ode
        database
        materials
        interpolation
        integration
        circuit
        iwg
        homology
        kernel
        postproc
        maxwell
        thermal
         )


# the old executables are obsolete and have been retired. belfem replaces both:
# it selects magnetic-only or coupled h-phi/T from the deck. The sources are kept
# because the documentation still cites them as reference drivers.
#set( EXECNAME hphirun )
#set( MAIN hphirun.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

#set( EXECNAME hphiTrun )
#set( MAIN hphiTrun.cpp)
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

set( EXECNAME belfem )
set( MAIN belfem.cpp)
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
