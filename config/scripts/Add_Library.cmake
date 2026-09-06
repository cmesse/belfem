include( ${BELFEM_CONFIG_DIR}/globals.cmake )

# add default includes
include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/containers )
include_directories( ${BELFEM_SOURCE_DIR}/linalg )
include_directories( ${BELFEM_SOURCE_DIR}/linalg/lapack )
if ( USE_MATRIX_ARMADILLO )
    include_directories( ${BELFEM_SOURCE_DIR}/linalg/armadillo )
elseif( USE_MATRIX_BLAZE )
    include_directories( ${BELFEM_SOURCE_DIR}/linalg/blaze )
endif()
include_directories( ${BELFEM_SOURCE_DIR}/linalg/operators )
include_directories( ${BELFEM_SOURCE_DIR}/io )

include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/sparse )

# cl_Material.hpp exposes Bezier in its interface, so every target that sees a
# material needs this on the include path
include_directories( ${BELFEM_SOURCE_DIR}/numerics/bezier )

# generated belfem_version.hpp (see src/core/CMakeLists.txt)
include_directories( ${CMAKE_BINARY_DIR}/generated )

# each module is an OBJECT library: it compiles on its own ( `make belfem_<name>`
# is the per-module compile gate ) but is never linked on its own. The objects
# of every module are collected into the single project library `belfem`, see
# Add_BelfemLibrary.cmake. One library instead of one archive per module,
# because the modules depend on each other circularly ( core <-> comm,
# core <-> io, kernel <-> iwg, ... ) and a cycle between shared libraries does
# not link on Darwin.
set( BELFEM_MODULE_TARGET ${LIBPREFIX}_${LIBNAME} )

add_library( ${BELFEM_MODULE_TARGET} OBJECT ${SOURCES} )

set_property( GLOBAL APPEND PROPERTY BELFEM_MODULE_TARGETS ${BELFEM_MODULE_TARGET} )

# make sure the version header is generated before this library compiles. the
# belfem_version target is defined in core (the first subdirectory); the guard
# skips core itself, which wires up the dependency explicitly after this point.
if( TARGET belfem_version )
    add_dependencies( ${BELFEM_MODULE_TARGET} belfem_version )
endif()
