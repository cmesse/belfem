set( LIBNAME comm )

include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/containers )
include_directories( ${BELFEM_SOURCE_DIR}/linalg )

set( SOURCES
        cl_Communicator.cpp
        commtools.cpp)

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

# where `make install` puts share/, so an installed binary finds its data
# without $BELFEM_DATA. Only this module sees it; the value is a prefix,
# never a path into the source tree ( see Communicator::set_globals() )
target_compile_definitions( ${BELFEM_MODULE_TARGET} PRIVATE
        BELFEM_INSTALL_DATADIR="${CMAKE_INSTALL_FULL_DATADIR}/${LIBPREFIX}" )
