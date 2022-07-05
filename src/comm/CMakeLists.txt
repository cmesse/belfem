set( LIBNAME comm )

include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/containers )
include_directories( ${BELFEM_SOURCE_DIR}/linalg )

set( SOURCES
        cl_Communicator.cpp
        commtools.cpp)

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )