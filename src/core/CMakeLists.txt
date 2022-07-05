set( LIBNAME core )

set( SOURCES
        stringtools.cpp
        banner.cpp
        cl_Logger.cpp
        cl_Profiler.cpp
        cl_Arguments.cpp
        cl_Progressbar.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/containers )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )