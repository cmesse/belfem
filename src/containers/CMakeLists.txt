set( LIBNAME containers )

set( SOURCES
        cl_StringList.cpp
        cl_DynamicBitset.cpp
        cl_Cell.cpp
)

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )