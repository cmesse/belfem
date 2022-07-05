set( LIBNAME integration )

set( SOURCES
        intpoints.f90
        en_IntegrationScheme.cpp
        fn_intpoints.cpp
        fn_intpoints_auto_integration_order.cpp )

include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )

set ( LIBLIST
        graph
        mesh )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )