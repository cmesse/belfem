if( USE_NLOPT )
    belfem_find_package(
            NLOPT
            REQUIRED
            HEADERS nlopt.h
            LIBRARIES nlopt )

    list( FIND BELFEM_DEFS "BELFEM_NLOPT" ITEM_INDEX )
    if( ITEM_INDEX EQUAL -1 )
        list( APPEND BELFEM_DEFS "BELFEM_NLOPT" )
    endif()

    list( APPEND BELFEM_INCLUDES "${NLOPT_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${NLOPT_LIB_DIR}" )
    belfem_link_libraries( 2 ${NLOPT_LIBRARIES} )
endif ()
