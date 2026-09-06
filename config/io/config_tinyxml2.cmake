if( USE_TINYXML2 )
    belfem_find_package(
            TINYXML2
            REQUIRED
            HEADERS tinyxml2.h
            LIBRARIES tinyxml2 )

    list( APPEND BELFEM_DEFS "BELFEM_XML" )
    list( APPEND BELFEM_INCLUDES "${TINYXML2_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${TINYXML2_LIB_DIR}" )
    belfem_link_libraries( 2 ${TINYXML2_LIBRARIES} )
endif()
