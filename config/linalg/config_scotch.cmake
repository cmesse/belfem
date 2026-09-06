if( USE_SCOTCH )
    list(FIND BELFEM_DEFS "BELFEM_SCOTCH" ITEM_INDEX )
    if(ITEM_INDEX EQUAL -1)
        list(APPEND BELFEM_DEFS "BELFEM_SCOTCH")
    endif()
    list(FIND BELFEM_DEFS "BELFEM_PTSCOTCH" ITEM_INDEX )
    if(ITEM_INDEX EQUAL -1)
        list(APPEND BELFEM_DEFS "BELFEM_PTSCOTCH")
    endif()

    belfem_find_package(
            SCOTCH
            REQUIRED
            HEADERS scotch.h ptscotch.h
            LIBRARIES ptscotch ptscotcherr scotch scotcherr
            INCLUDE_SUFFIXES . scotch )

    list( APPEND BELFEM_INCLUDES "${SCOTCH_INCLUDE_DIR}" )
    list( APPEND BELFEM_RPATH "${SCOTCH_LIB_DIR}" )

    belfem_link_libraries( 5 ${SCOTCH_LIBRARIES} )
endif()
