# -------------------------------------------------------------------------
# CXX FLAGS
# -------------------------------------------------------------------------
include( ${BELFEM_CONFIG_DIR}/globals.cmake )

# add compiler flags
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BELFEM_CXXFLAGS}" )

# add includes
foreach( ITEM ${BELFEM_INCLUDES} )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${ITEM}" )
endforeach()

# tidy up
string(STRIP "${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS)

# -------------------------------------------------------------------------
# Definitions
# -------------------------------------------------------------------------

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${BELFEM_FCFLAGS}" )

foreach( ITEM ${BELFEM_DEFS} )
    add_definitions( -D${ITEM} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D${ITEM}" )
endforeach()

foreach( ITEM ${BELFEM_INCLUDES} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${ITEM}" )
endforeach()


# tidy up
string(STRIP "${CMAKE_Fortran_FLAGS}" CMAKE_Fortran_FLAGS)

# -------------------------------------------------------------------------
# Assemble third-party link list by descending dependency rank
# -------------------------------------------------------------------------

# see belfem_link_libraries() in belfem_find_package.cmake: dependents must
# link before their dependencies, so higher ranks come first
set( BELFEM_TPL_LIBS )
foreach( _BELFEM_RANK 10 9 8 7 6 5 4 3 2 1 0 )
    if( DEFINED BELFEM_TPL_RANK_${_BELFEM_RANK} )
        list( APPEND BELFEM_TPL_LIBS ${BELFEM_TPL_RANK_${_BELFEM_RANK}} )
    endif()
endforeach()

# -------------------------------------------------------------------------
# Link List Hygiene
# -------------------------------------------------------------------------

foreach( _BELFEM_LINK_LIST
         BELFEM_FORTRANLIBS
         BELFEM_OPENMPLIBS )
    if( ${_BELFEM_LINK_LIST} )
        list( REMOVE_DUPLICATES ${_BELFEM_LINK_LIST} )
    endif()
endforeach()

if( BELFEM_RPATH )
    list( REMOVE_DUPLICATES BELFEM_RPATH )
endif()

# -------------------------------------------------------------------------
#  LD FLAGS
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# ASSEMBLE LIBLIST
# -------------------------------------------------------------------------

# add project and third-party library paths
get_property( _BELFEM_LINK_DIRECTORIES DIRECTORY PROPERTY LINK_DIRECTORIES )
list( APPEND _BELFEM_LINK_DIRECTORIES ${CMAKE_BINARY_DIR}/${LIBDIR} ${BELFEM_RPATH} )
list( REMOVE_DUPLICATES _BELFEM_LINK_DIRECTORIES )
set_property( DIRECTORY PROPERTY LINK_DIRECTORIES ${_BELFEM_LINK_DIRECTORIES} )
