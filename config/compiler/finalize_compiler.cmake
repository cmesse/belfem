# -------------------------------------------------------------------------
# CXX FLAGS
# -------------------------------------------------------------------------
include( ${BELFEM_CONFIG_DIR}/globals.cmake )

# add compiler flags
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BELFEM_CXXFLAGS}" )

# Third-party include directories go in as system headers, so a warning raised inside a
# vendor header is suppressed and cannot trip -Werror. This is preventative: Ubuntu's GCC
# driver injects -Wformat-security, which under -Wall -Werror broke the 0.9.0 build in two
# BELFEM headers (fixed in 0.9.1); the same mechanism would turn the next vendor-header
# warning into a build failure. Only directories in this list are covered -- mpi.h, for one,
# comes through the compiler wrapper. BELFEM's own directories never pass through this list:
# they are added with include_directories() in config/scripts/Add_*.cmake and stay ordinary
# -I, so their diagnostics are unaffected.
# A directory the compiler already searches by default keeps -I (which it then ignores):
# -isystem on such a directory would move it ahead of the standard library headers.
# normalize first (several TPL configs append the same SCLS dir, some with a trailing
# slash), so the deduplicated list is what both the C++ and the Fortran loop see
set( _BELFEM_INCS )
foreach( ITEM ${BELFEM_INCLUDES} )
    string( REGEX REPLACE "/+$" "" ITEM "${ITEM}" )
    list( APPEND _BELFEM_INCS "${ITEM}" )
endforeach()
if( _BELFEM_INCS )
    list( REMOVE_DUPLICATES _BELFEM_INCS )
endif()
set( BELFEM_INCLUDES ${_BELFEM_INCS} )
unset( _BELFEM_INCS )
foreach( ITEM ${BELFEM_INCLUDES} )
    if( ITEM IN_LIST CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES )
        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${ITEM}" )
    else()
        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${ITEM}" )
    endif()
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
    # whoever appended a system libdir, it does not go on the rpath
    belfem_prune_system_libdirs( BELFEM_RPATH )
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
