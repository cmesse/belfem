# -------------------------------------------------------------------------
# CXX FLAGS
# -------------------------------------------------------------------------
include( ${BELFEM_CONFIG_DIR}/globals.cmake )

# add compiler flags
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BELFEM_CXXFLAGS} " )

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
#  LD FLAGS
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# ASSEMBLE LIBLIST
# -------------------------------------------------------------------------

# add project path
link_directories(${CMAKE_BINARY_DIR}/${LIBDIR})

# add other paths
foreach( ITEM ${BELFEM_RPATH} )
    link_directories(${ITEM})
endforeach()