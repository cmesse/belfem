set( LIBNAME homology )

set( SOURCES
        cl_Chain.cpp
        cl_Cochain.cpp
        cl_SimplicialComplex.cpp
        cl_Homology.cpp
        cl_Cohomology.cpp
        fn_Smith.cpp
)

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )