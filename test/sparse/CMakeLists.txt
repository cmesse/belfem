# List source files
set( TESTNAME sparse )

set( SOURCES
        cl_SpMatrix_CSR.cpp
        cl_SpMatrix_CSC.cpp
        )

if ( USE_PETSC )
    list( APPEND SOURCES cl_SolverPETSC.cpp )
endif()

set ( LIBLIST
        sparse )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )