set( LIBNAME sparse )

set( SOURCES
        cl_SpMatrix.cpp
        cl_SolverWrapper.cpp
        cl_SolverUMFPACK.cpp
        cl_SolverMUMPS.cpp
        cl_SolverPARDISO.cpp
        cl_SolverPETSC.cpp
        cl_Solver.cpp
        cl_SolverPetscDistributor.cpp
        cl_SolverSTRUMPACK.cpp
        cl_SolverStrumpackDistributor.cpp
        en_SolverEnums.cpp
        fn_rcond.cpp
        petsctools.cpp )

if ( USE_MUMPS )
    list( APPEND SOURCES mumpstools.f90)
endif()

if ( USE_PARDISO )
    list( APPEND SOURCES pardisotools.f90)
endif()



include_directories(  ${BELFEM_SOURCE_DIR}/math/graph )
include_directories(  ${BELFEM_SOURCE_DIR}/io )
include_directories(  ${BELFEM_SOURCE_DIR}/sparse )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )
if( USE_EXAMPLES )
    set( EXECNAME solvertest)
    set( MAIN     solvertest.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()
