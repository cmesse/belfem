set( LIBNAME sparse )

set( SOURCES
        cl_SpMatrix.cpp
        cl_SolverWrapper.cpp
        cl_SolverUMFPACK.cpp
        cl_SolverSUPERLU.cpp
        cl_SolverMUMPS.cpp
        cl_SolverPARDISO.cpp
        cl_SolverPETSC.cpp
        cl_Solver.cpp
        cl_SolverSTRUMPACK.cpp
        cl_SolverDistMatrix.cpp
        cl_SolverParameters.cpp
        en_SolverEnums.cpp
        fn_rcond.cpp
        fn_compute_permutation.cpp
        fn_create_graph_from_matrix.cpp
        strumpacktools.cpp
        petsctools.cpp
        splinalg.f90 )

if ( USE_MUMPS )
    list( APPEND SOURCES mumpstools.f90)
endif()

if ( USE_PARDISO )
    list( APPEND SOURCES pardisotools.f90)
endif()

if ( USE_ARPACK )
    # the error decoders in the cpp are only meaningful together with the
    # fortran driver, so the pair is gated as a unit
    list( APPEND SOURCES arpacktools.cpp arpacktools.f90 )
endif()

if ( USE_PARPACK )
    list( APPEND SOURCES parpacktools.f90 )
endif()

include_directories(  ${BELFEM_SOURCE_DIR}/math/graph )
include_directories(  ${BELFEM_SOURCE_DIR}/io )
include_directories(  ${BELFEM_SOURCE_DIR}/sparse )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )
