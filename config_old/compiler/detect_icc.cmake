#-------------------------------------------------------------------------------
# check the other executables
#-------------------------------------------------------------------------------
# get the version of the c++ compiler
execute_process( COMMAND bash "-c" "$CXX -dumpversion" OUTPUT_VARIABLE BELFEM_CXX_VERSION )

# tidy up string
string(STRIP "${BELFEM_CXX_VERSION}" BELFEM_CXX_VERSION)

# test if there is a c compiler
if( DEFINED ENV{CC} )

    # get version of c-compiler
    execute_process( COMMAND bash "-c" "$CC -dumpversion" OUTPUT_VARIABLE BELFEM_CC_VERSION )

    # tidy up string
    string(STRIP "${BELFEM_CC_VERSION}" BELFEM_CC_VERSION)

    # test if versions match
    if( NOT BELFEM_CXX_VERSION VERSION_EQUAL BELFEM_CC_VERSION )
        message( FATAL_ERROR "The version numbers of icc and icpc do not match" )
    endif()

else()
    message( FATAL_ERROR "The CC envoronment variable is not set" )
endif()

# test if there is a fortran compiler
if( DEFINED ENV{FC} )

    # get version of fortran compiler
    execute_process( COMMAND bash "-c" "$FC --version | head -1 | cut -d' ' -f 3" OUTPUT_VARIABLE BELFEM_FC_VERSION )

    # tidy up string
    string(STRIP "${BELFEM_FC_VERSION}" BELFEM_FC_VERSION)

    # test if versions match
    if( NOT BELFEM_CXX_VERSION VERSION_EQUAL BELFEM_FC_VERSION )
        message( FATAL_ERROR "The version numbers of icc and ifort do not match" )
    endif()
else()
    message( FATAL_ERROR "The FC envoronment variable is not set" )
endif()

#-------------------------------------------------------------------------------
# Find Library Directory for the Intel compiler
#-------------------------------------------------------------------------------

# step 1: resolve any symlink ( can't use readlink -f, because this must work on mac also )
execute_process( COMMAND bash "-c" "TPATH=$( which $CXX ) && while [ -L $TPATH ]; do TPATH=$( readlink $TPATH ) ; done && echo $TPATH" OUTPUT_VARIABLE BELFEM_TPATH )

# tidy up string
string(STRIP "${BELFEM_TPATH}" BELFEM_TPATH)

# step 2: get the home directory
execute_process( COMMAND bash "-c" "dirname $( dirname ${BELFEM_TPATH} )" OUTPUT_VARIABLE BELFEM_TPATH )
string(STRIP "${BELFEM_TPATH}" BELFEM_TPATH)

if ( APPLE )
    message( FATAL_ERROR "Intel on Mac currently not supported" )
endif()

# step 3: find library directory
if( IS_DIRECTORY ${INTELROOT}/oneapi/compiler/latest/linux/lib )
    set( BELFEM_INTEL_LIBS ${INTELROOT}/oneapi/compiler/latest/linux/lib )
else()
    message( FATAL_ERROR "Could not find Intel library directory" )
endif()

# add this directory to the rpath
list( APPEND BELFEM_RPATH ${BELFEM_INTEL_LIBS} )

# test if we want to use mpi

if( USE_MPI )
    list( APPEND BELFEM_DEFS "BELFEM_MPI" )

    # check of mpicc exists
    execute_process( COMMAND bash "-c" "which mpicc" OUTPUT_VARIABLE BELFEM_MPICC )
    string( STRIP "${BELFEM_MPICC}" BELFEM_MPICC )

    if( NOT EXISTS ${BELFEM_MPICC} )
        message( FATAL_ERROR "Could not find mpicc executable" )
    endif()

    # get the mpicc version number
    execute_process( COMMAND bash "-c" "mpicc -dumpversion" OUTPUT_VARIABLE BELFEM_MPICC_VERSION )

    # tidy up string
    string( STRIP "${BELFEM_MPICC_VERSION}" BELFEM_MPICC_VERSION )

    # make sure that version number is the same
    if( NOT BELFEM_CC_VERSION VERSION_EQUAL BELFEM_MPICC_VERSION )
        message( FATAL_ERROR "The version numbers of the icc and mpicc do not match" )
    endif()

    # check of mpic++ exists
    execute_process( COMMAND bash "-c" "which mpicxx" OUTPUT_VARIABLE BELFEM_MPICXX )
    string( STRIP "${BELFEM_MPICXX}" BELFEM_MPICXX )

    if( NOT EXISTS ${BELFEM_MPICXX} )
        message( FATAL_ERROR "Could not find mpicxx executable" )
    endif()

    # get the mpicxx version number
    execute_process( COMMAND bash "-c" "mpicxx -dumpversion" OUTPUT_VARIABLE BELFEM_MPICXX_VERSION )

    # tidy up string
    string( STRIP "${BELFEM_MPICXX_VERSION}" BELFEM_MPICXX_VERSION )

    # make sure that version number is the same
    if( NOT BELFEM_CXX_VERSION VERSION_EQUAL BELFEM_MPICXX_VERSION )
        message( FATAL_ERROR "The version numbers of icpc and mpicxx do not match" )
    endif()

    # check of mpifortran exists
    execute_process( COMMAND bash "-c" "which mpifort" OUTPUT_VARIABLE BELFEM_MPIFORT )
    string( STRIP "${BELFEM_MPIFORT}" BELFEM_MPIFORT )

    if( NOT EXISTS ${BELFEM_MPIFORT} )
        message( FATAL_ERROR "Could not find mpifort executable" )
    endif()

    # get the mpiifortversion number
    execute_process( COMMAND bash "-c" "mpifort --version |  head -1 | cut -d' ' -f 3" OUTPUT_VARIABLE BELFEM_MPIFORT_VERSION )

    # tidy up string
    string( STRIP "${BELFEM_MPIFORT_VERSION}" BELFEM_MPIFORT_VERSION )

    # make sure that version number is the same
    if( NOT BELFEM_FC_VERSION VERSION_EQUAL BELFEM_MPIFORT_VERSION )
        message( FATAL_ERROR "The version numbers of the ifort and mpifort do not match" )
    endif()

    # now set the compiler variables
    set( CMAKE_C_COMPILER  mpicc )
    set( CMAKE_CXX_COMPILER mpicxx )
    set( CMAKE_Fortran_COMPILER mpifort )

else()
    # now set the compiler variables
    set( CMAKE_C_COMPILER  "$ENV{CC}" )
    set( CMAKE_CXX_COMPILER "$ENV{CXX}" )
    set( CMAKE_Fortran_COMPILER "$ENV{FC}" )
endif()

set( BELFEM_COMPILER_VERSION ${BELFEM_CC_VERSION} )