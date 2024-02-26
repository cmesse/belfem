#-------------------------------------------------------------------------------
# check the other executables
#-------------------------------------------------------------------------------
# get the version of the c++ compiler
if( DEFINED ENV{CXX} )
    execute_process( COMMAND bash "-c" "$CXX -dumpversion" OUTPUT_VARIABLE BELFEM_CXX_VERSION )
else()
    execute_process( COMMAND bash "-c" "g++ -dumpversion" OUTPUT_VARIABLE BELFEM_CXX_VERSION )
endif()

# tidy up string
string(STRIP "${BELFEM_CXX_VERSION}" BELFEM_CXX_VERSION)

# test if there is a c compiler
if( DEFINED ENV{CC} )
    # get version of c-compiler
    execute_process( COMMAND bash "-c" "$CC -dumpversion" OUTPUT_VARIABLE BELFEM_CC_VERSION )
    execute_process( COMMAND bash "-c" "which $CC" OUTPUT_VARIABLE BELFEM_CC_PATH )
else()
    execute_process( COMMAND bash "-c" "gcc -dumpversion" OUTPUT_VARIABLE BELFEM_CC_VERSION )
    execute_process( COMMAND bash "-c" "which gcc" OUTPUT_VARIABLE BELFEM_CC_PATH )
endif()

# tidy up string
string(STRIP "${BELFEM_CC_VERSION}" BELFEM_CC_VERSION)

# test if versions match
if( NOT BELFEM_CXX_VERSION VERSION_EQUAL BELFEM_CC_VERSION AND NOT BELFEM_USE_CLANG )
    message( FATAL_ERROR "The version numbers of gcc and g++ do not match" )
endif()


# test if there is a fortran compiler
if( DEFINED ENV{FC} )
    # get version of fortran compiler
    execute_process( COMMAND bash "-c" "$FC -dumpversion" OUTPUT_VARIABLE BELFEM_FC_VERSION )
else()
    execute_process( COMMAND bash "-c" "gfortran -dumpversion" OUTPUT_VARIABLE BELFEM_FC_VERSION )
endif()

# tidy up string
string(STRIP "${BELFEM_FC_VERSION}" BELFEM_FC_VERSION)

# test if versions match
if( NOT BELFEM_CXX_VERSION VERSION_EQUAL BELFEM_FC_VERSION AND NOT BELFEM_USE_CLANG )
    message( FATAL_ERROR "The version numbers of g++ and gfortran do not match" )
endif()


#-------------------------------------------------------------------------------
# Find Library Directory for the GNU compiler
#-------------------------------------------------------------------------------

# step 1: resolve any symlink ( can't use readlink -f, because this must work on mac also )
execute_process( COMMAND bash "-c" "TPATH=$( which $CXX ) && while [ -L $TPATH ]; do TPATH=$( readlink $TPATH ) ; done && echo $TPATH" OUTPUT_VARIABLE BELFEM_TPATH )

# tidy up string
string(STRIP "${BELFEM_TPATH}" BELFEM_TPATH)

# step 2: get the home directory
execute_process( COMMAND bash "-c" "dirname $( dirname ${BELFEM_TPATH} )" OUTPUT_VARIABLE BELFEM_TPATH )
string(STRIP "${BELFEM_TPATH}" BELFEM_TPATH)

# step 3: find library directory
if( IS_DIRECTORY ${BELFEM_TPATH}/lib64 )
    set( BELFEM_GCC_LIBS ${BELFEM_TPATH}/lib64 )
    list( APPEND BELFEM_RPATH ${BELFEM_GCC_LIBS} )
elseif( IS_DIRECTORY ${BELFEM_TPATH}/lib )
    set( BELFEM_GCC_LIBS ${BELFEM_TPATH}/lib )
    list( APPEND BELFEM_RPATH ${BELFEM_GCC_LIBS} )

endif()

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
        message( FATAL_ERROR "The version numbers of gcc and mpicc do not match gcc: ${BELFEM_CC_VERSION} mpicc: ${BELFEM_MPICC_VERSION}" )
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
        message( FATAL_ERROR "The version numbers of g++ and mpicxx do not match" )
    endif()

    # check of mpifortran exists
    execute_process( COMMAND bash "-c" "which mpifort" OUTPUT_VARIABLE BELFEM_MPIFORT )
    string( STRIP "${BELFEM_MPIFORT}" BELFEM_MPIFORT )

    if( NOT EXISTS ${BELFEM_MPIFORT} )
        message( FATAL_ERROR "Could not find mpifort executable" )
    endif()

    # get the mpicxx version number
    execute_process( COMMAND bash "-c" "mpifort -dumpversion" OUTPUT_VARIABLE BELFEM_MPIFORT_VERSION )

    # tidy up string
    string( STRIP "${BELFEM_MPIFORT_VERSION}" BELFEM_MPIFORT_VERSION )

    # make sure that version number is the same
    if( NOT BELFEM_FC_VERSION VERSION_EQUAL BELFEM_MPIFORT_VERSION )
        message( FATAL_ERROR "The version numbers of gfortran and mpifort do not match" )
    endif()

    # now set the compiler variables
    set( CMAKE_C_COMPILER  mpicc )
    set( CMAKE_CXX_COMPILER mpicxx )
    set( CMAKE_Fortran_COMPILER mpifort )

elseif((DEFINED ENV{CXX}))
    set( CMAKE_C_COMPILER  "$ENV{CC}" )
    set( CMAKE_CXX_COMPILER "$ENV{CXX}" )
    set(b CMAKE_Fortran_COMPILER "$ENV{FC}" )
elseif((DEFINED ENV{CXX}))
    if( BELFEM_USE_CLANG )
        set( CMAKE_C_COMPILER  "clang" )
        set( CMAKE_CXX_COMPILER "clang++" )
    else()
        set( CMAKE_C_COMPILER  "gcc" )
        set( CMAKE_CXX_COMPILER "g++" )
    endif()
    set( CMAKE_Fortran_COMPILER "gfortran" )
endif()

set( BELFEM_COMPILER_VERSION ${BELFEM_CC_VERSION} )