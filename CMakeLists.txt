cmake_minimum_required(VERSION 3.11)
project(belfem)
enable_language(Fortran)


# -------------------------------------------------------------------------
# User Defined Parameters
# -------------------------------------------------------------------------

#set( MKLROOT "/opt/intel/oneapi/mkl/current" )
set( INTELROOT "/opt/intel" )
set (CMAKE_CXX_STANDARD 14)

# -------------------------------------------------------------------------
# User Settings
# -------------------------------------------------------------------------

option( USE_MPI "Use Message Parsing Library" OFF )

option( USE_OPENMP "Use OpenMP" off )
option( USE_MKL "Use Intel Math Kernel Library" OFF )
option( USE_MKL_64BIT_API "Use the 64 bit api instead of the 32 bit one" OFF )
if( Apple )
    option( USE_MATRIX_ARMADILLO "Use Armadillo Matrix Library." OFF )
    option( USE_MATRIX_BLAZE "Use Blaze Matrix Library." ON )
else()
    option( USE_MATRIX_ARMADILLO "Use Armadillo Matrix Library." ON )
    option( USE_MATRIX_BLAZE "Use Blaze Matrix Library." OFF )
endif()
option( USE_SUITESPARSE "Use Suitesparse Library." ON )
option( USE_MUMPS "Use Mumps Library." OFF )
option( USE_PARDISO "Use Pardiso Library." OFF )
option( USE_PETSC "Use PETSc Library." OFF )
# option( USE_NLOPT "Use NLOPT library" OFF )
option( USE_HDF5 "Use HDF5 library" ON )
option( USE_EXODUS "Use Exodus library" ON )

option( USE_DEBUG "Compile with debug flags" ON )
option( USE_WARNINGS "Use pedantic warnings" ON )
option( USE_PROFILER "Use Google Profiling Tools" OFF )
option( USE_TEST "Build Tests" OFF )
option( USE_EXAMPLES "Build Examples" ON )
option( USE_MAXWELL "Use Maxwell Modules" ON )
option( USE_GASMODELS "Build Gas Library" OFF )
option( USE_STRUMPACK "Link against STRUMPACK" OFF )
option( USE_STRUMPACK_CUDA "Link against STRUMPACK with CUDA" OFF )

# -------------------------------------------------------------------------
# Global Parameters
# -------------------------------------------------------------------------

# Path to cmake configiration
set( BELFEM_CONFIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/config )

# Path to main sources
set( BELFEM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )

# make sure that this is a unix system
include( ${BELFEM_CONFIG_DIR}/system/find_unix.cmake )

# compiler setup
include( ${BELFEM_CONFIG_DIR}/compiler/detect_compiler.cmake )

# third party libraries
include( ${BELFEM_CONFIG_DIR}/system/find_tplibs.cmake )

# message parsing interface
include( ${BELFEM_CONFIG_DIR}/system/find_mpi.cmake )

# basic linear algebra
include( ${BELFEM_CONFIG_DIR}/linalg/config_mkl.cmake )

# matrix interface
include( ${BELFEM_CONFIG_DIR}/linalg/config_matrix.cmake )

# solver interface
include( ${BELFEM_CONFIG_DIR}/linalg/config_suitesparse.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_mumps.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_strumpack.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_petsc.cmake )

# input output
# xml is always on
set( BELFEM_IO_LIBS "-ltinyxml2" )
include( ${BELFEM_CONFIG_DIR}/io/config_hdf5.cmake )
include( ${BELFEM_CONFIG_DIR}/io/config_exodus.cmake )

if( USE_GASMODELS )
    list( APPEND BELFEM_DEFS "BELFEM_GASMODELS" )
endif()

# -------------------------------------------------------------------------
# profiler
# -------------------------------------------------------------------------

# make sure that this is a unix system
include( ${BELFEM_CONFIG_DIR}/system/config_gperftools.cmake )

# -------------------------------------------------------------------------
# finalize compiler
# -------------------------------------------------------------------------

# finish up compiler flags
include( ${BELFEM_CONFIG_DIR}/compiler/finalize_compiler.cmake )

# -------------------------------------------------------------------------
# guard against in-source builds
# -------------------------------------------------------------------------
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# -------------------------------------------------------------------------
# MY LIBS
# -------------------------------------------------------------------------

include( ${BELFEM_CONFIG_DIR}/globals.cmake )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${BINDIR} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${LIBDIR} )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}" )
if( TPLSLIBDIR  )
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${TPLSLIBDIR}" )
endif()

add_subdirectory( src )

# -------------------------------------------------------------------------
# INCLUDE LIB DIRECTORIES
# -------------------------------------------------------------------------

include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/containers )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/math )
include_directories( ${BELFEM_SOURCE_DIR}/sparse )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )

list( APPEND BELFEM_LIBLIST core )
list( APPEND BELFEM_LIBLIST comm )
list( APPEND BELFEM_LIBLIST io )

list( APPEND BELFEM_LIBLIST graph )
list( APPEND BELFEM_LIBLIST sparse )
list( APPEND BELFEM_LIBLIST spline )

# -------------------------------------------------------------------------
# ASSEMBLE LIBLIST
# -------------------------------------------------------------------------

foreach( LIBITEM ${BELFEM_LIBLIST} )
    set( BELFEM_LIBS "-l${LIBPREFIX}_${LIBITEM} ${BELFEM_LIBS}" )
endforeach()

include( ${BELFEM_CONFIG_DIR}/summary.cmake )
# -------------------------------------------------------------------------
# EXECUTABLE
# -------------------------------------------------------------------------

set( LIBLIST )

set( EXECNAME banner )
set( MAIN     main.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

# -------------------------------------------------------------------------
# CLEAN SCRIPT
# -------------------------------------------------------------------------
add_custom_target(reset
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND bash -c "if [ -d ./bin  ]; then rm -rvf ./bin && mkdir bin ; fi ; if [ -d ./lib  ]; then rm -rvf ./lib && mkdir lib ; fi "
        COMMENT "Removing binaries before rebuild"
        VERBATIM
         )

# -------------------------------------------------------------------------
# TEST
# -------------------------------------------------------------------------
if( USE_TEST AND NOT BELFEM_USE_CLANG )
    enable_testing()
    add_subdirectory( test )

add_custom_target(tests
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND bash -c "pushd  test && make && popd && ctest "
        VERBATIM
        )
endif()

# -------------------------------------------------------------------------
# DOC
# -------------------------------------------------------------------------
# Find Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
    # Input proper paths #> need to replace CMAKE_SOURCE_DIR in there
    configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

    # Create a doxygen make target
    add_custom_target(doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating belfem documentation with Doxygen."
            VERBATIM )
endif()
