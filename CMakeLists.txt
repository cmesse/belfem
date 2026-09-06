cmake_minimum_required(VERSION 3.11)

# Use the SCLS toolchain when no compiler was explicitly requested.
#
# This MUST stay above project(): by the time project() returns, the compilers
# are already detected and written to the cache, and setting these variables
# afterwards has no effect on the detection. The cache entries carry no FORCE,
# so an explicit -DCMAKE_C_COMPILER=... on the command line still wins.
#
# Why this block is needed at all: /opt/scls/bin ships gcc, g++ and c++ but no
# "cc", and CMake's default C-compiler search order is cc;gcc;... . It therefore
# finds /usr/bin/cc (AppleClang) even when the SCLS directory comes first in
# PATH, then uses /usr/bin as the hint for the C++ search - leaving a mismatched
# AppleClang / gfortran pair.
if( NOT DEFINED CMAKE_CXX_COMPILER AND NOT DEFINED ENV{CXX} AND DEFINED ENV{SCLS} )
    if( EXISTS "$ENV{SCLS}/bin/g++" )
        set( CMAKE_C_COMPILER       "$ENV{SCLS}/bin/gcc"      CACHE FILEPATH "C compiler" )
        set( CMAKE_CXX_COMPILER     "$ENV{SCLS}/bin/g++"      CACHE FILEPATH "C++ compiler" )
        set( CMAKE_Fortran_COMPILER "$ENV{SCLS}/bin/gfortran" CACHE FILEPATH "Fortran compiler" )
    endif()
endif()


# BELFEM semantic version (semver.org: MAJOR.MINOR.PATCH). This is the single
# canonical source of the version number; bump it here on release. CMake exposes
# it as PROJECT_VERSION and PROJECT_VERSION_MAJOR/MINOR/PATCH, which are fed into
# the generated belfem_version.hpp header (see src/core/CMakeLists.txt).
project(belfem VERSION 0.9.0)
enable_language(Fortran)


if( POLICY CMP0156 )
    cmake_policy( SET CMP0156 NEW )
endif()

if( POLICY CMP0179 )
    cmake_policy( SET CMP0179 NEW )
endif()

# match the deployment target to the OS so linker warnings about
# mismatched dylib versions are avoided
if( APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET )
    set( CMAKE_OSX_DEPLOYMENT_TARGET "$ENV{MACOSX_DEPLOYMENT_TARGET}" CACHE STRING "" )
    if( NOT CMAKE_OSX_DEPLOYMENT_TARGET )
        execute_process( COMMAND sw_vers -productVersion
                         OUTPUT_VARIABLE _macos_ver OUTPUT_STRIP_TRAILING_WHITESPACE )
        set( CMAKE_OSX_DEPLOYMENT_TARGET "${_macos_ver}" CACHE STRING "" FORCE )
    endif()
endif()

# always export compile_commands.json so clangd / IDE tooling stays in sync with
# the active backend (e.g. Blaze vs Armadillo). FORCE re-emits it on every
# configure, so touching cmake/ccmake can never leave a stale command database.
set( CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json for tooling" FORCE )

# -------------------------------------------------------------------------
# User Defined Parameters
# -------------------------------------------------------------------------

#set( MKLROOT "/opt/intel/oneapi/mkl/current" )
set( INTELROOT "/opt/intel" )
set (CMAKE_CXX_STANDARD 17)

# -------------------------------------------------------------------------
# User Settings
# -------------------------------------------------------------------------

# Path to cmake configuration
set( BELFEM_CONFIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/config )

# The SCLS flavor supplies the defaults for USE_MKL, USE_PARDISO and USE_DEBUG,
# so it has to be read before those options are declared: option() sets a
# default only when the cache entry does not exist yet. Presets nothing when
# $SCLS is unset; an explicit -D or a ccmake edit always wins.
include( ${BELFEM_CONFIG_DIR}/system/find_scls_flavor.cmake )

option( USE_MPI "Use Message Parsing Library" ON )

option( USE_OPENMP "Use OpenMP" ON )

# BELFEM's OWN kernels threading is a SEPARATE switch from OpenMP availability.
# USE_OPENMP stays ON: the third-party solvers ( STRUMPACK, MKL, threaded BLAS,
# PARDISO ) need -fopenmp and their thread budget is queried under the OMP define.
# USE_BELFEM_OPENMP governs only the !$omp directives in BELFEM's three Fortran
# kernels, and it is OFF because they never earned their keep: those matvecs are
# master-only and are dwarfed by assembly and factorization, while the array
# reduction in matvec_csc put an n * 8 byte private copy on a 2 MiB OpenMP worker
# stack and crashed above 262,144 rows. See doc/parallel_execution.md.
option( USE_BELFEM_OPENMP "Let BELFEM's own Fortran kernels use OpenMP (see doc/parallel_execution.md)" OFF )
option( USE_MKL "Use Intel Math Kernel Library" ${BELFEM_DEFAULT_USE_MKL} )
option( USE_MKL_64BIT_API "Use the 64 bit api instead of the 32 bit one" OFF )
if( APPLE )
    option( USE_MATRIX_ARMADILLO "Use Armadillo Matrix Library." OFF )
    option( USE_MATRIX_BLAZE "Use Blaze Matrix Library." ON )
else()
    option( USE_MATRIX_ARMADILLO "Use Armadillo Matrix Library." ON )
    option( USE_MATRIX_BLAZE "Use Blaze Matrix Library." OFF )
endif()
option( USE_SUITESPARSE "Use Suitesparse Library (not for binary redistribution)" OFF )
option( USE_SUPERLU "Use SuperLU Library." ON )
option( USE_MUMPS "Use MUMPS Library." ON )
option( USE_PARDISO "Use PARDISO Library." ${BELFEM_DEFAULT_USE_PARDISO} )
option( USE_PETSC "Use PETSc Library." ON )
option( USE_NLOPT "Use NLOPT library" ON )
option( USE_HDF5 "Use HDF5 library" ON )
option( USE_EXODUS "Use Exodus library" ON )
option( USE_METIS "Use METIS library" ON )
option( USE_SCOTCH "Use SCOTCH library" ON )
option( USE_DEBUG "Compile with debug flags" ${BELFEM_DEFAULT_USE_DEBUG} )
option( USE_WARNINGS "Use pedantic warnings" ON )
option( USE_PROFILER "Use Google Profiling Tools" OFF )
option( USE_TEST "Build Tests" ON )
option( USE_SHARED_LIBS "Build the project library libbelfem as a shared library" OFF )
option( USE_EXAMPLES "Build Coding Examples" OFF )
option( USE_MAXWELL "Use Maxwell Modules" ON )
option( USE_GASMODELS "Build Gas Library" OFF )
option( USE_COMBUSTION "Build finite-rate combustion (nonfree; needs USE_NONFREE and USE_GASMODELS)" OFF )
option( USE_STRUMPACK "Link against STRUMPACK" ON )
option( USE_STRUMPACK_CUDA "Link against STRUMPACK with CUDA" OFF )
option( USE_STRUMPACK_ZFP "Link against STRUMPACK with ZFP" ON )
option( USE_ARPACK "Link against ARPACK-ng" ON )
option( USE_PARPACK "Link against PARPACK ( needs USE_ARPACK and USE_MPI )" ON )
option( USE_VTK "Use Visualization Toolkit" OFF )
option( USE_NONFREE "Compile nonfree modules" OFF )
option( USE_TINYXML2 "Use TinyXML2" OFF )

# -------------------------------------------------------------------------
# Global Parameters
# -------------------------------------------------------------------------

# USE_DEBUG is the single mode switch; CMAKE_BUILD_TYPE follows it so IDEs
# and tools that key off the standard variable agree with the build. A
# conflicting value passed on the command line is overridden with a notice.
if( USE_DEBUG )
    set( BELFEM_BUILD_TYPE "Debug" )
else()
    set( BELFEM_BUILD_TYPE "Release" )
endif()
if( CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE STREQUAL BELFEM_BUILD_TYPE )
    message( STATUS "CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}' overridden to '${BELFEM_BUILD_TYPE}' ( derived from USE_DEBUG )" )
endif()
set( CMAKE_BUILD_TYPE ${BELFEM_BUILD_TYPE} CACHE STRING "derived from USE_DEBUG" FORCE )

# the optimization policy is hand-set in config/compiler ( -Og -g debug,
# -O2 release ); blank the per-config flag sets so CMake cannot append
# -O3 or a stray -g behind the policy's back
set( CMAKE_C_FLAGS_DEBUG "" )
set( CMAKE_CXX_FLAGS_DEBUG "" )
set( CMAKE_Fortran_FLAGS_DEBUG "" )
set( CMAKE_C_FLAGS_RELEASE "" )
set( CMAKE_CXX_FLAGS_RELEASE "" )
set( CMAKE_Fortran_FLAGS_RELEASE "" )

if( USE_MKL_64BIT_API )
    list( APPEND BELFEM_DEFS "BELFEM_INT64" )
endif ()

# Path to main sources
set( BELFEM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )

# Path to nonfree sources
set( BELFEM_NONFREE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/nonfree )

# make sure that this is a unix system
include( ${BELFEM_CONFIG_DIR}/system/find_unix.cmake )

# package finder helper
include( ${BELFEM_CONFIG_DIR}/scripts/belfem_find_package.cmake )

# compiler setup
include( ${BELFEM_CONFIG_DIR}/compiler/detect_compiler.cmake )

# BELFEM_OMP gates BELFEM's OWN OpenMP directives ( src/sparse/splinalg.f90,
# arpacktools.f90, parpacktools.f90 ) and NOTHING else. It requires OpenMP to be
# available at all, so it is the AND of the two options. The plain OMP define,
# appended by the compiler config above, stays independent of this and keeps
# guarding the third-party thread queries -- see the comments at
# src/sparse/cl_SolverWrapper.cpp ( hatch_turtle ), src/sparse/pardisotools.f90
# and src/core/banner.cpp, all of which MUST remain on OMP.
if( USE_OPENMP AND USE_BELFEM_OPENMP )
    list( APPEND BELFEM_DEFS "BELFEM_OMP" )
endif()

# third party libraries
include( ${BELFEM_CONFIG_DIR}/system/find_scls.cmake )

# message parsing interface
include( ${BELFEM_CONFIG_DIR}/system/find_mpi.cmake )

# basic linear algebra
include( ${BELFEM_CONFIG_DIR}/linalg/config_mkl.cmake )

# matrix interface
include( ${BELFEM_CONFIG_DIR}/linalg/config_matrix.cmake )

# arpack interface
include( ${BELFEM_CONFIG_DIR}/linalg/config_arpack.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_parpack.cmake )

# solver interface
include( ${BELFEM_CONFIG_DIR}/linalg/config_suitesparse.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_superlu.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_mumps.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_strumpack.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_petsc.cmake )

# graph interface
include( ${BELFEM_CONFIG_DIR}/linalg/config_metis.cmake )
include( ${BELFEM_CONFIG_DIR}/linalg/config_scotch.cmake )

# input output
include( ${BELFEM_CONFIG_DIR}/io/config_hdf5.cmake )
include( ${BELFEM_CONFIG_DIR}/io/config_exodus.cmake )
include( ${BELFEM_CONFIG_DIR}/io/config_tinyxml2.cmake )

if( USE_GASMODELS )
    list( APPEND BELFEM_DEFS "BELFEM_GASMODELS" )
endif()

# Finite-rate combustion is nonfree and builds on the gas models, so it needs
# both trees. The gas models themselves no longer do.
if( USE_COMBUSTION )
    if( NOT USE_NONFREE )
        message( FATAL_ERROR
                 "USE_COMBUSTION requires USE_NONFREE: combustion lives in nonfree/physics" )
    endif()
    if( NOT USE_GASMODELS )
        message( FATAL_ERROR
                 "USE_COMBUSTION requires USE_GASMODELS: combustion links gastables and gasmodels" )
    endif()
    list( APPEND BELFEM_DEFS "BELFEM_COMBUSTION" )
endif()


if( USE_VTK )
    # calls find_package( VTK ) with the proper hints
    include( ${BELFEM_CONFIG_DIR}/system/find_vtk.cmake )
    list( APPEND BELFEM_DEFS "BELFEM_VTK" )
endif ()

# numerics libraries
include( ${BELFEM_CONFIG_DIR}/numerics/find_nlopt.cmake )

# -------------------------------------------------------------------------
# profiler
# -------------------------------------------------------------------------

# profiler configuration
include( ${BELFEM_CONFIG_DIR}/system/config_gperftools.cmake )

# -------------------------------------------------------------------------
# finalize compiler
# -------------------------------------------------------------------------

# finish up compiler flags
include( ${BELFEM_CONFIG_DIR}/compiler/finalize_compiler.cmake )

# -------------------------------------------------------------------------
# Install layout and RPATH
# -------------------------------------------------------------------------

# lib vs lib64 follows the platform ( el9 installs to lib64 ); the install
# rules and the RPATH below both derive from CMAKE_INSTALL_LIBDIR so the two
# cannot disagree
include( GNUInstallDirs )

# build tree: binaries run without LD_LIBRARY_PATH / DYLD_LIBRARY_PATH, the
# third-party directories collected in BELFEM_RPATH are recorded in the
# binary. install tree: the project library is found relative to the
# executable, so a staged install ( DESTDIR ) stays relocatable, and the
# third-party directories follow. Set through the CMake variables, never as a
# raw -Wl,-rpath flag: $ORIGIN would be expanded by the shell.
set( CMAKE_SKIP_BUILD_RPATH OFF )
set( CMAKE_BUILD_WITH_INSTALL_RPATH OFF )
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH ON )
set( CMAKE_BUILD_RPATH ${BELFEM_RPATH} )
if( APPLE )
    set( BELFEM_RPATH_ORIGIN "@loader_path" )
    set( CMAKE_MACOSX_RPATH ON )
    set( CMAKE_INSTALL_NAME_DIR "@rpath" )
else()
    set( BELFEM_RPATH_ORIGIN "$ORIGIN" )
endif()
set( CMAKE_INSTALL_RPATH "${BELFEM_RPATH_ORIGIN}/../${CMAKE_INSTALL_LIBDIR}" ${BELFEM_RPATH} )

# -------------------------------------------------------------------------
# guard against in-source builds
# -------------------------------------------------------------------------
if( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" )
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# -------------------------------------------------------------------------
# pin the generator
#
# BELFEM is built with make. A Ninja tree looks fine in the IDE but leaves the
# terminal without a working Makefile, so reject it before anything is written.
# -------------------------------------------------------------------------
option( ALLOW_NINJA "permit generators other than Unix Makefiles" OFF )

if( NOT CMAKE_GENERATOR MATCHES "Makefiles" AND NOT ALLOW_NINJA )
    message( FATAL_ERROR
            "BELFEM expects the 'Unix Makefiles' generator, but this tree was "
            "configured with '${CMAKE_GENERATOR}'.\n"
            "  CLion: Settings -> Build, Execution, Deployment -> CMake -> "
            "Debug -> Generator: Unix Makefiles, then Reset Cache and Reload.\n"
            "  Shell: cmake -G \"Unix Makefiles\" -S <source> -B <build>\n"
            "Pass -DALLOW_NINJA=ON to override." )
endif()

# -------------------------------------------------------------------------
# guard against generator collisions
#
# Resetting the cache from an IDE wipes CMakeCache.txt and CMakeFiles/, but
# leaves the top-level Makefile behind. If the tree is then regenerated with
# Ninja, a terminal 'make' still finds that orphan, re-runs cmake (which
# writes Ninja files), and dies on the missing CMakeFiles/Makefile2. Drop the
# entry points of whichever generator is not the active one.
# -------------------------------------------------------------------------
if( CMAKE_GENERATOR MATCHES "Ninja" )
    set( STALE_BUILDFILES Makefile CMakeFiles/Makefile2 CMakeFiles/Makefile.cmake )
elseif( CMAKE_GENERATOR MATCHES "Makefiles" )
    set( STALE_BUILDFILES build.ninja .ninja_deps .ninja_log CMakeFiles/rules.ninja )
else()
    set( STALE_BUILDFILES )
endif()

foreach( STALEFILE ${STALE_BUILDFILES} )
    if( EXISTS ${CMAKE_BINARY_DIR}/${STALEFILE} )
        message( STATUS "Removing ${STALEFILE}: leftover from a different generator" )
        file( REMOVE ${CMAKE_BINARY_DIR}/${STALEFILE} )
    endif()
endforeach()

# -------------------------------------------------------------------------
# Project library paths
# -------------------------------------------------------------------------

include( ${BELFEM_CONFIG_DIR}/globals.cmake )

# every object is position independent, in all three languages: the project
# library may be built shared, and the Fortran kernels in sparse/ end up in it
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${BINDIR} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${LIBDIR} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${LIBDIR} )
if( SCLSLIBDIR )
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${SCLSLIBDIR}" )
endif()

add_subdirectory( src )

if ( USE_NONFREE )
    add_subdirectory( nonfree )
endif()

# the project library, assembled from every module registered above
include( ${BELFEM_CONFIG_DIR}/scripts/Add_BelfemLibrary.cmake )

# -------------------------------------------------------------------------
# INSTALL: headers and data. The library and the executables install
# themselves where they are defined ( Add_BelfemLibrary.cmake,
# Add_Executable.cmake ), the rest of the tree is deployed here.
# -------------------------------------------------------------------------

# public headers, module layout preserved: <prefix>/include/belfem/<module>/...
install( DIRECTORY ${BELFEM_SOURCE_DIR}/
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${LIBPREFIX}
         FILES_MATCHING
         PATTERN "*.hpp"
         PATTERN "doc" EXCLUDE
         PATTERN "CMakeFiles" EXCLUDE )

# generated version header
install( FILES ${CMAKE_BINARY_DIR}/generated/belfem_version.hpp
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${LIBPREFIX} )

# data files and example decks. install( DIRECTORY ) walks the tree at
# install time, so a file added to share/ or examples/ after the configure is
# picked up by the next `make install` with no re-run of cmake. The trailing
# slash on the source keeps the directory name out of the destination.
# Deleted files are not removed from an existing prefix.
# python bytecode caches appear beside share/python/database/howtoread.py as
# soon as anyone imports it from the source tree; they must not be deployed
install( DIRECTORY ${CMAKE_SOURCE_DIR}/share/
         DESTINATION ${CMAKE_INSTALL_DATADIR}/${LIBPREFIX}
         USE_SOURCE_PERMISSIONS
         PATTERN "__pycache__"         EXCLUDE
         PATTERN "*.pyc"               EXCLUDE )

# run output is excluded; the pattern list mirrors examples/scripts/Allclean —
# keep the two in step
install( DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
         DESTINATION ${CMAKE_INSTALL_DATADIR}/${LIBPREFIX}/examples
         USE_SOURCE_PERMISSIONS
         PATTERN "*.exo"               EXCLUDE
         PATTERN "*.e-s.*"             EXCLUDE
         PATTERN "*.bfm"               EXCLUDE
         PATTERN "*.csv"               EXCLUDE
         PATTERN "memdump.hdf5"        EXCLUDE
         PATTERN "belfem.log"          EXCLUDE
         PATTERN "out.txt"             EXCLUDE
         PATTERN "err.txt"             EXCLUDE
         PATTERN "CircuitResults.txt"  EXCLUDE
         PATTERN "circuitAnalysis.out" EXCLUDE
         PATTERN "slurm-*"             EXCLUDE )

# plugin build templates and the annotated example plugin sources. They live
# beside the code they document, so the header rule above -- which matches
# "*.hpp" -- does not pick them up, and a prefix would otherwise carry the
# plugin API with no recipe for building against it. The templates are the only
# build files that know the installed <prefix>/include/belfem/<module> layout.
install( FILES
         ${BELFEM_SOURCE_DIR}/physics/materials/UserMaterialTemplate.cmake
         ${BELFEM_SOURCE_DIR}/physics/materials/UserLibraryTemplate.cmake
         ${BELFEM_SOURCE_DIR}/physics/materials/example_user_material.cpp
         ${BELFEM_SOURCE_DIR}/physics/materials/example_user_defect.cpp
         ${BELFEM_SOURCE_DIR}/numerics/sources/example_user_source.cpp
         DESTINATION ${CMAKE_INSTALL_DATADIR}/${LIBPREFIX}/templates )

# every installed header refers the reader to "the top-level LICENSE file";
# ship it at the prefix root so that sentence is true on an installed copy
install( FILES ${CMAKE_SOURCE_DIR}/LICENSE
         DESTINATION . )

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
        COMMAND ${CMAKE_COMMAND} -E rm -rf bin lib test
        COMMAND ${CMAKE_COMMAND} -E make_directory bin lib test
        COMMENT "Removing binaries before rebuild"
        VERBATIM
         )

# -------------------------------------------------------------------------
# TEST
# -------------------------------------------------------------------------
if( USE_TEST )
    enable_testing()
    add_subdirectory( tests )

    # 'make check' builds all test executables, then runs ctest
    add_custom_target(check
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
            COMMENT "Building and running all tests"
            VERBATIM
            )

    # The TARGET guard makes optional suites safe to list: gastables and
    # gasmodels only exist under USE_GASMODELS (tests/physics/CMakeLists.txt),
    # and absent targets are skipped, not errors.
    # commmpi and sparsempi are Tier 2 ( multi-rank ) binaries rather than
    # directories: they live in tests/comm and tests/sparse beside the serial
    # suites. They are listed here and NOT in check-fast below, so `make check`
    # builds and runs them while `check-fast` stays serial. `make check` runs
    # ctest unfiltered, so the mpi-labelled tests are picked up without a second
    # invocation; a developer who wants to skip them on a small machine can run
    # `ctest -LE mpi` by choice.
    # A Tier 2 binary that is registered with ctest but missing from THIS list
    # fails as "Could not find executable" - measured 2026-09-04 with sparsempi.
    foreach( TESTDIR core containers linalg comm commmpi math sparse sparsempi mesh fem io ode
                     physics gastables gasmodels homology circuit )
        if( TARGET test_${TESTDIR} )
            add_dependencies( check test_${TESTDIR} )
        endif()
    endforeach()

    # 'make check-fast' builds and runs only the fast-labeled tests
    # (per-test wall time <~10 s, whole set <~3 min); includes the
    # interface/orientation regression battery when USE_MAXWELL is on
    add_custom_target(check-fast
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMAND ${CMAKE_CTEST_COMMAND} -L fast --output-on-failure
            COMMENT "Building and running fast tests"
            VERBATIM
            )

    foreach( TESTDIR core containers linalg math mesh fem gastables gasmodels
                     homology circuit )
        if( TARGET test_${TESTDIR} )
            add_dependencies( check-fast test_${TESTDIR} )
        endif()
    endforeach()

    # the backend-free material compile gate (tests/physics/backendfree) is an
    # OBJECT library: 'make' and 'make check' build it via ALL / test_physics,
    # but check-fast bypasses both, so hook it in explicitly
    if( TARGET test_material_backendfree )
        add_dependencies( check-fast test_material_backendfree )
    endif()

else()
    add_custom_target(check
            COMMAND ${CMAKE_COMMAND} -E echo "Error: tests are disabled. Re-run cmake with -DUSE_TEST=ON"
            COMMAND ${CMAKE_COMMAND} -E false
            )
    # override CMake's built-in 'test' target with a helpful error
    add_custom_target(test
            COMMAND ${CMAKE_COMMAND} -E echo "Error: use make check instead of make test."
            COMMAND ${CMAKE_COMMAND} -E false
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

    # The markdown page tree (anchors + doc/doxygen_nav.dox) is generated. Refresh it
    # before every run so a newly added document shows up in the navigation.
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        set(BELFEM_DOC_INDEX_CMD
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/update_doc_index.py)
    else()
        set(BELFEM_DOC_INDEX_CMD)
    endif()

    # The HTML header is generated from the INSTALLED Doxygen on every run and then
    # patched, rather than checked into the repository. Doxygen does not warn when a
    # custom header is stale, so a checked-in one would silently lose the navigation
    # tree, search or MathJax on a version upgrade. See config/doc/patch_doxygen_header.cmake.
    set(BELFEM_DOXYGEN_HEADER ${CMAKE_CURRENT_BINARY_DIR}/doxygen_header.html)
    set(BELFEM_PROJECT_URL "https://belfem.lbl.gov")

    # Doxyfile.in tracks a NEWER Doxygen than may be installed (it is maintained against
    # the version the site is developed on). An older Doxygen ignores every tag it does
    # not know -- 46 of them on 1.9.1 -- and silently drops the features they configure,
    # e.g. MATHJAX_VERSION, leaving the site on the MathJax 2 CDN with no error anywhere.
    #
    # So the config is normalised to the INSTALLED Doxygen with `doxygen -u`, which
    # rewrites a configuration to the running version's tag set while preserving every
    # value it still understands. Two separate files are required:
    #
    #   Doxyfile       configure_file output. PRISTINE -- never passed to `doxygen -u`.
    #   Doxyfile.work  the copy that is normalised and actually used.
    #
    # The split matters because `make doc` does not re-run CMake (a custom target is
    # always out of date, but that is not a reconfigure) and the Doxygen version is not
    # a configure dependency. Normalising the configure_file output in place would leave
    # a permanently downgraded file behind, and after a later Doxygen upgrade `-u` would
    # re-add the newer tags at their DEFAULTS rather than the values from Doxyfile.in.
    # Copying first makes the normalisation a pure function of (pristine, installed
    # Doxygen), so it is correct on every run and after any upgrade.
    set(BELFEM_DOXYFILE          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    set(BELFEM_DOXYFILE_WORK     ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.work)
    set(BELFEM_DOXYFILE_TEMPLATE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.template)

    # Create a doxygen make target
    add_custom_target(doc
            ${BELFEM_DOC_INDEX_CMD}

            # Normalise a COPY of the configured Doxyfile to the installed Doxygen.
            COMMAND ${CMAKE_COMMAND}
                    -DDOXYGEN=${DOXYGEN_EXECUTABLE}
                    -DPRISTINE=${BELFEM_DOXYFILE}
                    -DWORKING=${BELFEM_DOXYFILE_WORK}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/config/doc/normalize_doxyfile.cmake

            # Same treatment for the layout file, which is a full version-bound dump:
            # regenerate the installed Doxygen's default, re-apply BELFEM's custom nav
            # tabs and graph settings from doc/DoxygenLayout.xml, and repoint LAYOUT_FILE
            # at the result. The checked-in layout is only read, never written.
            COMMAND ${CMAKE_COMMAND}
                    -DDOXYGEN=${DOXYGEN_EXECUTABLE}
                    -DSOURCE_LAYOUT=${CMAKE_CURRENT_SOURCE_DIR}/doc/DoxygenLayout.xml
                    -DWORKING_LAYOUT=${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml
                    -DDOXYFILE=${BELFEM_DOXYFILE_WORK}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/config/doc/normalize_doxygen_layout.cmake

            # Generate the header template from a config copy whose HTML_HEADER is blank.
            # `doxygen -w html` checks that the config's HTML_HEADER EXISTS before writing
            # anything, and HTML_HEADER names the very file this command produces -- so
            # passing the real config here aborts the build on a tree that has no header
            # yet. Passing no config avoids the check but generates from Doxygen defaults,
            # and Doxygen documents the template as depending on the configuration (e.g.
            # GENERATE_TREEVIEW). Blanking the one self-referential tag gives both.
            COMMAND ${CMAKE_COMMAND}
                    -DIN=${BELFEM_DOXYFILE_WORK}
                    -DOUT=${BELFEM_DOXYFILE_TEMPLATE}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/config/doc/strip_doxygen_html_header.cmake
            COMMAND ${DOXYGEN_EXECUTABLE} -w html
                    ${BELFEM_DOXYGEN_HEADER}
                    ${CMAKE_CURRENT_BINARY_DIR}/doxygen_footer.html
                    ${CMAKE_CURRENT_BINARY_DIR}/doxygen_stylesheet.css
                    ${BELFEM_DOXYFILE_TEMPLATE}

            COMMAND ${CMAKE_COMMAND}
                    -DHEADER=${BELFEM_DOXYGEN_HEADER}
                    -DURL=${BELFEM_PROJECT_URL}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/config/doc/patch_doxygen_header.cmake
            COMMAND ${DOXYGEN_EXECUTABLE} ${BELFEM_DOXYFILE_WORK}
            COMMAND ${CMAKE_COMMAND}
                    -DLOG=${CMAKE_BINARY_DIR}/doc/doxygen_warnings.log
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/config/doc/report_doxygen_warnings.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating belfem documentation with Doxygen."
            VERBATIM )
endif()
