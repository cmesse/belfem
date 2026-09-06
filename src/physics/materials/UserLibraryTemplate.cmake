# =============================================================================
# BELFEM User Library Template
# =============================================================================
# This template helps you compile custom functions as shared libraries
# that can be dynamically loaded by BELFEM.
#
# USAGE:
# 1. Copy this file to your lib directory
# 2. Rename it to CMakeLists.txt
# 3. Edit the USER CONFIGURATION section below
# 4. Run: mkdir build && cd build && cmake .. && make
#
# Your library is compiled as a loadable module (libcustom.so) that BELFEM
# opens with dlopen() at run time.
#
# WHERE THE DECK'S `file :` IS LOOKED FOR. Since 2026-08-31 a source plugin is
# resolved exactly like a material plugin, in this order:
#
#   1. the path as written, relative to the run directory, or absolute
#   2. the same relative path below <data root>/material
#   3. for a path CONTAINING a slash, the file name alone below
#      <data root>/material  (a bare name was already tried at step 2)
#   4. otherwise the name is handed to dlopen unchanged
#
# The data root is BELFEM's gBelfemDataPath, NOT $BELFEM_DATA directly: an
# installed tree fills it from BELFEM_INSTALL_DATADIR even when the environment
# is silent, and an empty root skips steps 2 and 3 entirely.
#
# So the copy in your build directory works as-is (step 1), and a plugin
# installed into the shared data directory can be named without a path.
#
# One consequence, and it runs the OPPOSITE way to the obvious guess: dlopen
# consults the platform's dynamic-loader search path ($LD_LIBRARY_PATH on
# Linux, $DYLD_LIBRARY_PATH on macOS) only for a name carrying NO slash. Giving
# the path a directory component therefore DISABLES that search rather than
# enabling it. A bare name is the only spelling the loader path can serve, and
# even then the data root is now tried first.
#
# PLATFORMS: Linux and macOS. Windows is NOT supported -- the loader in
# cl_SourceFunction.cpp is POSIX dlopen/dlsym with no LoadLibrary path.
#
# THE PLUGIN API IS BACKEND-FREE. Everything <belfem_user_api> exposes --
# cl_Material.hpp, cl_SourceFunction.hpp, cl_JcFunction.hpp and the rest --
# is free of cl_Vector.hpp and cl_Matrix.hpp by contract, so an ordinary
# plugin compiles with neither BELFEM_ARMADILLO nor BELFEM_BLAZE defined and
# does not have to match the backend the host was built with. (For source
# functions this became true on 2026-08-29; cl_SourceFunction.hpp now carries
# an explicit "NEVER include cl_Vector.hpp here" rule. The contract is pinned
# by tests/physics/backendfree/, which compiles the shipped example with the
# backend macro stripped.)
#
# BELFEM_BACKEND below is therefore OPTIONAL, and OFF by default. Set it only
# if your own sources reach past the plugin API into src/linalg -- which is
# outside the supported surface, and means your library must then be rebuilt
# whenever the host's backend changes.
# =============================================================================

# 3.13 for target_link_options(), used in the APPLE branch below.
cmake_minimum_required(VERSION 3.13)

# =============================================================================
# USER CONFIGURATION - Edit these variables
# =============================================================================

# CMake project name. This is NOT the library name -- LIBRARY_NAME below is.
project(MyCustomLibrary)

# Path to the BELFEM installation or source directory.
#   Option 1: an installed BELFEM  -> the directory containing include/belfem
#   Option 2: a BELFEM source tree -> the directory containing src/
#
# Either edit the default below or pass -DBELFEM_DIR=... ; the command line
# wins and both survive a reconfigure. NOT a CACHE variable on purpose: a
# cached default is written on the first configure and then shadows every later
# edit to this file.
set(BELFEM_DIR_DEFAULT "/path/to/belfem")

if(NOT BELFEM_DIR)
    set(BELFEM_DIR "${BELFEM_DIR_DEFAULT}")
endif()

# OPT-IN, and normally leave it empty. Only a library that includes
# cl_Vector.hpp / cl_Matrix.hpp itself needs this: those headers select their
# implementation on this macro, and with neither defined they collapse to
# forward declarations. There is no way to detect the host's choice from the
# headers, so if you do need it, check the host's CMakeCache.txt for
# USE_MATRIX_ARMADILLO / USE_MATRIX_BLAZE -- and note that you are then
# outside the backend-free plugin contract.
set(BELFEM_BACKEND "" CACHE STRING "Host backend, only if you use linalg: ARMADILLO, BLAZE, or empty")
set_property(CACHE BELFEM_BACKEND PROPERTY STRINGS "" ARMADILLO BLAZE)

# Only needed alongside BELFEM_BACKEND above; leave empty otherwise.
# Where the Armadillo or Blaze headers themselves live. BELFEM's wrappers
# include <armadillo> / <blaze/Blaze.h> directly, so the backend's own include
# directory has to be on the path as well. Set it to whatever -I your host
# build used. Blaze additionally needs MKL's headers, because its
# config/BLAS.h includes mkl_cblas.h. Semicolon-separated list, for example:
#
#   -DBELFEM_TPL_INCLUDE_DIRS="/opt/scls/mkl/include;/opt/intel/oneapi/mkl/latest/include"
set(BELFEM_TPL_INCLUDE_DIRS "" CACHE STRING
        "Include directories for the backend headers (Armadillo or Blaze)")

# List your library source files here
set(LIB_SOURCES
        defect.cpp
        current.cpp
        # Add more .cpp files if needed
)

# Library name, without the "lib" prefix and without the extension.
# The output file is lib${LIBRARY_NAME}.so on both Linux and macOS.
set(LIBRARY_NAME "custom")

# =============================================================================
# BELFEM Configuration - Automatic (usually no need to edit)
# =============================================================================

if(BELFEM_DIR STREQUAL "/path/to/belfem")
    message(FATAL_ERROR
            "BELFEM_DIR still holds the placeholder path. Either edit "
            "BELFEM_DIR_DEFAULT near the top of this file, or re-run with:\n"
            "  cmake -DBELFEM_DIR=/path/to/belfem ..")
endif()

if(BELFEM_BACKEND AND NOT BELFEM_BACKEND MATCHES "^(ARMADILLO|BLAZE)$")
    message(FATAL_ERROR
            "BELFEM_BACKEND must be ARMADILLO, BLAZE, or empty, got '${BELFEM_BACKEND}'.")
endif()

# Probe for a header, not for a directory name: "${BELFEM_DIR}/include" exists
# after any `make install`, so a directory probe accepts the tree and the build
# then dies on a missing header instead of here. Headers install as
# <prefix>/include/belfem/<module>/..., so BOTH branches need the module list.
# fem/kernel carries belfem_user_api.hpp (the umbrella your sources include);
# the rest is its closure plus headroom.
set(BELFEM_MODULE_DIRS
        core containers math/graph io
        numerics numerics/sources physics/materials fem/kernel)

# Backend headers, only when BELFEM_BACKEND was set. The wrappers
# (cl_AR_*.hpp / cl_BZ_*.hpp) sit in their own subdirectory and cl_Vector.hpp
# includes them by bare name.
if(BELFEM_BACKEND)
    list(APPEND BELFEM_MODULE_DIRS linalg linalg/lapack linalg/operators sparse)
    if(BELFEM_BACKEND STREQUAL "ARMADILLO")
        list(APPEND BELFEM_MODULE_DIRS linalg/armadillo)
    else()
        list(APPEND BELFEM_MODULE_DIRS linalg/blaze)
    endif()
endif()

# Probe for a header, not a directory name.
if(EXISTS "${BELFEM_DIR}/src/numerics/sources/cl_SourceFunction.hpp")
    set(BELFEM_HEADER_ROOT "${BELFEM_DIR}/src")
    message(STATUS "Using BELFEM source tree: ${BELFEM_DIR}")
elseif(EXISTS "${BELFEM_DIR}/include/belfem/numerics/sources/cl_SourceFunction.hpp")
    set(BELFEM_HEADER_ROOT "${BELFEM_DIR}/include/belfem")
    message(STATUS "Using installed BELFEM: ${BELFEM_DIR}")
else()
    message(FATAL_ERROR
            "Cannot find cl_SourceFunction.hpp under BELFEM_DIR. Expected either\n"
            "  ${BELFEM_DIR}/src/numerics/sources/cl_SourceFunction.hpp   (source tree), or\n"
            "  ${BELFEM_DIR}/include/belfem/numerics/sources/cl_SourceFunction.hpp   (installed)\n"
            "Current BELFEM_DIR: ${BELFEM_DIR}")
endif()

set(BELFEM_INCLUDE_DIRS "${BELFEM_HEADER_ROOT}")
foreach(module IN LISTS BELFEM_MODULE_DIRS)
    list(APPEND BELFEM_INCLUDE_DIRS "${BELFEM_HEADER_ROOT}/${module}")
endforeach()

# With a backend opted in, the sources will include cl_Vector.hpp. Probe for it
# here rather than letting a partial installation configure cleanly and fail at
# the first compile.
if(BELFEM_BACKEND AND NOT EXISTS "${BELFEM_HEADER_ROOT}/linalg/cl_Vector.hpp")
    message(FATAL_ERROR
            "BELFEM_BACKEND=${BELFEM_BACKEND} was requested, but "
            "${BELFEM_HEADER_ROOT}/linalg/cl_Vector.hpp does not exist.")
endif()

if(BELFEM_TPL_INCLUDE_DIRS)
    list(APPEND BELFEM_INCLUDE_DIRS ${BELFEM_TPL_INCLUDE_DIRS})
endif()

# -----------------------------------------------------------------------------
# The <belfem_user_api> umbrella
# -----------------------------------------------------------------------------
# Your sources include the extension-less umbrella, which is what an INSTALLED
# BELFEM is meant to expose. Neither a source tree nor today's install ships
# that spelling -- both carry belfem_user_api.hpp under fem/kernel -- so
# generate a one-line forwarder and put it on the include path. Drop this block
# once the install ships the real header.
if(EXISTS "${BELFEM_HEADER_ROOT}/fem/kernel/belfem_user_api.hpp")
    set(BELFEM_USER_API_SHIM "${CMAKE_CURRENT_BINARY_DIR}/api_shim")
    file(MAKE_DIRECTORY "${BELFEM_USER_API_SHIM}")
    file(WRITE "${BELFEM_USER_API_SHIM}/belfem_user_api"
         "// generated by CMakeLists.txt -- forwards to the in-tree umbrella\n"
         "#include \"belfem_user_api.hpp\"\n")
    list(APPEND BELFEM_INCLUDE_DIRS "${BELFEM_USER_API_SHIM}")
    message(STATUS "Generated <belfem_user_api> forwarder")
elseif(NOT EXISTS "${BELFEM_HEADER_ROOT}/belfem_user_api")
    message(FATAL_ERROR
            "No <belfem_user_api> umbrella found under ${BELFEM_HEADER_ROOT}. "
            "The shipped example sources include it, so the compile cannot succeed.")
endif()

# -----------------------------------------------------------------------------
# Where your library reads its data tables from, if it reads any
# -----------------------------------------------------------------------------
# Compiled in, so the library works regardless of the directory the solver is
# launched from. Never locate data by deriving a path from __FILE__: that bakes
# the source directory into the binary and breaks as soon as anything moves.
set(USER_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data" CACHE PATH
        "Directory holding this library's data tables")

# =============================================================================
# Compiler Settings
# =============================================================================

# Match the host's build type. NDEBUG/DEBUG decide BELFEM_ASSERTIONS_ACTIVE
# (src/core/assert.hpp), so they change the body of inline BELFEM methods
# compiled into both your library and the host. With no build type set at all
# CMake optimizes nothing.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type for this library" FORCE)
endif()

# C++17 standard (required by BELFEM). CMAKE_CXX_EXTENSIONS is deliberately
# left alone: BELFEM does not set it either, so host and plugin both land on
# gnu++17.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Position Independent Code (required for shared libraries)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Compiler flags.
#
# -Wno-unused-parameter is deliberate: the callback signatures are FIXED by
# BELFEM ( a source function always takes the full argument list ), so a
# library that ignores one of them is writing correct code, not sloppy code.
# Without this, the common case warns on every build. It also silences -Wextra
# hits from BELFEM's own installed headers, which the author cannot fix.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wno-unused-parameter)
endif()

# =============================================================================
# Build Library
# =============================================================================

# MODULE, not SHARED: this artifact is only ever dlopen()ed, never linked
# against. MODULE also yields lib<name>.so on macOS as well as Linux, matching
# the name BELFEM's documentation uses.
add_library(${LIBRARY_NAME} MODULE ${LIB_SOURCES})

# Include BELFEM headers. Quoted: BELFEM_DIR may contain spaces.
target_include_directories(${LIBRARY_NAME} PRIVATE "${BELFEM_INCLUDE_DIRS}")

# Trailing slash is deliberate: the sources concatenate a bare file name onto it.
target_compile_definitions(${LIBRARY_NAME} PRIVATE
        BELFEM_USER_DATA_DIR="${USER_DATA_DIR}/")

# The backend macro cl_Vector.hpp / cl_Matrix.hpp switch on, defined ONLY if you
# opted in above. Without it those headers degrade to forward declarations and
# any use of Vector/Matrix fails to compile. When set, it MUST match the host.
if(BELFEM_BACKEND)
    target_compile_definitions(${LIBRARY_NAME} PRIVATE BELFEM_${BELFEM_BACKEND})
endif()

# Nothing is linked here on purpose: the host resolves the library's undefined
# BELFEM symbols when it dlopen()s the file. On Linux that works because
# BELFEM's executables are built with -rdynamic; an executable of your own
# embedding BELFEM needs -rdynamic too. Mach-O rejects undefined symbols at
# link time, so it needs the same policy libbelfem uses for gComm/gLog.
if(APPLE)
    target_link_options(${LIBRARY_NAME} PRIVATE "-Wl,-undefined,dynamic_lookup")
endif()

# =============================================================================
# Installation (optional)
# =============================================================================

# Installs under CMAKE_INSTALL_PREFIX, which CMake defaults to /usr/local --
# NOT the build directory. You do not have to install at all: BELFEM loads the
# library by path, so the copy in your build directory works as-is.
include(GNUInstallDirs)
install(TARGETS ${LIBRARY_NAME}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# =============================================================================
# Information
# =============================================================================

message(STATUS "==================================================")
message(STATUS "BELFEM User Library Configuration")
message(STATUS "==================================================")
message(STATUS "Library name:     ${LIBRARY_NAME}")
message(STATUS "Source files:      ${LIB_SOURCES}")
message(STATUS "BELFEM directory:  ${BELFEM_DIR}")
message(STATUS "Header root:       ${BELFEM_HEADER_ROOT}")
if(BELFEM_BACKEND)
    message(STATUS "Backend:           BELFEM_${BELFEM_BACKEND}")
else()
    message(STATUS "Backend:           none (backend-free plugin API)")
endif()
message(STATUS "Data directory:    ${USER_DATA_DIR}")
message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "Output library:    lib${LIBRARY_NAME}.so")
message(STATUS "==================================================")
message(STATUS "To build: cmake .. && make")
message(STATUS "To install: make install")
message(STATUS "==================================================")
