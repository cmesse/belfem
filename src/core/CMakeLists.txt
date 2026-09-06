set( LIBNAME core )

set( SOURCES
        assert.cpp
        stringtools.cpp
        banner.cpp
        cl_Logger.cpp
        cl_Profiler.cpp
        cl_Arguments.cpp
        cl_Progressbar.cpp
        fn_available_memory.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/containers )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

# -----------------------------------------------------------------------------
# version + git provenance header (consumed by banner.cpp)
#
# holds the semantic version (from project(... VERSION ...) in the top-level
# CMakeLists.txt) plus the git commit hash / branch / dirty flag. regenerated on
# EVERY build via an ALL custom target, so the hash and dirty flag stay correct
# after a plain `make`/`ninja` without re-running cmake. builds outside a git
# repo (e.g. a source tarball) fall back to "unknown" for the git fields.
# -----------------------------------------------------------------------------
set( BELFEM_VERSION_HEADER ${CMAKE_BINARY_DIR}/generated/belfem_version.hpp )

# run the generator as the custom TARGET's own COMMAND ( not via an
# add_custom_command OUTPUT, which is treated as up-to-date once the file
# exists and would never see a new commit / branch switch ): a custom target's
# COMMAND runs on EVERY build. The generator's configure_file only rewrites the
# header when the git fields actually change, so banner.cpp recompiles only then.
add_custom_target( belfem_version ALL
        COMMAND ${CMAKE_COMMAND}
            -DGIT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DOUTPUT_FILE=${BELFEM_VERSION_HEADER}
            -DVERSION_MAJOR=${PROJECT_VERSION_MAJOR}
            -DVERSION_MINOR=${PROJECT_VERSION_MINOR}
            -DVERSION_PATCH=${PROJECT_VERSION_PATCH}
            -P ${BELFEM_CONFIG_DIR}/scripts/GenerateVersion.cmake
        BYPRODUCTS ${BELFEM_VERSION_HEADER}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Updating BELFEM version header"
        VERBATIM )

add_dependencies( ${BELFEM_MODULE_TARGET} belfem_version )
target_include_directories( ${BELFEM_MODULE_TARGET}
        PRIVATE ${CMAKE_BINARY_DIR}/generated )