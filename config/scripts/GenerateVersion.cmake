# Generates the auto-generated belfem_version.hpp header holding both the
# semantic version (assigned in the top-level CMakeLists.txt) and the git commit
# provenance of the build. Runs on every build (not just at configure time), so
# the hash and dirty flag stay correct after a plain `make`/`ninja`.
#
# Inputs (passed via -D on the command line):
#   GIT_SOURCE_DIR : repo root (where .git lives)
#   OUTPUT_FILE    : full path of the header to write
#   VERSION_MAJOR  : PROJECT_VERSION_MAJOR
#   VERSION_MINOR  : PROJECT_VERSION_MINOR
#   VERSION_PATCH  : PROJECT_VERSION_PATCH

# Defaults for a non-git build (e.g. released source tarball)
set(BELFEM_GIT_HASH        "unknown")
set(BELFEM_GIT_HASH_SHORT  "unknown")
set(BELFEM_GIT_BRANCH      "unknown")
set(BELFEM_GIT_DIRTY       0)

find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${GIT_SOURCE_DIR}/.git")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${GIT_SOURCE_DIR}
        OUTPUT_VARIABLE BELFEM_GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${GIT_SOURCE_DIR}
        OUTPUT_VARIABLE BELFEM_GIT_HASH_SHORT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${GIT_SOURCE_DIR}
        OUTPUT_VARIABLE BELFEM_GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    # `--quiet HEAD` so staged-but-uncommitted changes also count as dirty
    execute_process(
        COMMAND ${GIT_EXECUTABLE} diff --quiet HEAD
        WORKING_DIRECTORY ${GIT_SOURCE_DIR}
        RESULT_VARIABLE BELFEM_GIT_DIRTY_RESULT
    )
    if(NOT BELFEM_GIT_DIRTY_RESULT EQUAL 0)
        set(BELFEM_GIT_DIRTY 1)
    endif()

    # if rev-parse produced nothing (e.g. empty repo with no commits), fall
    # back to the "unknown" defaults rather than emitting empty string literals
    if(BELFEM_GIT_HASH STREQUAL "")
        set(BELFEM_GIT_HASH "unknown")
        set(BELFEM_GIT_HASH_SHORT "unknown")
        set(BELFEM_GIT_BRANCH "unknown")
    endif()
endif()

set(BELFEM_VERSION_STRING "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# single decision point for "was this built from a git checkout": true unless the
# hash fell back to the "unknown" sentinel (no git, or an empty repo)
if(BELFEM_GIT_HASH STREQUAL "unknown")
    set(BELFEM_GIT_AVAILABLE 0)
else()
    set(BELFEM_GIT_AVAILABLE 1)
endif()

# fill in the header layout from the template. @ONLY so only @VAR@ placeholders
# are substituted (the header has no ${...} of its own to protect). the template
# sits next to this script, so CMAKE_CURRENT_LIST_DIR resolves it in `-P` mode.
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/belfem_version.hpp.in
    ${OUTPUT_FILE}
    @ONLY )
