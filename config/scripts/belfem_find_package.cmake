# register third-party link items at their dependency rank, following the
# SCLS build groups: a library must link BEFORE the libraries it depends on,
# so finalize_compiler.cmake assembles the link line from rank 10 down to 0.
# libraries of the same rank do not depend on each other; their relative
# order is free
macro( belfem_link_libraries RANK )
    list( APPEND BELFEM_TPL_RANK_${RANK} ${ARGN} )
endmacro()

function(belfem_find_package PACKAGE)
    set(options REQUIRED)
    set(oneValueArgs)
    set(multiValueArgs HEADERS LIBRARIES INCLUDE_SUFFIXES EXTRA_ROOTS)

    cmake_parse_arguments(BFP
            "${options}"
            "${oneValueArgs}"
            "${multiValueArgs}"
            ${ARGN}
    )

    if(NOT BFP_HEADERS)
        message(FATAL_ERROR "belfem_find_package(${PACKAGE}): no HEADERS given")
    endif()

    if(NOT BFP_LIBRARIES)
        message(FATAL_ERROR "belfem_find_package(${PACKAGE}): no LIBRARIES given")
    endif()

    set(ROOT_VAR "${PACKAGE}_DIR")

    set(ROOT_DIRS)

    # 1. -D<PACKAGE>_DIR=/path
    if(DEFINED ${ROOT_VAR} AND NOT "${${ROOT_VAR}}" STREQUAL "")
        list(APPEND ROOT_DIRS "${${ROOT_VAR}}")
    endif()

    # 2. Environment variable <PACKAGE>_DIR
    if(DEFINED ENV{${ROOT_VAR}} AND NOT "$ENV{${ROOT_VAR}}" STREQUAL "")
        list(APPEND ROOT_DIRS "$ENV{${ROOT_VAR}}")
    endif()

    # 3. SCLS root, if present
    if(DEFINED ENV{SCLS} AND NOT "$ENV{SCLS}" STREQUAL "")
        list(APPEND ROOT_DIRS "$ENV{SCLS}")
    endif()

    # 4. Extra roots passed by caller
    list(APPEND ROOT_DIRS ${BFP_EXTRA_ROOTS})

    # 5. System fallbacks
    list(APPEND ROOT_DIRS "/usr" "/usr/local")

    list(REMOVE_DUPLICATES ROOT_DIRS)

    # Include suffixes:
    #   ""      -> <root>/include
    #   nlopt   -> <root>/include/nlopt
    if(NOT BFP_INCLUDE_SUFFIXES)
        set(BFP_INCLUDE_SUFFIXES ".")
    endif()

    # Prefer lib64 over lib; on Debian-family systems also try the
    # multiarch directory, e.g. lib/x86_64-linux-gnu
    set(BFP_LIB_SUFFIXES lib64 lib)
    if(CMAKE_LIBRARY_ARCHITECTURE)
        list(APPEND BFP_LIB_SUFFIXES "lib/${CMAKE_LIBRARY_ARCHITECTURE}")
    endif()

    foreach(ROOT IN LISTS ROOT_DIRS)

        if(ROOT STREQUAL "")
            continue()
        endif()

        foreach(INC_SUFFIX IN LISTS BFP_INCLUDE_SUFFIXES)

            if(INC_SUFFIX STREQUAL ".")
                set(CANDIDATE_INCLUDE_DIR "${ROOT}/include")
            else()
                set(CANDIDATE_INCLUDE_DIR "${ROOT}/include/${INC_SUFFIX}")
            endif()

            if(NOT EXISTS "${CANDIDATE_INCLUDE_DIR}")
                continue()
            endif()

            set(HEADERS_FOUND TRUE)

            foreach(HEADER IN LISTS BFP_HEADERS)
                if(NOT EXISTS "${CANDIDATE_INCLUDE_DIR}/${HEADER}")
                    set(HEADERS_FOUND FALSE)
                    break()
                endif()
            endforeach()

            if(NOT HEADERS_FOUND)
                continue()
            endif()

            foreach(LIB_SUFFIX IN LISTS BFP_LIB_SUFFIXES)

                set(CANDIDATE_LIB_DIR "${ROOT}/${LIB_SUFFIX}")

                if(NOT EXISTS "${CANDIDATE_LIB_DIR}")
                    continue()
                endif()

                set(LIBRARIES_FOUND TRUE)
                set(FOUND_LIBRARIES)

                foreach(LIBRARY_NAME IN LISTS BFP_LIBRARIES)

                    unset(FOUND_LIBRARY CACHE)

                    find_library(
                            FOUND_LIBRARY
                            NAMES ${LIBRARY_NAME}
                            PATHS "${CANDIDATE_LIB_DIR}"
                            NO_DEFAULT_PATH
                    )

                    if(NOT FOUND_LIBRARY)
                        set(LIBRARIES_FOUND FALSE)
                        break()
                    endif()

                    list(APPEND FOUND_LIBRARIES "${FOUND_LIBRARY}")

                endforeach()

                if(LIBRARIES_FOUND)

                    set(${PACKAGE}_FOUND TRUE PARENT_SCOPE)
                    set(${PACKAGE}_ROOT_DIR "${ROOT}" PARENT_SCOPE)
                    set(${PACKAGE}_INCLUDE_DIR "${CANDIDATE_INCLUDE_DIR}" PARENT_SCOPE)
                    set(${PACKAGE}_LIB_DIR "${CANDIDATE_LIB_DIR}" PARENT_SCOPE)
                    set(${PACKAGE}_LIBRARIES "${FOUND_LIBRARIES}" PARENT_SCOPE)

                    message(STATUS "Found ${PACKAGE}:")
                    message(STATUS "  root:    ${ROOT}")
                    message(STATUS "  include: ${CANDIDATE_INCLUDE_DIR}")
                    message(STATUS "  lib:     ${CANDIDATE_LIB_DIR}")
                    message(STATUS "  libs:    ${FOUND_LIBRARIES}")

                    return()

                endif()

            endforeach()

        endforeach()

    endforeach()

    set(${PACKAGE}_FOUND FALSE PARENT_SCOPE)

    if(BFP_REQUIRED)
        message(FATAL_ERROR
                "${PACKAGE} not found. Tried roots: ${ROOT_DIRS}. "
                "Set -D${PACKAGE}_DIR=/path/to/${PACKAGE} or export ${PACKAGE}_DIR."
        )
    endif()

endfunction()
