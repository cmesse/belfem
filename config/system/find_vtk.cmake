set(VTK_FOUND FALSE)

# List of potential VTK directories to check
set(VTK_DIRS ${VTK_DIR} $ENV{VTK_DIR} $ENV{SCLS} "/usr" "/usr/local")

# Initialize variables to store the found VTK directory
set(VTK_FOUND FALSE)
set(VTK_VERSION "")

# Iterate over the directories
foreach(DIR IN LISTS VTK_DIRS)

    set(VTK_INCLUDE_DIR "${DIR}/include")
    if(NOT VTK_FOUND AND EXISTS "${VTK_INCLUDE_DIR}")
        file(GLOB VTK_VERSION_SUBDIRS "${VTK_INCLUDE_DIR}/vtk-*")

        if(VTK_VERSION_SUBDIRS)

            # Find the highest version
            foreach(SUBDIR IN LISTS VTK_VERSION_SUBDIRS)
                string(REGEX MATCH "vtk-([0-9]+\\.[0-9]+)" _ ${SUBDIR})
                set(CURRENT_VTK_VERSION ${CMAKE_MATCH_1})

                if (CURRENT_VTK_VERSION VERSION_GREATER VTK_VERSION OR VTK_VERSION STREQUAL "")
                    set(VTK_VERSION ${CURRENT_VTK_VERSION})
                    set(FINAL_VTK_DIR ${DIR})
                    set(FINAL_VTK_INCLUDES "${SUBDIR}" )
                endif()
            endforeach()

            # Mark as found to stop further searching
            set(VTK_FOUND TRUE)
        endif() # This endif corresponds to if(VTK_VERSION_SUBDIRS)
    endif() # This endif corresponds to if(NOT VTK_FOUND AND EXISTS "${VTK_INCLUDE_DIR}")
endforeach()

if(NOT VTK_FOUND)
    message(FATAL_ERROR "VTK not found in specified directories.")
else()
    # Set the VTK_DIR to the version found
    set(VTK_DIR "${FINAL_VTK_DIR}/lib/cmake/vtk-${VTK_VERSION}")

    # Now find the VTK package
    find_package(VTK REQUIRED HINTS ${VTK_DIR})
endif()

# Check if lib64 directory exists
if(EXISTS "${FINAL_VTK_DIR}/lib64")
    set(VTK_LIBDIR "${FINAL_VTK_DIR}/lib64")
else()
    set(VTK_LIBDIR "${FINAL_VTK_DIR}/lib")
endif()

# If VTK is found, configure it for use
if(VTK_FOUND)
    list(APPEND BELFEM_INCLUDES "${FINAL_VTK_INCLUDES}" )
    list(APPEND BELFEM_RPATH "${VTK_LIBDIR}")
    link_directories("${VTK_LIBDIR}")
else()
    message(FATAL_ERROR "VTK not found")
endif()

