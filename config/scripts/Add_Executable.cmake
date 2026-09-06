include( ${BELFEM_CONFIG_DIR}/globals.cmake )

# add default includes
include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/containers )
include_directories( ${BELFEM_SOURCE_DIR}/linalg )
include_directories( ${BELFEM_SOURCE_DIR}/linalg/lapack )
if ( USE_MATRIX_ARMADILLO )
    include_directories( ${BELFEM_SOURCE_DIR}/linalg/armadillo )
elseif( USE_MATRIX_BLAZE )
    include_directories( ${BELFEM_SOURCE_DIR}/linalg/blaze )
endif()
include_directories( ${BELFEM_SOURCE_DIR}/linalg/operators )
include_directories( ${BELFEM_SOURCE_DIR}/io )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/sparse )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
# cl_Material.hpp exposes Bezier in its interface, so every target that sees a
# material needs this on the include path
include_directories( ${BELFEM_SOURCE_DIR}/numerics/bezier )

# generated belfem_version.hpp (see src/core/CMakeLists.txt)
include_directories( ${CMAKE_BINARY_DIR}/generated )

add_executable( ${EXECNAME} ${MAIN})

# generate the version header before the executable's own sources compile
if( TARGET belfem_version )
    add_dependencies( ${EXECNAME} belfem_version )
endif()

if( APPLE )
    if( COMMAND target_link_options )
        target_link_options( ${EXECNAME} PRIVATE "-Wl,-no_warn_duplicate_libraries" )
    else()
        set_target_properties( ${EXECNAME} PROPERTIES LINK_FLAGS "-Wl,-no_warn_duplicate_libraries" )
    endif()
endif()

# the whole framework is one library target ( Add_BelfemLibrary.cmake ), which
# carries the third-party, Fortran and OpenMP link interface with it. LIBLIST,
# set by the callers, is no longer consulted: a module's objects are in the
# library whether or not the executable names it.
target_link_libraries( ${EXECNAME} ${LIBTARGET} )

# ad-hoc signing keeps the debugger happy on macOS; the linker already signs
# on arm64, this makes x86_64 match
if( APPLE )
    if( NOT DEFINED BELFEM_CODESIGN_IDENTITY )
        set( BELFEM_CODESIGN_IDENTITY "-" CACHE STRING "macOS codesign identity for executables; '-' uses ad-hoc signing" )
    endif()

    add_custom_command(TARGET ${EXECNAME}
            POST_BUILD
            COMMAND codesign -f -s "${BELFEM_CODESIGN_IDENTITY}" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXECNAME}
            )
endif()

# install the executables named in BELFEM_INSTALL_EXECUTABLES ( globals.cmake ).
# The rule lives here, next to the target, because install( TARGETS ) on a
# target from another directory needs CMake 3.13 and the tree admits 3.11.
if( ${EXECNAME} IN_LIST BELFEM_INSTALL_EXECUTABLES )
    install( TARGETS ${EXECNAME}
             EXPORT ${LIBPREFIX}Targets
             RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )

    # installing rewrites the load commands ( build RPATH -> install RPATH ),
    # which invalidates the signature; on arm64 the kernel refuses to exec
    # such a binary, so it is signed again in place. DESTDIR-aware.
    if( APPLE )
        install( CODE "execute_process( COMMAND codesign -f -s \"${BELFEM_CODESIGN_IDENTITY}\" \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/${EXECNAME}\" )" )
    endif()
endif()
