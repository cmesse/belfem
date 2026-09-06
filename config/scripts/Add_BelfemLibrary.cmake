# -----------------------------------------------------------------------------
# The project library: every module's objects in one target.
#
# Included from the top-level CMakeLists.txt once all module subdirectories
# ( src/, and nonfree/ when present ) have registered themselves in the
# BELFEM_MODULE_TARGETS global property ( Add_Library.cmake ). Executables and
# tests link the ${LIBTARGET} target and nothing else from the project. The
# target is named belfem_lib so the plain name `belfem` stays free for the
# solver executable; the archive on disk is still lib${LIBPREFIX}.a .
# -----------------------------------------------------------------------------
include( ${BELFEM_CONFIG_DIR}/globals.cmake )

get_property( BELFEM_MODULE_TARGETS GLOBAL PROPERTY BELFEM_MODULE_TARGETS )

set( BELFEM_MODULE_OBJECTS )
foreach( _BELFEM_MODULE ${BELFEM_MODULE_TARGETS} )
    list( APPEND BELFEM_MODULE_OBJECTS $<TARGET_OBJECTS:${_BELFEM_MODULE}> )
endforeach()

if( USE_SHARED_LIBS )
    add_library( ${LIBTARGET} SHARED ${BELFEM_MODULE_OBJECTS} )
else()
    add_library( ${LIBTARGET} STATIC ${BELFEM_MODULE_OBJECTS} )
endif()

# libbelfem.so.MAJOR.MINOR.PATCH -> libbelfem.so.MAJOR -> libbelfem.so
set_target_properties( ${LIBTARGET} PROPERTIES
        OUTPUT_NAME ${LIBPREFIX}
        VERSION     ${PROJECT_VERSION}
        SOVERSION   ${PROJECT_VERSION_MAJOR} )

# the link interface: everything an executable needed to name before, it now
# inherits from the library. third-party libraries are pre-sorted by descending
# dependency rank, see finalize_compiler.cmake
target_link_libraries( ${LIBTARGET} PUBLIC
        ${BELFEM_TPL_LIBS}
        ${BELFEM_FORTRANLIBS}
        ${BELFEM_OPENMPLIBS} )

if( USE_VTK )
    target_link_libraries( ${LIBTARGET} PUBLIC ${VTK_LIBRARIES} )
endif()

if( APPLE )
    if( COMMAND target_link_options )
        target_link_options( ${LIBTARGET} PRIVATE "-Wl,-no_warn_duplicate_libraries" )
    endif()

    # gComm and gLog are declared extern in the library and DEFINED in each
    # executable's main ( ~50 sites in src/, tests/ and nonfree/, so a main can
    # pick its own verbosity ). ELF tolerates the dangling reference in a .so;
    # Mach-O's two-level namespace does not, so these two symbols are resolved
    # at load time against the executable. The cost is that ANY undefined
    # symbol in the dylib is now reported at load rather than at link — see
    # O7 in todo/shared_library_and_install_plan.md for the alternative.
    if( USE_SHARED_LIBS AND COMMAND target_link_options )
        target_link_options( ${LIBTARGET} PRIVATE "-Wl,-undefined,dynamic_lookup" )
    endif()
endif()

install( TARGETS ${LIBTARGET}
         EXPORT ${LIBPREFIX}Targets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )

# see Add_Executable.cmake: the install-time RPATH rewrite voids the signature
if( APPLE AND USE_SHARED_LIBS )
    install( CODE "execute_process( COMMAND codesign -f -s - \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/lib${LIBPREFIX}.${PROJECT_VERSION}.dylib\" )" )
endif()
