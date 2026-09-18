if( USE_MPI )
    # test if MPI_HOME is set
    if( DEFINED ENV{MPI_HOME} )
        # for summary output
        set( BELFEM_MPIHOME $ENV{MPI_HOME} )

        # using an MPI that does not belong to SCLS is possible, but mixing
        # toolchains is asking for trouble
        if( DEFINED ENV{SCLS} AND NOT "$ENV{MPI_HOME}" STREQUAL "$ENV{SCLS}" )
            message( WARNING "MPI_HOME ($ENV{MPI_HOME}) points outside SCLS; using an MPI that was not built with the SCLS toolchain is not recommended." )
        endif()

        # find library directory; duplicate rpath entries are
        # removed in finalize_compiler.cmake
        if( IS_DIRECTORY "$ENV{MPI_HOME}/lib64" )
            list( APPEND BELFEM_RPATH "$ENV{MPI_HOME}/lib64" )
        elseif( IS_DIRECTORY "$ENV{MPI_HOME}/lib" )
            list( APPEND BELFEM_RPATH "$ENV{MPI_HOME}/lib" )
        else()
            message( FATAL_ERROR "Could not find MPI library directory in $ENV{MPI_HOME}" )
        endif()
    endif()

    # -------------------------------------------------------------------------
    # which MPI is this?
    #
    # BELFEM is built and validated against Open MPI only. MPICH - and Intel
    # MPI, which is MPICH-derived - is untested, and PETSc has been observed to
    # crash when called on that path. A build that merely links there is not
    # evidence that it works, so refuse the configure rather than let the
    # failure surface later inside a solver. See doc/mpi_support.md.
    #
    # The probe preprocesses mpi.h and reads the implementation macros the
    # header itself defines; nothing here depends on wrapper flags or on the
    # name of the compiler. CMAKE_CXX_COMPILER is already mpicxx at this point
    # ( detect_compiler.cmake runs before this file ), so mpi.h is on its
    # default include path.
    # -------------------------------------------------------------------------
    option( ALLOW_UNTESTED_MPI "permit MPI implementations other than Open MPI" OFF )

    set( _BELFEM_MPI_PROBE "${CMAKE_CURRENT_BINARY_DIR}/belfem_mpi_probe.cpp" )

    file( WRITE "${_BELFEM_MPI_PROBE}"
            "#include <mpi.h>\n"
            "#if defined( OPEN_MPI )\n"
            "BELFEM_MPI_FLAVOR OPENMPI\n"
            "#elif defined( MPICH ) || defined( MPICH_VERSION ) || defined( I_MPI_VERSION )\n"
            "BELFEM_MPI_FLAVOR MPICH\n"
            "#else\n"
            "BELFEM_MPI_FLAVOR UNKNOWN\n"
            "#endif\n" )

    # -E only: the marker lines are not valid C++, and are never compiled
    execute_process(
            COMMAND ${CMAKE_CXX_COMPILER} -E "${_BELFEM_MPI_PROBE}"
            OUTPUT_VARIABLE _BELFEM_MPI_PROBE_OUT
            ERROR_VARIABLE  _BELFEM_MPI_PROBE_ERR
            RESULT_VARIABLE _BELFEM_MPI_PROBE_RC )

    if( NOT _BELFEM_MPI_PROBE_RC EQUAL 0 )
        # could not read mpi.h. That is not evidence of an unsupported MPI, so
        # it must not stop the build - say so and carry on
        message( WARNING
                "Could not determine the MPI implementation: preprocessing mpi.h with "
                "${CMAKE_CXX_COMPILER} failed.\n"
                "BELFEM is validated against Open MPI only ( see doc/mpi_support.md ); "
                "this build is proceeding unchecked." )
    elseif( _BELFEM_MPI_PROBE_OUT MATCHES "BELFEM_MPI_FLAVOR[ \t]+OPENMPI" )
        set( BELFEM_MPI_FLAVOR "Open MPI" )
    elseif( _BELFEM_MPI_PROBE_OUT MATCHES "BELFEM_MPI_FLAVOR[ \t]+MPICH" )
        set( BELFEM_MPI_FLAVOR "MPICH family ( MPICH or Intel MPI )" )

        if( NOT ALLOW_UNTESTED_MPI )
            message( FATAL_ERROR
                    "This tree is configured against the MPICH family ( MPICH or Intel MPI ), "
                    "which BELFEM does not support.\n"
                    "  BELFEM is built and validated against Open MPI only. On MPICH, PETSc "
                    "has been observed to crash when called, so a build that links is not "
                    "evidence that it runs.\n"
                    "  Using the Intel compilers is fine - build Open MPI with them rather "
                    "than substituting Intel MPI.\n"
                    "  Point MPI_HOME / PATH at an Open MPI installation, or pass "
                    "-DALLOW_UNTESTED_MPI=ON to configure anyway.\n"
                    "  Background: doc/mpi_support.md" )
        endif()

        message( WARNING
                "Building against the MPICH family with ALLOW_UNTESTED_MPI=ON. This "
                "configuration is untested and PETSc is the known casualty; the MUMPS and "
                "MKL BLACS link flags also name Open MPI libraries explicitly." )
    else()
        set( BELFEM_MPI_FLAVOR "unrecognized" )

        message( WARNING
                "Could not recognize the MPI implementation from mpi.h: neither OPEN_MPI "
                "nor the MPICH macros are defined.\n"
                "BELFEM is validated against Open MPI only ( see doc/mpi_support.md ); "
                "this build is proceeding unchecked." )
    endif()

    # -------------------------------------------------------------------------
    # can the wrapper link an MPI program?
    #
    # CMake's own compiler test ran at project() with the plain C++ compiler; the
    # swap to mpicxx happens afterwards, so nothing so far has linked against
    # libmpi. A wrapper that cannot link these two calls cannot link belfem, and
    # without this probe the failure surfaces at the first executable, after the
    # whole library has compiled, as a wall of undefined references out of
    # libmpi.so. The probe carries the rpath entries collected so far ( the
    # toolchain prefix; the third-party directories come later ), so ld resolves
    # libmpi's own dependencies from the prefix first, as the real link does.
    # MPI_Init must be called: under --as-needed an empty main drops -lmpi and
    # the probe proves nothing. The binary is never executed.
    # -------------------------------------------------------------------------
    set( _BELFEM_MPI_LINK_PROBE_SRC "${CMAKE_BINARY_DIR}/CMakeFiles/belfem_mpi_link_probe.cpp" )
    set( _BELFEM_MPI_LINK_PROBE_BIN "${CMAKE_BINARY_DIR}/CMakeFiles/belfem_mpi_link_probe" )

    file( WRITE "${_BELFEM_MPI_LINK_PROBE_SRC}"
            "#include <mpi.h>\n"
            "int main( int argc, char** argv )\n"
            "{\n"
            "    MPI_Init( &argc, &argv );\n"
            "    MPI_Finalize();\n"
            "    return 0;\n"
            "}\n" )

    # same rpath the targets get: finalize_compiler.cmake prunes the list the
    # same way, but runs after this file. One -rpath per directory: GNU ld also
    # takes a colon-separated list, Apple's ld64 does not
    set( _BELFEM_MPI_LINK_PROBE_RPATH ${BELFEM_RPATH} )
    belfem_prune_system_libdirs( _BELFEM_MPI_LINK_PROBE_RPATH )
    set( _BELFEM_MPI_LINK_PROBE_FLAGS )
    foreach( _BELFEM_MPI_LINK_PROBE_DIR ${_BELFEM_MPI_LINK_PROBE_RPATH} )
        list( APPEND _BELFEM_MPI_LINK_PROBE_FLAGS "-Wl,-rpath,${_BELFEM_MPI_LINK_PROBE_DIR}" )
    endforeach()

    execute_process(
            COMMAND ${CMAKE_CXX_COMPILER} ${_BELFEM_MPI_LINK_PROBE_FLAGS}
                    "${_BELFEM_MPI_LINK_PROBE_SRC}" -o "${_BELFEM_MPI_LINK_PROBE_BIN}"
            OUTPUT_VARIABLE _BELFEM_MPI_LINK_PROBE_OUT
            ERROR_VARIABLE  _BELFEM_MPI_LINK_PROBE_ERR
            RESULT_VARIABLE _BELFEM_MPI_LINK_PROBE_RC )

    if( NOT _BELFEM_MPI_LINK_PROBE_RC EQUAL 0 )
        set( _BELFEM_MPI_LINK_PROBE_TEXT "${_BELFEM_MPI_LINK_PROBE_OUT}${_BELFEM_MPI_LINK_PROBE_ERR}" )
        set( _BELFEM_MPI_LINK_PROBE_HINT "" )

        # Open MPI 5 needs the PMIx it was built against; an older or absent
        # libpmix.so.2 leaves its PMIx_* imports unresolved. Only Open MPI has
        # this dependency, so the hint is not offered for another flavor
        if( BELFEM_MPI_FLAVOR STREQUAL "Open MPI"
            AND ( _BELFEM_MPI_LINK_PROBE_TEXT MATCHES "PMIx_"
                  OR _BELFEM_MPI_LINK_PROBE_TEXT MATCHES "libpmix\\.so\\.2[^\n]*not found"
                  OR _BELFEM_MPI_LINK_PROBE_TEXT MATCHES "cannot find libpmix\\.so\\.2" ) )

            if( BELFEM_SCLS_FLAVOR )
                set( _BELFEM_MPI_PMIX_SENTENCE "  Under SCLS, install scls-${BELFEM_SCLS_FLAVOR}-pmix.\n" )
            else()
                set( _BELFEM_MPI_PMIX_SENTENCE "  Install the PMIx this Open MPI was built against.\n" )
            endif()

            # the wrapper knows where libmpi lives; the hint names the first
            # directory that holds it, and no hint is better than a wrong one
            set( _BELFEM_MPI_LDD_HINT "" )
            execute_process(
                    COMMAND ${CMAKE_CXX_COMPILER} --showme:libdirs
                    OUTPUT_VARIABLE _BELFEM_MPI_LIBDIRS
                    ERROR_QUIET
                    RESULT_VARIABLE _BELFEM_MPI_LIBDIRS_RC )
            if( _BELFEM_MPI_LIBDIRS_RC EQUAL 0 )
                string( STRIP "${_BELFEM_MPI_LIBDIRS}" _BELFEM_MPI_LIBDIRS )
                separate_arguments( _BELFEM_MPI_LIBDIRS )
                foreach( _BELFEM_MPI_LIBDIR ${_BELFEM_MPI_LIBDIRS} )
                    if( NOT _BELFEM_MPI_LDD_HINT AND EXISTS "${_BELFEM_MPI_LIBDIR}/libmpi.so" )
                        set( _BELFEM_MPI_LDD_HINT "  Check: ldd ${_BELFEM_MPI_LIBDIR}/libmpi.so | grep pmix\n" )
                    endif()
                endforeach()
            endif()

            string( CONCAT _BELFEM_MPI_LINK_PROBE_HINT
                    "  Open MPI's PMIx dependency is missing or incompatible with this Open MPI.\n"
                    "${_BELFEM_MPI_PMIX_SENTENCE}"
                    "${_BELFEM_MPI_LDD_HINT}" )
        endif()

        if( NOT BELFEM_MPI_FLAVOR )
            set( BELFEM_MPI_FLAVOR "unknown MPI" )
        endif()
        message( FATAL_ERROR
                "Could not link an MPI program with ${CMAKE_CXX_COMPILER} ( ${BELFEM_MPI_FLAVOR} ):\n"
                "${_BELFEM_MPI_LINK_PROBE_TEXT}\n"
                "${_BELFEM_MPI_LINK_PROBE_HINT}"
                "  Background: doc/mpi_support.md" )
    endif()

    # -------------------------------------------------------------------------
    # which PMIx does the loader pick?
    #
    # The link above resolved libpmix.so.2 through -rpath, where the toolchain
    # prefix comes first. The loader reads LD_LIBRARY_PATH before RUNPATH, so an
    # environment that carries another PMIx runs the program against that one.
    # Only meaningful when the prefix ships a PMIx at all; ldd is Linux-only, and
    # a failure here is never an error - the link already succeeded.
    # -------------------------------------------------------------------------
    if( CMAKE_SYSTEM_NAME STREQUAL "Linux" AND DEFINED ENV{SCLS}
        AND SCLSLIBDIR AND EXISTS "${SCLSLIBDIR}/libpmix.so.2" )

        execute_process(
                COMMAND ldd "${_BELFEM_MPI_LINK_PROBE_BIN}"
                OUTPUT_VARIABLE _BELFEM_MPI_LDD_OUT
                ERROR_QUIET
                RESULT_VARIABLE _BELFEM_MPI_LDD_RC )

        if( _BELFEM_MPI_LDD_RC EQUAL 0 )
            if( _BELFEM_MPI_LDD_OUT MATCHES "libpmix\\.so\\.2 => not found" )
                message( WARNING
                        "The loader cannot resolve libpmix.so.2 for an MPI program in this "
                        "environment although ${SCLSLIBDIR}/libpmix.so.2 exists. The link "
                        "succeeded; running may not." )
            elseif( _BELFEM_MPI_LDD_OUT MATCHES "libpmix\\.so\\.2 => ([^ \n]+)" )
                set( BELFEM_PMIX_LIBRARY "${CMAKE_MATCH_1}" )
                get_filename_component( _BELFEM_PMIX_REAL "${BELFEM_PMIX_LIBRARY}" REALPATH )
                get_filename_component( _BELFEM_SCLS_REAL "$ENV{SCLS}" REALPATH )
                string( FIND "${_BELFEM_PMIX_REAL}" "${_BELFEM_SCLS_REAL}/" _BELFEM_PMIX_POS )
                if( NOT _BELFEM_PMIX_POS EQUAL 0 )
                    message( WARNING
                            "ldd resolves libpmix.so.2 to ${BELFEM_PMIX_LIBRARY}, not to the copy in "
                            "${_BELFEM_SCLS_REAL}. The loader reads LD_LIBRARY_PATH before RUNPATH, "
                            "so a run in this environment may load that PMIx." )
                endif()
            endif()
        endif()
    endif()
endif()
