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
endif()