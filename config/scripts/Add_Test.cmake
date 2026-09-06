include( ${BELFEM_CONFIG_DIR}/globals.cmake )

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
include_directories( ${BELFEM_SOURCE_DIR}/numerics/bezier )
include_directories( ${CMAKE_BINARY_DIR}/generated )
# shared test-only headers ( the Tier 2 launcher sentinel )
include_directories( ${CMAKE_SOURCE_DIR}/tests/common )
# The list above mirrors Add_Executable.cmake so a test tree is self-sufficient
# wherever it is added. include_directories() is directory-scoped, and the only
# root-scope call sites are the Add_Executable.cmake include( s ) in the
# top-level CMakeLists.txt, which run after add_subdirectory( nonfree ) and
# before add_subdirectory( tests ). The open-source tests/ tree therefore
# inherits those paths and nonfree/tests does not — and if the last root-level
# executable is ever retired, neither tree does. This list, not inheritance, is
# what a test tree is entitled to rely on.

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test )

add_executable( test_${TESTNAME} test_${TESTNAME}_main.cpp ${SOURCES} )

# the data lookup below is pointed at the source tree's share/. Tests run with
# CWD = their build directory, from where the CWD-relative fallbacks in the
# library cannot reach it, and the shell's BELFEM_DATA must not decide whether a
# data-backed test runs or skips. The property overrides an inherited value on
# purpose: the suite is hermetic. Kept out of any compiled define, so the
# configuring machine's source path never ends up inside a distributed binary.
set( BELFEM_TEST_ENV "BELFEM_DATA=${CMAKE_SOURCE_DIR}/share" )

# optional labels, e.g. set( TESTLABELS fast ) before the include; the
# criterion for "fast" is per-test wall time (<~10 s each), not module
# membership — `make check-fast` runs ctest -L fast

# -----------------------------------------------------------------------------
# Tier 1 / Tier 2 registration
#
# Tier 1 is the default and is unchanged: one serial ctest entry named after the
# suite.
#
# Tier 2 is opt-in, for the suites whose SUBJECT is communication. Such a suite
# writes
#
#     set( TESTRANKS 2 4 )
#
# before including this file, and is then registered once per rank count as
# ${TESTNAME}_np<N> under the MPI launcher — and NOT registered serially as
# well. The serial twin is deliberately absent: a Tier 2 binary run at one rank
# compiles its MPI bodies out under #ifdef BELFEM_MPI, so the remainder can pass
# vacuously and report the parallel capability green while nothing ever ran in
# parallel. Put 1 in TESTRANKS if a one-rank case is wanted; it is then an
# _np1 registration like any other, carrying the same sentinel.
#
# Every Tier 2 binary must carry the launcher sentinel and the verdict fold
# (both in tests/common/tier2_launcher_sentinel.hpp; tests/comm/test_commmpi_main.cpp
# and tests/sparse/test_sparsempi_main.cpp instantiate them). A GTEST_SKIP returns 0 from
# RUN_ALL_TESTS() and ctest reads only the process exit code, so a launcher that
# silently degrades to a single rank is otherwise indistinguishable from a run
# that passed.
# -----------------------------------------------------------------------------
if( USE_MPI AND TESTRANKS )

    # The launcher is resolved once per configure and cached. find_mpi.cmake has
    # already identified the implementation by preprocessing mpi.h
    # ( BELFEM_MPI_FLAVOR ) — a stronger check than a launcher's name or its
    # --version banner — so that verdict is reused here instead of probed again.
    if( NOT BELFEM_MPIEXEC )
        find_program( BELFEM_MPIEXEC
                NAMES mpirun mpiexec
                HINTS ${BELFEM_MPIHOME}/bin $ENV{MPI_HOME}/bin
                DOC "MPI launcher used to run Tier 2 ( multi-rank ) tests" )
    endif()

    if( NOT BELFEM_MPIEXEC )
        message( FATAL_ERROR
                "Test suite '${TESTNAME}' asks for Tier 2 ( multi-rank ) tests at ranks "
                "'${TESTRANKS}', but no MPI launcher could be found.\n"
                "  Searched ${BELFEM_MPIHOME}/bin, $ENV{MPI_HOME}/bin and PATH for "
                "mpirun / mpiexec.\n"
                "  Set MPI_HOME to the Open MPI installation, or configure with "
                "-DUSE_MPI=OFF for a serial tree ( the Tier 2 suites are then skipped )." )
    endif()

    # --oversubscribe is an Open MPI spelling. This tree refuses to configure
    # against anything else unless ALLOW_UNTESTED_MPI was passed
    # ( config/system/find_mpi.cmake ), so the flag is safe in the supported
    # configuration and is withheld — with a warning — in the escape-hatch one.
    #
    # It is passed unconditionally rather than only when ranks exceed the core
    # count: these are correctness runs, not performance runs, and Open MPI
    # refuses np > slots outright. Making the flag conditional on a detected
    # core count would make a Tier 2 suite green on a workstation and a configure
    # error on a two-core CI runner.
    # set() with no value, NOT unset(): unset() removes the normal variable and
    # thereby exposes any cache entry of the same name underneath, which would
    # then leak into the COMMAND below. An empty normal variable masks the cache
    # and expands to no argument at all.
    # The flag is withheld only when the implementation is KNOWN to be something
    # else — never merely because the probe could not read mpi.h.
    # find_mpi.cmake leaves BELFEM_MPI_FLAVOR unset when preprocessing fails and
    # only warns ( it is not evidence of an unsupported MPI ), so testing
    # STREQUAL "Open MPI" alone would silently drop --oversubscribe on that
    # edge. A two-core runner would then fail to LAUNCH np=4 rather than
    # time-share it: fail-closed, but a CI break on a configure edge that has
    # nothing to do with the tests.
    set( BELFEM_MPIEXEC_FLAGS )
    if( BELFEM_MPI_FLAVOR AND NOT BELFEM_MPI_FLAVOR STREQUAL "Open MPI" )
        message( WARNING
                "Tier 2 tests are being registered against '${BELFEM_MPI_FLAVOR}' rather "
                "than Open MPI. --oversubscribe is not passed, so any rank count above "
                "the available slots will fail to launch rather than time-share." )
    else()
        set( BELFEM_MPIEXEC_FLAGS --oversubscribe )
    endif()

    # a rank that fails inside a collective leaves its partners blocked in the
    # matching send/receive. A hung nightly is not a verdict, so every Tier 2
    # test carries a wall-clock ceiling; override with TESTTIMEOUT if a suite
    # legitimately needs longer.
    if( NOT TESTTIMEOUT )
        set( TESTTIMEOUT 300 )
    endif()

    # declared so `cmake -LH` lists it; empty means uncapped, which is the
    # normal case. See the announcement below - a cap is never silent.
    set( BELFEM_TEST_MAX_RANKS "" CACHE STRING
            "Cap on Tier 2 test rank counts; empty = uncapped. Ranks above the cap are NOT registered." )

    unset( BELFEM_TIER2_REGISTERED )
    unset( BELFEM_TIER2_SKIPPED )

    foreach( NP ${TESTRANKS} )

        # A capped rank count must be announced, never silently dropped: a rank
        # count that quietly fails to register is indistinguishable from one
        # that ran and passed, which is the failure this whole tier exists to
        # close.
        if( BELFEM_TEST_MAX_RANKS AND NP GREATER BELFEM_TEST_MAX_RANKS )
            list( APPEND BELFEM_TIER2_SKIPPED ${NP} )
        else()
            add_test( NAME ${TESTNAME}_np${NP}
                    COMMAND ${BELFEM_MPIEXEC}
                            -np ${NP}
                            ${BELFEM_MPIEXEC_FLAGS}
                            $<TARGET_FILE:test_${TESTNAME}> )

            # BELFEM_TESTRANKS is what the launcher sentinel reads. Every
            # property below must be set per generated name — the pre-Tier-2
            # code applied them to a single test name, and a per-rank loop that
            # forgets them yields tests that cannot find share/ and are
            # invisible to ctest -L.
            # 'mpi' always; any suite labels after it. Built conditionally so an
            # unset TESTLABELS does not leave a trailing empty label on the test.
            set( BELFEM_TIER2_LABELS "mpi" )
            if( TESTLABELS )
                list( APPEND BELFEM_TIER2_LABELS ${TESTLABELS} )
            endif()

            set_tests_properties( ${TESTNAME}_np${NP} PROPERTIES
                    ENVIRONMENT "${BELFEM_TEST_ENV};BELFEM_TESTRANKS=${NP}"
                    LABELS      "${BELFEM_TIER2_LABELS}"
                    PROCESSORS  ${NP}
                    TIMEOUT     ${TESTTIMEOUT} )

            list( APPEND BELFEM_TIER2_REGISTERED ${NP} )
        endif()
    endforeach()

    message( STATUS
            "Tier 2 test '${TESTNAME}': registered at ranks [${BELFEM_TIER2_REGISTERED}]" )
    if( BELFEM_TIER2_SKIPPED )
        message( WARNING
                "Tier 2 test '${TESTNAME}': ranks [${BELFEM_TIER2_SKIPPED}] NOT registered, "
                "capped by BELFEM_TEST_MAX_RANKS=${BELFEM_TEST_MAX_RANKS}. Those rank counts "
                "are not being tested by this build." )
    endif()

else()

    if( TESTRANKS AND NOT USE_MPI )
        message( STATUS
                "Test suite '${TESTNAME}' requests Tier 2 ranks [${TESTRANKS}] but this is a "
                "USE_MPI=OFF tree; the suite is not registered." )
    else()
        add_test( NAME ${TESTNAME} COMMAND test_${TESTNAME} )

        set_tests_properties( ${TESTNAME} PROPERTIES ENVIRONMENT "${BELFEM_TEST_ENV}" )

        if( TESTLABELS )
            set_tests_properties( ${TESTNAME} PROPERTIES LABELS "${TESTLABELS}" )
        endif()
    endif()

endif()

target_include_directories( test_${TESTNAME} BEFORE PRIVATE $ENV{SCLS}/include )

# the whole framework is one library target ( Add_BelfemLibrary.cmake ) that
# carries its third-party link interface; LIBLIST is no longer consulted
target_link_libraries( test_${TESTNAME} ${LIBTARGET} -lgtest -lgtest_main )
if( NOT APPLE )
    target_link_libraries( test_${TESTNAME} -pthread )
endif()

if( APPLE )
    add_custom_command(TARGET test_${TESTNAME}
            POST_BUILD
            COMMAND codesign -f -s - ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${TESTNAME}
            )
endif()


# -----------------------------------------------------------------------------
# Clear the per-suite inputs.
#
# None of these were reset before, which was harmless only because every
# tests/*/CMakeLists.txt included this file exactly once. Tier 2 breaks that
# assumption: a directory now registers a serial suite and an MPI suite side by
# side, and an inherited TESTLABELS ( tests/mesh and tests/fem both set it to
# 'fast' ) would put minutes of multi-rank work inside `make check-fast`.
# -----------------------------------------------------------------------------
unset( TESTNAME )
unset( SOURCES )
unset( TESTLABELS )
unset( TESTRANKS )
unset( TESTTIMEOUT )
