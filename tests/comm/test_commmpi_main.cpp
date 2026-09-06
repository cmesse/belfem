/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Test main for the Tier 2 ( MPI ) comm tests.
 *
 * This file owns gComm and gLog for the Tier 2 binary, following the same
 * convention as every other test_*_main.cpp. test_CommMPI.cpp declares them
 * extern; if both files defined them the split of the old in-source main()
 * would only have moved the duplicate symbol from main to the globals.
 */

#include <cstdlib>

#include <gtest/gtest.h>

#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "assert.hpp"
#include "tier2_launcher_sentinel.hpp"

belfem::Communicator gComm;
belfem::Logger       gLog( 5 );

// =============================================================================
// The launcher sentinel
//
// A GTEST_SKIP returns 0 from RUN_ALL_TESTS(), and ctest reads nothing but the
// process exit code ( config/scripts/Add_Test.cmake ; the tree has no
// gtest_discover_tests ). Almost every test in test_CommMPI.cpp skips below two
// ranks. So a binary launched WITHOUT mpirun - or with fewer ranks than it was
// registered for - reports a green ctest line while skipping everything that
// matters: the parallel capability would claim success having never run in
// parallel.
//
// This test therefore carries no skip guard of any kind. It is the one
// assertion that fails when the launcher does. Add_Test.cmake's Tier 2 path
// sets BELFEM_TESTRANKS to the rank count each test was registered for.
//
// The second Tier 2 suite ( tests/sparse, sparsempi ) appeared on 2026-09-04,
// so the sentinel body and the verdict fold now live in the shared header
// tests/common/tier2_launcher_sentinel.hpp; this file only instantiates them.
// =============================================================================

BELFEM_TIER2_LAUNCHER_SENTINEL( 2 )   // tests/comm/CMakeLists.txt registers 2 and 4

// =============================================================================

int
main( int    argc,
      char * argv[] )
{
    gComm.init( argc, argv );

    // NOTE: the two assert reactions are treated differently here, on purpose.
    //
    // test_comm_main.cpp sets BOTH set_throw_on_error( true ) and
    // set_syslog_on_error( false ), because its tests exercise error paths
    // deliberately with EXPECT_THROW. This file has no throw expectation at all
    // ( grep EXPECT_THROW/ASSERT_THROW over test_CommMPI.cpp: none; the one
    // throw-related line is an EXPECT_NO_THROW ), and test_Communicator.cpp,
    // which does exercise them, compiles into test_comm rather than into this
    // binary. So neither call could be copied over on that justification.
    //
    // throw_on_error is therefore NOT set: forcing it true would override the
    // RELEASE reaction, where a failed check must call MPI_Abort rather than
    // throw a exception GoogleTest then catches.
    //
    // syslog IS suppressed, for a reason of its own rather than by analogy.
    // The system log exists to recover what the terminal lost - the message and
    // rank of a run that died after hours, whose output is gone. A test never
    // has that problem: ctest captures it locally and the nightly job archives
    // it. The entry would be redundant, and worse than redundant: the nightly
    // is the highest-frequency writer that would ever reach this identity, and
    // its failures would dilute `journalctl -t belfem` for the long production
    // crash the hook is actually there to diagnose.
    //
    // Leaving throw_on_error alone keeps the framework contract exactly as
    // ruled: gThrowOnError = BELFEM_ASSERTIONS_ACTIVE, so debug throws at any
    // rank count and stays inspectable under one debugger per rank, while
    // release aborts the job. A rank-local throw does strand its peers in a
    // matching send/receive; that is the accepted price of a debuggable
    // failure, and the CTest TIMEOUT set by Add_Test.cmake is what stops the
    // nightly sitting on it.
    belfem::assert::set_syslog_on_error( false );

    testing::InitGoogleTest( &argc, argv );

    int aResult = RUN_ALL_TESTS();

    // Fold every rank's verdict into every rank's exit code.
    //
    // The launcher-folding assumption this removes was MEASURED on 2026-08-30:
    // a deliberate rank-3-only failure turned commmpi_np4 red while np2 stayed
    // green, so the reduction is what carries the verdict. belfem::allreduce
    // is a MAX-reduction over gComm.world(), which at this point in main() is
    // MPI_COMM_WORLD ( init has run, finalize has not ); its serial branch is
    // the identity copy, so no #ifdef is needed here.
    aResult = belfem::test::fold_verdict( aResult );

    gComm.finalize();

    return aResult;
}
