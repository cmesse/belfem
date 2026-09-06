/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Test main for the Tier 2 ( MPI ) sparse tests: parallel nested dissection.
 *
 * Owns gComm and gLog for the binary; test_SparseMPI.cpp declares them extern.
 * The launcher sentinel and the verdict fold come from the shared header -
 * see tests/common/tier2_launcher_sentinel.hpp for why both are mandatory.
 */

#include <gtest/gtest.h>

#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "assert.hpp"
#include "tier2_launcher_sentinel.hpp"

belfem::Communicator gComm;
belfem::Logger       gLog( 5 );

BELFEM_TIER2_LAUNCHER_SENTINEL( 2 )   // tests/sparse/CMakeLists.txt registers 2 and 4

int
main( int    argc,
      char * argv[] )
{
    gComm.init( argc, argv );

    // throw_on_error is deliberately NOT set: a rank-local throw would strand
    // the peers in the very collectives these tests exist to exercise, and the
    // release reaction ( MPI_Abort ) is what a hang here must turn into. The
    // CTest TIMEOUT set by Add_Test.cmake bounds the debug case. Syslog is
    // suppressed for the reason test_commmpi_main.cpp gives.
    belfem::assert::set_syslog_on_error( false );

    testing::InitGoogleTest( &argc, argv );

    int aResult = RUN_ALL_TESTS();

    aResult = belfem::test::fold_verdict( aResult );

    gComm.finalize();
    return aResult;
}
