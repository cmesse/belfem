/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Shared by every Tier 2 ( multi-rank ) test main. Two things live here, and
 * both are load-bearing:
 *
 * 1. The launcher sentinel. A GTEST_SKIP returns 0 from RUN_ALL_TESTS(), and
 *    ctest reads nothing but the process exit code, so a Tier 2 binary
 *    launched without mpirun - or with fewer ranks than it was registered
 *    for - would skip every rank-guarded test and still report green. The
 *    sentinel is the one test that fails when the launcher does.
 *    Add_Test.cmake's Tier 2 path sets BELFEM_TESTRANKS to the rank count each
 *    registration asked for.
 *
 * 2. The verdict fold. Measured 2026-08-30: a deliberate rank-3-only failure
 *    turned commmpi_np4 red while np2 stayed green - the launcher does not
 *    fold exit codes, the MAX-allreduce below does.
 *
 * A sentinel that only some Tier 2 binaries carry is worse than none, because
 * it makes the guarantee look uniform when it is not. Include this header;
 * do not copy it.
 */

#ifndef BELFEM_TIER2_LAUNCHER_SENTINEL_HPP
#define BELFEM_TIER2_LAUNCHER_SENTINEL_HPP

#include <cstdlib>

#include <gtest/gtest.h>

#include "commtools.hpp"

namespace belfem
{
    namespace test
    {
        /**
         * body of the sentinel test. aRankFloor is the smallest rank count the
         * binary's own guards assume ( the smallest entry of its TESTRANKS );
         * the environment variable alone is not enough, because a hand-run
         * `BELFEM_TESTRANKS=1 ./test_x` would satisfy the equality, skip every
         * rank-guarded test and still exit green.
         */
        inline void
        check_launcher( const long aRankFloor )
        {
            const char * tExpected = std::getenv( "BELFEM_TESTRANKS" );

            ASSERT_NE( tExpected, nullptr )
                << "BELFEM_TESTRANKS is not set. This is a Tier 2 binary and must be "
                   "registered through the TESTRANKS path in Add_Test.cmake, which sets it. "
                   "Running it by hand proves nothing about the parallel paths.";

            // strtol with an end-pointer check, not atoi: atoi( "2junk" ) is 2 and
            // atoi( "junk" ) is 0, both silently
            char *     tEnd       = nullptr ;
            const long tRequested = std::strtol( tExpected, &tEnd, 10 );

            ASSERT_TRUE( tEnd != tExpected && tEnd != nullptr && *tEnd == '\0' )
                << "BELFEM_TESTRANKS is not a plain integer: '" << tExpected << "'";

            EXPECT_EQ( static_cast< long >( belfem::comm_size() ), tRequested )
                << "launched at " << belfem::comm_size() << " rank(s), registered for "
                << tRequested << ". The MPI launcher did not do its job, and every "
                   "rank-guarded test in this binary just skipped.";

            EXPECT_GE( static_cast< long >( belfem::comm_size() ), aRankFloor )
                << "this binary is meaningless below " << aRankFloor
                << " rank(s): its tests guard on comm_size() and would skip.";
        }

        /**
         * fold every rank's RUN_ALL_TESTS() verdict into every rank's exit
         * code. belfem::allreduce is a MAX-reduction over gComm.world(),
         * which between init and finalize is MPI_COMM_WORLD; its serial
         * branch is the identity copy, so no #ifdef is needed.
         */
        inline int
        fold_verdict( int aResult )
        {
            int tGlobal = aResult;
            belfem::allreduce( &aResult, &tGlobal, 1 );
            return tGlobal;
        }
    }
}

//! defines the sentinel test in the including main; aFloor = smallest TESTRANKS
#define BELFEM_TIER2_LAUNCHER_SENTINEL( aFloor )                       \
    TEST( Tier2Launcher, RankCountMatchesRegistration )                \
    {                                                                  \
        belfem::test::check_launcher( aFloor );                        \
    }

#endif // BELFEM_TIER2_LAUNCHER_SENTINEL_HPP
