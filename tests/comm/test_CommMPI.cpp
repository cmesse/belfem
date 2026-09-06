/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Tier 2 MPI tests: send/receive, broadcast, distribute/collect,
 * chunking protocol, comm_tag.
 * See: tests_06_comm.md §1.3, §3, §4, §5
 *
 * Requires: mpirun -np N  (minimum 2 ranks, some tests require 4). Registered
 * through the TESTRANKS path in config/scripts/Add_Test.cmake, never by hand.
 *
 * main(), gComm/gLog and the MPI lifecycle live in test_commmpi_main.cpp.
 *
 * Tests that need more ranks than available use GTEST_SKIP(). That is safe ONLY
 * because test_commmpi_main.cpp carries an unskippable rank sentinel: a skip
 * returns 0 from RUN_ALL_TESTS(), so without it a binary launched at one rank
 * would skip this entire file and still report success.
 */

#include <gtest/gtest.h>
#include <string>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

// defined in test_commmpi_main.cpp, which owns the MPI lifecycle for this
// binary ( see cl_Communicator.hpp for the canonical extern declarations )
extern belfem::Communicator gComm;
extern belfem::Logger       gLog;

namespace
{
    const belfem::real tEps = 1e-12;
}

// =============================================================================
// §1.3  comm_tag  [semantic] — MPI-only (BUG-C1: % 0 in non-MPI)
// =============================================================================

#ifdef BELFEM_MPI

TEST( CommTag, CommTagSymmetric )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    EXPECT_EQ( belfem::comm_tag( 0, 1 ), belfem::comm_tag( 1, 0 ) );
    if( belfem::comm_size() >= 3 )
    {
        EXPECT_EQ( belfem::comm_tag( 0, 2 ), belfem::comm_tag( 2, 0 ) );
        EXPECT_EQ( belfem::comm_tag( 1, 2 ), belfem::comm_tag( 2, 1 ) );
    }
}

TEST( CommTag, CommTagEven )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    EXPECT_EQ( belfem::comm_tag( 0, 1 ) % 2, 0 );
    if( belfem::comm_size() >= 3 )
    {
        EXPECT_EQ( belfem::comm_tag( 0, 2 ) % 2, 0 );
    }
}

TEST( CommTag, CommTagDistinctPairs )
{
    if( belfem::comm_size() < 3 ) GTEST_SKIP() << "Requires 3 ranks";

    int tTag01 = belfem::comm_tag( 0, 1 );
    int tTag02 = belfem::comm_tag( 0, 2 );
    int tTag12 = belfem::comm_tag( 1, 2 );

    EXPECT_NE( tTag01, tTag02 );
    EXPECT_NE( tTag01, tTag12 );
    EXPECT_NE( tTag02, tTag12 );
}

TEST( CommTag, CommTagNonNegative )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    EXPECT_GE( belfem::comm_tag( 0, 1 ), 0 );
}

TEST( CommTag, CommTagWithinMaxTag )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    EXPECT_LT( belfem::comm_tag( 0, 1 ), gComm.max_tag() );
}

#endif // BELFEM_MPI

// =============================================================================
// §3.1  Scalar Send/Receive  [semantic]
// =============================================================================

#ifdef BELFEM_MPI

// =============================================================================
// allreduce_min — the MIN sibling of the MAX wrapper the test main uses
// =============================================================================

TEST( CommReduce, AllreduceMinTakesSmallestRank )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    // rank r contributes 100 + r, so the minimum is rank 0's value
    belfem::int_t tLocal  = 100 + ( belfem::int_t ) belfem::comm_rank();
    belfem::int_t tGlobal = -1 ;
    belfem::allreduce_min( &tLocal, &tGlobal, 1 );
    EXPECT_EQ( tGlobal, 100 );

    // and the value every rank sees is the same one
    belfem::int_t tCheck = tGlobal ;
    belfem::int_t tMax   = 0 ;
    belfem::allreduce( &tCheck, &tMax, 1 );
    EXPECT_EQ( tMax, tGlobal );
}

TEST( CommReduce, AllreduceMinZeroWins )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    // the memory-budget contract: one rank that could not measure ( 0 )
    // pulls every rank to 0
    belfem::int_t tLocal  = belfem::comm_rank() == 1 ? 0 : 5000 ;
    belfem::int_t tGlobal = -1 ;
    belfem::allreduce_min( &tLocal, &tGlobal, 1 );
    EXPECT_EQ( tGlobal, 0 );
}

TEST( CommSendRecv, SendReceiveInt )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        int tVal = 42;
        belfem::send( tVal, 1 );
    }
    else if( tRank == 1 )
    {
        int tVal = 0;
        belfem::receive( tVal, 0 );
        EXPECT_EQ( tVal, 42 );
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendReceiveReal )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        belfem::real tVal = 3.14159;
        belfem::send( tVal, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::real tVal = 0.0;
        belfem::receive( tVal, 0 );
        EXPECT_NEAR( tVal, 3.14159, tEps );
    }

    belfem::comm_barrier();
}

// =============================================================================
// §3.2  Vector Send/Receive  [semantic]
// =============================================================================

TEST( CommSendRecv, SendReceiveVectorSmall )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        belfem::Vector< belfem::real > tVec( 10 );
        for( belfem::index_t i = 0; i < 10; ++i )
        {
            tVec( i ) = static_cast< belfem::real >( i * i );
        }
        belfem::send( tVec, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Vector< belfem::real > tVec;
        belfem::receive( tVec, 0 );

        ASSERT_EQ( tVec.length(), 10u );
        for( belfem::index_t i = 0; i < 10; ++i )
        {
            EXPECT_NEAR( tVec( i ), static_cast< belfem::real >( i * i ), tEps );
        }
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendReceiveVectorEmpty )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        belfem::Vector< belfem::real > tVec;
        belfem::send( tVec, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Vector< belfem::real > tVec;
        belfem::receive( tVec, 0 );
        EXPECT_EQ( tVec.length(), 0u );
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendReceiveVectorLargeChunked )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();
    belfem::index_t tLen = 100000;   // exceeds 64K chunk

    if( tRank == 0 )
    {
        belfem::Vector< belfem::real > tVec( tLen );
        for( belfem::index_t i = 0; i < tLen; ++i )
        {
            tVec( i ) = std::sin( 0.001 * i );
        }
        belfem::send( tVec, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Vector< belfem::real > tVec;
        belfem::receive( tVec, 0 );

        ASSERT_EQ( tVec.length(), tLen );
        for( belfem::index_t i = 0; i < tLen; ++i )
        {
            EXPECT_NEAR( tVec( i ), std::sin( 0.001 * i ), tEps );
        }
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendReceiveVectorExactChunkBoundary )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();
    belfem::index_t tLen = belfem::gMaxCommChunkLength;

    if( tRank == 0 )
    {
        belfem::Vector< belfem::real > tVec( tLen, 7.0 );
        belfem::send( tVec, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Vector< belfem::real > tVec;
        belfem::receive( tVec, 0 );

        ASSERT_EQ( tVec.length(), tLen );
        EXPECT_NEAR( tVec( 0 ), 7.0, tEps );
        EXPECT_NEAR( tVec( tLen - 1 ), 7.0, tEps );
    }

    belfem::comm_barrier();
}

// =============================================================================
// §3.4  Matrix Send/Receive  [semantic]
// =============================================================================

TEST( CommSendRecv, SendReceiveMatrixSmall )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        belfem::Matrix< belfem::real > tMat( 3, 4 );
        for( size_t i = 0; i < 3; ++i )
        for( size_t j = 0; j < 4; ++j )
        {
            tMat( i, j ) = static_cast< belfem::real >( 10 * i + j );
        }
        belfem::send( tMat, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Matrix< belfem::real > tMat;
        belfem::receive( tMat, 0 );

        EXPECT_EQ( tMat.n_rows(), 3u );
        EXPECT_EQ( tMat.n_cols(), 4u );
        EXPECT_NEAR( tMat( 2, 3 ), 23.0, tEps );
    }

    belfem::comm_barrier();
}

// Regression: a persistent matrix that once held a larger shape keeps its old
// allocation ( Blaze resize(...,false) never shrinks capacity_ ). The transfer
// must send the padded footprint of the CURRENT shape, not capacity() — the
// stale-large count overflowed the receiver's exact-fit buffer.
TEST( CommSendRecv, SendReceiveMatrixAfterShrink )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        // grow, then shrink: capacity() now exceeds the current footprint
        belfem::Matrix< belfem::real > tMat( 20, 20 );
        tMat.set_size( 3, 4 );
        for( size_t i = 0; i < 3; ++i )
        for( size_t j = 0; j < 4; ++j )
        {
            tMat( i, j ) = static_cast< belfem::real >( 10 * i + j );
        }
        belfem::send( tMat, 1 );
    }
    else if( tRank == 1 )
    {
        // fresh receiver: exact-fit allocation for ( 3, 4 )
        belfem::Matrix< belfem::real > tMat;
        belfem::receive( tMat, 0 );

        ASSERT_EQ( tMat.n_rows(), 3u );
        ASSERT_EQ( tMat.n_cols(), 4u );
        for( size_t i = 0; i < 3; ++i )
        for( size_t j = 0; j < 4; ++j )
        {
            EXPECT_NEAR( tMat( i, j ),
                static_cast< belfem::real >( 10 * i + j ), tEps );
        }
    }

    belfem::comm_barrier();
}

// =============================================================================
// §3.4b  Cell Send/Receive  [semantic]
// =============================================================================

TEST( CommSendRecv, SendReceiveCellInt )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        belfem::Cell< int > tData{ 10, 20, 30, 40, 50 };
        belfem::send( tData, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Cell< int > tData;
        belfem::receive( tData, 0 );

        ASSERT_EQ( tData.size(), 5u );
        EXPECT_EQ( tData( 0 ), 10 );
        EXPECT_EQ( tData( 4 ), 50 );
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendReceiveCellEmpty )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        belfem::Cell< int > tData;
        belfem::send( tData, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Cell< int > tData;
        belfem::receive( tData, 0 );
        EXPECT_EQ( tData.size(), 0u );
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendReceiveMatrixEmpty )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        belfem::Matrix< belfem::real > tMat( 0, 0 );
        belfem::send( tMat, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Matrix< belfem::real > tMat;
        belfem::receive( tMat, 0 );
        EXPECT_EQ( tMat.n_rows(), 0u );
        EXPECT_EQ( tMat.n_cols(), 0u );
    }

    belfem::comm_barrier();
}

// =============================================================================
// §3.5  String Send/Receive  [semantic]
// =============================================================================

TEST( CommSendRecv, SendReceiveString )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        std::string tMsg = "Hello BELFEM";
        belfem::send( tMsg, 1 );
    }
    else if( tRank == 1 )
    {
        std::string tMsg;
        belfem::receive( tMsg, 0 );
        EXPECT_EQ( tMsg, "Hello BELFEM" );
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendReceiveEmptyString )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        std::string tMsg;
        belfem::send( tMsg, 1 );
    }
    else if( tRank == 1 )
    {
        std::string tMsg;
        belfem::receive( tMsg, 0 );
        EXPECT_TRUE( tMsg.empty() );
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendReceiveStringSpecialChars )
{
    // (idea from ChatGPT)
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        std::string tMsg = "Hello BELFEM\nwith spaces\tand tabs";
        belfem::send( tMsg, 1 );
    }
    else if( tRank == 1 )
    {
        std::string tMsg;
        belfem::receive( tMsg, 0 );
        EXPECT_EQ( tMsg, "Hello BELFEM\nwith spaces\tand tabs" );
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, SendToSelfIsNoOp )
{
    // (idea from Gemini) — should not deadlock or corrupt
    int tVal = 7;
    EXPECT_NO_THROW( belfem::send( tVal, belfem::comm_rank() ) );
}

// =============================================================================
// §4.1  Broadcast  [semantic]
// =============================================================================

TEST( CommBroadcast, BroadcastScalar )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::real tVal = ( belfem::comm_rank() == 0 ) ? 99.0 : 0.0;
    belfem::broadcast( tVal, 0 );
    EXPECT_NEAR( tVal, 99.0, tEps );
}

TEST( CommBroadcast, BroadcastVector )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::Vector< belfem::real > tVec;
    if( belfem::comm_rank() == 0 )
    {
        tVec.set_size( 100 );
        for( belfem::index_t i = 0; i < 100; ++i )
        {
            tVec( i ) = static_cast< belfem::real >( i );
        }
    }

    belfem::broadcast( tVec, 0 );

    ASSERT_EQ( tVec.length(), 100u );
    EXPECT_NEAR( tVec( 0 ), 0.0, tEps );
    EXPECT_NEAR( tVec( 99 ), 99.0, tEps );
}

TEST( CommBroadcast, BroadcastVectorEmpty )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::Vector< belfem::real > tVec;
    belfem::broadcast( tVec, 0 );
    EXPECT_EQ( tVec.length(), 0u );
}

TEST( CommBroadcast, BroadcastCellString )
{
    // (idea from ChatGPT/Gemini)
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::Cell< std::string > tData;
    if( belfem::comm_rank() == 0 )
    {
        tData.push( "foo" );
        tData.push( "bar" );
        tData.push( "baz qux" );
    }

    belfem::broadcast( tData, 0 );

    ASSERT_EQ( tData.size(), 3u );
    EXPECT_EQ( tData( 0 ), "foo" );
    EXPECT_EQ( tData( 1 ), "bar" );
    EXPECT_EQ( tData( 2 ), "baz qux" );
}

// =============================================================================
// §4.2  Distribute / Collect  [semantic]
// =============================================================================

TEST( CommDistributeCollect, DistributeCollectScalar )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();
    belfem::proc_t tSize = belfem::comm_size();

    // each rank distributes its rank number
    belfem::Cell< int > tSendData( tSize, 0 );
    for( belfem::proc_t p = 0; p < tSize; ++p )
    {
        tSendData( p ) = tRank * 100 + p;
    }
    belfem::distribute( tSendData );

    // collect from all ranks
    belfem::Cell< int > tRecvData;
    belfem::collect( tRecvData, tRank * 100 + tRank );

    ASSERT_EQ( static_cast< belfem::proc_t >( tRecvData.size() ), tSize );
    for( belfem::proc_t p = 0; p < tSize; ++p )
    {
        EXPECT_EQ( tRecvData( p ), p * 100 + tRank );
    }

    belfem::comm_barrier();
}

TEST( CommSendRecv, RingPassVector4Ranks )
{
    if( belfem::comm_size() < 4 ) GTEST_SKIP() << "Requires 4 ranks";

    belfem::proc_t tRank = belfem::comm_rank();
    belfem::proc_t tSize = belfem::comm_size();
    belfem::proc_t tNext = ( tRank + 1 ) % tSize;
    belfem::proc_t tPrev = ( tRank + tSize - 1 ) % tSize;

    // each rank sends a vector stamped with its rank to the next rank
    belfem::Vector< belfem::real > tSend( 10 );
    for( belfem::index_t i = 0; i < 10; ++i )
    {
        tSend( i ) = 100.0 * tRank + static_cast< belfem::real >( i );
    }

    // non-blocking: send to next, receive from previous
    belfem::send( tSend, tNext );

    belfem::Vector< belfem::real > tRecv;
    belfem::receive( tRecv, tPrev );

    // verify: received vector should be stamped with previous rank's number
    ASSERT_EQ( tRecv.length(), 10u );
    for( belfem::index_t i = 0; i < 10; ++i )
    {
        EXPECT_NEAR( tRecv( i ), 100.0 * tPrev + static_cast< belfem::real >( i ), tEps );
    }

    belfem::comm_barrier();
}

// =============================================================================
// §4.3  Barrier  [semantic]
// =============================================================================

TEST( CommBarrier, BarrierDoesNotDeadlock )
{
    // trivial — just verify it returns on all ranks
    belfem::comm_barrier();
}

// =============================================================================
// §5.1  Chunking Protocol  [semantic]
// =============================================================================

TEST( CommChunking, VectorChunkBoundaryMinus1 )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();
    belfem::index_t tLen = belfem::gMaxCommChunkLength - 1;

    if( tRank == 0 )
    {
        belfem::Vector< belfem::real > tVec( tLen, 1.0 );
        belfem::send( tVec, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Vector< belfem::real > tVec;
        belfem::receive( tVec, 0 );
        ASSERT_EQ( tVec.length(), tLen );
        EXPECT_NEAR( tVec( tLen - 1 ), 1.0, tEps );
    }

    belfem::comm_barrier();
}

TEST( CommChunking, VectorChunkBoundaryPlus1 )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();
    belfem::index_t tLen = belfem::gMaxCommChunkLength + 1;

    if( tRank == 0 )
    {
        belfem::Vector< belfem::real > tVec( tLen, 2.0 );
        belfem::send( tVec, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Vector< belfem::real > tVec;
        belfem::receive( tVec, 0 );
        ASSERT_EQ( tVec.length(), tLen );
        EXPECT_NEAR( tVec( 0 ), 2.0, tEps );
        EXPECT_NEAR( tVec( tLen - 1 ), 2.0, tEps );
    }

    belfem::comm_barrier();
}

TEST( CommChunking, MultipleConsecutiveTransfers )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    belfem::proc_t tRank = belfem::comm_rank();

    if( tRank == 0 )
    {
        belfem::Vector< belfem::real > tA( 10, 1.0 );
        belfem::Vector< belfem::real > tB( 20, 2.0 );
        belfem::send( tA, 1 );
        belfem::send( tB, 1 );
    }
    else if( tRank == 1 )
    {
        belfem::Vector< belfem::real > tA;
        belfem::Vector< belfem::real > tB;
        belfem::receive( tA, 0 );
        belfem::receive( tB, 0 );

        ASSERT_EQ( tA.length(), 10u );
        ASSERT_EQ( tB.length(), 20u );
        EXPECT_NEAR( tA( 0 ), 1.0, tEps );
        EXPECT_NEAR( tB( 0 ), 2.0, tEps );
    }

    belfem::comm_barrier();
}

// =============================================================================
// §6  share / receive  [semantic] — the chunked, ASYMMETRIC transfer
//
// Added 2026-08-30. `belfem::share` had no test call site anywhere under
// tests/ ( searched tests/**/*.cpp for `share(` ), despite being the path
// doc/coding_philosophy.md MANDATES for large Vector/Cell payloads, and despite
// its share/receive pair being non-collective: calling either on the wrong rank
// is a deadlock, not a wrong answer.
//
// Payloads deliberately cross the 64 K chunk boundary with a partial tail
// ( 2 * gMaxCommChunkLength + 7 → three chunks ), because the whole reason
// share() exists rather than broadcast() is that it chunks. A payload that fits
// in one chunk would exercise none of that.
// =============================================================================

TEST( ShareReceive, ShareVectorCrossesChunkBoundary )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    const belfem::index_t tLength =
            2 * static_cast< belfem::index_t >( belfem::gMaxCommChunkLength ) + 7 ;

    belfem::Vector< belfem::real > tData;

    if( belfem::comm_rank() == 0 )
    {
        // root fills a reproducible, position-dependent pattern: a constant
        // payload would pass even if a chunk landed at the wrong offset
        tData.set_size( tLength );
        for( belfem::index_t k = 0; k < tLength; ++k )
        {
            tData( k ) = static_cast< belfem::real >( k ) * 0.5 - 1.0 ;
        }

        belfem::share( tData );
    }
    else
    {
        belfem::receive( tData );

        ASSERT_EQ( tData.length(), tLength )
            << "receive() must size the destination from the wire, not from the caller";

        for( belfem::index_t k = 0; k < tLength; ++k )
        {
            ASSERT_NEAR( tData( k ),
                         static_cast< belfem::real >( k ) * 0.5 - 1.0,
                         tEps ) << "payload differs at index " << k;
        }
    }
}

TEST( ShareReceive, ShareCellCrossesChunkBoundary )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    const belfem::index_t tSize =
            2 * static_cast< belfem::index_t >( belfem::gMaxCommChunkLength ) + 7 ;

    belfem::Cell< belfem::index_t > tData;

    if( belfem::comm_rank() == 0 )
    {
        tData.set_size( tSize, 0 );
        for( belfem::index_t k = 0; k < tSize; ++k )
        {
            tData( k ) = k * 3 + 1 ;
        }

        belfem::share( tData );
    }
    else
    {
        belfem::receive( tData );

        ASSERT_EQ( tData.size(), tSize );

        for( belfem::index_t k = 0; k < tSize; ++k )
        {
            ASSERT_EQ( tData( k ), k * 3 + 1 ) << "payload differs at index " << k;
        }
    }
}

TEST( ShareReceive, ShareEmptyVectorIsWellDefined )
{
    if( belfem::comm_size() < 2 ) GTEST_SKIP() << "Requires 2 ranks";

    // data() on an empty container is well-defined but possibly null, and the
    // zero-length transfer must still agree on the size rather than hang
    belfem::Vector< belfem::real > tData;

    if( belfem::comm_rank() == 0 )
    {
        tData.set_size( 0 );
        belfem::share( tData );
    }
    else
    {
        belfem::receive( tData );
        EXPECT_EQ( tData.length(), 0u );
    }
}

#endif // BELFEM_MPI

// =============================================================================
// The MPI lifecycle main() lives in test_commmpi_main.cpp, together with the
// launcher sentinel and the Allreduce that folds every rank's verdict into
// every rank's exit code. It is not here because Add_Test.cmake always compiles
// test_<name>_main.cpp into the binary, so a main() in this file would be a
// duplicate symbol.
// =============================================================================
