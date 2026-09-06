/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Tier 1 tests for comm utility functions: comm_split, comm_splitcount,
 * comm_type mapping. These run without MPI in a single process.
 * See: tests_06_comm.md §1
 *
 * NOTE: comm_tag() tests are in Tier 2 (test_CommMPI.cpp) because
 * comm_tag() uses % gComm.max_tag(), which is 0 in non-MPI mode (BUG-C1).
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "commtypes.hpp"
#include "cl_Cell.hpp"

// =============================================================================
// §1.1  comm_split  [semantic]
// =============================================================================

TEST( CommSplit, CommSplitZeroLength )
{
    belfem::Cell< int > tChunks = belfem::comm_split( 0 );
    EXPECT_EQ( tChunks.size(), 0u );
}

TEST( CommSplit, CommSplitSmallMessage )
{
    // less than 64K → 1 chunk of that length
    belfem::Cell< int > tChunks = belfem::comm_split( 100 );
    ASSERT_EQ( tChunks.size(), 1u );
    EXPECT_EQ( tChunks( 0 ), 100 );
}

TEST( CommSplit, CommSplitExactBoundary )
{
    // exactly 64K → 1 chunk of 64K
    belfem::index_t tLen = belfem::gMaxCommChunkLength;
    belfem::Cell< int > tChunks = belfem::comm_split( tLen );
    ASSERT_EQ( tChunks.size(), 1u );
    EXPECT_EQ( tChunks( 0 ), static_cast< int >( tLen ) );
}

TEST( CommSplit, CommSplitBoundaryPlusOne )
{
    // 64K + 1 → 2 chunks: 64K and 1
    belfem::index_t tLen = belfem::gMaxCommChunkLength + 1;
    belfem::Cell< int > tChunks = belfem::comm_split( tLen );
    ASSERT_EQ( tChunks.size(), 2u );
    EXPECT_EQ( tChunks( 0 ), belfem::gMaxCommChunkLength );
    EXPECT_EQ( tChunks( 1 ), 1 );
}

TEST( CommSplit, CommSplitLargeMessage )
{
    // 3×64K + 100 → 4 chunks: three 64K and one 100
    belfem::index_t tLen = 3 * belfem::gMaxCommChunkLength + 100;
    belfem::Cell< int > tChunks = belfem::comm_split( tLen );
    ASSERT_EQ( tChunks.size(), 4u );
    EXPECT_EQ( tChunks( 0 ), belfem::gMaxCommChunkLength );
    EXPECT_EQ( tChunks( 1 ), belfem::gMaxCommChunkLength );
    EXPECT_EQ( tChunks( 2 ), belfem::gMaxCommChunkLength );
    EXPECT_EQ( tChunks( 3 ), 100 );
}

TEST( CommSplit, CommSplitTotalEqualsInput )
{
    // verify for several sizes that sum of chunks == input
    belfem::index_t tSizes[] = { 0, 1, 100, 65535, 65536, 65537, 200000 };
    for( auto tLen : tSizes )
    {
        belfem::Cell< int > tChunks = belfem::comm_split( tLen );
        belfem::index_t tTotal = 0;
        for( belfem::index_t i = 0; i < tChunks.size(); ++i )
        {
            tTotal += tChunks( i );
        }
        EXPECT_EQ( tTotal, tLen );
    }
}

TEST( CommSplit, CommSplitNoChunkExceedsMax )
{
    belfem::Cell< int > tChunks = belfem::comm_split( 200000 );
    for( belfem::index_t i = 0; i < tChunks.size(); ++i )
    {
        EXPECT_LE( tChunks( i ), belfem::gMaxCommChunkLength );
    }
}

// =============================================================================
// §1.2  comm_splitcount  [semantic]
// =============================================================================

TEST( CommSplitcount, CommSplitcountZeroLength )
{
    EXPECT_EQ( belfem::comm_splitcount( 0 ), 0u );
}

TEST( CommSplitcount, CommSplitcountConsistentWithSplit )
{
    // comm_splitcount(length) == comm_split(length).size() * (comm_size()-1)
    belfem::index_t tLen = 200000;
    belfem::Cell< int > tChunks = belfem::comm_split( tLen );
    belfem::index_t tExpected = tChunks.size() * ( belfem::comm_size() - 1 );
    EXPECT_EQ( belfem::comm_splitcount( tLen ), tExpected );
}

TEST( CommSplitcount, CommSplitcountVectorConsistentWithSplit )
{
    // (idea from ChatGPT) — vector overload with root exclusion
    belfem::Vector< belfem::index_t > tLengths = {
        0,
        static_cast< belfem::index_t >( belfem::gMaxCommChunkLength + 3 ),
        5,
        static_cast< belfem::index_t >( 2 * belfem::gMaxCommChunkLength ) };
    belfem::proc_t tRoot = 2;

    belfem::index_t tExpected = 0;
    for( belfem::proc_t p = 0;
         p < static_cast< belfem::proc_t >( tLengths.length() ); ++p )
    {
        if( p == tRoot ) continue;
        tExpected += belfem::comm_split( tLengths( p ) ).size();
    }

    EXPECT_EQ( belfem::comm_splitcount( tLengths, tRoot ), tExpected );
}

// =============================================================================
// §1.4  comm_type Mapping  [semantic]
// =============================================================================

#ifndef BELFEM_MPI

// In non-MPI builds, all comm_type<T>() return 0
TEST( CommType, CommTypeNonMpiReturnsZero )
{
    EXPECT_EQ( belfem::comm_type< int >(), 0 );
    EXPECT_EQ( belfem::comm_type< double >(), 0 );
    EXPECT_EQ( belfem::comm_type< float >(), 0 );
    EXPECT_EQ( belfem::comm_type< char >(), 0 );
    EXPECT_EQ( belfem::comm_type< unsigned long int >(), 0 );
    EXPECT_EQ( belfem::comm_type< bool >(), 0 );
}

#else

// In MPI builds, verify against known MPI constants
TEST( CommType, CommTypeInt )
{
    EXPECT_EQ( belfem::comm_type< int >(), MPI_INT );
}

TEST( CommType, CommTypeDouble )
{
    EXPECT_EQ( belfem::comm_type< double >(), MPI_DOUBLE );
}

TEST( CommType, CommTypeFloat )
{
    EXPECT_EQ( belfem::comm_type< float >(), MPI_FLOAT );
}

TEST( CommType, CommTypeChar )
{
    EXPECT_EQ( belfem::comm_type< char >(), MPI_CHAR );
}

TEST( CommType, CommTypeUnsignedLong )
{
    EXPECT_EQ( belfem::comm_type< unsigned long int >(), MPI_UNSIGNED_LONG );
}

TEST( CommType, CommTypeBool )
{
    EXPECT_EQ( belfem::comm_type< bool >(), MPI_CXX_BOOL );
}

#endif // BELFEM_MPI
