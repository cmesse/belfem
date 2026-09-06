/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for DynamicBitset container
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §2 (test matrix)
 *
 * All tests in this file are [valgrind] candidates (manual malloc/free).
 */

#include <gtest/gtest.h>
#include <cstdint>
#include <cstring>
#include <string>

#include <algorithm>
#include <random>
#include <set>
#include <vector>

#include "typedefs.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Cell.hpp"

namespace
{
    //! Reads the data pointer through the const overload. The mutable overload
    //! is private on purpose: a write through it bypasses the summary bitmaps,
    //! after which where() and reset() silently miss the affected words.
    const uint64_t *
    data_ptr( const belfem::DynamicBitset & aBitset )
    {
        return aBitset.data();
    }

    //! Naive reference for where(): tests every bit individually.
    std::vector< belfem::index_t >
    reference_where( const belfem::DynamicBitset & aBitset )
    {
        std::vector< belfem::index_t > aIndices;

        for( belfem::index_t k = 0; k < aBitset.size(); ++k )
        {
            if( aBitset.test( k ) )
            {
                aIndices.push_back( k );
            }
        }
        return aIndices;
    }

    //! Asserts that where() matches the naive reference, is strictly ascending,
    //! agrees between both aAssumeSparse settings, and left the summaries tight.
    void
    expect_where_matches_reference( const belfem::DynamicBitset & aBitset )
    {
        belfem::Cell< belfem::index_t > tSparse;
        aBitset.where( tSparse );

        const std::vector< belfem::index_t > tReference = reference_where( aBitset );

        ASSERT_EQ( tSparse.size(), tReference.size() );

        for( belfem::index_t k = 0; k < tSparse.size(); ++k )
        {
            EXPECT_EQ( tSparse( k ), tReference[ k ] );
        }

        for( belfem::index_t k = 1; k < tSparse.size(); ++k )
        {
            EXPECT_GT( tSparse( k ), tSparse( k - 1 ) );
        }

        belfem::Cell< belfem::index_t > tDense;
        aBitset.where( tDense, false );

        ASSERT_EQ( tDense.size(), tSparse.size() );

        for( belfem::index_t k = 0; k < tDense.size(); ++k )
        {
            EXPECT_EQ( tDense( k ), tSparse( k ) );
        }

        EXPECT_TRUE( aBitset.summaries_are_tight() );
    }
}

// =============================================================================
// Value-parameterized test infrastructure
// =============================================================================
// Sizes chosen to exercise block boundaries and edge cases:
//   1       — single bit
//   7       — less than one byte
//   63      — one bit short of a full block
//   64      — exactly one block
//   65      — one bit into second block
//   127     — one bit short of two full blocks
//   128     — exactly two blocks
//   1024    — large, multi-block

class DynamicBitset : public ::testing::TestWithParam< belfem::index_t > {};

INSTANTIATE_TEST_SUITE_P(
    Sizes,
    DynamicBitset,
    ::testing::Values( 1, 7, 63, 64, 65, 127, 128, 1024 )
);

// =============================================================================
// §2.1  Construction & Destruction  [semantic]
// =============================================================================

TEST_P( DynamicBitset, ConstructorInitializesAllZero )
{
    belfem::DynamicBitset tBs( GetParam() );
    EXPECT_EQ( tBs.count(), 0u );
}

TEST_P( DynamicBitset, SizeMatchesConstructorArg )
{
    belfem::DynamicBitset tBs( GetParam() );
    EXPECT_EQ( tBs.size(), GetParam() );
}

TEST_P( DynamicBitset, MemorySizeCorrect )
{
    belfem::index_t tN = GetParam();
    belfem::DynamicBitset tBs( tN );
    EXPECT_EQ( tBs.memory(), ( tN + 63 ) / 64 );
}

TEST_P( DynamicBitset, CopyConstructorDeepCopies )
{
    belfem::index_t tN = GetParam();
    belfem::DynamicBitset tOriginal( tN );

    tOriginal.set( 0 );
    if( tN > 1 )
    {
        tOriginal.set( tN - 1 );
    }

    belfem::DynamicBitset tCopy( tOriginal );

    EXPECT_EQ( tCopy.size(), tOriginal.size() );
    EXPECT_EQ( tCopy.count(), tOriginal.count() );
    EXPECT_TRUE( tCopy.test( 0 ) );

    // modify copy, verify original unchanged
    tCopy.reset( 0 );
    EXPECT_TRUE( tOriginal.test( 0 ) );
}

TEST( DynamicBitset, CopyPreservesLockState )
{
    belfem::DynamicBitset tOriginal( 64 );
    tOriginal.set( 0 );
    tOriginal.set( 63 );
    tOriginal.lock();

    belfem::DynamicBitset tCopy( tOriginal );

    EXPECT_TRUE( tCopy.is_locked() );
    EXPECT_EQ( tCopy.hash(), tOriginal.hash() );
}

TEST_P( DynamicBitset, MoveConstructorTransfers )
{
    belfem::index_t tN = GetParam();
    belfem::DynamicBitset tSource( tN );
    tSource.set( 0 );

    const uint64_t * tOldData = data_ptr( tSource );

    belfem::DynamicBitset tTarget( std::move( tSource ) );

    EXPECT_EQ( tTarget.size(), tN );
    EXPECT_TRUE( tTarget.test( 0 ) );
    EXPECT_EQ( data_ptr( tTarget ), tOldData );

    // the summary bitmaps must have been stolen along with the data, not left
    // behind pointing at the source
    EXPECT_TRUE( tTarget.summaries_are_tight() );
    expect_where_matches_reference( tTarget );

    // source is nullified (DynamicBitset guarantees this, unlike STL wrappers)
    EXPECT_EQ( tSource.size(), 0u );
    EXPECT_EQ( data_ptr( tSource ), nullptr );
    EXPECT_TRUE( tSource.summaries_are_tight() );
}

TEST( DynamicBitset, MoveConstructorPreservesHash )
{
    belfem::DynamicBitset tSource( 128 );
    tSource.set( 0 );
    tSource.set( 100 );
    tSource.lock();
    size_t tHash = tSource.hash();

    belfem::DynamicBitset tTarget( std::move( tSource ) );

    EXPECT_TRUE( tTarget.is_locked() );
    EXPECT_EQ( tTarget.hash(), tHash );
}

TEST( DynamicBitset, CopyAssignmentDifferentSizes )
{
    belfem::DynamicBitset tSmall( 32 );
    belfem::DynamicBitset tLarge( 128 );
    tLarge.set( 100 );

    tSmall = tLarge;

    EXPECT_EQ( tSmall.size(), 128u );
    EXPECT_EQ( tSmall.memory(), tLarge.memory() );
    EXPECT_TRUE( tSmall.test( 100 ) );
}

TEST( DynamicBitset, CopyAssignmentSameSize )
{
    belfem::DynamicBitset tA( 64 );
    belfem::DynamicBitset tB( 64 );
    tB.set( 10 );

    belfem::index_t tMemBefore = tA.memory();

    tA = tB;

    EXPECT_EQ( tA.memory(), tMemBefore );
    EXPECT_TRUE( tA.test( 10 ) );
}

TEST( DynamicBitset, MoveAssignment )
{
    belfem::DynamicBitset tSource( 128 );
    tSource.set( 50 );

    const uint64_t * tOldData = data_ptr( tSource );

    belfem::DynamicBitset tTarget( 64 );
    tTarget = std::move( tSource );

    EXPECT_EQ( tTarget.size(), 128u );
    EXPECT_TRUE( tTarget.test( 50 ) );
    EXPECT_EQ( data_ptr( tTarget ), tOldData );

    // the summaries must follow the data through the move assignment
    EXPECT_TRUE( tTarget.summaries_are_tight() );
    expect_where_matches_reference( tTarget );

    EXPECT_EQ( tSource.size(), 0u );
    EXPECT_EQ( data_ptr( tSource ), nullptr );
    EXPECT_TRUE( tSource.summaries_are_tight() );
}

TEST( DynamicBitset, SelfCopyAssignment )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set( 10 );

    tBs = tBs;

    EXPECT_EQ( tBs.size(), 64u );
    EXPECT_TRUE( tBs.test( 10 ) );
}

TEST( DynamicBitset, SelfMoveAssignment )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set( 10 );

#pragma GCC diagnostic push
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif
    tBs = std::move( tBs );
#pragma GCC diagnostic pop

    EXPECT_EQ( tBs.size(), 64u );
    EXPECT_TRUE( tBs.test( 10 ) );
}

// =============================================================================
// §2.2  Bit Manipulation  [semantic]
// =============================================================================

TEST_P( DynamicBitset, SetAndTest )
{
    belfem::index_t tN = GetParam();
    belfem::DynamicBitset tBs( tN );

    for( belfem::index_t i = 0; i < tN; ++i )
    {
        tBs.set( i );
        EXPECT_TRUE( tBs.test( i ) );
    }
}

TEST_P( DynamicBitset, ResetSingleBit )
{
    belfem::index_t tN = GetParam();
    belfem::DynamicBitset tBs( tN );

    tBs.set( 0 );
    EXPECT_TRUE( tBs.test( 0 ) );

    tBs.reset( 0 );
    EXPECT_FALSE( tBs.test( 0 ) );
}

TEST_P( DynamicBitset, ResetAll )
{
    belfem::index_t tN = GetParam();
    belfem::DynamicBitset tBs( tN );

    for( belfem::index_t i = 0; i < tN; ++i )
    {
        tBs.set( i );
    }
    EXPECT_EQ( tBs.count(), tN );

    tBs.reset();
    EXPECT_EQ( tBs.count(), 0u );
}

TEST( DynamicBitset, ResetAllUnlocks )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set( 0 );
    tBs.lock();
    EXPECT_TRUE( tBs.is_locked() );

    tBs.reset();

    EXPECT_FALSE( tBs.is_locked() );
    EXPECT_EQ( tBs.count(), 0u );
}

TEST_P( DynamicBitset, FlipSingleBit )
{
    belfem::DynamicBitset tBs( GetParam() );

    tBs.flip( 0 );
    EXPECT_TRUE( tBs.test( 0 ) );

    tBs.flip( 0 );
    EXPECT_FALSE( tBs.test( 0 ) );
}

TEST( DynamicBitset, FlipAll )
{
    belfem::DynamicBitset tBs( 128 );
    for( belfem::index_t i = 0; i < 128; i += 2 )
    {
        tBs.set( i );
    }
    EXPECT_EQ( tBs.count(), 64u );

    tBs.flip();

    EXPECT_EQ( tBs.count(), 64u );
    EXPECT_FALSE( tBs.test( 0 ) );
    EXPECT_TRUE( tBs.test( 1 ) );
}

TEST_P( DynamicBitset, SetWithBoolValue )
{
    belfem::DynamicBitset tBs( GetParam() );

    tBs.set( 0, true );
    EXPECT_TRUE( tBs.test( 0 ) );

    tBs.set( 0, false );
    EXPECT_FALSE( tBs.test( 0 ) );
}

TEST( DynamicBitset, BlockBoundaryBits )
{
    belfem::DynamicBitset tBs( 128 );

    belfem::index_t tBoundaryBits[] = { 62, 63, 64, 65 };
    for( auto tBit : tBoundaryBits )
    {
        tBs.set( tBit );
        EXPECT_TRUE( tBs.test( tBit ) );
    }

    EXPECT_EQ( tBs.count(), 4u );

    tBs.reset( 63 );
    EXPECT_FALSE( tBs.test( 63 ) );
    EXPECT_TRUE( tBs.test( 62 ) );
    EXPECT_TRUE( tBs.test( 64 ) );
}

TEST( DynamicBitset, LastBitInOddSize )
{
    belfem::DynamicBitset tBs( 65 );
    tBs.set( 64 );

    EXPECT_TRUE( tBs.test( 64 ) );
    EXPECT_EQ( tBs.count(), 1u );
}

// =============================================================================
// §2.3  Counting  [semantic]
// =============================================================================

TEST_P( DynamicBitset, CountAllZeros )
{
    belfem::DynamicBitset tBs( GetParam() );
    EXPECT_EQ( tBs.count(), 0u );
}

TEST_P( DynamicBitset, CountAllOnes )
{
    belfem::index_t tN = GetParam();
    belfem::DynamicBitset tBs( tN );

    for( belfem::index_t i = 0; i < tN; ++i )
    {
        tBs.set( i );
    }
    EXPECT_EQ( tBs.count(), tN );
}

TEST_P( DynamicBitset, CountSingleBit )
{
    belfem::DynamicBitset tBs( GetParam() );
    tBs.set( 0 );
    EXPECT_EQ( tBs.count(), 1u );
}

TEST( DynamicBitset, CountSparse )
{
    belfem::DynamicBitset tBs( 128 );
    tBs.set( 0 );
    tBs.set( 32 );
    tBs.set( 64 );
    EXPECT_EQ( tBs.count(), 3u );
}

TEST( DynamicBitset, CountIgnoresTrailingBits )
{
    belfem::DynamicBitset tBs( 65 );

    for( belfem::index_t i = 0; i < 65; ++i )
    {
        tBs.set( i );
    }
    EXPECT_EQ( tBs.count(), 65u );

    // This test used to reach through the mutable data() to set the trailing
    // bits of the last block by hand. That accessor is now private, because a
    // write through it bypasses the summary bitmaps. The invariant it was
    // probing is asserted directly instead: no public mutator can leave a bit
    // set at or beyond size(), which is what lets where() skip its per-bit
    // bounds check.
    tBs.flip();
    EXPECT_EQ( tBs.count(), 0u );
    EXPECT_TRUE( tBs.summaries_are_tight() );

    tBs.flip();
    EXPECT_EQ( tBs.count(), 65u );
    EXPECT_TRUE( tBs.summaries_are_tight() );

    // a full flip on an all-zero bitset must set exactly size() bits, never
    // the padding of the last block
    belfem::DynamicBitset tPadded( 65 );
    tPadded.flip();
    EXPECT_EQ( tPadded.count(), 65u );
    expect_where_matches_reference( tPadded );
}

// =============================================================================
// §2.4  Lock / Unlock / Hash  [semantic]
// =============================================================================

TEST( DynamicBitset, LockSetsLockedState )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set( 0 );
    tBs.lock();
    EXPECT_TRUE( tBs.is_locked() );
}

TEST( DynamicBitset, UnlockClearsLockedState )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set( 0 );
    tBs.lock();
    tBs.unlock();
    EXPECT_FALSE( tBs.is_locked() );
}

TEST( DynamicBitset, HashDeterministic )
{
    belfem::DynamicBitset tBs( 128 );
    tBs.set( 10 );
    tBs.set( 100 );

    tBs.lock();
    size_t tHash1 = tBs.hash();

    tBs.unlock();
    tBs.lock();
    size_t tHash2 = tBs.hash();

    EXPECT_EQ( tHash1, tHash2 );
}

TEST( DynamicBitset, EqualBitsetsEqualHash )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );

    tA.set( 10 );
    tA.set( 100 );
    tB.set( 10 );
    tB.set( 100 );

    tA.lock();
    tB.lock();

    EXPECT_EQ( tA.hash(), tB.hash() );
}

TEST( DynamicBitset, DifferentBitsetsDifferentHash )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );

    tA.set( 10 );
    tB.set( 11 );

    tA.lock();
    tB.lock();

    EXPECT_NE( tA.hash(), tB.hash() );
}

TEST( DynamicBitset, IndexSurvivesLock )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set_index( 42 );
    tBs.lock();
    EXPECT_EQ( tBs.index(), 42u );
}

// =============================================================================
// §2.5  Lock / Unlock / Hash  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( DynamicBitsetDebug, SetOnLockedThrows )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.lock();
    EXPECT_THROW( tBs.set( 0 ), std::runtime_error );
}

TEST( DynamicBitsetDebug, ResetBitOnLockedThrows )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set( 0 );
    tBs.lock();
    EXPECT_THROW( tBs.reset( 0 ), std::runtime_error );
}

TEST( DynamicBitsetDebug, FlipOnLockedThrows )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.lock();
    EXPECT_THROW( tBs.flip( 0 ), std::runtime_error );
}

TEST( DynamicBitsetDebug, FlipAllOnLockedThrows )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.lock();
    EXPECT_THROW( tBs.flip(), std::runtime_error );
}

TEST( DynamicBitsetDebug, HashOnUnlockedThrows )
{
    belfem::DynamicBitset tBs( 64 );
    EXPECT_THROW( tBs.hash(), std::runtime_error );
}

TEST( DynamicBitsetDebug, CompoundOrOnLockedThrows )
{
    belfem::DynamicBitset tA( 64 );
    belfem::DynamicBitset tB( 64 );
    tA.lock();
    EXPECT_THROW( tA |= tB, std::runtime_error );
}

TEST( DynamicBitsetDebug, CompoundXorOnLockedThrows )
{
    belfem::DynamicBitset tA( 64 );
    belfem::DynamicBitset tB( 64 );
    tA.lock();
    EXPECT_THROW( tA ^= tB, std::runtime_error );
}

TEST( DynamicBitsetDebug, CompoundAndOnLockedThrows )
{
    belfem::DynamicBitset tA( 64 );
    belfem::DynamicBitset tB( 64 );
    tA.lock();
    EXPECT_THROW( tA &= tB, std::runtime_error );
}

TEST( DynamicBitsetDebug, OutOfBoundsSetThrows )
{
    belfem::DynamicBitset tBs( 64 );
    EXPECT_THROW( tBs.set( 64 ), std::runtime_error );
}

TEST( DynamicBitsetDebug, OutOfBoundsResetThrows )
{
    belfem::DynamicBitset tBs( 64 );
    EXPECT_THROW( tBs.reset( 64 ), std::runtime_error );
}

TEST( DynamicBitsetDebug, OutOfBoundsFlipThrows )
{
    belfem::DynamicBitset tBs( 64 );
    EXPECT_THROW( tBs.flip( 64 ), std::runtime_error );
}

TEST( DynamicBitsetDebug, OutOfBoundsTestThrows )
{
    belfem::DynamicBitset tBs( 64 );
    EXPECT_THROW( tBs.test( 64 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §2.6  Bitwise Operators  [semantic]
// =============================================================================

TEST( DynamicBitset, OrOperator )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );

    tA.set( 0 );
    tA.set( 10 );
    tB.set( 10 );
    tB.set( 100 );

    belfem::DynamicBitset tResult = tA | tB;

    EXPECT_TRUE( tResult.test( 0 ) );
    EXPECT_TRUE( tResult.test( 10 ) );
    EXPECT_TRUE( tResult.test( 100 ) );
    EXPECT_EQ( tResult.count(), 3u );
}

TEST( DynamicBitset, OrAssignOperator )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );

    tA.set( 0 );
    tB.set( 100 );

    tA |= tB;

    EXPECT_TRUE( tA.test( 0 ) );
    EXPECT_TRUE( tA.test( 100 ) );
    EXPECT_EQ( tA.count(), 2u );
}

TEST( DynamicBitset, AndOperator )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );

    tA.set( 0 );
    tA.set( 10 );
    tB.set( 10 );
    tB.set( 100 );

    belfem::DynamicBitset tResult = tA & tB;

    EXPECT_FALSE( tResult.test( 0 ) );
    EXPECT_TRUE( tResult.test( 10 ) );
    EXPECT_FALSE( tResult.test( 100 ) );
    EXPECT_EQ( tResult.count(), 1u );
}

TEST( DynamicBitset, AndAssignOperator )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );

    tA.set( 0 );
    tA.set( 10 );
    tB.set( 10 );

    tA &= tB;

    EXPECT_FALSE( tA.test( 0 ) );
    EXPECT_TRUE( tA.test( 10 ) );
    EXPECT_EQ( tA.count(), 1u );
}

TEST( DynamicBitset, XorOperator )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );

    tA.set( 0 );
    tA.set( 10 );
    tB.set( 10 );
    tB.set( 100 );

    belfem::DynamicBitset tResult = tA ^ tB;

    EXPECT_TRUE( tResult.test( 0 ) );
    EXPECT_FALSE( tResult.test( 10 ) );
    EXPECT_TRUE( tResult.test( 100 ) );
    EXPECT_EQ( tResult.count(), 2u );
}

TEST( DynamicBitset, XorAssignOperator )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );

    tA.set( 0 );
    tA.set( 10 );
    tB.set( 10 );

    tA ^= tB;

    EXPECT_TRUE( tA.test( 0 ) );
    EXPECT_FALSE( tA.test( 10 ) );
    EXPECT_EQ( tA.count(), 1u );
}

TEST( DynamicBitset, OperatorWithSelf )
{
    belfem::DynamicBitset tA( 128 );
    tA.set( 0 );
    tA.set( 50 );
    tA.set( 100 );

    belfem::DynamicBitset tOr = tA | tA;
    EXPECT_EQ( tOr.count(), 3u );

    belfem::DynamicBitset tAnd = tA & tA;
    EXPECT_EQ( tAnd.count(), 3u );

    belfem::DynamicBitset tXor = tA ^ tA;
    EXPECT_EQ( tXor.count(), 0u );
}

TEST( DynamicBitset, OperatorWithAllZeros )
{
    belfem::DynamicBitset tA( 128 );
    tA.set( 10 );
    tA.set( 100 );

    belfem::DynamicBitset tZeros( 128 );

    belfem::DynamicBitset tOr = tA | tZeros;
    EXPECT_EQ( tOr.count(), 2u );

    belfem::DynamicBitset tAnd = tA & tZeros;
    EXPECT_EQ( tAnd.count(), 0u );
}

TEST( DynamicBitset, OperatorWithAllOnes )
{
    belfem::DynamicBitset tA( 64 );
    tA.set( 10 );
    tA.set( 50 );

    belfem::DynamicBitset tOnes( 64 );
    for( belfem::index_t i = 0; i < 64; ++i )
    {
        tOnes.set( i );
    }

    belfem::DynamicBitset tAnd = tA & tOnes;
    EXPECT_EQ( tAnd.count(), 2u );

    belfem::DynamicBitset tOr = tA | tOnes;
    EXPECT_EQ( tOr.count(), 64u );
}

// Equality / Inequality: operator== requires both operands locked.

TEST( DynamicBitset, EqualityOnIdenticalBitsets )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );
    tA.set( 10 );
    tA.set( 100 );
    tB.set( 10 );
    tB.set( 100 );

    tA.lock();
    tB.lock();

    EXPECT_TRUE( tA == tB );
    EXPECT_FALSE( tA != tB );
}

TEST( DynamicBitset, InequalityOnDifferentBitsets )
{
    belfem::DynamicBitset tA( 128 );
    belfem::DynamicBitset tB( 128 );
    tA.set( 10 );
    tB.set( 11 );

    tA.lock();
    tB.lock();

    EXPECT_TRUE( tA != tB );
    EXPECT_FALSE( tA == tB );
}

// =============================================================================
// §2.7  Bitwise Operators  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( DynamicBitsetDebug, OrMismatchedSizeThrows )
{
    belfem::DynamicBitset tA( 32 );
    belfem::DynamicBitset tB( 64 );
    EXPECT_THROW( tA | tB, std::runtime_error );
}

TEST( DynamicBitsetDebug, AndMismatchedSizeThrows )
{
    belfem::DynamicBitset tA( 32 );
    belfem::DynamicBitset tB( 64 );
    EXPECT_THROW( tA & tB, std::runtime_error );
}

TEST( DynamicBitsetDebug, XorMismatchedSizeThrows )
{
    belfem::DynamicBitset tA( 32 );
    belfem::DynamicBitset tB( 64 );
    EXPECT_THROW( tA ^ tB, std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §2.8  Conversions  [semantic]
// =============================================================================

TEST( DynamicBitset, ToStringBinaryRepresentation )
{
    belfem::DynamicBitset tBs( 4 );
    tBs.set( 0 );
    tBs.set( 2 );

    belfem::string tStr = tBs.to_string();

    EXPECT_EQ( tStr.size(), 4u );
    EXPECT_EQ( tStr, "0101" );
}

TEST( DynamicBitset, ToHexRoundTrip )
{
    belfem::DynamicBitset tBs( 128 );
    tBs.set( 0 );
    tBs.set( 7 );
    tBs.set( 64 );
    tBs.set( 127 );

    belfem::string tHex = tBs.to_hex();

    belfem::DynamicBitset tBs2( 128 );
    tBs2.set_from_hex( tHex );

    EXPECT_TRUE( tBs2.test( 0 ) );
    EXPECT_TRUE( tBs2.test( 7 ) );
    EXPECT_TRUE( tBs2.test( 64 ) );
    EXPECT_TRUE( tBs2.test( 127 ) );
    EXPECT_EQ( tBs2.count(), 4u );
}

TEST( DynamicBitset, ToHexEmptyBitset )
{
    // Codex finding: tests the zero-size branch in to_hex()
    belfem::DynamicBitset tBs( 0 );
    belfem::string tHex = tBs.to_hex();
    EXPECT_TRUE( tHex.empty() );
}

TEST( DynamicBitset, SetFromHexKnownPattern )
{
    // "A5" = 1010 0101 in binary (MSB first)
    // LSB-first bit ordering: bit 0=1, 1=0, 2=1, 3=0, 4=0, 5=1, 6=0, 7=1
    belfem::DynamicBitset tBs( 8 );
    tBs.set_from_hex( "A5" );

    EXPECT_TRUE( tBs.test( 0 ) );
    EXPECT_FALSE( tBs.test( 1 ) );
    EXPECT_TRUE( tBs.test( 2 ) );
    EXPECT_FALSE( tBs.test( 3 ) );
    EXPECT_FALSE( tBs.test( 4 ) );
    EXPECT_TRUE( tBs.test( 5 ) );
    EXPECT_FALSE( tBs.test( 6 ) );
    EXPECT_TRUE( tBs.test( 7 ) );
    EXPECT_EQ( tBs.count(), 4u );
}

TEST( DynamicBitset, ToIntPartial )
{
    belfem::DynamicBitset tBs( 8 );
    tBs.set( 0 );
    tBs.set( 2 );
    tBs.set( 4 );
    // 1 + 4 + 16 = 21

    EXPECT_EQ( tBs.to_int(), 21u );
}

TEST( DynamicBitset, ToIntFull )
{
    belfem::index_t tBits = sizeof( belfem::index_t ) * 8;
    belfem::DynamicBitset tBs( tBits );
    tBs.set( 0 );
    tBs.set( 1 );

    EXPECT_EQ( tBs.to_int(), 3u );
}

TEST( DynamicBitset, ToIntAllZeros )
{
    belfem::DynamicBitset tBs( 16 );
    EXPECT_EQ( tBs.to_int(), 0u );
}

TEST( DynamicBitset, ToRawString )
{
    belfem::DynamicBitset tBs( 65 );
    tBs.set( 0 );
    tBs.set( 64 );

    belfem::string tRaw = tBs.to_raw_string();

    ASSERT_EQ( tRaw.size(), 2 * sizeof( uint64_t ) );

    uint64_t tBlock0 = 0;
    uint64_t tBlock1 = 0;
    std::memcpy( &tBlock0, tRaw.data(), sizeof( uint64_t ) );
    std::memcpy( &tBlock1, tRaw.data() + sizeof( uint64_t ), sizeof( uint64_t ) );

    EXPECT_EQ( tBlock0, uint64_t( 1 ) );
    EXPECT_EQ( tBlock1, uint64_t( 1 ) );
}

// =============================================================================
// §2.9  Conversions  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( DynamicBitsetDebug, ToIntTooLargeThrows )
{
    belfem::DynamicBitset tBs( 128 );
    EXPECT_THROW( tBs.to_int(), std::runtime_error );
}

TEST( DynamicBitsetDebug, SetFromHexInvalidCharacterThrows )
{
    belfem::DynamicBitset tBs( 8 );
    EXPECT_THROW( tBs.set_from_hex( "G1" ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §2.10  Where (Index Extraction)  [semantic]
// =============================================================================

TEST( DynamicBitset, WhereEmptyBitset )
{
    belfem::DynamicBitset tBs( 64 );

    belfem::Cell< belfem::index_t > tBits{ 99, 100 };

    tBs.where( tBits );
    EXPECT_TRUE( tBits.empty() );
}

TEST( DynamicBitset, WhereAllSet )
{
    belfem::index_t tN = 65;
    belfem::DynamicBitset tBs( tN );
    for( belfem::index_t i = 0; i < tN; ++i )
    {
        tBs.set( i );
    }

    belfem::Cell< belfem::index_t > tBits;
    tBs.where( tBits, false );

    EXPECT_EQ( tBits.size(), tN );

    for( belfem::index_t i = 0; i < tN; ++i )
    {
        EXPECT_EQ( tBits( i ), i );
    }
}

TEST( DynamicBitset, WhereSparseMode )
{
    belfem::DynamicBitset tBs( 128 );
    tBs.set( 0 );
    tBs.set( 63 );
    tBs.set( 64 );

    belfem::Cell< belfem::index_t > tBits;
    tBs.where( tBits, true );

    ASSERT_EQ( tBits.size(), 3u );
    EXPECT_EQ( tBits( 0 ), 0u );
    EXPECT_EQ( tBits( 1 ), 63u );
    EXPECT_EQ( tBits( 2 ), 64u );
}

TEST( DynamicBitset, WhereDenseMode )
{
    belfem::DynamicBitset tBs( 128 );
    tBs.set( 0 );
    tBs.set( 63 );
    tBs.set( 64 );

    belfem::Cell< belfem::index_t > tBits;
    tBs.where( tBits, false );

    ASSERT_EQ( tBits.size(), 3u );
    EXPECT_EQ( tBits( 0 ), 0u );
    EXPECT_EQ( tBits( 1 ), 63u );
    EXPECT_EQ( tBits( 2 ), 64u );
}

TEST( DynamicBitset, WhereSparseVsDenseConsistency )
{
    belfem::DynamicBitset tBs( 128 );
    tBs.set( 5 );
    tBs.set( 63 );
    tBs.set( 64 );
    tBs.set( 127 );

    belfem::Cell< belfem::index_t > tSparse;
    belfem::Cell< belfem::index_t > tDense;

    tBs.where( tSparse, true );
    tBs.where( tDense, false );

    ASSERT_EQ( tSparse.size(), tDense.size() );

    for( size_t i = 0; i < tSparse.size(); ++i )
    {
        EXPECT_EQ( tSparse( i ), tDense( i ) );
    }
}

TEST( DynamicBitset, WhereZeroSizeBitset )
{
    belfem::DynamicBitset tBs( 0 );

    belfem::Cell< belfem::index_t > tBits{ 1, 2, 3 };
    tBs.where( tBits );
    EXPECT_TRUE( tBits.empty() );
}

// =============================================================================
// §2.11  Index  [semantic]
// =============================================================================

TEST( DynamicBitset, DefaultIndexIsNoIndex )
{
    belfem::DynamicBitset tBs( 64 );
    EXPECT_EQ( tBs.index(), belfem::gNoIndex );
}

TEST( DynamicBitset, SetIndexAndRetrieve )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set_index( 42 );
    EXPECT_EQ( tBs.index(), 42u );
}

TEST( DynamicBitset, IndexSurvivesCopy )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set_index( 99 );

    belfem::DynamicBitset tCopy( tBs );
    EXPECT_EQ( tCopy.index(), 99u );
}

TEST( DynamicBitset, IndexSurvivesMove )
{
    belfem::DynamicBitset tBs( 64 );
    tBs.set_index( 77 );

    belfem::DynamicBitset tTarget( std::move( tBs ) );

    EXPECT_EQ( tTarget.index(), 77u );
    EXPECT_EQ( tBs.index(), belfem::gNoIndex );
}

// =============================================================================
// §2.12  Summary Bitmaps  [semantic]
// =============================================================================
//
// where() and reset() walk a two-level summary bitmap instead of scanning the
// whole data array. The tests below pin the two properties that make that
// safe: the summaries stay tight under every mutator, and where() still agrees
// with a naive per-bit reference.
//
// Sizes here deliberately exceed the parameterized set at the top of this file:
// with 64 bits per data word and 64 data words per level-1 word, a bitset must
// pass 4096 bits before it uses more than one level-1 word, and 262144 before
// it uses more than one level-2 word. Every size in the Sizes/ suite is below
// 4096, so none of them exercise the walk at all.

TEST( DynamicBitset, SummariesTightAcrossAllMutators )
{
    const belfem::index_t tSizes[] =
        { 1, 63, 64, 65, 4095, 4096, 4097, 262144, 262145, 300000 };

    std::mt19937_64 tRandom( 20260818 );

    for( belfem::index_t tN : tSizes )
    {
        belfem::DynamicBitset tBs( tN );
        EXPECT_TRUE( tBs.summaries_are_tight() );

        std::uniform_int_distribution< belfem::index_t > tPos( 0, tN - 1 );

        for( int tRound = 0; tRound < 300; ++tRound )
        {
            switch( tRandom() % 6 )
            {
                case 0 : tBs.set( tPos( tRandom ) );                        break;
                case 1 : tBs.reset( tPos( tRandom ) );                      break;
                case 2 : tBs.set( tPos( tRandom ), ( tRandom() & 1 ) != 0 );break;
                case 3 : tBs.flip( tPos( tRandom ) );                       break;
                case 4 : tBs.reset();                                       break;
                default: tBs.set( tPos( tRandom ) );                        break;
            }

            ASSERT_TRUE( tBs.summaries_are_tight() )
                << "size " << tN << ", round " << tRound;
        }

        expect_where_matches_reference( tBs );
    }
}

TEST( DynamicBitset, SummariesTightAcrossBitwiseOperators )
{
    const belfem::index_t tN = 5000;   // spans several level-1 words

    std::mt19937_64 tRandom( 4711 );
    std::uniform_int_distribution< belfem::index_t > tPos( 0, tN - 1 );

    belfem::DynamicBitset tA( tN );
    belfem::DynamicBitset tB( tN );

    for( int k = 0; k < 200; ++k )
    {
        tA.set( tPos( tRandom ) );
        tB.set( tPos( tRandom ) );
    }

    // XOR and AND can zero a word that was nonzero in both operands, so their
    // summaries cannot be derived from the inputs and must be rebuilt
    expect_where_matches_reference( tA | tB );
    expect_where_matches_reference( tA ^ tB );
    expect_where_matches_reference( tA & tB );

    belfem::DynamicBitset tOr( tA );  tOr  |= tB;  expect_where_matches_reference( tOr );
    belfem::DynamicBitset tXor( tA ); tXor ^= tB;  expect_where_matches_reference( tXor );
    belfem::DynamicBitset tAnd( tA ); tAnd &= tB;  expect_where_matches_reference( tAnd );

    // self-XOR empties the bitset: every summary bit must come back down
    belfem::DynamicBitset tSelf( tA );
    tSelf ^= tA;
    EXPECT_EQ( tSelf.count(), 0u );
    EXPECT_TRUE( tSelf.summaries_are_tight() );
}

TEST( DynamicBitset, SummariesSurviveCopyAndAssignment )
{
    const belfem::index_t tN = 300000;

    belfem::DynamicBitset tBs( tN );
    for( belfem::index_t k = 0; k < tN; k += 997 )
    {
        tBs.set( k );
    }

    belfem::DynamicBitset tCopy( tBs );
    expect_where_matches_reference( tCopy );

    // assignment from a differently sized bitset must reallocate the summaries
    belfem::DynamicBitset tSmall( 64 );
    tSmall = tBs;
    expect_where_matches_reference( tSmall );

    // and assignment onto an equally sized one must overwrite them
    belfem::DynamicBitset tSame( tN );
    tSame.set( 12345 );
    tSame = tBs;
    expect_where_matches_reference( tSame );
}

TEST( DynamicBitset, FullResetCyclingMatchesReference )
{
    // the shape used by the sparsity-pattern builders: reset the whole
    // workspace, mark a handful of bits, extract them, repeat
    const belfem::index_t tN = 1u << 21;

    belfem::DynamicBitset tBs( tN );
    belfem::Cell< belfem::index_t > tIndices;

    std::mt19937_64 tRandom( 90210 );
    std::uniform_int_distribution< belfem::index_t > tPos( 0, tN - 1 );

    for( int tRound = 0; tRound < 200; ++tRound )
    {
        tBs.reset();

        std::set< belfem::index_t > tMarked;
        for( int k = 0; k < 64; ++k )
        {
            const belfem::index_t tBit = tPos( tRandom );
            tBs.set( tBit );
            tMarked.insert( tBit );
        }

        // the builders drop the self-link before extracting
        const belfem::index_t tDropped = *tMarked.begin();
        tBs.reset( tDropped );
        tMarked.erase( tDropped );

        tBs.where( tIndices );

        ASSERT_EQ( tIndices.size(), tMarked.size() ) << "round " << tRound;

        belfem::index_t k = 0;
        for( belfem::index_t tBit : tMarked )
        {
            EXPECT_EQ( tIndices( k++ ), tBit );
        }

        ASSERT_TRUE( tBs.summaries_are_tight() ) << "round " << tRound;
    }

    tBs.reset();
    EXPECT_EQ( tBs.count(), 0u );
    tBs.where( tIndices );
    EXPECT_EQ( tIndices.size(), 0u );
}

TEST( DynamicBitset, PerBitClearingKeepsSummariesTight )
{
    // Regression tripwire. A summary design that lets reset( aPos ) leave a
    // stale bit behind is still *correct*, so no equivalence test catches it —
    // but a caller that clears bits one by one instead of calling reset() then
    // saturates the summaries and where() silently decays to a full scan. That
    // is exactly what the sparsity-pattern builders and the shortest-path
    // solver do, so it is pinned here with an exact predicate rather than a
    // wall-clock threshold.
    const belfem::index_t tN = 1u << 20;

    belfem::DynamicBitset tBs( tN );
    belfem::Cell< belfem::index_t > tIndices;

    std::mt19937_64 tRandom( 31337 );
    std::uniform_int_distribution< belfem::index_t > tPos( 0, tN - 1 );

    for( int tRound = 0; tRound < 300; ++tRound )
    {
        for( int k = 0; k < 32; ++k )
        {
            tBs.set( tPos( tRandom ) );
        }

        tBs.where( tIndices );

        const std::vector< belfem::index_t > tReference = reference_where( tBs );
        ASSERT_EQ( tIndices.size(), tReference.size() ) << "round " << tRound;

        for( belfem::index_t k = 0; k < tIndices.size(); ++k )
        {
            EXPECT_EQ( tIndices( k ), tReference[ k ] );
        }

        // clear only the bits we just listed — never a full reset()
        for( belfem::index_t k = 0; k < tIndices.size(); ++k )
        {
            tBs.reset( tIndices( k ) );
        }

        ASSERT_TRUE( tBs.summaries_are_tight() )
            << "summaries went stale in round " << tRound;
        ASSERT_EQ( tBs.count(), 0u ) << "round " << tRound;
    }
}

TEST( DynamicBitset, FlipInterleavedWithSetAndReset )
{
    // the shape used by the Maxwell factory: flip() / reset() / set() / where()
    // interleaved on a size whose last data word AND last summary words are
    // partial
    const belfem::index_t tN = 262145;   // one bit past a full level-2 word

    belfem::DynamicBitset tBs( tN );

    tBs.flip();
    EXPECT_EQ( tBs.count(), tN );
    EXPECT_TRUE( tBs.summaries_are_tight() );

    tBs.reset( 0 );
    tBs.reset( tN - 1 );
    EXPECT_EQ( tBs.count(), tN - 2 );
    EXPECT_TRUE( tBs.summaries_are_tight() );

    tBs.flip();
    EXPECT_EQ( tBs.count(), 2u );
    expect_where_matches_reference( tBs );

    tBs.reset();
    tBs.set( 262144 );
    tBs.flip( 262144 );
    EXPECT_EQ( tBs.count(), 0u );
    EXPECT_TRUE( tBs.summaries_are_tight() );

    // flipping a single bit ON in a word whose summary bit is clear must be
    // visible to where() — this is a correctness requirement, not tidiness
    tBs.flip( 200000 );

    belfem::Cell< belfem::index_t > tIndices;
    tBs.where( tIndices );

    ASSERT_EQ( tIndices.size(), 1u );
    EXPECT_EQ( tIndices( 0 ), 200000u );
    EXPECT_TRUE( tBs.summaries_are_tight() );
}

TEST( DynamicBitset, ZeroSizeBitsetSurvivesEveryOperation )
{
    belfem::DynamicBitset tBs( 0 );

    // flip() used to underflow its loop bound here and walk a null pointer
    tBs.flip();
    tBs.reset();

    belfem::Cell< belfem::index_t > tIndices;
    tBs.where( tIndices );

    EXPECT_EQ( tIndices.size(), 0u );
    EXPECT_EQ( tBs.count(), 0u );
    EXPECT_EQ( tBs.size(), 0u );
    EXPECT_TRUE( tBs.summaries_are_tight() );

    // to_int() used to dispatch to the partial conversion here, which reads
    // mData[ 0 ] — null for an empty bitset
    EXPECT_EQ( tBs.to_int(), 0u );

    EXPECT_EQ( tBs.to_hex(), belfem::string( "" ) );
    EXPECT_EQ( tBs.to_string(), belfem::string( "" ) );
}
