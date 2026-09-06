/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Bitset<N> container
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §4 (test matrix)
 *
 * Bitset<N> is a thin wrapper around std::bitset<N>.  All methods are
 * direct STL pass-throughs, so these are smoke tests per the strategy doc.
 * No BELFEM_ASSERT/BELFEM_ERROR in this class → no debug-only tests.
 */

#include <gtest/gtest.h>
#include <type_traits>
#include <bitset>

#include "typedefs.hpp"
#include "cl_Bitset.hpp"

// =============================================================================
// Typed test infrastructure
// =============================================================================
// Test plan §4: typed tests over N = 8, 64, 128, 256.
// GTest typed tests require types, so we use std::integral_constant
// to encode the non-type template parameter as a type.

template< typename T >
class Bitset : public ::testing::Test {};

using BitsetSizes = ::testing::Types<
    std::integral_constant< belfem::index_t, 8 >,
    std::integral_constant< belfem::index_t, 64 >,
    std::integral_constant< belfem::index_t, 128 >,
    std::integral_constant< belfem::index_t, 256 >
>;

TYPED_TEST_SUITE( Bitset, BitsetSizes );

// =============================================================================
// §4.1  Tests  [semantic]
// =============================================================================

TYPED_TEST( Bitset, DefaultConstructorAllZero )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;

    EXPECT_EQ( tBs.count(), 0u );
}

TYPED_TEST( Bitset, SizeEqualsN )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;

    EXPECT_EQ( tBs.size(), N );
}

TYPED_TEST( Bitset, SetAndTest )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;

    // test first, last, and a middle bit
    tBs.set( 0 );
    EXPECT_TRUE( tBs.test( 0 ) );

    tBs.set( N - 1 );
    EXPECT_TRUE( tBs.test( N - 1 ) );

    belfem::index_t tMid = N / 2;
    tBs.set( tMid );
    EXPECT_TRUE( tBs.test( tMid ) );

    EXPECT_EQ( tBs.count(), 3u );
}

TYPED_TEST( Bitset, ResetSingleBit )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;

    tBs.set( 0 );
    EXPECT_TRUE( tBs.test( 0 ) );

    tBs.reset( 0 );
    EXPECT_FALSE( tBs.test( 0 ) );
    EXPECT_EQ( tBs.count(), 0u );
}

TYPED_TEST( Bitset, ResetAll )
{
    // NOTE: the test plan says "verify this exists; if absent, skip."
    // Source inspection confirms reset() (no arg) exists at line 112.
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;

    // set several bits
    tBs.set( 0 );
    tBs.set( N / 2 );
    tBs.set( N - 1 );
    ASSERT_EQ( tBs.count(), 3u );   // verify setup (idea from ChatGPT)

    tBs.reset();
    EXPECT_EQ( tBs.count(), 0u );

    // spot-check that specific bits are cleared
    EXPECT_FALSE( tBs.test( 0 ) );
    EXPECT_FALSE( tBs.test( N / 2 ) );
    EXPECT_FALSE( tBs.test( N - 1 ) );
}

TYPED_TEST( Bitset, FlipBit )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;

    // 0 → 1
    tBs.flip( 0 );
    EXPECT_TRUE( tBs.test( 0 ) );

    // 1 → 0
    tBs.flip( 0 );
    EXPECT_FALSE( tBs.test( 0 ) );
}

TYPED_TEST( Bitset, CountMatchesManual )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;

    // set bits at indices 0, 2, 4 (up to N)
    belfem::index_t tExpected = 0;
    for( belfem::index_t i = 0; i < N; i += 2 )
    {
        tBs.set( i );
        ++tExpected;
    }

    EXPECT_EQ( tBs.count(), tExpected );
}

TYPED_TEST( Bitset, CopyConstructorDeepCopies )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tOriginal;
    tOriginal.set( 0 );
    tOriginal.set( N - 1 );

    belfem::Bitset< N > tCopy( tOriginal );

    EXPECT_EQ( tCopy.count(), 2u );
    EXPECT_TRUE( tCopy.test( 0 ) );
    EXPECT_TRUE( tCopy.test( N - 1 ) );

    // modify copy, verify original unchanged
    tCopy.reset( 0 );
    EXPECT_TRUE( tOriginal.test( 0 ) );
}

TYPED_TEST( Bitset, MoveConstructor )
{
    // std::bitset is stack-allocated, so move is effectively copy.
    // We test for API completeness, not observable move semantics.
    // Do NOT assert source is empty — it retains its data after move.
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tSource;
    tSource.set( 0 );
    tSource.set( N - 1 );

    belfem::Bitset< N > tTarget( std::move( tSource ) );

    EXPECT_EQ( tTarget.count(), 2u );
    EXPECT_TRUE( tTarget.test( 0 ) );
    EXPECT_TRUE( tTarget.test( N - 1 ) );
}

TYPED_TEST( Bitset, CopyAssignment )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tA;
    tA.set( 0 );
    tA.set( N / 2 );

    belfem::Bitset< N > tB;
    tB = tA;

    EXPECT_EQ( tB.count(), 2u );
    EXPECT_TRUE( tB.test( 0 ) );
    EXPECT_TRUE( tB.test( N / 2 ) );

    // verify independence
    tB.reset( 0 );
    EXPECT_TRUE( tA.test( 0 ) );
}

TYPED_TEST( Bitset, CopyAssignmentSelfSafe )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;
    tBs.set( 0 );

    tBs = tBs;

    EXPECT_EQ( tBs.count(), 1u );
    EXPECT_TRUE( tBs.test( 0 ) );
}

TYPED_TEST( Bitset, MoveAssignment )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tSource;
    tSource.set( 0 );
    tSource.set( N - 1 );

    belfem::Bitset< N > tTarget;
    tTarget = std::move( tSource );

    EXPECT_EQ( tTarget.count(), 2u );
    EXPECT_TRUE( tTarget.test( 0 ) );
    EXPECT_TRUE( tTarget.test( N - 1 ) );
}

TYPED_TEST( Bitset, MoveAssignmentSelfSafe )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;
    tBs.set( 0 );

#pragma GCC diagnostic push
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif
    tBs = std::move( tBs );
#pragma GCC diagnostic pop

    EXPECT_EQ( tBs.count(), 1u );
    EXPECT_TRUE( tBs.test( 0 ) );
}

TYPED_TEST( Bitset, EqualityOperator )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tA;
    belfem::Bitset< N > tB;

    // both empty → equal
    EXPECT_TRUE( tA == tB );

    // same bits set → equal
    tA.set( 0 );
    tA.set( N - 1 );
    tB.set( 0 );
    tB.set( N - 1 );
    EXPECT_TRUE( tA == tB );

    // different bits → not equal
    tB.set( 1 );
    EXPECT_FALSE( tA == tB );
}

TYPED_TEST( Bitset, InequalityOperator )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tA;
    belfem::Bitset< N > tB;

    tA.set( 0 );
    EXPECT_TRUE( tA != tB );

    tB.set( 0 );
    EXPECT_FALSE( tA != tB );
}

TYPED_TEST( Bitset, DataExposesStdBitset )
{
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;
    tBs.set( 0 );

    std::bitset< N > & tRef = tBs.data();

    // compile-time check: data() returns a reference, not a copy (idea from Grok)
    static_assert(
        std::is_same< decltype( tBs.data() ), std::bitset< N > & >::value,
        "data() must return reference to internal std::bitset<N>" );

    // verify the reference points to the real internal data
    EXPECT_TRUE( tRef.test( 0 ) );

    // modify through the reference, verify Bitset sees the change
    tRef.set( N - 1 );
    EXPECT_TRUE( tBs.test( N - 1 ) );
}

TYPED_TEST( Bitset, DataMutationMatchesWrapperMutation )
{
    // Reverse-direction aliasing test: write through wrapper, manipulate
    // through data() ref, verify coherence (idea from ChatGPT)
    constexpr belfem::index_t N = TypeParam::value;
    belfem::Bitset< N > tBs;

    tBs.set( 0 );
    tBs.data().flip( 0 );
    tBs.data().set( N / 2 );

    EXPECT_FALSE( tBs.test( 0 ) );
    EXPECT_TRUE( tBs.test( N / 2 ) );
    EXPECT_EQ( tBs.count(), 1u );
}
