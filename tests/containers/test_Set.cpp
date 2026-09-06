/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Set<Key> container
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §6 (test matrix)
 *
 * Set<Key> wraps std::unordered_set<Key>.  Most member functions are thin
 * STL pass-throughs (smoke tests), but the set algebra operators (|, &, -, ^),
 * subset/superset checks, and equality are BELFEM-authored and get full coverage.
 *
 * No BELFEM_ASSERT/BELFEM_ERROR in this class → no debug-only tests.
 */

#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <unordered_set>

#include "typedefs.hpp"
#include "cl_Set.hpp"

// =============================================================================
// §6.1  Tests  [semantic] — int keys
// =============================================================================

TEST( Set, DefaultConstructorEmpty )
{
    belfem::Set< int > tSet;
    EXPECT_EQ( tSet.size(), 0u );
    EXPECT_TRUE( tSet.empty() );
}

TEST( Set, InitializerListConstructor )
{
    belfem::Set< int > tSet{ 1, 2, 3 };
    EXPECT_EQ( tSet.size(), 3u );
    EXPECT_TRUE( tSet.contains( 1 ) );
    EXPECT_TRUE( tSet.contains( 2 ) );
    EXPECT_TRUE( tSet.contains( 3 ) );
}

TEST( Set, InitializerListDeduplicates )
{
    // duplicates in initializer list should be removed (idea from Grok)
    belfem::Set< int > tSet{ 5, 2, 8, 2 };
    EXPECT_EQ( tSet.size(), 3u );
    EXPECT_TRUE( tSet.contains( 5 ) );
    EXPECT_TRUE( tSet.contains( 2 ) );
    EXPECT_TRUE( tSet.contains( 8 ) );
}

TEST( Set, IteratorRangeConstructor )
{
    // constructor from iterator pair (idea from ChatGPT)
    std::vector< int > tData{ 10, 20, 30, 20 };

    belfem::Set< int > tSet( tData.begin(), tData.end() );

    EXPECT_EQ( tSet.size(), 3u );
    EXPECT_TRUE( tSet.contains( 10 ) );
    EXPECT_TRUE( tSet.contains( 20 ) );
    EXPECT_TRUE( tSet.contains( 30 ) );
}

TEST( Set, InsertNewElement )
{
    belfem::Set< int > tSet;
    auto tResult = tSet.insert( 42 );

    EXPECT_TRUE( tResult.second );     // inserted
    EXPECT_EQ( tSet.size(), 1u );
    EXPECT_TRUE( tSet.contains( 42 ) );
}

TEST( Set, InsertDuplicate )
{
    belfem::Set< int > tSet;
    tSet.insert( 42 );
    auto tResult = tSet.insert( 42 );

    EXPECT_FALSE( tResult.second );    // not inserted
    EXPECT_EQ( tSet.size(), 1u );
}

TEST( Set, InsertRange )
{
    // range-insert overload (idea from ChatGPT)
    std::vector< int > tData{ 10, 20, 30, 20 };
    belfem::Set< int > tSet;

    tSet.insert( tData.begin(), tData.end() );

    EXPECT_EQ( tSet.size(), 3u );
    EXPECT_TRUE( tSet.contains( 10 ) );
    EXPECT_TRUE( tSet.contains( 20 ) );
    EXPECT_TRUE( tSet.contains( 30 ) );
}

TEST( Set, Emplace )
{
    // emplace method (idea from ChatGPT)
    belfem::Set< int > tSet;
    auto tResult = tSet.emplace( 42 );

    EXPECT_TRUE( tResult.second );
    EXPECT_EQ( tSet.size(), 1u );
    EXPECT_TRUE( tSet.contains( 42 ) );
}

TEST( Set, ContainsAndKeyExists )
{
    belfem::Set< int > tSet{ 10, 20 };

    // both methods should agree
    EXPECT_TRUE( tSet.contains( 10 ) );
    EXPECT_TRUE( tSet.key_exists( 10 ) );

    EXPECT_FALSE( tSet.contains( 99 ) );
    EXPECT_FALSE( tSet.key_exists( 99 ) );
}

TEST( Set, Count )
{
    belfem::Set< int > tSet{ 10, 20 };

    EXPECT_EQ( tSet.count( 10 ), 1u );
    EXPECT_EQ( tSet.count( 99 ), 0u );
}

TEST( Set, EraseByKey )
{
    belfem::Set< int > tSet{ 10, 20, 30 };

    size_t tRemoved = tSet.erase( 20 );

    EXPECT_EQ( tRemoved, 1u );
    EXPECT_EQ( tSet.size(), 2u );
    EXPECT_FALSE( tSet.contains( 20 ) );
}

TEST( Set, EraseAbsentKey )
{
    belfem::Set< int > tSet{ 10, 20 };

    size_t tRemoved = tSet.erase( 99 );

    EXPECT_EQ( tRemoved, 0u );
    EXPECT_EQ( tSet.size(), 2u );
}

TEST( Set, EraseByIterator )
{
    belfem::Set< int > tSet{ 10, 20, 30 };

    auto tIt = tSet.find( 20 );
    ASSERT_NE( tIt, tSet.end() );

    tSet.erase( tIt );

    EXPECT_EQ( tSet.size(), 2u );
    EXPECT_FALSE( tSet.contains( 20 ) );
}

TEST( Set, EraseRange )
{
    // range-erase overload (idea from ChatGPT)
    belfem::Set< int > tSet{ 10, 20, 30, 40 };

    tSet.erase( tSet.begin(), tSet.end() );

    EXPECT_TRUE( tSet.empty() );
    EXPECT_EQ( tSet.size(), 0u );
}

TEST( Set, Clear )
{
    belfem::Set< int > tSet{ 1, 2, 3 };
    tSet.clear();

    EXPECT_EQ( tSet.size(), 0u );
    EXPECT_TRUE( tSet.empty() );
}

TEST( Set, Reserve )
{
    belfem::Set< int > tSet{ 1, 2 };
    tSet.reserve( 100 );

    // reserve does not change size
    EXPECT_EQ( tSet.size(), 2u );
}

TEST( Set, Swap )
{
    belfem::Set< int > tA{ 1, 2 };
    belfem::Set< int > tB{ 3, 4, 5 };

    tA.swap( tB );

    EXPECT_EQ( tA.size(), 3u );
    EXPECT_TRUE( tA.contains( 3 ) );
    EXPECT_EQ( tB.size(), 2u );
    EXPECT_TRUE( tB.contains( 1 ) );
}

TEST( Set, CopyConstructor )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tB( tA );

    EXPECT_EQ( tB.size(), 3u );
    EXPECT_TRUE( tB.contains( 1 ) );

    // verify independence
    tB.insert( 99 );
    EXPECT_FALSE( tA.contains( 99 ) );
}

TEST( Set, CopyAssignment )
{
    // separate from copy constructor (idea from ChatGPT)
    belfem::Set< int > tA{ 1, 2 };
    belfem::Set< int > tB{ 99 };

    tB = tA;
    tB.insert( 3 );

    EXPECT_TRUE( tA.contains( 1 ) );
    EXPECT_TRUE( tA.contains( 2 ) );
    EXPECT_FALSE( tA.contains( 3 ) );

    EXPECT_TRUE( tB.contains( 1 ) );
    EXPECT_TRUE( tB.contains( 2 ) );
    EXPECT_TRUE( tB.contains( 3 ) );
}

TEST( Set, MoveConstructor )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tB( std::move( tA ) );

    EXPECT_EQ( tB.size(), 3u );
    EXPECT_TRUE( tB.contains( 1 ) );

    // source remains valid (do NOT assert empty — STL-backed)
    EXPECT_NO_THROW(
    {
        ( void ) tA.size();
        tA.clear();
        tA.insert( 42 );
    });
}

TEST( Set, MoveAssignment )
{
    // separate from move constructor (idea from ChatGPT)
    belfem::Set< int > tA{ 1, 2 };
    belfem::Set< int > tB{ 99 };

    tB = std::move( tA );

    EXPECT_EQ( tB.size(), 2u );
    EXPECT_TRUE( tB.contains( 1 ) );
    EXPECT_TRUE( tB.contains( 2 ) );

    // source remains valid (Codex finding: was missing source validity check)
    EXPECT_NO_THROW(
    {
        ( void ) tA.size();
        tA.clear();
        tA.insert( 42 );
    });
}

TEST( Set, Iteration )
{
    belfem::Set< int > tSet{ 10, 20, 30 };

    size_t tCount = 0;
    bool tSaw10 = false;
    bool tSaw20 = false;
    bool tSaw30 = false;

    for( const auto & tKey : tSet )
    {
        ++tCount;
        if( tKey == 10 ) tSaw10 = true;
        if( tKey == 20 ) tSaw20 = true;
        if( tKey == 30 ) tSaw30 = true;
    }

    EXPECT_EQ( tCount, 3u );
    EXPECT_TRUE( tSaw10 );
    EXPECT_TRUE( tSaw20 );
    EXPECT_TRUE( tSaw30 );
}

TEST( Set, SetDataExposesUnderlying )
{
    belfem::Set< int > tSet{ 1, 2 };

    std::unordered_set< int > & tRef = tSet.set_data();

    // compile-time check (idea from Grok)
    static_assert(
        std::is_same< decltype( tSet.set_data() ), std::unordered_set< int > & >::value,
        "set_data() must return reference to internal unordered_set" );

    // insert through underlying, verify wrapper sees it
    tRef.insert( 99 );
    EXPECT_EQ( tSet.size(), 3u );
    EXPECT_TRUE( tSet.contains( 99 ) );
}

// =============================================================================
// §6.1  Tests  [semantic] — string keys
// =============================================================================

TEST( Set, InitializerListConstructorString )
{
    belfem::Set< belfem::string > tSet{ "alpha", "beta", "gamma" };
    EXPECT_EQ( tSet.size(), 3u );
    EXPECT_TRUE( tSet.contains( "alpha" ) );
}

TEST( Set, InsertMoveSemantics )
{
    belfem::Set< belfem::string > tSet;
    belfem::string tVal = "moveable";
    tSet.insert( std::move( tVal ) );

    EXPECT_EQ( tSet.size(), 1u );
    EXPECT_TRUE( tSet.contains( "moveable" ) );
    // tVal is in valid-but-unspecified state
}

TEST( Set, CopySemantics )
{
    belfem::Set< belfem::string > tA{ "x", "y" };
    belfem::Set< belfem::string > tB( tA );

    tB.insert( "z" );
    EXPECT_FALSE( tA.contains( "z" ) );
}

TEST( Set, EraseByKeyString )
{
    belfem::Set< belfem::string > tSet{ "a", "b", "c" };
    tSet.erase( "b" );

    EXPECT_EQ( tSet.size(), 2u );
    EXPECT_FALSE( tSet.contains( "b" ) );
    EXPECT_TRUE( tSet.contains( "a" ) );
    EXPECT_TRUE( tSet.contains( "c" ) );
}

// =============================================================================
// §6.2  Set Algebra  [semantic]
// =============================================================================

TEST( Set, UnionOperator )
{
    belfem::Set< int > tA{ 1, 2 };
    belfem::Set< int > tB{ 2, 3 };

    belfem::Set< int > tResult = tA | tB;

    EXPECT_EQ( tResult.size(), 3u );
    EXPECT_TRUE( tResult.contains( 1 ) );
    EXPECT_TRUE( tResult.contains( 2 ) );
    EXPECT_TRUE( tResult.contains( 3 ) );
}

TEST( Set, IntersectionOperator )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tB{ 2, 3, 4 };

    belfem::Set< int > tResult = tA & tB;

    EXPECT_EQ( tResult.size(), 2u );
    EXPECT_TRUE( tResult.contains( 2 ) );
    EXPECT_TRUE( tResult.contains( 3 ) );
}

TEST( Set, DifferenceOperator )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tB{ 2, 4 };

    belfem::Set< int > tResult = tA - tB;

    EXPECT_EQ( tResult.size(), 2u );
    EXPECT_TRUE( tResult.contains( 1 ) );
    EXPECT_TRUE( tResult.contains( 3 ) );
}

TEST( Set, SymmetricDifferenceOperator )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tB{ 2, 3, 4 };

    belfem::Set< int > tResult = tA ^ tB;

    EXPECT_EQ( tResult.size(), 2u );
    EXPECT_TRUE( tResult.contains( 1 ) );
    EXPECT_TRUE( tResult.contains( 4 ) );
}

TEST( Set, UnionWithEmpty )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tEmpty;

    belfem::Set< int > tResult = tA | tEmpty;

    // use operator== for identity check (idea from ChatGPT)
    EXPECT_TRUE( tResult == tA );
}

TEST( Set, IntersectionWithEmpty )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tEmpty;

    belfem::Set< int > tResult = tA & tEmpty;

    EXPECT_EQ( tResult.size(), 0u );
}

TEST( Set, UnionWithSelf )
{
    belfem::Set< int > tA{ 1, 2, 3 };

    belfem::Set< int > tResult = tA | tA;

    EXPECT_TRUE( tResult == tA );
}

TEST( Set, IntersectionWithSelf )
{
    belfem::Set< int > tA{ 1, 2, 3 };

    belfem::Set< int > tResult = tA & tA;

    EXPECT_TRUE( tResult == tA );
}

TEST( Set, DifferenceWithSelf )
{
    belfem::Set< int > tA{ 1, 2, 3 };

    belfem::Set< int > tResult = tA - tA;

    EXPECT_EQ( tResult.size(), 0u );
}

TEST( Set, SymmetricDifferenceWithSelf )
{
    belfem::Set< int > tA{ 1, 2, 3 };

    belfem::Set< int > tResult = tA ^ tA;

    EXPECT_EQ( tResult.size(), 0u );
}

TEST( Set, DisjointSetsIntersection )
{
    belfem::Set< int > tA{ 1, 2 };
    belfem::Set< int > tB{ 3, 4 };

    belfem::Set< int > tResult = tA & tB;

    EXPECT_EQ( tResult.size(), 0u );
}

TEST( Set, EqualityOperator )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tB{ 3, 1, 2 };   // same elements, different insertion order

    EXPECT_TRUE( tA == tB );
    EXPECT_FALSE( tA != tB );
}

TEST( Set, InequalityOperator )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    belfem::Set< int > tB{ 1, 2, 4 };

    EXPECT_TRUE( tA != tB );
    EXPECT_FALSE( tA == tB );
}

TEST( Set, EqualityDifferentSizes )
{
    belfem::Set< int > tA{ 1, 2 };
    belfem::Set< int > tB{ 1, 2, 3 };

    EXPECT_FALSE( tA == tB );
    EXPECT_TRUE( tA != tB );
}

TEST( Set, SubsetCheck )
{
    // test plan: "NOTE: verify this method exists; if not, skip."
    // Source confirms is_subset_of() exists.
    belfem::Set< int > tSmall{ 1, 2 };
    belfem::Set< int > tLarge{ 1, 2, 3 };

    EXPECT_TRUE( tSmall.is_subset_of( tLarge ) );
    EXPECT_FALSE( tLarge.is_subset_of( tSmall ) );
}

TEST( Set, SubsetOfSelf )
{
    belfem::Set< int > tA{ 1, 2, 3 };
    EXPECT_TRUE( tA.is_subset_of( tA ) );
}

TEST( Set, EmptyIsSubsetOfAnything )
{
    belfem::Set< int > tEmpty;
    belfem::Set< int > tA{ 1, 2 };

    EXPECT_TRUE( tEmpty.is_subset_of( tA ) );
    EXPECT_TRUE( tEmpty.is_subset_of( tEmpty ) );
}

TEST( Set, SupersetCheck )
{
    // test plan: "Same caveat" — is_superset_of() also exists.
    belfem::Set< int > tSmall{ 1, 2 };
    belfem::Set< int > tLarge{ 1, 2, 3 };

    EXPECT_TRUE( tLarge.is_superset_of( tSmall ) );
    EXPECT_FALSE( tSmall.is_superset_of( tLarge ) );
}

// =============================================================================
// §6.2  Set Algebra with string keys  [semantic]
// =============================================================================
// Verify set algebra also works with non-int keys.

TEST( Set, UnionOperatorString )
{
    belfem::Set< belfem::string > tA{ "a", "b" };
    belfem::Set< belfem::string > tB{ "b", "c" };

    belfem::Set< belfem::string > tResult = tA | tB;

    EXPECT_EQ( tResult.size(), 3u );
    EXPECT_TRUE( tResult.contains( "a" ) );
    EXPECT_TRUE( tResult.contains( "b" ) );
    EXPECT_TRUE( tResult.contains( "c" ) );
}

TEST( Set, IntersectionOperatorString )
{
    belfem::Set< belfem::string > tA{ "a", "b", "c" };
    belfem::Set< belfem::string > tB{ "b", "c", "d" };

    belfem::Set< belfem::string > tResult = tA & tB;

    EXPECT_EQ( tResult.size(), 2u );
    EXPECT_TRUE( tResult.contains( "b" ) );
    EXPECT_TRUE( tResult.contains( "c" ) );
}

TEST( Set, DifferenceOperatorString )
{
    belfem::Set< belfem::string > tA{ "a", "b", "c" };
    belfem::Set< belfem::string > tB{ "b" };

    belfem::Set< belfem::string > tResult = tA - tB;

    EXPECT_EQ( tResult.size(), 2u );
    EXPECT_TRUE( tResult.contains( "a" ) );
    EXPECT_TRUE( tResult.contains( "c" ) );
}
