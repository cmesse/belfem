/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Map<K,V> and OrderedMap<K,V> containers
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §5 (test matrix)
 *
 * Both containers are thin STL wrappers (unordered_map / map).
 * Shared tests are typed over 4 combinations: {Map, OrderedMap} × {<string,int>, <index_t,real>}.
 * Map-specific tests (get_entry) and OrderedMap-specific tests (sorted iteration)
 * are in separate groups.
 */

#include <gtest/gtest.h>
#include <string>

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_OrderedMap.hpp"
#include "cl_Set.hpp"

// =============================================================================
// Typed test infrastructure
// =============================================================================
// Traits encode the map type, key type, value type, and factory functions
// so that typed tests work across all 4 combinations without string-literal
// keys breaking index_t instantiations.

template< typename T >
struct MapTraits;

// --- Map<string, int> ---
struct MapStringInt_Tag {};
template<> struct MapTraits< MapStringInt_Tag >
{
    using MapType   = belfem::Map< belfem::string, int >;
    using KeyType   = belfem::string;
    using ValueType = int;

    static belfem::string key( int i )   { return "key_" + std::to_string( i ); }
    static int    val( int i )   { return i * 10; }
    static belfem::string missingKey()   { return "nonexistent"; }
};

// --- Map<index_t, real> ---
struct MapIndexReal_Tag {};
template<> struct MapTraits< MapIndexReal_Tag >
{
    using MapType   = belfem::Map< belfem::index_t, belfem::real >;
    using KeyType   = belfem::index_t;
    using ValueType = belfem::real;

    static belfem::index_t key( int i )  { return static_cast< belfem::index_t >( i + 100 ); }
    static belfem::real    val( int i )  { return 3.14 * i; }
    static belfem::index_t missingKey()  { return 999999; }
};

// --- OrderedMap<string, int> ---
struct OMapStringInt_Tag {};
template<> struct MapTraits< OMapStringInt_Tag >
{
    using MapType   = belfem::OrderedMap< belfem::string, int >;
    using KeyType   = belfem::string;
    using ValueType = int;

    static belfem::string key( int i )   { return "key_" + std::to_string( i ); }
    static int    val( int i )   { return i * 10; }
    static belfem::string missingKey()   { return "nonexistent"; }
};

// --- OrderedMap<index_t, real> ---
struct OMapIndexReal_Tag {};
template<> struct MapTraits< OMapIndexReal_Tag >
{
    using MapType   = belfem::OrderedMap< belfem::index_t, belfem::real >;
    using KeyType   = belfem::index_t;
    using ValueType = belfem::real;

    static belfem::index_t key( int i )  { return static_cast< belfem::index_t >( i + 100 ); }
    static belfem::real    val( int i )  { return 3.14 * i; }
    static belfem::index_t missingKey()  { return 999999; }
};

// Type list for shared tests: all 4 combinations
template< typename T >
class Map : public ::testing::Test {};

using AllMapTypes = ::testing::Types<
    MapStringInt_Tag,
    MapIndexReal_Tag,
    OMapStringInt_Tag,
    OMapIndexReal_Tag
>;

TYPED_TEST_SUITE( Map, AllMapTypes );

// =============================================================================
// §5.1  Shared Tests  [semantic]
// =============================================================================

TYPED_TEST( Map, DefaultConstructorEmpty )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    EXPECT_EQ( tMap.size(), 0u );
    EXPECT_TRUE( tMap.empty() );
}

TYPED_TEST( Map, InsertViaSubscript )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    tMap[ Traits::key( 1 ) ] = Traits::val( 1 );

    EXPECT_EQ( tMap.size(), 1u );
    EXPECT_FALSE( tMap.empty() );
}

TYPED_TEST( Map, LookupViaParenthesis )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    auto tKey = Traits::key( 1 );
    auto tVal = Traits::val( 1 );
    tMap[ tKey ] = tVal;

    EXPECT_EQ( tMap( tKey ), tVal );
}

TYPED_TEST( Map, ConstLookupViaParenthesis )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    auto tKey = Traits::key( 1 );
    auto tVal = Traits::val( 1 );
    tMap[ tKey ] = tVal;

    const typename Traits::MapType & tRef = tMap;
    EXPECT_EQ( tRef( tKey ), tVal );
}

TYPED_TEST( Map, KeyExists )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    auto tKey = Traits::key( 1 );
    tMap[ tKey ] = Traits::val( 1 );

    EXPECT_TRUE( tMap.key_exists( tKey ) );
    EXPECT_FALSE( tMap.key_exists( Traits::missingKey() ) );
}

TYPED_TEST( Map, EraseKey )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    auto tKey = Traits::key( 1 );
    tMap[ tKey ] = Traits::val( 1 );
    EXPECT_EQ( tMap.size(), 1u );

    tMap.erase_key( tKey );

    EXPECT_EQ( tMap.size(), 0u );
    EXPECT_FALSE( tMap.key_exists( tKey ) );
}

TYPED_TEST( Map, Clear )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    tMap[ Traits::key( 1 ) ] = Traits::val( 1 );
    tMap[ Traits::key( 2 ) ] = Traits::val( 2 );

    tMap.clear();

    EXPECT_EQ( tMap.size(), 0u );
    EXPECT_TRUE( tMap.empty() );
}

TYPED_TEST( Map, CopySemantics )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tA;

    auto tKey = Traits::key( 1 );
    auto tVal = Traits::val( 1 );
    tA[ tKey ] = tVal;

    typename Traits::MapType tB( tA );

    EXPECT_EQ( tB.size(), 1u );
    EXPECT_EQ( tB( tKey ), tVal );

    // verify independence: modify copy, check original
    tB[ tKey ] = Traits::val( 99 );
    EXPECT_EQ( tA( tKey ), tVal );
}

TYPED_TEST( Map, MoveSemantics )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tSource;

    auto tKey = Traits::key( 1 );
    auto tVal = Traits::val( 1 );
    tSource[ tKey ] = tVal;

    typename Traits::MapType tTarget( std::move( tSource ) );

    EXPECT_EQ( tTarget.size(), 1u );
    EXPECT_EQ( tTarget( tKey ), tVal );

    // source remains valid (do NOT assert empty — STL-backed)
    // verify multiple operations are safe on moved-from source (idea from Grok)
    EXPECT_NO_THROW(
    {
        ( void ) tSource.size();
        tSource.clear();
        tSource[ Traits::key( 2 ) ] = Traits::val( 2 );
    });
}

TYPED_TEST( Map, CopyAssignment )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tA;

    auto tKey = Traits::key( 1 );
    auto tVal = Traits::val( 1 );
    tA[ tKey ] = tVal;

    typename Traits::MapType tB;
    tB[ Traits::key( 99 ) ] = Traits::val( 99 );

    tB = tA;

    EXPECT_EQ( tB.size(), 1u );
    EXPECT_EQ( tB( tKey ), tVal );

    // verify independence
    tB[ tKey ] = Traits::val( 77 );
    EXPECT_EQ( tA( tKey ), tVal );
}

TYPED_TEST( Map, MoveAssignment )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tSource;

    auto tKey = Traits::key( 1 );
    auto tVal = Traits::val( 1 );
    tSource[ tKey ] = tVal;

    typename Traits::MapType tTarget;
    tTarget[ Traits::key( 99 ) ] = Traits::val( 99 );

    tTarget = std::move( tSource );

    EXPECT_EQ( tTarget.size(), 1u );
    EXPECT_EQ( tTarget( tKey ), tVal );

    // source remains valid (do NOT assert empty — STL-backed)
    EXPECT_NO_THROW(
    {
        ( void ) tSource.size();
        tSource.clear();
        tSource[ Traits::key( 2 ) ] = Traits::val( 2 );
    });
}

TYPED_TEST( Map, OverwriteExistingKey )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    auto tKey = Traits::key( 1 );
    tMap[ tKey ] = Traits::val( 1 );
    tMap[ tKey ] = Traits::val( 2 );

    EXPECT_EQ( tMap.size(), 1u );
    EXPECT_EQ( tMap( tKey ), Traits::val( 2 ) );
}

TYPED_TEST( Map, FindReturnsIterator )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    auto tKey = Traits::key( 1 );
    tMap[ tKey ] = Traits::val( 1 );

    // existing key → valid iterator
    auto tIt = tMap.find( tKey );
    EXPECT_NE( tIt, tMap.end() );
    EXPECT_EQ( tIt->second, Traits::val( 1 ) );

    // missing key → end()
    auto tMissing = tMap.find( Traits::missingKey() );
    EXPECT_EQ( tMissing, tMap.end() );
}

TYPED_TEST( Map, Iteration )
{
    using Traits = MapTraits< TypeParam >;
    typename Traits::MapType tMap;

    auto tK1 = Traits::key( 1 );
    auto tK2 = Traits::key( 2 );
    auto tK3 = Traits::key( 3 );
    tMap[ tK1 ] = Traits::val( 1 );
    tMap[ tK2 ] = Traits::val( 2 );
    tMap[ tK3 ] = Traits::val( 3 );

    // track which specific keys were visited, not just the count
    // (idea from ChatGPT — stronger for unordered maps)
    size_t tCount = 0;
    bool tSaw1 = false;
    bool tSaw2 = false;
    bool tSaw3 = false;

    for( const auto & tEntry : tMap )
    {
        ++tCount;
        if( tEntry.first == tK1 ) tSaw1 = true;
        if( tEntry.first == tK2 ) tSaw2 = true;
        if( tEntry.first == tK3 ) tSaw3 = true;
    }

    EXPECT_EQ( tCount, 3u );
    EXPECT_TRUE( tSaw1 );
    EXPECT_TRUE( tSaw2 );
    EXPECT_TRUE( tSaw3 );
}

// =============================================================================
// §5.2  Shared Tests  [debug]
// =============================================================================
//
// operator() uses BELFEM_ASSERT in debug (throws std::runtime_error)
// and BELFEM_ERROR in release (aborts).

#ifndef NDEBUG

TEST( MapDebug, LookupMissingKeyThrows_StringInt )
{
    belfem::Map< belfem::string, int > tMap;
    tMap[ "exists" ] = 42;

    EXPECT_THROW( tMap( "nonexistent" ), std::runtime_error );
}

TEST( MapDebug, LookupMissingKeyThrows_IndexReal )
{
    belfem::Map< belfem::index_t, belfem::real > tMap;
    tMap[ 1 ] = 3.14;

    EXPECT_THROW( tMap( 999u ), std::runtime_error );
}

TEST( OrderedMapDebug, LookupMissingKeyThrows_StringInt )
{
    belfem::OrderedMap< belfem::string, int > tMap;
    tMap[ "exists" ] = 42;

    EXPECT_THROW( tMap( "nonexistent" ), std::runtime_error );
}

TEST( OrderedMapDebug, LookupMissingKeyThrows_IndexReal )
{
    belfem::OrderedMap< belfem::index_t, belfem::real > tMap;
    tMap[ 1 ] = 3.14;

    EXPECT_THROW( tMap( 999u ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §5.1 supplement: Map-only tests  [semantic]
// =============================================================================
// get_entry() exists on Map only, NOT on OrderedMap.

TEST( Map, GetEntryByIndex )
{
    belfem::Map< belfem::string, int > tMap;
    tMap[ "alpha" ] = 1;
    tMap[ "beta" ]  = 2;
    tMap[ "gamma" ] = 3;

    // get_entry(0) should return a valid key-value pair.
    // Iteration order is unspecified for unordered_map, so we only
    // verify the returned entry actually exists in the map.
    auto tEntry = tMap.get_entry( 0 );
    EXPECT_TRUE( tMap.key_exists( tEntry.first ) );
    EXPECT_EQ( tEntry.second, tMap( tEntry.first ) );

    // verify all 3 entries are reachable via get_entry, with no duplicates
    // (Codex finding: old test would pass even if same entry returned for all indices)
    belfem::Set< belfem::string > tSeen;
    for( belfem::index_t i = 0; i < 3; ++i )
    {
        auto tE = tMap.get_entry( i );
        EXPECT_TRUE( tMap.key_exists( tE.first ) );
        tSeen.insert( tE.first );
    }
    EXPECT_EQ( tSeen.size(), 3u );
}

// --- Map-only debug ---

#ifndef NDEBUG

TEST( MapDebug, GetEntryOutOfBoundsThrows )
{
    belfem::Map< belfem::string, int > tMap;
    tMap[ "a" ] = 1;

    EXPECT_THROW( tMap.get_entry( 1 ), std::runtime_error );
    EXPECT_THROW( tMap.get_entry( 100 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §5.1 supplement: map_data() exposure  [semantic]
// =============================================================================

TEST( Map, MapDataExposesUnderlying )
{
    belfem::Map< belfem::string, int > tMap;
    tMap[ "a" ] = 1;

    auto & tInternal = tMap.map_data();
    tInternal[ "b" ] = 2;

    EXPECT_EQ( tMap.size(), 2u );
    EXPECT_TRUE( tMap.key_exists( "b" ) );
}

TEST( OrderedMap, MapDataExposesUnderlying )
{
    belfem::OrderedMap< belfem::string, int > tMap;
    tMap[ "a" ] = 1;

    auto & tInternal = tMap.map_data();
    tInternal[ "b" ] = 2;

    EXPECT_EQ( tMap.size(), 2u );
    EXPECT_TRUE( tMap.key_exists( "b" ) );
}

// =============================================================================
// §5.3  OrderedMap-Specific  [semantic]
// =============================================================================

TEST( OrderedMap, IterationOrderIsSorted_IntKeys )
{
    belfem::OrderedMap< int, int > tMap;

    // insert in non-sorted order
    tMap[ 3 ] = 30;
    tMap[ 1 ] = 10;
    tMap[ 2 ] = 20;

    // iteration should yield keys in sorted order: 1, 2, 3
    // use iterator-based checking with ASSERT_NE guards (idea from ChatGPT)
    auto tIt = tMap.begin();

    ASSERT_NE( tIt, tMap.end() );
    EXPECT_EQ( tIt->first, 1 );
    ++tIt;

    ASSERT_NE( tIt, tMap.end() );
    EXPECT_EQ( tIt->first, 2 );
    ++tIt;

    ASSERT_NE( tIt, tMap.end() );
    EXPECT_EQ( tIt->first, 3 );
    ++tIt;

    EXPECT_EQ( tIt, tMap.end() );
}

TEST( OrderedMap, IterationOrderIsSorted_StringKeys )
{
    belfem::OrderedMap< belfem::string, int > tMap;

    tMap[ "cherry" ] = 3;
    tMap[ "apple" ]  = 1;
    tMap[ "banana" ] = 2;

    // iteration should yield keys in lexicographic order
    auto tIt = tMap.begin();

    ASSERT_NE( tIt, tMap.end() );
    EXPECT_EQ( tIt->first, "apple" );
    ++tIt;

    ASSERT_NE( tIt, tMap.end() );
    EXPECT_EQ( tIt->first, "banana" );
    ++tIt;

    ASSERT_NE( tIt, tMap.end() );
    EXPECT_EQ( tIt->first, "cherry" );
    ++tIt;

    EXPECT_EQ( tIt, tMap.end() );
}

// =============================================================================
// §5.4  KeyToString  [semantic]
// =============================================================================

TEST( Map, StringKeyToString )
{
    EXPECT_EQ( belfem::map::KeyToString( belfem::string( "hello" ) ), "hello" );
}

TEST( Map, IntKeyToString )
{
    EXPECT_EQ( belfem::map::KeyToString( 42 ), "42" );
}

TEST( Map, UnsignedKeyToString )
{
    EXPECT_EQ( belfem::map::KeyToString( 42u ), "42" );
}

TEST( Map, LongUnsignedKeyToString )
{
    EXPECT_EQ( belfem::map::KeyToString( 42LU ), "42" );
}

TEST( Map, UnknownTypeReturnsUnknown )
{
    EXPECT_EQ( belfem::map::KeyToString( 3.14 ), "unknown" );
}
