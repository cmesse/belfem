/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for StringList container
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §8 (test matrix)
 *
 * All tests in this file are [valgrind] candidates (manual malloc/free).
 *
 * StringList is explicitly non-copyable and non-movable (copy/move
 * constructors and assignments are = delete).  Do NOT write copy/move tests.
 */

#include <gtest/gtest.h>
#include <cstring>    // strcmp, strlen
#include <string>
#include <type_traits>

#include "typedefs.hpp"
#include "cl_StringList.hpp"

// =============================================================================
// §8.1  Tests  [semantic]
// =============================================================================

TEST( StringList, TypeTraitsNonCopyableNonMovable )
{
    // compile-time proof that copy/move are deleted (idea from ChatGPT)
    static_assert( !std::is_copy_constructible< belfem::StringList >::value,
            "StringList must not be copy-constructible" );
    static_assert( !std::is_copy_assignable< belfem::StringList >::value,
            "StringList must not be copy-assignable" );
    static_assert( !std::is_move_constructible< belfem::StringList >::value,
            "StringList must not be move-constructible" );
    static_assert( !std::is_move_assignable< belfem::StringList >::value,
            "StringList must not be move-assignable" );

    SUCCEED();
}

TEST( StringList, ConstructAndPush )
{
    belfem::StringList tList( 5 );

    tList.push( "alpha" );
    tList.push( "beta" );
    tList.push( "gamma" );

    EXPECT_STREQ( tList.item( 0 ), "alpha" );
    EXPECT_STREQ( tList.item( 1 ), "beta" );
    EXPECT_STREQ( tList.item( 2 ), "gamma" );
}

TEST( StringList, ItemReturnsCString )
{
    belfem::StringList tList( 3 );
    tList.push( "hello" );

    // item() returns a const char*, verify with strcmp
    EXPECT_EQ( std::strcmp( tList.item( 0 ), "hello" ), 0 );
}

TEST( StringList, DataPointerValid )
{
    belfem::StringList tList( 3 );
    tList.push( "first" );
    tList.push( "second" );

    char ** tData = tList.data();

    EXPECT_NE( tData, nullptr );

    // data()[i] should match item(i) — both content and pointer identity
    EXPECT_STREQ( tData[ 0 ], tList.item( 0 ) );
    EXPECT_STREQ( tData[ 1 ], tList.item( 1 ) );
    EXPECT_EQ( tData[ 0 ], tList.item( 0 ) );
    EXPECT_EQ( tData[ 1 ], tList.item( 1 ) );
}

TEST( StringList, PushEmptyString )
{
    belfem::StringList tList( 3 );
    tList.push( "" );

    EXPECT_STREQ( tList.item( 0 ), "" );
    EXPECT_EQ( std::strlen( tList.item( 0 ) ), 0u );
}

TEST( StringList, PushLongString )
{
    // 1000-char string should survive the malloc + strcpy round-trip
    std::string tLong( 1000, 'x' );

    belfem::StringList tList( 1 );
    tList.push( tLong );

    EXPECT_STREQ( tList.item( 0 ), tLong.c_str() );
    EXPECT_EQ( std::strlen( tList.item( 0 ) ), 1000u );
}

TEST( StringList, DistinctStoragePerEntry )
{
    // push the same string twice — each push does its own malloc,
    // so the pointers must differ even though content is identical
    // (idea from ChatGPT)
    belfem::StringList tList( 2 );

    tList.push( "alpha" );
    tList.push( "alpha" );

    EXPECT_STREQ( tList.data()[ 0 ], "alpha" );
    EXPECT_STREQ( tList.data()[ 1 ], "alpha" );

    // different malloc allocations → different pointers
    EXPECT_NE( tList.data()[ 0 ], tList.data()[ 1 ] );
}

TEST( StringList, EarlierEntriesRemainStableAfterLaterPush )
{
    // verify that later pushes don't invalidate earlier pointers
    // (catches realloc-style bugs in the pointer array)
    // (idea from ChatGPT)
    belfem::StringList tList( 3 );

    tList.push( "first" );
    const char * tFirstBefore = tList.item( 0 );

    tList.push( "second" );
    tList.push( "third" );

    // pointer to first entry should not have changed
    EXPECT_EQ( tList.item( 0 ), tFirstBefore );
    EXPECT_STREQ( tFirstBefore, "first" );

    // all entries still correct
    EXPECT_STREQ( tList.item( 0 ), "first" );
    EXPECT_STREQ( tList.item( 1 ), "second" );
    EXPECT_STREQ( tList.item( 2 ), "third" );
}

// =============================================================================
// §8.2  Tests  [debug]
// =============================================================================
//
// push() uses BELFEM_ASSERT(mCount < mMemory) — debug only.
// item() uses BELFEM_ASSERT(aIndex < mCount) — debug only.
// Both throw std::runtime_error in debug builds.

#ifndef NDEBUG

TEST( StringListDebug, PushBeyondCapacityThrows )
{
    belfem::StringList tList( 2 );
    tList.push( "one" );
    tList.push( "two" );

    // third push exceeds capacity of 2
    EXPECT_THROW( tList.push( "three" ), std::runtime_error );
}

TEST( StringListDebug, ItemOutOfBoundsThrows )
{
    belfem::StringList tList( 5 );
    tList.push( "only" );

    // mCount is 1, so index 1 is out of bounds
    EXPECT_THROW( tList.item( 1 ), std::runtime_error );

    // index 0 should still work
    EXPECT_NO_THROW( tList.item( 0 ) );
}

TEST( StringListDebug, ItemOnEmptyListThrows )
{
    // edge case: zero pushes → mCount == 0 → even index 0 is OOB
    // (idea from ChatGPT)
    belfem::StringList tList( 3 );

    EXPECT_THROW( tList.item( 0 ), std::runtime_error );
}

#endif // NDEBUG
