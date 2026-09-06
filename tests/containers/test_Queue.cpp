/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Queue<T> container
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §7 (test matrix)
 *
 * Queue<T> is a thin wrapper around std::queue<T>.  Low priority, Phase 3.
 * The only BELFEM-authored method is pop(), which combines front()+pop()
 * into one call returning the value (STL's pop() returns void).
 *
 * No BELFEM_ASSERT/BELFEM_ERROR in this class → no debug-only tests.
 *
 * NOTE: Tests use plain TEST() rather than typed tests because
 * TypeParam(int_literal) won't compile for std::string — there is no
 * std::string(int) constructor.  String-specific tests use string literals.
 */

#include <gtest/gtest.h>
#include <string>

#include "typedefs.hpp"
#include "cl_Queue.hpp"
#include "cl_Cell.hpp"

// =============================================================================
// §7.1  Tests  [semantic]
// =============================================================================

TEST( Queue, DefaultConstructorEmpty )
{
    belfem::Queue< int > tQueue;

    EXPECT_EQ( tQueue.size(), 0u );
    EXPECT_TRUE( tQueue.empty() );
}

TEST( Queue, PushIncreasesSize )
{
    belfem::Queue< int > tQueue;

    tQueue.push( 10 );
    EXPECT_EQ( tQueue.size(), 1u );
    EXPECT_FALSE( tQueue.empty() );

    tQueue.push( 20 );
    EXPECT_EQ( tQueue.size(), 2u );

    tQueue.push( 30 );
    EXPECT_EQ( tQueue.size(), 3u );
}

TEST( Queue, PopReturnsFront )
{
    // FIFO: push 1, 2, 3 → pop returns 1, then 2, then 3
    belfem::Queue< int > tQueue;
    tQueue.push( 1 );
    tQueue.push( 2 );
    tQueue.push( 3 );

    EXPECT_EQ( tQueue.pop(), 1 );
    EXPECT_EQ( tQueue.pop(), 2 );
    EXPECT_EQ( tQueue.pop(), 3 );
}

TEST( Queue, PopDecreasesSize )
{
    belfem::Queue< int > tQueue;
    tQueue.push( 10 );
    tQueue.push( 20 );
    tQueue.push( 30 );
    EXPECT_EQ( tQueue.size(), 3u );

    tQueue.pop();
    EXPECT_EQ( tQueue.size(), 2u );

    tQueue.pop();
    EXPECT_EQ( tQueue.size(), 1u );

    tQueue.pop();
    EXPECT_EQ( tQueue.size(), 0u );
    EXPECT_TRUE( tQueue.empty() );
}

TEST( Queue, ConstructFromCell )
{
    // Queue(Cell<T>&) converts from Cell — note non-const reference
    belfem::Cell< int > tCell{ 1, 2, 3 };

    belfem::Queue< int > tQueue( tCell );

    EXPECT_EQ( tQueue.size(), 3u );

    // FIFO: should pop in Cell order: 1, 2, 3
    EXPECT_EQ( tQueue.pop(), 1 );
    EXPECT_EQ( tQueue.pop(), 2 );
    EXPECT_EQ( tQueue.pop(), 3 );
    EXPECT_TRUE( tQueue.empty() );

    // verify source Cell is unchanged — constructor iterates by value,
    // not by move (idea from ChatGPT)
    EXPECT_EQ( tCell.size(), 3u );
    EXPECT_EQ( tCell( 0 ), 1 );
    EXPECT_EQ( tCell( 1 ), 2 );
    EXPECT_EQ( tCell( 2 ), 3 );
}

TEST( Queue, ConstructFromEmptyCell )
{
    // edge case: empty Cell → empty Queue (idea from ChatGPT)
    belfem::Cell< int > tCell;

    belfem::Queue< int > tQueue( tCell );

    EXPECT_TRUE( tQueue.empty() );
    EXPECT_EQ( tQueue.size(), 0u );
}

TEST( Queue, CopyConstructor )
{
    belfem::Queue< int > tA;
    tA.push( 10 );
    tA.push( 20 );

    belfem::Queue< int > tB( tA );

    // verify copy has same contents
    EXPECT_EQ( tB.size(), 2u );
    EXPECT_EQ( tB.pop(), 10 );
    EXPECT_EQ( tB.pop(), 20 );

    // verify original is unaffected by popping from copy
    EXPECT_EQ( tA.size(), 2u );
    EXPECT_EQ( tA.pop(), 10 );
}

TEST( Queue, CopyAssignment )
{
    // separate from copy constructor (idea from ChatGPT)
    belfem::Queue< int > tSource;
    tSource.push( 5 );
    tSource.push( 6 );

    belfem::Queue< int > tTarget;
    tTarget.push( 99 );

    tTarget = tSource;

    // target has source's data
    EXPECT_EQ( tTarget.size(), 2u );
    EXPECT_EQ( tTarget.pop(), 5 );
    EXPECT_EQ( tTarget.pop(), 6 );

    // source still has its data
    EXPECT_EQ( tSource.size(), 2u );
    EXPECT_EQ( tSource.pop(), 5 );
    EXPECT_EQ( tSource.pop(), 6 );
}

TEST( Queue, MoveConstructor )
{
    belfem::Queue< int > tA;
    tA.push( 10 );
    tA.push( 20 );
    tA.push( 30 );

    belfem::Queue< int > tB( std::move( tA ) );

    // target has the data
    EXPECT_EQ( tB.size(), 3u );
    EXPECT_EQ( tB.pop(), 10 );
    EXPECT_EQ( tB.pop(), 20 );
    EXPECT_EQ( tB.pop(), 30 );

    // source remains valid after move (do NOT assert empty — STL-backed,
    // Codex finding: empty() is not portably guaranteed for moved-from std::queue)
    EXPECT_NO_THROW(
    {
        ( void ) tA.size();
        tA.push( 99 );
    });
}

TEST( Queue, MoveAssignment )
{
    // separate from move constructor (idea from ChatGPT)
    belfem::Queue< int > tSource;
    tSource.push( 11 );
    tSource.push( 12 );

    belfem::Queue< int > tTarget;
    tTarget.push( 99 );

    tTarget = std::move( tSource );

    // target has source's data
    EXPECT_EQ( tTarget.size(), 2u );
    EXPECT_EQ( tTarget.pop(), 11 );
    EXPECT_EQ( tTarget.pop(), 12 );
    EXPECT_TRUE( tTarget.empty() );
}

TEST( Queue, SelfCopyAssignment )
{
    // (idea from ChatGPT)
    belfem::Queue< int > tQueue;
    tQueue.push( 1 );
    tQueue.push( 2 );

    tQueue = tQueue;

    EXPECT_EQ( tQueue.size(), 2u );
    EXPECT_EQ( tQueue.pop(), 1 );
    EXPECT_EQ( tQueue.pop(), 2 );
}

TEST( Queue, SelfMoveAssignment )
{
    // Queue has no self-move guard (defaulted std::queue move assignment).
    // Self-move is UB per the standard and may empty the queue.
    // We only verify the object remains usable afterwards.
    belfem::Queue< int > tQueue;
    tQueue.push( 1 );
    tQueue.push( 2 );

#pragma GCC diagnostic push
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif
    tQueue = std::move( tQueue );
#pragma GCC diagnostic pop

    // must remain in a valid state (can push/pop without crash)
    EXPECT_NO_THROW( tQueue.push( 99 ) );
}

// --- string type smoke tests ---

TEST( Queue, PopReturnsFront_String )
{
    belfem::Queue< belfem::string > tQueue;
    tQueue.push( "first" );
    tQueue.push( "second" );
    tQueue.push( "third" );

    EXPECT_EQ( tQueue.pop(), "first" );
    EXPECT_EQ( tQueue.pop(), "second" );
    EXPECT_EQ( tQueue.pop(), "third" );
}

TEST( Queue, ConstructFromCell_String )
{
    belfem::Cell< belfem::string > tCell{ "alpha", "beta", "gamma" };

    belfem::Queue< belfem::string > tQueue( tCell );

    EXPECT_EQ( tQueue.size(), 3u );
    EXPECT_EQ( tQueue.pop(), "alpha" );
    EXPECT_EQ( tQueue.pop(), "beta" );
    EXPECT_EQ( tQueue.pop(), "gamma" );

    // source Cell unchanged
    EXPECT_EQ( tCell.size(), 3u );
}
