/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Cell<T> container
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §1 (test matrix)
 */

#include <gtest/gtest.h>
#include <string>
#include <cstring>

#include "typedefs.hpp"
#include "cl_Cell.hpp"

// =============================================================================
// Typed test infrastructure
// =============================================================================

namespace
{
    // Per-type value factories for typed tests.
    template< typename T >
    struct CellTestTraits;

    template<>
    struct CellTestTraits< int >
    {
        static int val( int i ) { return i * 7; }
        static int defaultVal() { return 0; }
    };

    template<>
    struct CellTestTraits< belfem::real >
    {
        static belfem::real val( int i ) { return 3.14 * i; }
        static belfem::real defaultVal() { return 0.0; }
    };

    template<>
    struct CellTestTraits< belfem::index_t >
    {
        static belfem::index_t val( int i )
        {
            return static_cast< belfem::index_t >( i + 100 );
        }
        static belfem::index_t defaultVal() { return 0; }
    };

    template<>
    struct CellTestTraits< std::string >
    {
        static std::string val( int i ) { return "str_" + std::to_string( i ); }
        static std::string defaultVal() { return std::string(); }
    };
}

// Type list: matches the test plan §3.3
using CellTypes = ::testing::Types< int, belfem::real, belfem::index_t, std::string >;

// Single typed test suite for all Cell<T> semantic tests
template< typename T >
class Cell : public ::testing::Test {};
TYPED_TEST_SUITE( Cell, CellTypes );

// =============================================================================
// §1.1  Construction & Destruction  [semantic]
// =============================================================================

TYPED_TEST( Cell, DefaultConstructorIsEmpty )
{
    belfem::Cell< TypeParam > tCell;
    EXPECT_EQ( tCell.size(), 0u );
    EXPECT_TRUE( tCell.empty() );
}

TYPED_TEST( Cell, SizedConstructorWithValue )
{
    auto tVal = CellTestTraits< TypeParam >::val( 42 );
    belfem::Cell< TypeParam > tCell( 5, tVal );
    EXPECT_EQ( tCell.size(), 5u );
    for( size_t i = 0; i < 5; ++i )
    {
        EXPECT_EQ( tCell( i ), tVal );
    }
}

TYPED_TEST( Cell, InitializerListConstructor )
{
    auto tV0 = CellTestTraits< TypeParam >::val( 1 );
    auto tV1 = CellTestTraits< TypeParam >::val( 2 );
    auto tV2 = CellTestTraits< TypeParam >::val( 3 );

    belfem::Cell< TypeParam > tCell{ tV0, tV1, tV2 };
    EXPECT_EQ( tCell.size(), 3u );
    EXPECT_EQ( tCell( 0 ), tV0 );
    EXPECT_EQ( tCell( 1 ), tV1 );
    EXPECT_EQ( tCell( 2 ), tV2 );
}

TYPED_TEST( Cell, CopyConstructorDeepCopies )
{
    auto tVal = CellTestTraits< TypeParam >::val( 10 );
    belfem::Cell< TypeParam > tOriginal( 3, tVal );

    belfem::Cell< TypeParam > tCopy( tOriginal );

    EXPECT_EQ( tCopy.size(), tOriginal.size() );
    for( size_t i = 0; i < tCopy.size(); ++i )
    {
        EXPECT_EQ( tCopy( i ), tOriginal( i ) );
    }

    // modify copy, verify original unchanged
    tCopy( 0 ) = CellTestTraits< TypeParam >::val( 99 );
    EXPECT_EQ( tOriginal( 0 ), tVal );
}

TYPED_TEST( Cell, MoveConstructorTransfersOwnership )
{
    auto tV0 = CellTestTraits< TypeParam >::val( 1 );
    auto tV1 = CellTestTraits< TypeParam >::val( 2 );

    belfem::Cell< TypeParam > tSource{ tV0, tV1 };
    belfem::Cell< TypeParam > tTarget( std::move( tSource ) );

    EXPECT_EQ( tTarget.size(), 2u );
    EXPECT_EQ( tTarget( 0 ), tV0 );
    EXPECT_EQ( tTarget( 1 ), tV1 );

    // moved-from state: must remain valid (do NOT assert empty — STL-backed)
    EXPECT_NO_THROW(
    {
        ( void ) tSource.size();
        ( void ) tSource.empty();
        tSource.clear();
    });
}

TYPED_TEST( Cell, CopyAssignmentDeepCopies )
{
    auto tValA = CellTestTraits< TypeParam >::val( 5 );
    auto tValB = CellTestTraits< TypeParam >::val( 9 );

    belfem::Cell< TypeParam > tA( 3, tValA );
    belfem::Cell< TypeParam > tB( 2, tValB );

    tB = tA;

    EXPECT_EQ( tB.size(), 3u );
    for( size_t i = 0; i < 3; ++i )
    {
        EXPECT_EQ( tB( i ), tValA );
    }

    // verify independence
    tB( 0 ) = CellTestTraits< TypeParam >::val( 77 );
    EXPECT_EQ( tA( 0 ), tValA );
}

TYPED_TEST( Cell, MoveAssignmentTransfersOwnership )
{
    auto tV0 = CellTestTraits< TypeParam >::val( 1 );
    auto tV1 = CellTestTraits< TypeParam >::val( 2 );

    belfem::Cell< TypeParam > tSource{ tV0, tV1 };
    belfem::Cell< TypeParam > tTarget;

    tTarget = std::move( tSource );

    EXPECT_EQ( tTarget.size(), 2u );
    EXPECT_EQ( tTarget( 0 ), tV0 );
    EXPECT_EQ( tTarget( 1 ), tV1 );

    EXPECT_NO_THROW(
    {
        ( void ) tSource.size();
        ( void ) tSource.empty();
        tSource.clear();
    });
}

TYPED_TEST( Cell, SelfCopyAssignment )
{
    auto tVal = CellTestTraits< TypeParam >::val( 42 );
    belfem::Cell< TypeParam > tCell( 3, tVal );

    tCell = tCell;

    EXPECT_EQ( tCell.size(), 3u );
    EXPECT_EQ( tCell( 0 ), tVal );
}

TYPED_TEST( Cell, SelfMoveAssignment )
{
    auto tVal = CellTestTraits< TypeParam >::val( 42 );
    belfem::Cell< TypeParam > tCell( 3, tVal );

    // Self-move is technically UB per the standard, but Cell's self-assignment
    // guard (if this != &other) makes this safe.
#pragma GCC diagnostic push
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif
    tCell = std::move( tCell );
#pragma GCC diagnostic pop

    EXPECT_EQ( tCell.size(), 3u );
    EXPECT_EQ( tCell( 0 ), tVal );
}

// =============================================================================
// §1.2  Element Access  [semantic]
// =============================================================================

TYPED_TEST( Cell, ParenthesisOperatorReadWrite )
{
    belfem::Cell< TypeParam > tCell( 3, CellTestTraits< TypeParam >::val( 0 ) );

    auto tNew = CellTestTraits< TypeParam >::val( 55 );
    tCell( 1 ) = tNew;
    EXPECT_EQ( tCell( 1 ), tNew );
}

TYPED_TEST( Cell, DataPointerMatchesFirstElement )
{
    belfem::Cell< TypeParam > tCell( 3, CellTestTraits< TypeParam >::val( 1 ) );
    EXPECT_EQ( tCell.data(), &tCell( 0 ) );
}

TYPED_TEST( Cell, VectorDataReturnsInternalVector )
{
    belfem::Cell< TypeParam > tCell( 4, CellTestTraits< TypeParam >::val( 1 ) );
    EXPECT_EQ( tCell.vector_data().size(), tCell.size() );
}

TYPED_TEST( Cell, FirstAndLast )
{
    auto tV0 = CellTestTraits< TypeParam >::val( 10 );
    auto tV1 = CellTestTraits< TypeParam >::val( 20 );
    auto tV2 = CellTestTraits< TypeParam >::val( 30 );

    belfem::Cell< TypeParam > tCell{ tV0, tV1, tV2 };
    EXPECT_EQ( tCell.first(), tV0 );
    EXPECT_EQ( tCell.last(),  tV2 );
}

TYPED_TEST( Cell, FirstAndLastAreMutable )
{
    auto tV0 = CellTestTraits< TypeParam >::val( 10 );
    auto tV1 = CellTestTraits< TypeParam >::val( 20 );
    auto tV2 = CellTestTraits< TypeParam >::val( 30 );

    belfem::Cell< TypeParam > tCell{ tV0, tV1, tV2 };

    auto tNewFirst = CellTestTraits< TypeParam >::val( 88 );
    auto tNewLast  = CellTestTraits< TypeParam >::val( 99 );

    tCell.first() = tNewFirst;
    tCell.last()  = tNewLast;

    EXPECT_EQ( tCell( 0 ), tNewFirst );
    EXPECT_EQ( tCell( 2 ), tNewLast );
    EXPECT_EQ( tCell( 1 ), tV1 );
}

TYPED_TEST( Cell, ConstAccess )
{
    auto tVal = CellTestTraits< TypeParam >::val( 7 );
    belfem::Cell< TypeParam > tCell( 3, tVal );

    const belfem::Cell< TypeParam > & tRef = tCell;
    EXPECT_EQ( tRef( 0 ), tVal );
    EXPECT_EQ( tRef.size(), 3u );
}

// =============================================================================
// §1.3  Element Access  [debug]
// =============================================================================
//
// NOTE ON EXCEPTION TYPES:
//   operator()  uses BELFEM_ASSERT  → throws std::runtime_error
//   first()/last() use vector::at() → throws std::out_of_range

#ifndef NDEBUG

TEST( CellDebug, OutOfBoundsThrows )
{
    belfem::Cell< int > tCell( 5, 0 );
    EXPECT_THROW( tCell( 5 ), std::runtime_error );
    EXPECT_THROW( tCell( 100 ), std::runtime_error );
}

TEST( CellDebug, FirstOnEmptyThrows )
{
    belfem::Cell< int > tCell;
    EXPECT_THROW( tCell.first(), std::out_of_range );
}

TEST( CellDebug, LastOnEmptyThrows )
{
    belfem::Cell< int > tCell;
    EXPECT_THROW( tCell.last(), std::out_of_range );
}

#endif // NDEBUG

// =============================================================================
// §1.4  Mutation  [semantic]
// =============================================================================

TYPED_TEST( Cell, PushIncreasesSize )
{
    belfem::Cell< TypeParam > tCell;
    auto tVal = CellTestTraits< TypeParam >::val( 1 );
    tCell.push( tVal );

    EXPECT_EQ( tCell.size(), 1u );
    EXPECT_EQ( tCell.last(), tVal );

    tCell.push( CellTestTraits< TypeParam >::val( 2 ) );
    EXPECT_EQ( tCell.size(), 2u );
}

TYPED_TEST( Cell, PopReturnsLastAndShrinks )
{
    auto tV0 = CellTestTraits< TypeParam >::val( 1 );
    auto tV1 = CellTestTraits< TypeParam >::val( 2 );
    belfem::Cell< TypeParam > tCell{ tV0, tV1 };

    auto tPopped = tCell.pop();
    EXPECT_EQ( tPopped, tV1 );
    EXPECT_EQ( tCell.size(), 1u );
}

TYPED_TEST( Cell, PopFromSingleElement )
{
    belfem::Cell< TypeParam > tCell{ CellTestTraits< TypeParam >::val( 1 ) };
    tCell.pop();
    EXPECT_TRUE( tCell.empty() );
}

TYPED_TEST( Cell, SetSizeGrows )
{
    auto tV0 = CellTestTraits< TypeParam >::val( 1 );
    auto tV1 = CellTestTraits< TypeParam >::val( 2 );
    belfem::Cell< TypeParam > tCell{ tV0, tV1 };

    tCell.set_size( 5 );
    EXPECT_EQ( tCell.size(), 5u );

    // existing elements survive
    EXPECT_EQ( tCell( 0 ), tV0 );
    EXPECT_EQ( tCell( 1 ), tV1 );

    // appended slots are default-initialized (Codex finding)
    auto tDefault = CellTestTraits< TypeParam >::defaultVal();
    EXPECT_EQ( tCell( 2 ), tDefault );
    EXPECT_EQ( tCell( 3 ), tDefault );
    EXPECT_EQ( tCell( 4 ), tDefault );
}

TYPED_TEST( Cell, SetSizeShrinks )
{
    belfem::Cell< TypeParam > tCell( 5, CellTestTraits< TypeParam >::val( 0 ) );
    tCell.set_size( 2 );
    EXPECT_EQ( tCell.size(), 2u );
}

TYPED_TEST( Cell, SetSizeWithValue )
{
    auto tVal = CellTestTraits< TypeParam >::val( 42 );
    belfem::Cell< TypeParam > tCell( 2, CellTestTraits< TypeParam >::val( 0 ) );
    tCell.set_size( 5, tVal );

    EXPECT_EQ( tCell( 2 ), tVal );
    EXPECT_EQ( tCell( 3 ), tVal );
    EXPECT_EQ( tCell( 4 ), tVal );
}

TYPED_TEST( Cell, ClearMakesEmpty )
{
    belfem::Cell< TypeParam > tCell( 5, CellTestTraits< TypeParam >::val( 0 ) );
    tCell.clear();
    EXPECT_EQ( tCell.size(), 0u );
    EXPECT_TRUE( tCell.empty() );
}

TYPED_TEST( Cell, ReserveDoesNotChangeSize )
{
    belfem::Cell< TypeParam > tCell( 3, CellTestTraits< TypeParam >::val( 0 ) );
    tCell.reserve( 100 );
    EXPECT_EQ( tCell.size(), 3u );
    EXPECT_GE( tCell.capacity(), 100u );
}

TYPED_TEST( Cell, SwapExchangesContents )
{
    auto tValA = CellTestTraits< TypeParam >::val( 1 );
    auto tValB = CellTestTraits< TypeParam >::val( 2 );
    belfem::Cell< TypeParam > tA( 3, tValA );
    belfem::Cell< TypeParam > tB( 5, tValB );

    tA.swap( tB );

    EXPECT_EQ( tA.size(), 5u );
    EXPECT_EQ( tB.size(), 3u );
    EXPECT_EQ( tA( 0 ), tValB );
    EXPECT_EQ( tB( 0 ), tValA );
}

// =============================================================================
// Non-typed semantic tests
// =============================================================================

// Push-move and emplace: tested with std::string for observable move semantics.

TEST( Cell, PushMoveVersion )
{
    belfem::Cell< std::string > tCell;
    std::string tStr = "hello_move";
    tCell.push( std::move( tStr ) );

    EXPECT_EQ( tCell.size(), 1u );
    EXPECT_EQ( tCell( 0 ), "hello_move" );
}

TEST( Cell, EmplaceConstructsInPlace )
{
    belfem::Cell< std::string > tCell;
    tCell.emplace( 5, 'x' );

    EXPECT_EQ( tCell.size(), 1u );
    EXPECT_EQ( tCell( 0 ), "xxxxx" );
}

TEST( Cell, InsertAtBegin )
{
    belfem::Cell< int > tCell{ 10, 20, 30 };
    tCell.insert( tCell.begin(), 5 );

    EXPECT_EQ( tCell.size(), 4u );
    EXPECT_EQ( tCell( 0 ), 5 );
    EXPECT_EQ( tCell( 1 ), 10 );
}

TEST( Cell, InsertAtMiddle )
{
    belfem::Cell< int > tCell{ 10, 20, 30 };
    tCell.insert( tCell.begin() + 1, 15 );

    EXPECT_EQ( tCell.size(), 4u );
    EXPECT_EQ( tCell( 0 ), 10 );
    EXPECT_EQ( tCell( 1 ), 15 );
    EXPECT_EQ( tCell( 2 ), 20 );
}

TEST( Cell, InsertMoveVersion )
{
    belfem::Cell< std::string > tCell;
    tCell.push( "aaa" );
    tCell.push( "ccc" );

    std::string tVal = "bbb";
    tCell.insert( tCell.begin() + 1, std::move( tVal ) );

    EXPECT_EQ( tCell.size(), 3u );
    EXPECT_EQ( tCell( 1 ), "bbb" );
}

TEST( Cell, EraseAtPosition )
{
    belfem::Cell< int > tCell{ 10, 20, 30 };
    tCell.erase( tCell.begin() + 1 );

    EXPECT_EQ( tCell.size(), 2u );
    EXPECT_EQ( tCell( 0 ), 10 );
    EXPECT_EQ( tCell( 1 ), 30 );
}

TEST( Cell, EraseRange )
{
    belfem::Cell< int > tCell{ 10, 20, 30, 40, 50 };
    tCell.erase( tCell.begin() + 1, tCell.begin() + 4 );

    EXPECT_EQ( tCell.size(), 2u );
    EXPECT_EQ( tCell( 0 ), 10 );
    EXPECT_EQ( tCell( 1 ), 50 );
}

TEST( Cell, ShrinkToFit )
{
    belfem::Cell< int > tCell;
    tCell.reserve( 100 );
    tCell.push( 1 );
    tCell.push( 2 );

    size_t tCapBefore = tCell.capacity();
    tCell.shrink_to_fit();
    // shrink_to_fit is non-binding, but capacity should not increase
    EXPECT_LE( tCell.capacity(), tCapBefore );
    EXPECT_EQ( tCell.size(), 2u );
}

// =============================================================================
// §1.5  Free Functions  [semantic]
// =============================================================================

TEST( Cell, SortAscending )
{
    belfem::Cell< int > tCell{ 5, 2, 8, 1, 9 };
    belfem::sort( tCell );

    EXPECT_EQ( tCell( 0 ), 1 );
    EXPECT_EQ( tCell( 1 ), 2 );
    EXPECT_EQ( tCell( 2 ), 5 );
    EXPECT_EQ( tCell( 3 ), 8 );
    EXPECT_EQ( tCell( 4 ), 9 );
}

TEST( Cell, SortWithComparator )
{
    belfem::Cell< int > tCell{ 5, 2, 8, 1, 9 };
    auto tComp = []( int a, int b ){ return a > b; };
    belfem::sort( tCell, tComp );

    EXPECT_EQ( tCell( 0 ), 9 );
    EXPECT_EQ( tCell( 1 ), 8 );
    EXPECT_EQ( tCell( 2 ), 5 );
    EXPECT_EQ( tCell( 3 ), 2 );
    EXPECT_EQ( tCell( 4 ), 1 );
}

TEST( Cell, SortPartial )
{
    belfem::Cell< int > tCell{ 5, 2, 8, 1, 9 };
    auto tComp = []( int a, int b ){ return a < b; };
    belfem::sort( tCell, tComp, 3 );

    // first 3 elements sorted: 2, 5, 8
    EXPECT_EQ( tCell( 0 ), 2 );
    EXPECT_EQ( tCell( 1 ), 5 );
    EXPECT_EQ( tCell( 2 ), 8 );
    // remaining elements untouched
    EXPECT_EQ( tCell( 3 ), 1 );
    EXPECT_EQ( tCell( 4 ), 9 );
}

TEST( Cell, UniqueRemovesDuplicates )
{
    belfem::Cell< int > tCell{ 5, 2, 5, 1, 2 };
    belfem::unique( tCell );

    EXPECT_EQ( tCell.size(), 3u );
    EXPECT_EQ( tCell( 0 ), 1 );
    EXPECT_EQ( tCell( 1 ), 2 );
    EXPECT_EQ( tCell( 2 ), 5 );
}

TEST( Cell, UniqueOnAlreadyUnique )
{
    belfem::Cell< int > tCell{ 1, 2, 3 };
    belfem::unique( tCell );
    EXPECT_EQ( tCell.size(), 3u );
}

TEST( Cell, ReverseFlipsOrder )
{
    belfem::Cell< int > tCell{ 1, 2, 3 };
    belfem::reverse( tCell );

    EXPECT_EQ( tCell( 0 ), 3 );
    EXPECT_EQ( tCell( 1 ), 2 );
    EXPECT_EQ( tCell( 2 ), 1 );
}

TEST( Cell, AppendConcatenates )
{
    belfem::Cell< int > tA{ 1, 2, 3 };
    belfem::Cell< int > tB{ 4, 5 };

    belfem::append( tA, tB );

    EXPECT_EQ( tA.size(), 5u );
    EXPECT_EQ( tA( 3 ), 4 );
    EXPECT_EQ( tA( 4 ), 5 );
    EXPECT_EQ( tB.size(), 2u );
}

TEST( Cell, AppendMoveClearsSource )
{
    belfem::Cell< int > tA{ 1, 2 };
    belfem::Cell< int > tB{ 3, 4 };

    belfem::append_move( tA, tB );

    EXPECT_EQ( tA.size(), 4u );
    EXPECT_EQ( tA( 2 ), 3 );
    EXPECT_EQ( tA( 3 ), 4 );
    EXPECT_EQ( tB.size(), 0u );
}

TEST( Cell, SwapFreeFunction )
{
    belfem::Cell< int > tA{ 1, 2, 3 };
    belfem::Cell< int > tB{ 4, 5 };

    belfem::swap( tA, tB );

    EXPECT_EQ( tA.size(), 2u );
    EXPECT_EQ( tA( 0 ), 4 );
    EXPECT_EQ( tB.size(), 3u );
    EXPECT_EQ( tB( 0 ), 1 );
}

// =============================================================================
// §1.6  Pointer Semantics  [semantic]
// =============================================================================

TEST( Cell, CellOfPointersNonOwning )
{
    int tA = 10;
    int tB = 20;
    int tC = 30;

    {
        belfem::Cell< int* > tCell;
        tCell.push( &tA );
        tCell.push( &tB );
        tCell.push( &tC );

        EXPECT_EQ( *tCell( 0 ), 10 );
        EXPECT_EQ( *tCell( 1 ), 20 );
        EXPECT_EQ( *tCell( 2 ), 30 );

        tCell.clear();
    }

    // pointed-to values survive Cell destruction and clear
    EXPECT_EQ( tA, 10 );
    EXPECT_EQ( tB, 20 );
    EXPECT_EQ( tC, 30 );
}

// =============================================================================
// §1.7  Print  [semantic] (smoke test only)
// =============================================================================

TEST( Cell, PrintIntProducesOutput )
{
    belfem::Cell< int > tCell{ 1, 2, 3 };

    testing::internal::CaptureStdout();
    tCell.print( "test" );
    std::string tOutput = testing::internal::GetCapturedStdout();

    EXPECT_FALSE( tOutput.empty() );
    EXPECT_NE( tOutput.find( "test" ), std::string::npos );
}

TEST( Cell, PrintRealSpecialization )
{
    belfem::Cell< belfem::real > tCell{ 1.0, 2.0 };

    testing::internal::CaptureStdout();
    tCell.print( "reals" );
    std::string tOutput = testing::internal::GetCapturedStdout();

    EXPECT_FALSE( tOutput.empty() );
    // verify scientific notation from %+.15e format
    EXPECT_NE( tOutput.find( "e+" ), std::string::npos );
}

// =============================================================================
// §1.8  Iteration  [semantic]
// =============================================================================

TEST( Cell, RangeBasedForLoop )
{
    belfem::Cell< int > tCell{ 10, 20, 30 };
    int tSum = 0;
    for( auto & tVal : tCell )
    {
        tSum += tVal;
    }
    EXPECT_EQ( tSum, 60 );
}

TEST( Cell, ConstIteration )
{
    belfem::Cell< int > tCell{ 10, 20, 30 };
    const belfem::Cell< int > & tRef = tCell;

    int tSum = 0;
    for( const auto & tVal : tRef )
    {
        tSum += tVal;
    }
    EXPECT_EQ( tSum, 60 );
}
