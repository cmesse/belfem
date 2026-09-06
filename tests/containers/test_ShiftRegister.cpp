/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for ShiftRegister<T> container
 * See: tests_00_strategy.md (conventions), tests_01_containers.md §3 (test matrix)
 *
 * All tests in this file are [valgrind] candidates (manual malloc/free).
 *
 * NOTE: ShiftRegister uses malloc/free internally.  Trivially-copyable types
 * (int, real, index_t, …) are stored raw.  Owning types are supported only when
 * whitelisted via belfem::is_shift_register_safe (Vector<T>, Matrix<T>): for
 * those the register manages element lifetime (placement-new on reserve,
 * destroy before free).  A non-whitelisted owning type such as std::string is
 * (correctly) rejected by the static_assert.  The typed suite below covers the
 * trivially-copyable path; the OwningType section covers Vector<real>.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_ShiftRegister.hpp"

// =============================================================================
// Typed test infrastructure
// =============================================================================
// Test plan §3.3: typed tests over int, real.

template< typename T >
class ShiftRegister : public ::testing::Test {};

using ShiftTypes = ::testing::Types< int, belfem::real >;
TYPED_TEST_SUITE( ShiftRegister, ShiftTypes );

// Helper: cast int literal to the test type.
template< typename T >
T tval( int i ) { return static_cast< T >( i ); }

// =============================================================================
// §3.1  Construction & Destruction  [semantic]
// =============================================================================

TYPED_TEST( ShiftRegister, CapacityConstructor )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );

    EXPECT_EQ( tReg.capacity(), 5u );
    EXPECT_EQ( tReg.size(), 0u );
    EXPECT_TRUE( tReg.empty() );
}

TYPED_TEST( ShiftRegister, CapacityWithFillValue )
{
    auto tVal = tval< TypeParam >( 7 );
    belfem::ShiftRegister< TypeParam > tReg( 5, tVal );

    EXPECT_EQ( tReg.capacity(), 5u );
    EXPECT_EQ( tReg.size(), 5u );
    EXPECT_TRUE( tReg.full() );

    for( size_t i = 0; i < 5; ++i )
    {
        EXPECT_EQ( tReg( i ), tVal );
    }
}

TYPED_TEST( ShiftRegister, InitializerListConstructor )
{
    // initializer list stores in list order (NOT push order):
    // {1, 2, 3} → (0)==1, (1)==2, (2)==3
    belfem::ShiftRegister< TypeParam > tReg{
        tval< TypeParam >( 1 ),
        tval< TypeParam >( 2 ),
        tval< TypeParam >( 3 )
    };

    EXPECT_EQ( tReg.size(), 3u );
    EXPECT_EQ( tReg.capacity(), 3u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 1 ) );
    EXPECT_EQ( tReg( 1 ), tval< TypeParam >( 2 ) );
    EXPECT_EQ( tReg( 2 ), tval< TypeParam >( 3 ) );
}

TYPED_TEST( ShiftRegister, CopyConstructorDeepCopies )
{
    belfem::ShiftRegister< TypeParam > tA{
        tval< TypeParam >( 10 ),
        tval< TypeParam >( 20 )
    };

    belfem::ShiftRegister< TypeParam > tB( tA );

    EXPECT_EQ( tB.size(), 2u );
    EXPECT_EQ( tB( 0 ), tval< TypeParam >( 10 ) );
    EXPECT_EQ( tB( 1 ), tval< TypeParam >( 20 ) );

    // modify copy, verify original unchanged
    tB( 0 ) = tval< TypeParam >( 99 );
    EXPECT_EQ( tA( 0 ), tval< TypeParam >( 10 ) );
}

TYPED_TEST( ShiftRegister, CopyDoesNotPreserveRevertState )
{
    belfem::ShiftRegister< TypeParam > tA( 5 );
    tA.push( tval< TypeParam >( 1 ) );
    tA.push( tval< TypeParam >( 2 ) );
    // tA can now revert

    belfem::ShiftRegister< TypeParam > tB( tA );

    // copy explicitly resets revert state — revert should fail.
    // revert() uses BELFEM_ERROR, which throws std::runtime_error in debug.
#ifndef NDEBUG
    EXPECT_THROW( tB.revert(), std::runtime_error );
#endif
}

TYPED_TEST( ShiftRegister, MoveConstructor )
{
    belfem::ShiftRegister< TypeParam > tA{
        tval< TypeParam >( 1 ),
        tval< TypeParam >( 2 )
    };

    // capture raw pointer before move (idea from ChatGPT)
    TypeParam * tOldData = tA.data();

    belfem::ShiftRegister< TypeParam > tB( std::move( tA ) );

    EXPECT_EQ( tB.size(), 2u );
    EXPECT_EQ( tB( 0 ), tval< TypeParam >( 1 ) );
    EXPECT_EQ( tB.data(), tOldData );   // zero-copy pointer transfer

    // source nullified
    EXPECT_EQ( tA.size(), 0u );
    EXPECT_EQ( tA.capacity(), 0u );
    EXPECT_EQ( tA.data(), nullptr );
}

TYPED_TEST( ShiftRegister, CopyAssignment )
{
    belfem::ShiftRegister< TypeParam > tA{
        tval< TypeParam >( 10 ),
        tval< TypeParam >( 20 )
    };
    belfem::ShiftRegister< TypeParam > tB( 5 );

    tB = tA;

    EXPECT_EQ( tB.size(), 2u );
    // verify target capacity matches source (idea from ChatGPT)
    EXPECT_EQ( tB.capacity(), tA.capacity() );
    EXPECT_EQ( tB( 0 ), tval< TypeParam >( 10 ) );
    EXPECT_EQ( tB( 1 ), tval< TypeParam >( 20 ) );

    // verify independence
    tB( 0 ) = tval< TypeParam >( 99 );
    EXPECT_EQ( tA( 0 ), tval< TypeParam >( 10 ) );
}

TYPED_TEST( ShiftRegister, MoveAssignment )
{
    belfem::ShiftRegister< TypeParam > tA{
        tval< TypeParam >( 1 ),
        tval< TypeParam >( 2 )
    };

    // capture raw pointer (idea from ChatGPT)
    TypeParam * tOldData = tA.data();

    belfem::ShiftRegister< TypeParam > tB( 1 );
    tB = std::move( tA );

    EXPECT_EQ( tB.size(), 2u );
    EXPECT_EQ( tB( 0 ), tval< TypeParam >( 1 ) );
    EXPECT_EQ( tB.data(), tOldData );   // pointer transfer

    // source nullified
    EXPECT_EQ( tA.size(), 0u );
    EXPECT_EQ( tA.capacity(), 0u );
    EXPECT_EQ( tA.data(), nullptr );
}

TYPED_TEST( ShiftRegister, SelfCopyAssignment )
{
    belfem::ShiftRegister< TypeParam > tReg{
        tval< TypeParam >( 1 ),
        tval< TypeParam >( 2 )
    };

    tReg = tReg;

    EXPECT_EQ( tReg.size(), 2u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 1 ) );
}

TYPED_TEST( ShiftRegister, SelfMoveAssignment )
{
    belfem::ShiftRegister< TypeParam > tReg{
        tval< TypeParam >( 1 ),
        tval< TypeParam >( 2 )
    };

#pragma GCC diagnostic push
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif
    tReg = std::move( tReg );
#pragma GCC diagnostic pop

    // self-move guard: state preserved
    EXPECT_EQ( tReg.size(), 2u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 1 ) );
}

// =============================================================================
// §3.2  Push & Shift  [semantic]
// =============================================================================

TYPED_TEST( ShiftRegister, PushToEmpty )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 10 ) );

    EXPECT_EQ( tReg.size(), 1u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 10 ) );
}

TYPED_TEST( ShiftRegister, PushShiftsRight )
{
    // push inserts at index 0 and shifts existing elements right:
    // push(A), push(B), push(C) → (0)==C, (1)==B, (2)==A
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );   // A
    tReg.push( tval< TypeParam >( 2 ) );   // B
    tReg.push( tval< TypeParam >( 3 ) );   // C

    EXPECT_EQ( tReg.size(), 3u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 3 ) );  // newest
    EXPECT_EQ( tReg( 1 ), tval< TypeParam >( 2 ) );
    EXPECT_EQ( tReg( 2 ), tval< TypeParam >( 1 ) );  // oldest
}

TYPED_TEST( ShiftRegister, PushBeyondCapacity )
{
    belfem::ShiftRegister< TypeParam > tReg( 3 );
    tReg.push( tval< TypeParam >( 1 ) );
    tReg.push( tval< TypeParam >( 2 ) );
    tReg.push( tval< TypeParam >( 3 ) );
    tReg.push( tval< TypeParam >( 4 ) );   // drops oldest (1)

    EXPECT_EQ( tReg.size(), 3u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 4 ) );
    EXPECT_EQ( tReg( 1 ), tval< TypeParam >( 3 ) );
    EXPECT_EQ( tReg( 2 ), tval< TypeParam >( 2 ) );
}

TYPED_TEST( ShiftRegister, PushConstRef )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    const TypeParam tVal = tval< TypeParam >( 42 );
    tReg.push( tVal );

    EXPECT_EQ( tReg.size(), 1u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 42 ) );
}

TYPED_TEST( ShiftRegister, PushLvalueRef )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    TypeParam tVal = tval< TypeParam >( 42 );
    tReg.push( tVal );

    EXPECT_EQ( tReg.size(), 1u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 42 ) );
}

TYPED_TEST( ShiftRegister, SizeNeverExceedsCapacity )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    for( int i = 0; i < 100; ++i )
    {
        tReg.push( tval< TypeParam >( i ) );
    }
    EXPECT_EQ( tReg.size(), 5u );
}

TYPED_TEST( ShiftRegister, FullAndEmptyStates )
{
    belfem::ShiftRegister< TypeParam > tReg( 2 );

    EXPECT_TRUE( tReg.empty() );
    EXPECT_FALSE( tReg.full() );

    tReg.push( tval< TypeParam >( 1 ) );
    EXPECT_FALSE( tReg.empty() );
    EXPECT_FALSE( tReg.full() );

    tReg.push( tval< TypeParam >( 2 ) );
    EXPECT_FALSE( tReg.empty() );
    EXPECT_TRUE( tReg.full() );
}

// =============================================================================
// §3.3  Revert  [semantic]
// =============================================================================

TYPED_TEST( ShiftRegister, RevertUndoesPush )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );
    tReg.push( tval< TypeParam >( 2 ) );

    tReg.revert();

    EXPECT_EQ( tReg.size(), 1u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 1 ) );
}

TYPED_TEST( ShiftRegister, RevertWhenFull )
{
    // register was full before the last push → revert recovers the dropped value
    belfem::ShiftRegister< TypeParam > tReg( 2 );
    tReg.push( tval< TypeParam >( 1 ) );
    tReg.push( tval< TypeParam >( 2 ) );   // now full: (0)==2, (1)==1

    tReg.push( tval< TypeParam >( 3 ) );   // drops 1: (0)==3, (1)==2
    tReg.revert();                           // recovers: (0)==2, (1)==1

    EXPECT_EQ( tReg.size(), 2u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 2 ) );
    EXPECT_EQ( tReg( 1 ), tval< TypeParam >( 1 ) );
}

TYPED_TEST( ShiftRegister, RevertRestoresSize )
{
    // register was NOT full before push → revert decreases size
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );
    EXPECT_EQ( tReg.size(), 1u );

    tReg.push( tval< TypeParam >( 2 ) );
    EXPECT_EQ( tReg.size(), 2u );

    tReg.revert();
    EXPECT_EQ( tReg.size(), 1u );
}

TYPED_TEST( ShiftRegister, PushAfterRevertEnablesNewRevert )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );
    tReg.push( tval< TypeParam >( 2 ) );

    tReg.revert();
    tReg.push( tval< TypeParam >( 3 ) );

    // should be able to revert again
    EXPECT_NO_THROW( tReg.revert() );

    // verify final state: back to just the original push (idea from ChatGPT)
    EXPECT_EQ( tReg.size(), 1u );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 1 ) );
}

// =============================================================================
// §3.4  Revert  [debug]
// =============================================================================
//
// revert() uses BELFEM_ERROR (always active), which throws
// std::runtime_error in debug builds and aborts in release.
//
// RevertMakesSubsequentRevertFail from §3.3 is equivalent to
// DoubleRevertThrows — consolidated here in the debug section.

#ifndef NDEBUG

TYPED_TEST( ShiftRegister, RevertWithoutPushThrows )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    EXPECT_THROW( tReg.revert(), std::runtime_error );
}

TYPED_TEST( ShiftRegister, DoubleRevertThrows )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );
    tReg.push( tval< TypeParam >( 2 ) );

    tReg.revert();
    EXPECT_THROW( tReg.revert(), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §3.5  Access  [semantic]
// =============================================================================

TYPED_TEST( ShiftRegister, ParenthesisOperatorReadWrite )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 10 ) );

    tReg( 0 ) = tval< TypeParam >( 20 );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 20 ) );
}

TYPED_TEST( ShiftRegister, ConstAccess )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 10 ) );

    const belfem::ShiftRegister< TypeParam > & tRef = tReg;
    EXPECT_EQ( tRef( 0 ), tval< TypeParam >( 10 ) );
}

TYPED_TEST( ShiftRegister, DataPointerValid )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    EXPECT_NE( tReg.data(), nullptr );

    // also verify data pointer is usable after push
    tReg.push( tval< TypeParam >( 9 ) );
    EXPECT_EQ( tReg.data()[ 0 ], tval< TypeParam >( 9 ) );
}

// =============================================================================
// §3.6  Access  [debug]
// =============================================================================

#ifndef NDEBUG

TYPED_TEST( ShiftRegister, OutOfBoundsThrows )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );

    // size is 1, so index 1 is out of bounds
    EXPECT_THROW( tReg( 1 ), std::runtime_error );
}

TYPED_TEST( ShiftRegister, AccessOnEmptyThrows )
{
    // accessing index 0 on an empty register (idea from Grok)
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    EXPECT_THROW( tReg( 0 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §3.7  Other  [semantic]
// =============================================================================

TYPED_TEST( ShiftRegister, ClearResetsSize )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );
    tReg.push( tval< TypeParam >( 2 ) );

    tReg.clear();

    EXPECT_EQ( tReg.size(), 0u );
    EXPECT_TRUE( tReg.empty() );
}

TYPED_TEST( ShiftRegister, ClearResetsRevertState )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );
    tReg.push( tval< TypeParam >( 2 ) );

    tReg.clear();

#ifndef NDEBUG
    EXPECT_THROW( tReg.revert(), std::runtime_error );
#endif
}

TYPED_TEST( ShiftRegister, FillSetsAllElements )
{
    // CRITICAL: fill() operates on mSize, NOT mCapacity.
    // Must create with (capacity, value) to have size == capacity first,
    // then fill with a different value.
    belfem::ShiftRegister< TypeParam > tReg( 5, tval< TypeParam >( 0 ) );

    tReg.fill( tval< TypeParam >( 42 ) );

    EXPECT_EQ( tReg.size(), 5u );
    for( size_t i = 0; i < 5; ++i )
    {
        EXPECT_EQ( tReg( i ), tval< TypeParam >( 42 ) );
    }
}

TYPED_TEST( ShiftRegister, FillResetsRevertState )
{
    belfem::ShiftRegister< TypeParam > tReg( 5, tval< TypeParam >( 0 ) );
    tReg.push( tval< TypeParam >( 1 ) );
    // can now revert

    tReg.fill( tval< TypeParam >( 2 ) );

#ifndef NDEBUG
    EXPECT_THROW( tReg.revert(), std::runtime_error );
#endif
}

TYPED_TEST( ShiftRegister, ReserveDestroysPrevious )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );
    EXPECT_EQ( tReg.size(), 1u );

    tReg.reserve( 10 );

    EXPECT_EQ( tReg.capacity(), 10u );
    EXPECT_EQ( tReg.size(), 0u );
}

TYPED_TEST( ShiftRegister, ReserveSameCapacityNoOp )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );
    EXPECT_EQ( tReg.size(), 1u );

    // capture pointer to verify no reallocation (idea from ChatGPT)
    TypeParam * tOldData = tReg.data();

    tReg.reserve( 5 );

    // same capacity → early return, size and data preserved
    EXPECT_EQ( tReg.capacity(), 5u );
    EXPECT_EQ( tReg.size(), 1u );
    EXPECT_EQ( tReg.data(), tOldData );
    EXPECT_EQ( tReg( 0 ), tval< TypeParam >( 1 ) );
}

TYPED_TEST( ShiftRegister, FreeReleasesMemory )
{
    belfem::ShiftRegister< TypeParam > tReg( 5 );
    tReg.push( tval< TypeParam >( 1 ) );

    tReg.free();

    EXPECT_EQ( tReg.data(), nullptr );
    EXPECT_EQ( tReg.size(), 0u );
    EXPECT_EQ( tReg.capacity(), 0u );
}

TYPED_TEST( ShiftRegister, IteratorRange )
{
    belfem::ShiftRegister< TypeParam > tReg{
        tval< TypeParam >( 10 ),
        tval< TypeParam >( 20 ),
        tval< TypeParam >( 30 )
    };

    size_t tCount = 0;
    TypeParam tSum = tval< TypeParam >( 0 );
    for( auto & tX : tReg )
    {
        tSum += tX;
        ++tCount;
    }

    EXPECT_EQ( tCount, 3u );
    // 10 + 20 + 30 = 60
    EXPECT_EQ( tSum, tval< TypeParam >( 60 ) );
}

// =============================================================================
// §3.8  Owning element type  [semantic] [valgrind]
// =============================================================================
//
// Vector<real> owns a heap buffer (arma::Mat), so these tests exercise the
// lifetime-managed path: placement-new on reserve, element move/copy on
// push/shift/revert, and destroy-before-free.  Run under valgrind/ASan to
// confirm no leak, use-after-free, or double-free.

namespace
{
    using Vec = belfem::Vector< belfem::real >;

    // build a length-1 vector holding value aVal
    Vec
    vec1( belfem::real aVal )
    {
        Vec tV( 1 );
        tV( 0 ) = aVal;
        return tV;
    }

    // compare a stored vector against an expected length-1 payload
    ::testing::AssertionResult
    vec_is( const Vec & aVec, belfem::real aExpected )
    {
        if( aVec.length() != 1 )
        {
            return ::testing::AssertionFailure()
                   << "length " << aVec.length() << " != 1";
        }
        if( aVec( 0 ) != aExpected )
        {
            return ::testing::AssertionFailure()
                   << "value " << aVec( 0 ) << " != " << aExpected;
        }
        return ::testing::AssertionSuccess();
    }
}

TEST( ShiftRegisterOwning, PushShiftAndAccess )
{
    belfem::ShiftRegister< Vec > tReg( 5 );
    tReg.push( vec1( 1.0 ) );
    tReg.push( vec1( 2.0 ) );
    tReg.push( vec1( 3.0 ) );

    EXPECT_EQ( tReg.size(), 3u );
    EXPECT_TRUE( vec_is( tReg( 0 ), 3.0 ) );   // newest
    EXPECT_TRUE( vec_is( tReg( 1 ), 2.0 ) );
    EXPECT_TRUE( vec_is( tReg( 2 ), 1.0 ) );   // oldest
}

TEST( ShiftRegisterOwning, PushBeyondCapacityDropsOldest )
{
    // full-buffer shift exercises move-assignment through the backup slot
    belfem::ShiftRegister< Vec > tReg( 3 );
    tReg.push( vec1( 1.0 ) );
    tReg.push( vec1( 2.0 ) );
    tReg.push( vec1( 3.0 ) );
    tReg.push( vec1( 4.0 ) );   // drops 1.0

    EXPECT_EQ( tReg.size(), 3u );
    EXPECT_TRUE( vec_is( tReg( 0 ), 4.0 ) );
    EXPECT_TRUE( vec_is( tReg( 1 ), 3.0 ) );
    EXPECT_TRUE( vec_is( tReg( 2 ), 2.0 ) );
}

TEST( ShiftRegisterOwning, RevertRecoversDroppedValue )
{
    belfem::ShiftRegister< Vec > tReg( 2 );
    tReg.push( vec1( 1.0 ) );
    tReg.push( vec1( 2.0 ) );   // full: (0)=2, (1)=1
    tReg.push( vec1( 3.0 ) );   // drops 1: (0)=3, (1)=2
    tReg.revert();              // recovers: (0)=2, (1)=1

    EXPECT_EQ( tReg.size(), 2u );
    EXPECT_TRUE( vec_is( tReg( 0 ), 2.0 ) );
    EXPECT_TRUE( vec_is( tReg( 1 ), 1.0 ) );
}

TEST( ShiftRegisterOwning, CopyConstructorDeepCopies )
{
    belfem::ShiftRegister< Vec > tA( 5 );
    tA.push( vec1( 10.0 ) );
    tA.push( vec1( 20.0 ) );

    belfem::ShiftRegister< Vec > tB( tA );

    EXPECT_EQ( tB.size(), 2u );
    EXPECT_TRUE( vec_is( tB( 0 ), 20.0 ) );

    // mutate copy; original must keep its own heap buffer
    tB( 0 ) = vec1( 99.0 );
    EXPECT_TRUE( vec_is( tA( 0 ), 20.0 ) );
    EXPECT_TRUE( vec_is( tB( 0 ), 99.0 ) );
}

TEST( ShiftRegisterOwning, CopyAssignmentDeepCopies )
{
    belfem::ShiftRegister< Vec > tA( 5 );
    tA.push( vec1( 10.0 ) );
    tA.push( vec1( 20.0 ) );

    belfem::ShiftRegister< Vec > tB( 3 );
    tB.push( vec1( 7.0 ) );
    tB = tA;

    EXPECT_EQ( tB.size(), 2u );
    EXPECT_EQ( tB.capacity(), tA.capacity() );
    EXPECT_TRUE( vec_is( tB( 0 ), 20.0 ) );

    tB( 0 ) = vec1( 99.0 );
    EXPECT_TRUE( vec_is( tA( 0 ), 20.0 ) );
}

TEST( ShiftRegisterOwning, MoveTransfersBuffer )
{
    belfem::ShiftRegister< Vec > tA( 5 );
    tA.push( vec1( 1.0 ) );
    tA.push( vec1( 2.0 ) );

    belfem::ShiftRegister< Vec > tB( std::move( tA ) );

    EXPECT_EQ( tB.size(), 2u );
    EXPECT_TRUE( vec_is( tB( 0 ), 2.0 ) );
    EXPECT_EQ( tA.data(), nullptr );
    EXPECT_EQ( tA.size(), 0u );
}

TEST( ShiftRegisterOwning, ReserveDestroysAndRebuilds )
{
    // reserve() must destroy the old owning slots (no leak) and construct new
    belfem::ShiftRegister< Vec > tReg( 3 );
    tReg.push( vec1( 1.0 ) );
    tReg.push( vec1( 2.0 ) );

    tReg.reserve( 6 );
    EXPECT_EQ( tReg.capacity(), 6u );
    EXPECT_EQ( tReg.size(), 0u );

    tReg.push( vec1( 5.0 ) );
    EXPECT_TRUE( vec_is( tReg( 0 ), 5.0 ) );
}

TEST( ShiftRegisterOwning, FillAssignsOwningSlots )
{
    belfem::ShiftRegister< Vec > tReg( 3, vec1( 0.0 ) );
    EXPECT_EQ( tReg.size(), 3u );

    tReg.fill( vec1( 42.0 ) );
    for( size_t i = 0; i < 3; ++i )
    {
        EXPECT_TRUE( vec_is( tReg( i ), 42.0 ) );
    }
}

TEST( ShiftRegisterOwning, FreeReleasesSlots )
{
    belfem::ShiftRegister< Vec > tReg( 3 );
    tReg.push( vec1( 1.0 ) );

    tReg.free();   // destroys owning slots then frees buffer
    EXPECT_EQ( tReg.data(), nullptr );
    EXPECT_EQ( tReg.capacity(), 0u );
}
