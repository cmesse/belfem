/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Executed allocation parity of the linalg wrapper: for every expression shape the wrapper
 * must allocate exactly as often as the backend expression it forwards to. Both backends
 * allocate their buffers through posix_memalign ( armadillo_bits/memory.hpp, blaze/util/Memory.h );
 * this translation unit defines that symbol, forwards to the real one, and counts the events on
 * the calling thread. Armadillo on macOS allocates through std::malloc instead; the positive
 * control detects a platform the hook cannot see and skips the suite there.
 * See todo/linalg_wrapper_hardening_plan.md, gap #18.
 */

#include <gtest/gtest.h>
#include <cstddef>
#include <cstdlib>
#include <dlfcn.h>
#include <unistd.h>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_trans.hpp"

namespace
{
    thread_local std::size_t gAllocEvents = 0 ;
}

extern "C" int
posix_memalign( void ** aPtr, size_t aAlignment, size_t aSize )
{
    using Fn = int ( * )( void **, size_t, size_t );
    static Fn tReal = reinterpret_cast< Fn >( dlsym( RTLD_NEXT, "posix_memalign" ) );
    if ( tReal == nullptr )
    {
        // an allocator that cannot forward must not return; no allocation is possible here
        static const char tMessage[] = "test_LinalgAllocations: dlsym( RTLD_NEXT, posix_memalign ) failed\n" ;
        ( void ) ! write( 2, tMessage, sizeof( tMessage ) - 1 );
        std::abort();
    }
    ++gAllocEvents ;
    return tReal( aPtr, aAlignment, aSize );
}

namespace
{
    using belfem::real;
    using belfem::uint;
    using belfem::Matrix;
    using belfem::Vector;

    // counts the allocation events of one statement on this thread
    #define BELFEM_COUNT_ALLOCATIONS( aStatement ) \
        ( [&]() { const std::size_t tBefore = gAllocEvents; aStatement; return gAllocEvents - tBefore; }() )

    class AllocationParity : public ::testing::Test
    {
        protected:
            void
            SetUp() override
            {
                // positive control: a 48x48 backend temporary must be seen, or the hook does not
                // reach this backend's allocator on this platform
                Matrix< real > tP( 48, 48, 1.0 );
                const std::size_t tSeen = BELFEM_COUNT_ALLOCATIONS(
                    Matrix< real >::MatrixType tQ( tP.matrix_data() * 2.0 ); ( void ) tQ );
                if ( tSeen == 0 )
                {
                    GTEST_SKIP() << "the posix_memalign hook does not see this backend's allocations here "
                                    "( Armadillo on macOS uses std::malloc ); the parity gate runs on Linux and on Blaze" ;
                }
            }
    };

    // above Armadillo's 16-element local storage, and below Blaze's shared-memory thresholds:
    // a 64x64 product ( 4096 target elements ) exceeds BLAZE_SMP_DMATDMATMULT_THRESHOLD ( 3025 ) and
    // runs OpenMP-parallel under USE_OPENMP, where the runtime allocates on the calling thread
    // nondeterministically ( observed: counts of 1, 3 and 5 for the same statement ). 48x48 is 2304.
    const uint gMatrixSizes[ 2 ] = { 8, 48 };
    const uint gVectorSizes[ 2 ] = { 32, 256 };
}

TEST_F( AllocationParity, MatrixCompoundAssignmentFromExpression )
{
    for ( uint n : gMatrixSizes )
    {
        Matrix< real > tA( n, n, 1.0 ), tB( n, n, 2.0 );
        const real tDt = 0.5 ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tA += tB * tDt ),
                   BELFEM_COUNT_ALLOCATIONS( tA.matrix_data() += tB.matrix_data() * tDt ) ) << "+= n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tA -= tB * tDt ),
                   BELFEM_COUNT_ALLOCATIONS( tA.matrix_data() -= tB.matrix_data() * tDt ) ) << "-= n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tA += belfem::trans( tB ) ),
                   BELFEM_COUNT_ALLOCATIONS( tA.matrix_data() += belfem::trans( tB.matrix_data() ) ) ) << "+= trans n=" << n ;
        // the backend itself materializes an in-place product; parity, not zero
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tA *= tB * tDt ),
                   BELFEM_COUNT_ALLOCATIONS( tA.matrix_data() *= tB.matrix_data() * tDt ) ) << "*= n=" << n ;
    }
}

TEST_F( AllocationParity, MatrixAssignmentConstructionAndAliasing )
{
    for ( uint n : gMatrixSizes )
    {
        Matrix< real > tA( n, n, 1.0 ), tB( n, n, 2.0 ), tC( n, n, 0.0 );
        const real tDt = 0.5 ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tC = tB * tDt + tA ),
                   BELFEM_COUNT_ALLOCATIONS( tC.matrix_data() = tB.matrix_data() * tDt + tA.matrix_data() ) ) << "= n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( Matrix< real > tD( tB * tDt ); ( void ) tD ),
                   BELFEM_COUNT_ALLOCATIONS( Matrix< real >::MatrixType tD( tB.matrix_data() * tDt ); ( void ) tD ) ) << "ctor n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( Matrix< real > tE = tB * tDt; ( void ) tE ),
                   BELFEM_COUNT_ALLOCATIONS( Matrix< real >::MatrixType tE = tB.matrix_data() * tDt; ( void ) tE ) ) << "copy-init n=" << n ;
        // aliasing: the backend may allocate a temporary; the wrapper must not add one
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tA += tA * 2.0 ),
                   BELFEM_COUNT_ALLOCATIONS( tA.matrix_data() += tA.matrix_data() * 2.0 ) ) << "aliased += n=" << n ;
    }
}

TEST_F( AllocationParity, VectorCompoundAssignmentFromExpression )
{
    for ( uint n : gVectorSizes )
    {
        Vector< real > tY( n, 1.0 ), tX( n, 2.0 ), tZ( n, 3.0 );
        Matrix< real > tA( n, n, 1.0 );
        const real tDt = 0.5 ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY += tX * tDt ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() += tX.vector_data() * tDt ) ) << "+= x*dt n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY -= tX * tDt ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() -= tX.vector_data() * tDt ) ) << "-= x*dt n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY += tX + tZ ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() += tX.vector_data() + tZ.vector_data() ) ) << "+= x+z n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY += tA.col( 0 ) ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() += tA.col( 0 ) ) ) << "+= col n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY += tA * tX ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() += tA.matrix_data() * tX.vector_data() ) ) << "+= A*x n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY -= tA * tX ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() -= tA.matrix_data() * tX.vector_data() ) ) << "-= A*x n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( Vector< real > tW( tX * tDt ); ( void ) tW ),
                   BELFEM_COUNT_ALLOCATIONS( Vector< real >::VectorType tW( tX.vector_data() * tDt ); ( void ) tW ) ) << "ctor n=" << n ;
    }
}

TEST_F( AllocationParity, VectorRowOperandIsTransposedWithoutExtraAllocation )
{
    // the row operand is a contract the two backends implement differently, so the oracle is the
    // explicit backend form that transposes
    for ( uint n : gVectorSizes )
    {
        Vector< real > tY( n, 1.0 );
        Matrix< real > tA( n, n, 2.0 );
#ifdef BELFEM_ARMADILLO
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY += tA.row( 0 ) ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() += arma::strans( tA.row( 0 ) ) ) ) << "+= row n=" << n ;
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY = tA.row( 0 ) ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() = arma::strans( tA.row( 0 ) ) ) ) << "= row n=" << n ;
#else
        EXPECT_EQ( BELFEM_COUNT_ALLOCATIONS( tY += tA.row( 0 ) ),
                   BELFEM_COUNT_ALLOCATIONS( tY.vector_data() += blaze::trans( tA.row( 0 ) ) ) ) << "+= row n=" << n ;
#endif
    }
}
