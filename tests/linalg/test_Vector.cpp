/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Vector<T> wrapper and vector-centric free functions/operators
 * See: tests_00_strategy.md, tests_02_linalg.md §1, §3.1–3.2, §4.1–4.2, §4.10–4.14
 *
 * Vector<T> wraps arma::Mat<T> (Nx1) or blaze::DynamicVector<T>.
 * Tests go through the unified BELFEM API — no #ifdef for backend.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <string>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "fn_dot.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_sum.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"
#include "fn_linspace.hpp"
#include "fn_sort.hpp"
#include "fn_unique.hpp"
#include "fn_reverse.hpp"
#include "fn_append.hpp"
#include "fn_combine.hpp"
#include "op_VectorPlus.hpp"
#include "op_VectorMinus.hpp"
#include "op_VectorTimes.hpp"
#include "op_VectorDivide.hpp"
#include "op_VectorEqualEqual.hpp"
// operator% (elementwise multiply) included transitively via cl_Vector.hpp

namespace
{
    const belfem::real tEps = 1e-12;  // for exact-in-theory results
}

// =============================================================================
// §1.1  Construction & Destruction  [semantic]
// =============================================================================

TEST( Vector, DefaultConstructorEmpty )
{
    belfem::Vector< belfem::real > tVec;
    EXPECT_EQ( tVec.length(), 0u );
}

TEST( Vector, SizedConstructor )
{
    belfem::Vector< belfem::real > tVec( 5 );
    EXPECT_EQ( tVec.length(), 5u );
}

TEST( Vector, SizedWithFillValue )
{
    belfem::Vector< belfem::real > tVec( 5, 3.0 );
    EXPECT_EQ( tVec.length(), 5u );
    for( size_t i = 0; i < 5; ++i )
    {
        EXPECT_NEAR( tVec( i ), 3.0, tEps );
    }
}

TEST( Vector, InitializerListMultipleElements )
{
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
    EXPECT_EQ( tVec.length(), 3u );
    EXPECT_NEAR( tVec( 0 ), 1.0, tEps );
    EXPECT_NEAR( tVec( 1 ), 2.0, tEps );
    EXPECT_NEAR( tVec( 2 ), 3.0, tEps );
}

TEST( Vector, InitializerListSingleElement )
{
    belfem::Vector< belfem::real > tVec = { 42.0 };
    EXPECT_EQ( tVec.length(), 1u );
    EXPECT_NEAR( tVec( 0 ), 42.0, tEps );
}

TEST( Vector, CopyConstructorDeepCopies )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB( tA );

    tB( 0 ) = 99.0;
    EXPECT_NEAR( tA( 0 ), 1.0, tEps );   // original unchanged
    EXPECT_NEAR( tB( 0 ), 99.0, tEps );
}

TEST( Vector, MoveConstructor )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB( std::move( tA ) );

    EXPECT_EQ( tB.length(), 3u );
    EXPECT_NEAR( tB( 0 ), 1.0, tEps );

    // source is valid moved-from state
    EXPECT_NO_THROW( ( void ) tA.length() );
}

TEST( Vector, CopyAssignment )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0 };
    belfem::Vector< belfem::real > tB;
    tB = tA;
    tB( 0 ) = 99.0;
    EXPECT_NEAR( tA( 0 ), 1.0, tEps );
}

TEST( Vector, MoveAssignment )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0 };
    belfem::Vector< belfem::real > tB;
    tB = std::move( tA );
    EXPECT_EQ( tB.length(), 2u );
    EXPECT_NEAR( tB( 0 ), 1.0, tEps );
}

TEST( Vector, SelfCopyAssignment )
{
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
    tVec = tVec;
    EXPECT_EQ( tVec.length(), 3u );
    EXPECT_NEAR( tVec( 1 ), 2.0, tEps );
}

TEST( Vector, SelfMoveAssignment )
{
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
#pragma GCC diagnostic push
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif
    tVec = std::move( tVec );
#pragma GCC diagnostic pop
    EXPECT_EQ( tVec.length(), 3u );
    EXPECT_NEAR( tVec( 1 ), 2.0, tEps );
}

TEST( Vector, AssignFromScalar )
{
    belfem::Vector< belfem::real > tVec( 4 );
    tVec = 5.0;
    for( size_t i = 0; i < 4; ++i )
    {
        EXPECT_NEAR( tVec( i ), 5.0, tEps );
    }
}

TEST( Vector, AssignFromInitializerList )
{
    belfem::Vector< belfem::real > tVec( 10 );
    tVec = { 1.0, 2.0 };
    EXPECT_EQ( tVec.length(), 2u );
    EXPECT_NEAR( tVec( 0 ), 1.0, tEps );
    EXPECT_NEAR( tVec( 1 ), 2.0, tEps );
}

TEST( Vector, AssignFromSingleElementInitializerList )
{
    // tests the single-element branch of operator=(initializer_list)
    // which calls set_size(1, value) (idea from ChatGPT)
    belfem::Vector< belfem::real > tVec( 5, 0.0 );
    tVec = { 7.5 };
    EXPECT_EQ( tVec.length(), 1u );
    EXPECT_NEAR( tVec( 0 ), 7.5, tEps );
}

// =============================================================================
// §1.2  Memory & Access  [semantic]
// =============================================================================

TEST( Vector, DataPointerNonNull )
{
    belfem::Vector< belfem::real > tVec( 3 );
    EXPECT_NE( tVec.data(), nullptr );
}

TEST( Vector, DataPointerMatchesFirstElement )
{
    belfem::Vector< belfem::real > tVec = { 10.0, 20.0, 30.0 };
    for( size_t i = 0; i < 3; ++i )
    {
        EXPECT_NEAR( tVec.data()[ i ], tVec( i ), tEps );
    }
}

TEST( Vector, ColumnMajorLayout )
{
    // Vector is Nx1 — contiguous access via data()[i] must match operator()(i)
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    for( size_t i = 0; i < tVec.length(); ++i )
    {
        EXPECT_NEAR( tVec.data()[ i ], tVec( i ), tEps );
    }
}

TEST( Vector, ConstDataPointer )
{
    belfem::Vector< belfem::real > tVec = { 7.0, 8.0 };
    const belfem::Vector< belfem::real > & tRef = tVec;
    EXPECT_EQ( tRef.data(), tVec.data() );
}

TEST( Vector, VectorDataExposesBackend )
{
    // modify through data() pointer, verify wrapper sees changes
    // NOTE: use data() pointer, not backend operator, to stay backend-agnostic
    // (Armadillo uses operator(), Blaze uses operator[])
    belfem::Vector< belfem::real > tVec( 3, 0.0 );

    belfem::real * tRaw = tVec.data();
    tRaw[ 0 ] = 1.0;
    tRaw[ 1 ] = 2.0;
    tRaw[ 2 ] = 3.0;

    EXPECT_NEAR( tVec( 0 ), 1.0, tEps );
    EXPECT_NEAR( tVec( 1 ), 2.0, tEps );
    EXPECT_NEAR( tVec( 2 ), 3.0, tEps );
}

TEST( Vector, ParenthesisReadWrite )
{
    belfem::Vector< belfem::real > tVec( 3 );
    tVec( 0 ) = 42.0;
    tVec( 1 ) = -7.5;
    tVec( 2 ) = 0.0;
    EXPECT_NEAR( tVec( 0 ), 42.0, tEps );
    EXPECT_NEAR( tVec( 1 ), -7.5, tEps );
    EXPECT_NEAR( tVec( 2 ), 0.0,  tEps );
}

TEST( Vector, ConstParenthesisAccess )
{
    belfem::Vector< belfem::real > tVec = { 3.14 };
    const belfem::Vector< belfem::real > & tRef = tVec;
    EXPECT_NEAR( tRef( 0 ), 3.14, tEps );
}

TEST( Vector, LengthMatchesSized )
{
    belfem::Vector< belfem::real > tVec( 17 );
    EXPECT_EQ( tVec.length(), 17u );
}

TEST( Vector, Iteration )
{
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
    size_t tCount = 0;
    belfem::real tSum = 0.0;
    for( auto tVal : tVec )
    {
        tSum += tVal;
        ++tCount;
    }
    EXPECT_EQ( tCount, 3u );
    EXPECT_NEAR( tSum, 6.0, tEps );
}

// =============================================================================
// §1.3  Mutation  [semantic]
// =============================================================================

TEST( Vector, FillSetsAllElements )
{
    belfem::Vector< belfem::real > tVec( 5 );
    tVec.fill( 7.0 );
    for( size_t i = 0; i < 5; ++i )
    {
        EXPECT_NEAR( tVec( i ), 7.0, tEps );
    }
}

TEST( Vector, SetSizeChangesLength )
{
    belfem::Vector< belfem::real > tVec;
    tVec.set_size( 10 );
    EXPECT_EQ( tVec.length(), 10u );
}

TEST( Vector, SetSizeWithValue )
{
    belfem::Vector< belfem::real > tVec;
    tVec.set_size( 10, 3.0 );
    EXPECT_EQ( tVec.length(), 10u );
    for( size_t i = 0; i < 10; ++i )
    {
        EXPECT_NEAR( tVec( i ), 3.0, tEps );
    }
}

// =============================================================================
// §1.4  Compound Operators  [semantic]
// =============================================================================

TEST( Vector, PlusEqualsScalar )
{
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
    tVec += 2.0;
    EXPECT_NEAR( tVec( 0 ), 3.0, tEps );
    EXPECT_NEAR( tVec( 1 ), 4.0, tEps );
    EXPECT_NEAR( tVec( 2 ), 5.0, tEps );
}

TEST( Vector, PlusEqualsVector )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB = { 10.0, 20.0, 30.0 };
    tA += tB;
    EXPECT_NEAR( tA( 0 ), 11.0, tEps );
    EXPECT_NEAR( tA( 1 ), 22.0, tEps );
    EXPECT_NEAR( tA( 2 ), 33.0, tEps );
}

TEST( Vector, MinusEqualsScalar )
{
    belfem::Vector< belfem::real > tVec = { 5.0, 7.0, 9.0 };
    tVec -= 2.0;
    EXPECT_NEAR( tVec( 0 ), 3.0, tEps );
    EXPECT_NEAR( tVec( 1 ), 5.0, tEps );
    EXPECT_NEAR( tVec( 2 ), 7.0, tEps );
}

TEST( Vector, MinusEqualsVector )
{
    belfem::Vector< belfem::real > tA = { 10.0, 20.0, 30.0 };
    belfem::Vector< belfem::real > tB = { 1.0, 2.0, 3.0 };
    tA -= tB;
    EXPECT_NEAR( tA( 0 ), 9.0, tEps );
    EXPECT_NEAR( tA( 1 ), 18.0, tEps );
    EXPECT_NEAR( tA( 2 ), 27.0, tEps );
}

TEST( Vector, TimesEqualsScalar )
{
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
    tVec *= 3.0;
    EXPECT_NEAR( tVec( 0 ), 3.0, tEps );
    EXPECT_NEAR( tVec( 1 ), 6.0, tEps );
    EXPECT_NEAR( tVec( 2 ), 9.0, tEps );
}

TEST( Vector, DivideEqualsScalar )
{
    belfem::Vector< belfem::real > tVec = { 2.0, 4.0, 6.0 };
    tVec /= 2.0;
    EXPECT_NEAR( tVec( 0 ), 1.0, tEps );
    EXPECT_NEAR( tVec( 1 ), 2.0, tEps );
    EXPECT_NEAR( tVec( 2 ), 3.0, tEps );
}

TEST( Vector, ElementwiseMultiplyEquals )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB = { 4.0, 5.0, 6.0 };
    tA %= tB;
    EXPECT_NEAR( tA( 0 ), 4.0,  tEps );
    EXPECT_NEAR( tA( 1 ), 10.0, tEps );
    EXPECT_NEAR( tA( 2 ), 18.0, tEps );
}

// =============================================================================
// §1.5  Access  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( VectorDebug, OutOfBoundsThrows )
{
    belfem::Vector< belfem::real > tVec( 3 );
    EXPECT_THROW( tVec( 3 ), std::runtime_error );
}

TEST( VectorDebug, ElementwiseMultiplySizeMismatchThrows )
{
    belfem::Vector< belfem::real > tA( 3 );
    belfem::Vector< belfem::real > tB( 5 );
    EXPECT_THROW( tA %= tB, std::runtime_error );
}

TEST( VectorDebug, VectorPlusSizeMismatchThrows )
{
    belfem::Vector< belfem::real > tA( 3 );
    belfem::Vector< belfem::real > tB( 5 );
    EXPECT_THROW( tA + tB, std::runtime_error );
}

TEST( VectorDebug, VectorMinusSizeMismatchThrows )
{
    belfem::Vector< belfem::real > tA( 3 );
    belfem::Vector< belfem::real > tB( 5 );
    EXPECT_THROW( tA - tB, std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §1.6  Print  [semantic] (smoke test only)
// =============================================================================

TEST( Vector, PrintProducesOutput )
{
    // capture stdout, verify non-empty (idea from ChatGPT)
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0 };
    testing::internal::CaptureStdout();
    tVec.print( "tVec" );
    std::string tOut = testing::internal::GetCapturedStdout();
    EXPECT_FALSE( tOut.empty() );
}

// =============================================================================
// §3.1  Vector Binary Operators  [semantic]
// =============================================================================

TEST( Vector, VectorPlusVector )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB = { 4.0, 5.0, 6.0 };
    belfem::Vector< belfem::real > tC( tA + tB );
    EXPECT_NEAR( tC( 0 ), 5.0, tEps );
    EXPECT_NEAR( tC( 1 ), 7.0, tEps );
    EXPECT_NEAR( tC( 2 ), 9.0, tEps );
}

TEST( Vector, VectorPlusScalar )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tC( tA + 10.0 );
    EXPECT_NEAR( tC( 0 ), 11.0, tEps );
    EXPECT_NEAR( tC( 1 ), 12.0, tEps );
    EXPECT_NEAR( tC( 2 ), 13.0, tEps );
}

TEST( Vector, VectorMinusVector )
{
    belfem::Vector< belfem::real > tA = { 5.0, 7.0, 9.0 };
    belfem::Vector< belfem::real > tB = { 4.0, 5.0, 6.0 };
    belfem::Vector< belfem::real > tC( tA - tB );
    EXPECT_NEAR( tC( 0 ), 1.0, tEps );
    EXPECT_NEAR( tC( 1 ), 2.0, tEps );
    EXPECT_NEAR( tC( 2 ), 3.0, tEps );
}

TEST( Vector, VectorTimesScalar )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tC( tA * 2.0 );
    EXPECT_NEAR( tC( 0 ), 2.0, tEps );
    EXPECT_NEAR( tC( 1 ), 4.0, tEps );
    EXPECT_NEAR( tC( 2 ), 6.0, tEps );
}

TEST( Vector, ScalarTimesVector )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tC( 2.0 * tA );
    EXPECT_NEAR( tC( 0 ), 2.0, tEps );
    EXPECT_NEAR( tC( 1 ), 4.0, tEps );
    EXPECT_NEAR( tC( 2 ), 6.0, tEps );
}

TEST( Vector, VectorDivideScalar )
{
    belfem::Vector< belfem::real > tA = { 2.0, 4.0, 6.0 };
    belfem::Vector< belfem::real > tC( tA / 2.0 );
    EXPECT_NEAR( tC( 0 ), 1.0, tEps );
    EXPECT_NEAR( tC( 1 ), 2.0, tEps );
    EXPECT_NEAR( tC( 2 ), 3.0, tEps );
}

TEST( Vector, VectorEqualityTrue )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB = { 1.0, 2.0, 3.0 };
    EXPECT_TRUE( tA == tB );
}

TEST( Vector, VectorEqualityFalse )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB = { 1.0, 2.0, 4.0 };
    EXPECT_FALSE( tA == tB );
}

TEST( Vector, VectorEqualityScalar )
{
    belfem::Vector< belfem::real > tA = { 5.0, 5.0, 5.0 };
    EXPECT_TRUE( tA == 5.0 );
    belfem::Vector< belfem::real > tB = { 5.0, 5.0, 4.0 };
    EXPECT_FALSE( tB == 5.0 );
}

// =============================================================================
// §4.1  Dot Product  [semantic]
// =============================================================================

TEST( Vector, DotProductBasic )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB = { 4.0, 5.0, 6.0 };
    EXPECT_NEAR( belfem::dot( tA, tB ), 32.0, tEps );
}

TEST( Vector, DotProductOrthogonal )
{
    belfem::Vector< belfem::real > tA = { 1.0, 0.0 };
    belfem::Vector< belfem::real > tB = { 0.0, 1.0 };
    EXPECT_NEAR( belfem::dot( tA, tB ), 0.0, tEps );
}

TEST( Vector, DotProductSelf )
{
    belfem::Vector< belfem::real > tV = { 3.0, 4.0 };
    belfem::real tN = belfem::norm( tV );
    EXPECT_NEAR( belfem::dot( tV, tV ), tN * tN, tEps );
}

// =============================================================================
// §4.2  Cross Product  [semantic]
// =============================================================================

TEST( Vector, CrossProductBasis )
{
    belfem::Vector< belfem::real > tA = { 1.0, 0.0, 0.0 };
    belfem::Vector< belfem::real > tB = { 0.0, 1.0, 0.0 };
    belfem::Vector< belfem::real > tC( belfem::cross( tA, tB ) );
    EXPECT_NEAR( tC( 0 ), 0.0, tEps );
    EXPECT_NEAR( tC( 1 ), 0.0, tEps );
    EXPECT_NEAR( tC( 2 ), 1.0, tEps );
}

TEST( Vector, CrossProductAnticommutative )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tB = { 4.0, 5.0, 6.0 };
    belfem::Vector< belfem::real > tAB( belfem::cross( tA, tB ) );
    belfem::Vector< belfem::real > tBA( belfem::cross( tB, tA ) );

    for( size_t i = 0; i < 3; ++i )
    {
        EXPECT_NEAR( tAB( i ), -tBA( i ), tEps );
    }
}

TEST( Vector, CrossProductSelf )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tC( belfem::cross( tA, tA ) );
    for( size_t i = 0; i < 3; ++i )
    {
        EXPECT_NEAR( tC( i ), 0.0, tEps );
    }
}

// =============================================================================
// §4.2  Cross Product  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( VectorDebug, CrossProductWrongLengthThrows )
{
    belfem::Vector< belfem::real > tA( 2, 1.0 );
    belfem::Vector< belfem::real > tB( 2, 1.0 );
    EXPECT_THROW( belfem::cross( tA, tB ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §4.10  Norm  [semantic]
// =============================================================================

TEST( Vector, NormUnitVector )
{
    belfem::Vector< belfem::real > tV = { 1.0, 0.0, 0.0 };
    EXPECT_NEAR( belfem::norm( tV ), 1.0, tEps );
}

TEST( Vector, NormKnownVector )
{
    belfem::Vector< belfem::real > tV = { 3.0, 4.0 };
    EXPECT_NEAR( belfem::norm( tV ), 5.0, tEps );
}

TEST( Vector, NormZeroVector )
{
    belfem::Vector< belfem::real > tV = { 0.0, 0.0, 0.0 };
    EXPECT_NEAR( belfem::norm( tV ), 0.0, tEps );
}

// =============================================================================
// §4.11  Sum, Min, Max  [semantic]
// =============================================================================

TEST( Vector, SumVector )
{
    belfem::Vector< belfem::real > tV = { 1.0, 2.0, 3.0, 4.0 };
    EXPECT_NEAR( belfem::sum( tV ), 10.0, tEps );
}

TEST( Vector, MinVector )
{
    belfem::Vector< belfem::real > tV = { 3.0, 1.0, 4.0, 1.0, 5.0 };
    EXPECT_NEAR( belfem::min( tV ), 1.0, tEps );
}

TEST( Vector, MaxVector )
{
    belfem::Vector< belfem::real > tV = { 3.0, 1.0, 4.0, 1.0, 5.0 };
    EXPECT_NEAR( belfem::max( tV ), 5.0, tEps );
}

// =============================================================================
// §4.12  Linspace  [semantic]
// =============================================================================

TEST( Vector, LinspaceEndpoints )
{
    belfem::Vector< belfem::real > tV = belfem::linspace( 0.0, 1.0, 5 );
    EXPECT_NEAR( tV( 0 ), 0.0, tEps );
    EXPECT_NEAR( tV( tV.length() - 1 ), 1.0, tEps );
}

TEST( Vector, LinspaceLength )
{
    belfem::Vector< belfem::real > tV = belfem::linspace( 0.0, 10.0, 11 );
    EXPECT_EQ( tV.length(), 11u );
}

TEST( Vector, LinspaceUniformSpacing )
{
    belfem::Vector< belfem::real > tV = belfem::linspace( 0.0, 1.0, 5 );
    belfem::real tExpectedStep = 0.25;
    for( size_t i = 1; i < tV.length(); ++i )
    {
        EXPECT_NEAR( tV( i ) - tV( i - 1 ), tExpectedStep, tEps );
    }
}

// =============================================================================
// §4.13  Sort, Unique, Reverse  [semantic]
// =============================================================================

TEST( Vector, SortAscending )
{
    belfem::Vector< belfem::real > tV = { 3.0, 1.0, 4.0, 1.0, 5.0 };
    belfem::sort( tV );
    for( size_t i = 1; i < tV.length(); ++i )
    {
        EXPECT_LE( tV( i - 1 ), tV( i ) );
    }
}

TEST( Vector, UniqueRemovesDuplicates )
{
    belfem::Vector< belfem::real > tV = { 3.0, 1.0, 3.0, 2.0, 1.0 };
    belfem::unique( tV );

    // unique sorts and removes duplicates → {1, 2, 3}
    EXPECT_EQ( tV.length(), 3u );
    EXPECT_NEAR( tV( 0 ), 1.0, tEps );
    EXPECT_NEAR( tV( 1 ), 2.0, tEps );
    EXPECT_NEAR( tV( 2 ), 3.0, tEps );
}

TEST( Vector, ReverseFlipsOrder )
{
    belfem::Vector< belfem::real > tV = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tR( belfem::reverse( tV ) );
    EXPECT_NEAR( tR( 0 ), 3.0, tEps );
    EXPECT_NEAR( tR( 1 ), 2.0, tEps );
    EXPECT_NEAR( tR( 2 ), 1.0, tEps );
}

// =============================================================================
// §4.14  Append and Combine  [semantic]
// =============================================================================

TEST( Vector, AppendConcatenates )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0 };
    belfem::Vector< belfem::real > tB = { 3.0, 4.0 };
    belfem::append( tA, tB );
    EXPECT_EQ( tA.length(), 4u );
    EXPECT_NEAR( tA( 0 ), 1.0, tEps );
    EXPECT_NEAR( tA( 1 ), 2.0, tEps );
    EXPECT_NEAR( tA( 2 ), 3.0, tEps );
    EXPECT_NEAR( tA( 3 ), 4.0, tEps );
}

TEST( Vector, AppendConstOverload )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0 };
    const belfem::Vector< belfem::real > tB = { 3.0, 4.0 };
    belfem::append( tA, tB );
    EXPECT_EQ( tA.length(), 4u );
    EXPECT_NEAR( tA( 3 ), 4.0, tEps );
}

TEST( Vector, Combine2Vectors )
{
    belfem::Vector< belfem::real > tA = { 1.0, 2.0 };
    belfem::Vector< belfem::real > tB = { 3.0, 4.0 };
    belfem::Vector< belfem::real > tC;
    belfem::combine( tA, tB, tC );
    EXPECT_EQ( tC.length(), 4u );
    EXPECT_NEAR( tC( 0 ), 1.0, tEps );
    EXPECT_NEAR( tC( 3 ), 4.0, tEps );
}

TEST( Vector, Combine3Vectors )
{
    belfem::Vector< belfem::real > tA = { 1.0 };
    belfem::Vector< belfem::real > tB = { 2.0 };
    belfem::Vector< belfem::real > tC = { 3.0 };
    belfem::Vector< belfem::real > tD;
    belfem::combine( tA, tB, tC, tD );
    EXPECT_EQ( tD.length(), 3u );
    EXPECT_NEAR( tD( 0 ), 1.0, tEps );
    EXPECT_NEAR( tD( 1 ), 2.0, tEps );
    EXPECT_NEAR( tD( 2 ), 3.0, tEps );
}

TEST( Vector, Combine4Vectors )
{
    belfem::Vector< belfem::real > tA = { 1.0 };
    belfem::Vector< belfem::real > tB = { 2.0 };
    belfem::Vector< belfem::real > tC = { 3.0 };
    belfem::Vector< belfem::real > tD = { 4.0 };
    belfem::Vector< belfem::real > tE;
    belfem::combine( tA, tB, tC, tD, tE );
    EXPECT_EQ( tE.length(), 4u );
    EXPECT_NEAR( tE( 3 ), 4.0, tEps );
}

// =============================================================================
// Missing test rows (Codex finding #4)
// =============================================================================

TEST( Vector, ScalarPlusVector )
{
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
    // expression template requires explicit evaluation into Vector
    belfem::Vector< belfem::real > tResult( 10.0 + tVec );

    EXPECT_NEAR( tResult( 0 ), 11.0, tEps );
    EXPECT_NEAR( tResult( 1 ), 12.0, tEps );
    EXPECT_NEAR( tResult( 2 ), 13.0, tEps );
}

TEST( Vector, VectorMinusScalar )
{
    belfem::Vector< belfem::real > tVec = { 10.0, 20.0, 30.0 };
    belfem::Vector< belfem::real > tResult( tVec - 5.0 );

    EXPECT_NEAR( tResult( 0 ), 5.0,  tEps );
    EXPECT_NEAR( tResult( 1 ), 15.0, tEps );
    EXPECT_NEAR( tResult( 2 ), 25.0, tEps );
}

TEST( Vector, ScalarMinusVector )
{
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
    belfem::Vector< belfem::real > tResult( 10.0 - tVec );

    EXPECT_NEAR( tResult( 0 ), 9.0, tEps );
    EXPECT_NEAR( tResult( 1 ), 8.0, tEps );
    EXPECT_NEAR( tResult( 2 ), 7.0, tEps );
}

TEST( Vector, LinspaceOutputOverload )
{
    belfem::Vector< belfem::real > tOut;
    belfem::linspace( 0.0, 1.0, 5, tOut );

    EXPECT_EQ( tOut.length(), 5u );
    EXPECT_NEAR( tOut( 0 ), 0.0,  tEps );
    EXPECT_NEAR( tOut( 4 ), 1.0,  tEps );
    EXPECT_NEAR( tOut( 2 ), 0.5,  tEps );
}

TEST( Vector, LinearCombination )
{
    belfem::Vector< belfem::real > tA = { 1.0, 0.0 };
    belfem::Vector< belfem::real > tB = { 0.0, 1.0 };

    // 3a + 2b = {3, 2}
    belfem::Vector< belfem::real > tResult( 3.0 * tA + 2.0 * tB );

    EXPECT_NEAR( tResult( 0 ), 3.0, tEps );
    EXPECT_NEAR( tResult( 1 ), 2.0, tEps );
}
