/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Values and shapes of the expression overloads of Matrix and Vector: compound assignment
 * and construction from backend expressions, the row-operand orientation contract, aliasing,
 * complex operands. Operands are initialized independently and the expected result is
 * computed before anything is mutated. See todo/linalg_wrapper_hardening_plan.md, gap #33.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <complex>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace
{
    using belfem::real;
    using belfem::uint;
    using belfem::cplx;
    using belfem::Matrix;
    using belfem::Vector;

    constexpr real gTol = 1e-13 ;

    real a( const uint i, const uint j ) { return 1.0 + static_cast< real >( i ) + 2.0 * static_cast< real >( j ); }
    real b( const uint i, const uint j ) { return 3.0 + static_cast< real >( i ) - static_cast< real >( j ); }
    real x( const uint i ) { return 0.5 + static_cast< real >( i ); }
    real z( const uint i ) { return 2.0 - 0.25 * static_cast< real >( i ); }

    void
    fill( Matrix< real > & aM, real ( * aF )( const uint, const uint ) )
    {
        for ( uint j = 0; j < aM.n_cols(); ++j )
            for ( uint i = 0; i < aM.n_rows(); ++i )
                aM( i, j ) = aF( i, j );
    }
    void
    fill( Vector< real > & aV, real ( * aF )( const uint ) )
    {
        for ( uint i = 0; i < aV.length(); ++i ) aV( i ) = aF( i );
    }
}

// -----------------------------------------------------------------------------
// Matrix
// -----------------------------------------------------------------------------

TEST( OverloadSets, MatrixCompoundAssignmentFromExpression )
{
    const uint n = 5 ;
    Matrix< real > tA( n, n ), tB( n, n );
    fill( tA, a ); fill( tB, b );
    tA += tB * 0.5 ;
    for ( uint j = 0; j < n; ++j ) for ( uint i = 0; i < n; ++i )
        EXPECT_NEAR( tA( i, j ), a( i, j ) + 0.5 * b( i, j ), gTol );
    tA -= belfem::trans( tB );
    for ( uint j = 0; j < n; ++j ) for ( uint i = 0; i < n; ++i )
        EXPECT_NEAR( tA( i, j ), a( i, j ) + 0.5 * b( i, j ) - b( j, i ), gTol );
    ASSERT_EQ( tA.n_rows(), n ); ASSERT_EQ( tA.n_cols(), n );
}

TEST( OverloadSets, MatrixTimesEqualsExpressionIsAMatrixProduct )
{
    const uint n = 4 ;
    Matrix< real > tA( n, n ), tB( n, n ), tExpect( n, n, 0.0 );
    fill( tA, a ); fill( tB, b );
    for ( uint j = 0; j < n; ++j ) for ( uint i = 0; i < n; ++i ) for ( uint k = 0; k < n; ++k )
        tExpect( i, j ) += a( i, k ) * ( 2.0 * b( k, j ) );
    tA *= tB * 2.0 ;
    for ( uint j = 0; j < n; ++j ) for ( uint i = 0; i < n; ++i )
        EXPECT_NEAR( tA( i, j ), tExpect( i, j ), gTol * 100.0 );
}

TEST( OverloadSets, MatrixAliasedCompoundAssignment )
{
    const uint n = 6 ;
    Matrix< real > tA( n, n );
    fill( tA, a );
    tA += tA * 2.0 ;
    for ( uint j = 0; j < n; ++j ) for ( uint i = 0; i < n; ++i )
        EXPECT_NEAR( tA( i, j ), 3.0 * a( i, j ), gTol );
}

TEST( OverloadSets, MatrixConstructionFromExpression )
{
    const uint n = 3 ;
    Matrix< real > tB( n, n + 1 );
    fill( tB, b );
    Matrix< real > tD( tB * 0.5 );                    // direct-initialization
    Matrix< real > tE = tB * 0.5 + tB ;               // copy-initialization, legal on both backends
    Matrix< real > tT( belfem::trans( tB ) );
    ASSERT_EQ( tD.n_rows(), n ); ASSERT_EQ( tD.n_cols(), n + 1 );
    ASSERT_EQ( tT.n_rows(), n + 1 ); ASSERT_EQ( tT.n_cols(), n );
    for ( uint j = 0; j < n + 1; ++j ) for ( uint i = 0; i < n; ++i )
    {
        EXPECT_NEAR( tD( i, j ), 0.5 * b( i, j ), gTol );
        EXPECT_NEAR( tE( i, j ), 1.5 * b( i, j ), gTol );
        EXPECT_NEAR( tT( j, i ), b( i, j ), gTol );
    }
}

TEST( OverloadSets, IntegerMatrixCompoundAssignment )
{
    Matrix< int > tI( 2, 2, 1 ), tJ( 2, 2, 3 );
    tI += tJ * 2 ;
    EXPECT_EQ( tI( 0, 0 ), 7 ); EXPECT_EQ( tI( 1, 1 ), 7 );
}

// -----------------------------------------------------------------------------
// Vector
// -----------------------------------------------------------------------------

TEST( OverloadSets, VectorCompoundAssignmentFromExpression )
{
    const uint n = 7 ;
    Vector< real > tY( n, 1.0 ), tX( n ), tZ( n );
    fill( tX, x ); fill( tZ, z );
    tY += tX * 0.5 ;
    tY -= tX + tZ ;
    tY += tZ * 2.0 ;
    for ( uint i = 0; i < n; ++i )
        EXPECT_NEAR( tY( i ), 1.0 + 0.5 * x( i ) - x( i ) - z( i ) + 2.0 * z( i ), gTol );
    ASSERT_EQ( tY.length(), n );
}

TEST( OverloadSets, VectorCompoundAssignmentFromMatrixVectorProduct )
{
    const uint n = 4 ;
    Matrix< real > tA( n, n ); Vector< real > tX( n ), tY( n, 0.0 );
    fill( tA, a ); fill( tX, x );
    tY += tA * tX ;
    tY -= tA * tX * 0.5 ;
    for ( uint i = 0; i < n; ++i )
    {
        real tExpect = 0.0 ;
        for ( uint k = 0; k < n; ++k ) tExpect += 0.5 * a( i, k ) * x( k );
        EXPECT_NEAR( tY( i ), tExpect, gTol * 10.0 );
    }
}

TEST( OverloadSets, VectorColumnAndRowOperands )
{
    // a column view adds as is; a row view is transposed into the column on both backends
    const uint n = 5 ;
    Matrix< real > tA( n, n ); Vector< real > tY( n, 1.0 ), tW( n, 1.0 );
    fill( tA, a );
    tY += tA.col( 2 );
    tW += tA.row( 3 );
    for ( uint i = 0; i < n; ++i )
    {
        EXPECT_NEAR( tY( i ), 1.0 + a( i, 2 ), gTol );
        EXPECT_NEAR( tW( i ), 1.0 + a( 3, i ), gTol );
    }
    tW -= tA.row( 3 );
    for ( uint i = 0; i < n; ++i ) EXPECT_NEAR( tW( i ), 1.0, gTol );
    // a scaled row expression follows the same contract
    tW += tA.row( 1 ) * 2.0 ;
    for ( uint i = 0; i < n; ++i ) EXPECT_NEAR( tW( i ), 1.0 + 2.0 * a( 1, i ), gTol );
    // assignment from a row view keeps the column shape
    Vector< real > tR( 2, 0.0 );
    tR = tA.row( 0 );
    ASSERT_EQ( tR.length(), n );
    for ( uint i = 0; i < n; ++i ) EXPECT_NEAR( tR( i ), a( 0, i ), gTol );
}

TEST( OverloadSets, VectorTransposedColumnOperand )
{
    const uint n = 6 ;
    Vector< real > tY( n, 1.0 ), tX( n );
    fill( tX, x );
    tY += belfem::trans( tX.vector_data() );
    for ( uint i = 0; i < n; ++i ) EXPECT_NEAR( tY( i ), 1.0 + x( i ), gTol );
}

TEST( OverloadSets, VectorConstructionFromExpression )
{
    const uint n = 5 ;
    Vector< real > tX( n ), tZ( n );
    fill( tX, x ); fill( tZ, z );
    Vector< real > tW( tX * 0.5 );            // direct-initialization
    Vector< real > tV = tX + tZ ;             // copy-initialization
    ASSERT_EQ( tW.length(), n ); ASSERT_EQ( tV.length(), n );
    for ( uint i = 0; i < n; ++i )
    {
        EXPECT_NEAR( tW( i ), 0.5 * x( i ), gTol );
        EXPECT_NEAR( tV( i ), x( i ) + z( i ), gTol );
    }
}

TEST( OverloadSets, VectorConstructionFromTransposeAndRowViewIsAColumn )
{
    // the expression constructor orients: trans( x ) and A.row( k ) build n x 1 vectors
    const uint n = 5 ;
    Vector< real > tX( n ); Matrix< real > tA( n, n + 2 );
    fill( tX, x ); fill( tA, a );
    Vector< real > tT( belfem::trans( tX.vector_data() ) );
    Vector< real > tR( tA.row( 2 ) );
    Vector< real > tC( tA.col( 1 ) );
    ASSERT_EQ( tT.length(), n );
    ASSERT_EQ( tR.length(), n + 2 );
    ASSERT_EQ( tC.length(), n );
    for ( uint i = 0; i < n; ++i ) EXPECT_NEAR( tT( i ), x( i ), gTol );
    for ( uint j = 0; j < n + 2; ++j ) EXPECT_NEAR( tR( j ), a( 2, j ), gTol );
    for ( uint i = 0; i < n; ++i ) EXPECT_NEAR( tC( i ), a( i, 1 ), gTol );
}

TEST( OverloadSets, ScalarOperandsOfAnotherTypeReachTheScalarOverloads )
{
    // an int scalar must not be taken for an expression ( regression: Armadillo Proxy< int > )
    Vector< real > tY( 3, 1.0 );
    tY += 1 ;
    tY -= 2 ;
    tY *= 2 ;
    for ( uint i = 0; i < 3; ++i ) EXPECT_NEAR( tY( i ), 0.0, gTol );
    tY = 4 ;
    for ( uint i = 0; i < 3; ++i ) EXPECT_NEAR( tY( i ), 4.0, gTol );
}

TEST( OverloadSets, DotOfTwoBackendColumnViews )
{
    const uint n = 4 ;
    Matrix< real > tA( n, n );
    fill( tA, a );
    real tExpect = 0.0 ;
    for ( uint i = 0; i < n; ++i ) tExpect += a( i, 0 ) * a( i, 1 );
    EXPECT_NEAR( belfem::dot( tA.col( 0 ), tA.col( 1 ) ), tExpect, gTol * 10.0 );
}

TEST( OverloadSets, ComplexRowOperandIsNotConjugated )
{
    // the row-to-column transposition is a plain transpose on both backends
    Matrix< cplx > tZ( 2, 3 );
    for ( uint j = 0; j < 3; ++j ) for ( uint i = 0; i < 2; ++i )
        tZ( i, j ) = cplx( 1.0 + i, 2.0 + j );
    Vector< cplx > tW( 3, cplx( 0.0, 0.0 ) );
    tW += tZ.row( 1 );
    Vector< cplx > tC( tZ.row( 0 ) );
    for ( uint j = 0; j < 3; ++j )
    {
        EXPECT_NEAR( tW( j ).real(), 2.0, gTol );
        EXPECT_NEAR( tW( j ).imag(), 2.0 + j, gTol );       // +, not -
        EXPECT_NEAR( tC( j ).imag(), 2.0 + j, gTol );
    }
}

#ifdef BELFEM_BLAZE
TEST( OverloadSets, BlazeRowMajorMatrixScalarProductAddsItsFirstColumn )
{
    // the DMatScalarMultExpr< ET, T, rowMajor > overload adds column 0 of a row-major
    // matrix-scalar product; with BELFEM's column-major storage such an expression arises
    // from a transposed matrix, so column 0 of trans( A ) * s is s * A.row( 0 )
    const uint n = 4 ;
    Matrix< real > tA( n, n ); Vector< real > tY( n, 0.0 );
    fill( tA, a );
    tY += belfem::trans( tA ) * 2.0 ;
    for ( uint i = 0; i < n; ++i ) EXPECT_NEAR( tY( i ), 2.0 * a( 0, i ), gTol );
}
#endif
