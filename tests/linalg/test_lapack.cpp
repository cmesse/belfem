/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for the unified backend-agnostic LAPACK interface in
 * src/linalg/lapack: gesv, posv, gels, gemm, getrf/getri and gesvd.
 *
 * Every test runs as a typed suite over the four LAPACK datatypes
 * ( float, double, std::complex<float>, std::complex<double> ), so each
 * covers the s/d/c/z flavor of its routine. Solvers are verified against
 * manually computed references; under Blaze, odd matrix dimensions ensure
 * the SIMD-padded column stride ( spacing > n_rows ) is exercised, which
 * is what lapack::leading_dimension() exists for.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <complex>

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_gesv.hpp"
#include "fn_posv.hpp"
#include "fn_gels.hpp"
#include "fn_gemm.hpp"
#include "fn_getrf.hpp"
#include "fn_getri.hpp"
#include "fn_gesvd.hpp"
#include "fn_geev.hpp"
#include "fn_gees.hpp"

namespace
{
    using belfem::real ;
    using belfem::uint ;
    using belfem::int_t ;

    // scalar factory and per-type tolerance; real types drop the
    // imaginary part, so one set of test data serves all four flavors
    template< typename T >
    struct Scalar
    {
        static T    mk( const real aRe, const real ) { return static_cast< T >( aRe ); }
        static T    conj( const T & aVal )           { return aVal; }
        static real tol()                            { return 1e-11; }
    };

    template<>
    struct Scalar< float >
    {
        static float mk( const real aRe, const real ) { return static_cast< float >( aRe ); }
        static float conj( const float & aVal )       { return aVal; }
        static real  tol()                            { return 5e-4; }
    };

    template< typename R >
    struct Scalar< std::complex< R > >
    {
        static std::complex< R >
        mk( const real aRe, const real aIm )
        {
            return std::complex< R >( static_cast< R >( aRe ), static_cast< R >( aIm ) );
        }

        static std::complex< R >
        conj( const std::complex< R > & aVal )
        {
            return std::conj( aVal );
        }

        static real tol() { return sizeof( R ) == 4 ? 5e-4 : 1e-11; }
    };

    // well-conditioned 3x3 test matrix; hermitian positive definite when
    // aHermitian is set ( real diagonal, conjugate-symmetric off-diagonals )
    template< typename T >
    void
    fill_matrix( belfem::Matrix< T > & aA, const bool aHermitian )
    {
        using S = Scalar< T >;
        aA.set_size( 3, 3 );
        aA( 0, 0 ) = S::mk( 4.0, 0.0 );
        aA( 1, 1 ) = S::mk( 5.0, 0.0 );
        aA( 2, 2 ) = S::mk( 6.0, 0.0 );
        aA( 0, 1 ) = S::mk( 1.0,  0.5 );
        aA( 0, 2 ) = S::mk( 0.5, -0.25 );
        aA( 1, 2 ) = S::mk( 1.5,  0.75 );
        if ( aHermitian )
        {
            aA( 1, 0 ) = S::conj( aA( 0, 1 ) );
            aA( 2, 0 ) = S::conj( aA( 0, 2 ) );
            aA( 2, 1 ) = S::conj( aA( 1, 2 ) );
        }
        else
        {
            aA( 1, 0 ) = S::mk( 2.0, -1.0 );
            aA( 2, 0 ) = S::mk( 0.25, 0.5 );
            aA( 2, 1 ) = S::mk( 0.75, -0.5 );
        }
    }

    // reference solution vector
    template< typename T >
    void
    fill_solution( belfem::Vector< T > & aX )
    {
        using S = Scalar< T >;
        aX.set_size( 3 );
        aX( 0 ) = S::mk(  1.0, -1.0 );
        aX( 1 ) = S::mk( -2.0,  0.5 );
        aX( 2 ) = S::mk(  3.0,  2.0 );
    }

    // b = A * x, computed manually so the test does not depend on gemm
    template< typename T >
    void
    matvec( const belfem::Matrix< T > & aA, const belfem::Vector< T > & aX,
            belfem::Vector< T > & aB )
    {
        aB.set_size( aA.n_rows() );
        for ( uint i = 0; i < aA.n_rows(); ++i )
        {
            T tSum = Scalar< T >::mk( 0.0, 0.0 );
            for ( uint j = 0; j < aA.n_cols(); ++j )
            {
                tSum += aA( i, j ) * aX( j );
            }
            aB( i ) = tSum;
        }
    }

    // the Work length fn_gesvd.hpp considers sufficient, i.e. the length at
    // which the wrapper skips the workspace query and takes the
    // caller-provided branch: lwork entries of T in the head ( two reals per
    // entry for the complex flavors ) plus the complex-only 5*min(m,n) real
    // rwork tail.
    //
    // This mirrors the wrapper's own formula on purpose. It is not a second
    // opinion on what LAPACK needs — LAPACK gives that verdict itself: if the
    // formula understates the true minimum, ?gesvd returns info = -13 and the
    // tests below fail on the info check.
    //
    // It is SUFFICIENT, not NECESSARY, and the difference matters: this is the
    // reference-LAPACK minimum, and a vendor implementation may need less. MKL's
    // dgesvd workspace query returns 7 for a 4x3 'A','A' problem where this
    // formula gives 15, and then completes happily with lwork = 7. So a buffer
    // of this size always works, but a QUERIED size may legitimately be smaller
    // and must never be asserted against it.
    template< typename T >
    int_t
    gesvd_min_work( const int_t aM, const int_t aN )
    {
        constexpr int_t tRealsPerT =
            std::is_same< T, belfem::lapack::real_t< T > >::value ? 1 : 2 ;

        const int_t mn = std::min( aM, aN );

        const int_t lwork = std::max< int_t >( 1, tRealsPerT == 1 ?
            std::max( 3 * mn + std::max( aM, aN ), 5 * mn ) :
            2 * mn + std::max( aM, aN ) );

        return tRealsPerT * lwork + ( tRealsPerT == 1 ? 0 : 5 * mn );
    }

    // deterministic, well conditioned m x n test matrix
    template< typename T >
    void
    fill_svd_matrix( belfem::Matrix< T > & aA, const uint aM, const uint aN )
    {
        aA.set_size( aM, aN );
        for ( uint i = 0; i < aM; ++i )
            for ( uint j = 0; j < aN; ++j )
                aA( i, j ) = Scalar< T >::mk(
                    1.0 / ( 1.0 + i + j ) + ( i == j ? 2.0 : 0.0 ),
                    0.1 * i - 0.2 * j );
    }

    // max | A - U * diag(S) * VT | over the first min(m,n) modes
    template< typename T >
    real
    svd_error(
        const belfem::Matrix< T > & aA,
        const belfem::Vector< belfem::lapack::real_t< T > > & aS,
        const belfem::Matrix< T > & aU,
        const belfem::Matrix< T > & aVT )
    {
        real tErr = 0.0;
        for ( uint i = 0; i < aA.n_rows(); ++i )
            for ( uint j = 0; j < aA.n_cols(); ++j )
            {
                T tSum = Scalar< T >::mk( 0.0, 0.0 );
                for ( uint l = 0; l < aS.length(); ++l )
                {
                    tSum += aU( i, l ) * aS( l ) * aVT( l, j );
                }
                tErr = std::max( tErr, ( real ) std::abs( tSum - aA( i, j ) ) );
            }
        return tErr;
    }

    template< typename T >
    real
    max_error( const belfem::Vector< T > & aValue, const belfem::Vector< T > & aRef )
    {
        real tErr = 0.0;
        for ( uint i = 0; i < aRef.length(); ++i )
        {
            tErr = std::max( tErr, ( real ) std::abs( aValue( i ) - aRef( i ) ) );
        }
        return tErr;
    }
}

// =============================================================================
// typed suite over the four LAPACK datatypes
// =============================================================================

template< typename T >
class LapackTest : public ::testing::Test {};

using LapackTypes = ::testing::Types<
    float, double, std::complex< float >, std::complex< double > >;

TYPED_TEST_SUITE( LapackTest, LapackTypes );

// =============================================================================
// gesv
// =============================================================================

TYPED_TEST( LapackTest, GesvVectorRhs )
{
    belfem::Matrix< TypeParam > tA;
    belfem::Vector< TypeParam > tX, tB;
    fill_matrix( tA, false );
    fill_solution( tX );
    matvec( tA, tX, tB );

    belfem::Vector< belfem::int_t > tP( 3 );
    belfem::gesv( tA, tB, tP );   // tB now holds the solution

    EXPECT_LT( max_error( tB, tX ), Scalar< TypeParam >::tol() );
}

TYPED_TEST( LapackTest, GesvMatrixRhs )
{
    using S = Scalar< TypeParam >;

    belfem::Matrix< TypeParam > tA;
    fill_matrix( tA, false );

    // two right hand sides, columns are A * x1 and A * x2
    belfem::Vector< TypeParam > tX1, tX2, tB1, tB2;
    fill_solution( tX1 );
    tX2.set_size( 3 );
    tX2( 0 ) = S::mk( 0.5, 1.0 );
    tX2( 1 ) = S::mk( 4.0, -2.0 );
    tX2( 2 ) = S::mk( -1.0, 0.25 );
    matvec( tA, tX1, tB1 );
    matvec( tA, tX2, tB2 );

    belfem::Matrix< TypeParam > tB( 3, 2 );
    for ( uint i = 0; i < 3; ++i )
    {
        tB( i, 0 ) = tB1( i );
        tB( i, 1 ) = tB2( i );
    }

    belfem::Vector< belfem::int_t > tP( 3 );
    belfem::gesv( tA, tB, tP );

    real tErr = 0.0;
    for ( uint i = 0; i < 3; ++i )
    {
        tErr = std::max( tErr, ( real ) std::abs( tB( i, 0 ) - tX1( i ) ) );
        tErr = std::max( tErr, ( real ) std::abs( tB( i, 1 ) - tX2( i ) ) );
    }
    EXPECT_LT( tErr, S::tol() );
}

TYPED_TEST( LapackTest, GesvSingularNoAbort )
{
    using S = Scalar< TypeParam >;

    // singular system ( two identical rows ): with AbortOnError = false the
    // wrapper must report info != 0 instead of aborting, so that iterative
    // callers ( e.g. Anderson stabilization ) can react to the failure
    belfem::Matrix< TypeParam > tA( 3, 3 );
    fill_matrix( tA, false );
    for ( uint j = 0; j < 3; ++j )
    {
        tA( 2, j ) = tA( 1, j );
    }

    belfem::Vector< TypeParam > tB( 3 );
    tB( 0 ) = S::mk( 1.0, 0.0 );
    tB( 1 ) = S::mk( 2.0, 0.0 );
    tB( 2 ) = S::mk( 3.0, 0.0 );

    belfem::Vector< belfem::int_t > tP( 3 );
    belfem::int_t tInfo = belfem::gesv( tA, tB, tP, false );

    EXPECT_NE( tInfo, 0 );
}

// =============================================================================
// posv ( symmetric / hermitian positive definite )
// =============================================================================

TYPED_TEST( LapackTest, PosvVectorRhs )
{
    belfem::Matrix< TypeParam > tA;
    belfem::Vector< TypeParam > tX, tB;
    fill_matrix( tA, true );
    fill_solution( tX );
    matvec( tA, tX, tB );

    belfem::posv( tA, tB );

    EXPECT_LT( max_error( tB, tX ), Scalar< TypeParam >::tol() );
}

TYPED_TEST( LapackTest, PosvMatrixRhs )
{
    using S = Scalar< TypeParam >;

    belfem::Matrix< TypeParam > tA;
    fill_matrix( tA, true );

    belfem::Vector< TypeParam > tX1, tX2, tB1, tB2;
    fill_solution( tX1 );
    tX2.set_size( 3 );
    tX2( 0 ) = S::mk( -0.5, 0.5 );
    tX2( 1 ) = S::mk( 2.0, 1.5 );
    tX2( 2 ) = S::mk( 1.0, -1.0 );
    matvec( tA, tX1, tB1 );
    matvec( tA, tX2, tB2 );

    belfem::Matrix< TypeParam > tB( 3, 2 );
    for ( uint i = 0; i < 3; ++i )
    {
        tB( i, 0 ) = tB1( i );
        tB( i, 1 ) = tB2( i );
    }

    belfem::posv( tA, tB );

    real tErr = 0.0;
    for ( uint i = 0; i < 3; ++i )
    {
        tErr = std::max( tErr, ( real ) std::abs( tB( i, 0 ) - tX1( i ) ) );
        tErr = std::max( tErr, ( real ) std::abs( tB( i, 1 ) - tX2( i ) ) );
    }
    EXPECT_LT( tErr, S::tol() );
}

// =============================================================================
// gels ( least squares )
// =============================================================================

TYPED_TEST( LapackTest, GelsOverdetermined )
{
    using S = Scalar< TypeParam >;

    // exact linear model y = c0 + c1 * x sampled at 5 points: the least
    // squares solution must reproduce the coefficients
    const TypeParam tC0 = S::mk( 2.0, -1.0 );
    const TypeParam tC1 = S::mk( 3.0,  0.5 );

    belfem::Matrix< TypeParam > tA( 5, 2 );
    belfem::Vector< TypeParam > tB( 5 );
    for ( uint i = 0; i < 5; ++i )
    {
        const TypeParam tXi = S::mk( ( real ) i, 0.0 );
        tA( i, 0 ) = S::mk( 1.0, 0.0 );
        tA( i, 1 ) = tXi;
        tB( i ) = tC0 + tC1 * tXi;
    }

    belfem::Vector< TypeParam > tWork;
    belfem::gels( tA, tB, tWork );

    real tErr = std::max( ( real ) std::abs( tB( 0 ) - tC0 ),
                          ( real ) std::abs( tB( 1 ) - tC1 ) );
    EXPECT_LT( tErr, S::tol() );
}

TYPED_TEST( LapackTest, GelsUnderdeterminedMinNorm )
{
    using S = Scalar< TypeParam >;

    // 2 x 3 system: B must be allocated with max( m, n ) = 3 entries,
    // the first m hold the right hand side on input
    belfem::Matrix< TypeParam > tA( 2, 3 );
    tA( 0, 0 ) = S::mk( 1.0, 0.25 ); tA( 0, 1 ) = S::mk( 2.0, 0.0 ); tA( 0, 2 ) = S::mk( -1.0, 0.5 );
    tA( 1, 0 ) = S::mk( 0.5, 0.0 );  tA( 1, 1 ) = S::mk( -1.0, 1.0 ); tA( 1, 2 ) = S::mk( 3.0, 0.0 );

    belfem::Matrix< TypeParam > tAorig( tA );

    belfem::Vector< TypeParam > tB( 3 );
    tB( 0 ) = S::mk( 4.0, -1.0 );
    tB( 1 ) = S::mk( 1.0,  2.0 );
    tB( 2 ) = S::mk( 0.0,  0.0 );

    belfem::Vector< TypeParam > tRhs( 2 );
    tRhs( 0 ) = tB( 0 );
    tRhs( 1 ) = tB( 1 );

    belfem::Vector< TypeParam > tWork;
    belfem::gels( tA, tB, tWork );

    // the min-norm solution must satisfy A * x = rhs exactly
    real tErr = 0.0;
    for ( uint i = 0; i < 2; ++i )
    {
        TypeParam tSum = S::mk( 0.0, 0.0 );
        for ( uint j = 0; j < 3; ++j )
        {
            tSum += tAorig( i, j ) * tB( j );
        }
        tErr = std::max( tErr, ( real ) std::abs( tSum - tRhs( i ) ) );
    }
    EXPECT_LT( tErr, S::tol() );
}

// =============================================================================
// gemm
// =============================================================================

TYPED_TEST( LapackTest, GemmPlain )
{
    using S = Scalar< TypeParam >;

    // op(A) is 2x3, op(B) is 3x4 — non-square so that a wrong m/k or a
    // transposition-dependent leading dimension cannot cancel out
    belfem::Matrix< TypeParam > tA( 2, 3 );
    belfem::Matrix< TypeParam > tB( 3, 4 );
    for ( uint i = 0; i < 2; ++i )
        for ( uint j = 0; j < 3; ++j )
            tA( i, j ) = S::mk( 1.0 + i + 0.5 * j, 0.25 * j - 0.5 * i );
    for ( uint i = 0; i < 3; ++i )
        for ( uint j = 0; j < 4; ++j )
            tB( i, j ) = S::mk( 0.3 * ( i + 1.0 ) - 0.7 * j, 0.1 * i + 0.2 * j );

    belfem::Matrix< TypeParam > tR( 2, 4 );
    for ( uint i = 0; i < 2; ++i )
        for ( uint j = 0; j < 4; ++j )
        {
            TypeParam tSum = S::mk( 0.0, 0.0 );
            for ( uint l = 0; l < 3; ++l ) tSum += tA( i, l ) * tB( l, j );
            tR( i, j ) = tSum;
        }

    belfem::Matrix< TypeParam > tC;
    belfem::gemm( tA, tB, tC );

    real tErr = 0.0;
    for ( uint i = 0; i < 2; ++i )
        for ( uint j = 0; j < 4; ++j )
            tErr = std::max( tErr, ( real ) std::abs( tC( i, j ) - tR( i, j ) ) );
    EXPECT_LT( tErr, S::tol() );
}

TYPED_TEST( LapackTest, GemmTransposed )
{
    using S = Scalar< TypeParam >;

    // same product as GemmPlain, but both operands stored transposed
    belfem::Matrix< TypeParam > tAt( 3, 2 );
    belfem::Matrix< TypeParam > tBt( 4, 3 );
    belfem::Matrix< TypeParam > tR( 2, 4 );
    {
        belfem::Matrix< TypeParam > tA( 2, 3 );
        belfem::Matrix< TypeParam > tB( 3, 4 );
        for ( uint i = 0; i < 2; ++i )
            for ( uint j = 0; j < 3; ++j )
            {
                tA( i, j ) = S::mk( 1.0 + i + 0.5 * j, 0.25 * j - 0.5 * i );
                tAt( j, i ) = tA( i, j );
            }
        for ( uint i = 0; i < 3; ++i )
            for ( uint j = 0; j < 4; ++j )
            {
                tB( i, j ) = S::mk( 0.3 * ( i + 1.0 ) - 0.7 * j, 0.1 * i + 0.2 * j );
                tBt( j, i ) = tB( i, j );
            }
        for ( uint i = 0; i < 2; ++i )
            for ( uint j = 0; j < 4; ++j )
            {
                TypeParam tSum = S::mk( 0.0, 0.0 );
                for ( uint l = 0; l < 3; ++l ) tSum += tA( i, l ) * tB( l, j );
                tR( i, j ) = tSum;
            }
    }

    belfem::Matrix< TypeParam > tC;
    belfem::gemm( tAt, tBt, tC,
                  Scalar< TypeParam >::mk( 1.0, 0.0 ),
                  Scalar< TypeParam >::mk( 0.0, 0.0 ), 'T', 'T' );

    real tErr = 0.0;
    for ( uint i = 0; i < 2; ++i )
        for ( uint j = 0; j < 4; ++j )
            tErr = std::max( tErr, ( real ) std::abs( tC( i, j ) - tR( i, j ) ) );
    EXPECT_LT( tErr, S::tol() );
}

TYPED_TEST( LapackTest, GemmAccumulate )
{
    using S = Scalar< TypeParam >;

    // C := 2 * A * B + 3 * C with C preset to the plain product R,
    // so the result must be 5 * R
    belfem::Matrix< TypeParam > tA( 2, 3 );
    belfem::Matrix< TypeParam > tB( 3, 4 );
    for ( uint i = 0; i < 2; ++i )
        for ( uint j = 0; j < 3; ++j )
            tA( i, j ) = S::mk( 1.0 + i + 0.5 * j, 0.25 * j );
    for ( uint i = 0; i < 3; ++i )
        for ( uint j = 0; j < 4; ++j )
            tB( i, j ) = S::mk( 0.3 * ( i + 1.0 ) - 0.7 * j, 0.2 * j );

    belfem::Matrix< TypeParam > tR( 2, 4 );
    for ( uint i = 0; i < 2; ++i )
        for ( uint j = 0; j < 4; ++j )
        {
            TypeParam tSum = S::mk( 0.0, 0.0 );
            for ( uint l = 0; l < 3; ++l ) tSum += tA( i, l ) * tB( l, j );
            tR( i, j ) = tSum;
        }

    belfem::Matrix< TypeParam > tC( tR );
    belfem::gemm( tA, tB, tC,
                  S::mk( 2.0, 0.0 ), S::mk( 3.0, 0.0 ) );

    real tErr = 0.0;
    for ( uint i = 0; i < 2; ++i )
        for ( uint j = 0; j < 4; ++j )
            tErr = std::max( tErr, ( real ) std::abs(
                tC( i, j ) - S::mk( 5.0, 0.0 ) * tR( i, j ) ) );
    EXPECT_LT( tErr, S::tol() );
}

// =============================================================================
// getrf + getri ( inversion round trip )
// =============================================================================

TYPED_TEST( LapackTest, GetrfGetriRoundTrip )
{
    using S = Scalar< TypeParam >;

    belfem::Matrix< TypeParam > tA;
    fill_matrix( tA, false );
    belfem::Matrix< TypeParam > tAinv( tA );

    belfem::Vector< belfem::int_t > tP;
    belfem::Vector< TypeParam > tWork;
    belfem::getrf( tAinv, tP );
    belfem::getri( tAinv, tP, tWork );

    // A * inv(A) must be the identity
    real tErr = 0.0;
    for ( uint i = 0; i < 3; ++i )
        for ( uint j = 0; j < 3; ++j )
        {
            TypeParam tSum = S::mk( 0.0, 0.0 );
            for ( uint l = 0; l < 3; ++l ) tSum += tA( i, l ) * tAinv( l, j );
            tErr = std::max( tErr, ( real ) std::abs(
                tSum - S::mk( i == j ? 1.0 : 0.0, 0.0 ) ) );
        }
    EXPECT_LT( tErr, S::tol() );
}

// =============================================================================
// gesvd
// =============================================================================

TYPED_TEST( LapackTest, GesvdReconstruct )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;

    // full SVD of a 3x2 matrix; reconstruction uses only the first
    // min(m,n) columns of U and rows of VT
    belfem::Matrix< TypeParam > tA0( 3, 2 );
    tA0( 0, 0 ) = S::mk( 3.0, 0.5 );  tA0( 0, 1 ) = S::mk( 2.0, -1.0 );
    tA0( 1, 0 ) = S::mk( 2.0, 0.0 );  tA0( 1, 1 ) = S::mk( 3.0, 0.75 );
    tA0( 2, 0 ) = S::mk( 2.0, -0.5 ); tA0( 2, 1 ) = S::mk( -2.0, 0.25 );

    belfem::Matrix< TypeParam > tA( tA0 );
    belfem::Vector< Treal > tS;
    belfem::Matrix< TypeParam > tU, tVT;
    belfem::Vector< Treal > tWork;

    belfem::gesvd( tA, tS, tU, tVT, tWork );

    // shapes for jobu = jobvt = 'A'
    EXPECT_EQ( tU.n_rows(), 3u );
    EXPECT_EQ( tU.n_cols(), 3u );
    EXPECT_EQ( tVT.n_rows(), 2u );
    EXPECT_EQ( tVT.n_cols(), 2u );
    EXPECT_EQ( tS.length(), 2u );

    // singular values are sorted descending and positive
    EXPECT_GE( ( real ) tS( 0 ), ( real ) tS( 1 ) );
    EXPECT_GT( ( real ) tS( 1 ), 0.0 );

    // A == U * diag(S) * VT
    real tErr = 0.0;
    for ( uint i = 0; i < 3; ++i )
        for ( uint j = 0; j < 2; ++j )
        {
            TypeParam tSum = S::mk( 0.0, 0.0 );
            for ( uint l = 0; l < 2; ++l )
            {
                tSum += tU( i, l ) * tS( l ) * tVT( l, j );
            }
            tErr = std::max( tErr, ( real ) std::abs( tSum - tA0( i, j ) ) );
        }
    EXPECT_LT( tErr, S::tol() );
}

TYPED_TEST( LapackTest, GesvdEconomy )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;

    belfem::Matrix< TypeParam > tA0( 4, 3 );
    for ( uint i = 0; i < 4; ++i )
        for ( uint j = 0; j < 3; ++j )
            tA0( i, j ) = S::mk( 1.0 / ( 1.0 + i + j ) + ( i == j ? 2.0 : 0.0 ),
                                 0.1 * i - 0.2 * j );

    belfem::Matrix< TypeParam > tA( tA0 );
    belfem::Vector< Treal > tS;
    belfem::Matrix< TypeParam > tU, tVT;
    belfem::Vector< Treal > tWork;

    belfem::gesvd( tA, tS, tU, tVT, tWork, 'S', 'S' );

    // economy shapes: U is m x mn, VT is mn x n
    EXPECT_EQ( tU.n_rows(), 4u );
    EXPECT_EQ( tU.n_cols(), 3u );
    EXPECT_EQ( tVT.n_rows(), 3u );
    EXPECT_EQ( tVT.n_cols(), 3u );

    real tErr = 0.0;
    for ( uint i = 0; i < 4; ++i )
        for ( uint j = 0; j < 3; ++j )
        {
            TypeParam tSum = S::mk( 0.0, 0.0 );
            for ( uint l = 0; l < 3; ++l )
            {
                tSum += tU( i, l ) * tS( l ) * tVT( l, j );
            }
            tErr = std::max( tErr, ( real ) std::abs( tSum - tA0( i, j ) ) );
        }
    EXPECT_LT( tErr, S::tol() );
}

// -----------------------------------------------------------------------------
// gesvd — the caller-provided Work branch
//
// The two tests above enter gesvd with an empty Work, so every one of them
// takes the workspace-query branch. These cover the other half of the
// wrapper: the branch that trusts the buffer it is handed and derives lwork
// from its length ( fn_gesvd.hpp, the else of the Work.length() test ).
//
// That branch is where the head/tail split can go wrong, and it is worse for
// the complex flavors: work is reinterpreted in place as two reals per T
// entry and rwork sits behind it, so an off-by-one in either the length
// arithmetic or the rwork offset overruns the buffer instead of failing
// loudly. All tests pass AbortOnError = false and assert on info, so a
// rejected lwork surfaces as a red test rather than as an abort.
// -----------------------------------------------------------------------------

TYPED_TEST( LapackTest, GesvdProvidedWorkExactMinimum )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;

    belfem::Matrix< TypeParam > tA0;
    fill_svd_matrix( tA0, 4, 3 );

    belfem::Matrix< TypeParam > tA( tA0 );
    belfem::Vector< Treal > tS;
    belfem::Matrix< TypeParam > tU, tVT;

    // exactly the length the wrapper calls sufficient — one entry less and
    // it would query instead
    const int_t tMin = gesvd_min_work< TypeParam >( 4, 3 );
    belfem::Vector< Treal > tWork( tMin, Treal( 0 ) );

    const int_t tInfo =
        belfem::gesvd( tA, tS, tU, tVT, tWork, 'A', 'A', false );

    EXPECT_EQ( tInfo, 0 );

    // the query branch resizes Work; this one must not touch it
    EXPECT_EQ( ( int_t ) tWork.length(), tMin );

    EXPECT_LT( svd_error( tA0, tS, tU, tVT ), S::tol() );
}

TYPED_TEST( LapackTest, GesvdProvidedWorkOversized )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;

    belfem::Matrix< TypeParam > tA0;
    fill_svd_matrix( tA0, 4, 3 );

    belfem::Matrix< TypeParam > tA( tA0 );
    belfem::Vector< Treal > tS;
    belfem::Matrix< TypeParam > tU, tVT;

    // a generously oversized buffer: lwork is derived from the full length,
    // so the rwork tail moves far away from where the minimum would put it.
    // The odd length also puts the complex head/tail split on a boundary the
    // exact-minimum case cannot reach ( 2 * lwork < length - rwork ).
    const int_t tLen = gesvd_min_work< TypeParam >( 4, 3 ) + 137;
    belfem::Vector< Treal > tWork( tLen, Treal( 0 ) );

    const int_t tInfo =
        belfem::gesvd( tA, tS, tU, tVT, tWork, 'A', 'A', false );

    EXPECT_EQ( tInfo, 0 );
    EXPECT_EQ( ( int_t ) tWork.length(), tLen );
    EXPECT_LT( svd_error( tA0, tS, tU, tVT ), S::tol() );
}

TYPED_TEST( LapackTest, GesvdWorkBufferReuse )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;

    belfem::Vector< Treal > tS;
    belfem::Matrix< TypeParam > tU, tVT;

    // one buffer, three decompositions — the way a caller in a loop uses it

    // first call: empty Work, so the query branch grows it
    belfem::Vector< Treal > tWork;

    belfem::Matrix< TypeParam > tA0;
    fill_svd_matrix( tA0, 4, 3 );
    belfem::Matrix< TypeParam > tA( tA0 );

    EXPECT_EQ( belfem::gesvd( tA, tS, tU, tVT, tWork, 'A', 'A', false ), 0 );
    EXPECT_LT( svd_error( tA0, tS, tU, tVT ), S::tol() );

    // the query branch sized it, flooring the vendor query at the reference
    // minimum — so on any LAPACK the buffer now passes the wrapper's own
    // length test, and every call below takes the reuse branch
    const int_t tGrown = ( int_t ) tWork.length();
    EXPECT_GE( tGrown, gesvd_min_work< TypeParam >( 4, 3 ) );

    // second call, same shape: the buffer is accepted outright and comes
    // back untouched. tA0 is untouched by the call above — only the
    // working copy tA is destroyed
    tA = tA0;

    EXPECT_EQ( belfem::gesvd( tA, tS, tU, tVT, tWork, 'A', 'A', false ), 0 );
    EXPECT_EQ( ( int_t ) tWork.length(), tGrown );
    EXPECT_LT( svd_error( tA0, tS, tU, tVT ), S::tol() );

    // third call, a SMALLER problem, with the buffer from the larger one
    // still in hand — the reason the overload exists. The reuse branch
    // never resizes, so the length holds on any LAPACK
    fill_svd_matrix( tA0, 3, 2 );
    tA = tA0;

    EXPECT_EQ( belfem::gesvd( tA, tS, tU, tVT, tWork, 'A', 'A', false ), 0 );
    EXPECT_EQ( ( int_t ) tWork.length(), tGrown );
    EXPECT_EQ( tS.length(), 2u );
    EXPECT_LT( svd_error( tA0, tS, tU, tVT ), S::tol() );
}

// The acceptance boundary: a caller buffer of EXACTLY the reference minimum
// ( fn_gesvd.hpp, the length test before the query ) is taken outright and
// never regrown, on any LAPACK. The wrapper floors its own query sizing at
// that same minimum, so wrapper-sized buffers reach the reuse branch too —
// that path is covered by GesvdWorkBufferReuse above; this case pins the
// exact-minimum edge, which the query branch cannot produce when the vendor
// optimum exceeds it.
TYPED_TEST( LapackTest, GesvdWorkBufferReuseAtReferenceMinimum )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;

    belfem::Vector< Treal > tS;
    belfem::Matrix< TypeParam > tU, tVT;

    // sized so the wrapper accepts it outright, on any LAPACK
    const int_t tLen = gesvd_min_work< TypeParam >( 4, 3 );
    belfem::Vector< Treal > tWork( tLen, Treal( 0 ) );

    belfem::Matrix< TypeParam > tA0;
    fill_svd_matrix( tA0, 4, 3 );

    for ( uint k = 0; k < 3; ++k )
    {
        belfem::Matrix< TypeParam > tA( tA0 );

        EXPECT_EQ( belfem::gesvd( tA, tS, tU, tVT, tWork, 'A', 'A', false ), 0 );

        // the whole point: not touched, not regrown, on any iteration
        EXPECT_EQ( ( int_t ) tWork.length(), tLen );
        EXPECT_LT( svd_error( tA0, tS, tU, tVT ), S::tol() );
    }
}

TYPED_TEST( LapackTest, GesvdSingularValuesOnly )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;

    // jobu = jobvt = 'N': U and VT are never referenced and stay unsized,
    // so the wrapper passes ldu = ldvt = 1 and the data() of empty matrices
    belfem::Matrix< TypeParam > tA0;
    fill_svd_matrix( tA0, 4, 3 );

    belfem::Matrix< TypeParam > tA( tA0 );
    belfem::Vector< Treal > tS;
    belfem::Matrix< TypeParam > tU, tVT;
    belfem::Vector< Treal > tWork;

    EXPECT_EQ( belfem::gesvd( tA, tS, tU, tVT, tWork, 'N', 'N', false ), 0 );

    EXPECT_EQ( tS.length(), 3u );
    EXPECT_EQ( tU.n_rows(), 0u );
    EXPECT_EQ( tVT.n_rows(), 0u );

    // the singular values must be the ones the full decomposition gives
    belfem::Matrix< TypeParam > tAfull( tA0 );
    belfem::Vector< Treal > tSfull;
    belfem::Matrix< TypeParam > tUfull, tVTfull;
    belfem::Vector< Treal > tWorkFull;

    EXPECT_EQ( belfem::gesvd(
        tAfull, tSfull, tUfull, tVTfull, tWorkFull, 'A', 'A', false ), 0 );

    real tErr = 0.0;
    for ( uint k = 0; k < 3; ++k )
    {
        tErr = std::max( tErr, ( real ) std::abs( tS( k ) - tSfull( k ) ) );
    }
    EXPECT_LT( tErr, S::tol() );
}

// =============================================================================
// geev ( general eigenvalue problem )
// =============================================================================

TYPED_TEST( LapackTest, GeevRealSpectrum )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;
    typedef std::complex< Treal > Tcplx ;

    // upper triangular: eigenvalues are the diagonal { 1, 2, 3 }
    belfem::Matrix< TypeParam > tA( 3, 3 );
    tA( 0, 0 ) = S::mk( 1.0, 0.0 ); tA( 0, 1 ) = S::mk( 4.0, 0.0 ); tA( 0, 2 ) = S::mk( 5.0, 0.0 );
    tA( 1, 0 ) = S::mk( 0.0, 0.0 ); tA( 1, 1 ) = S::mk( 2.0, 0.0 ); tA( 1, 2 ) = S::mk( 6.0, 0.0 );
    tA( 2, 0 ) = S::mk( 0.0, 0.0 ); tA( 2, 1 ) = S::mk( 0.0, 0.0 ); tA( 2, 2 ) = S::mk( 3.0, 0.0 );

    belfem::Vector< Tcplx > tW;
    belfem::Vector< Treal > tWork;

    belfem::int_t tInfo = belfem::geev( tA, tW, tWork );
    EXPECT_EQ( tInfo, 0 );

    // eigenvalues come unordered: sort by real part
    Treal tRe[ 3 ] = { tW( 0 ).real(), tW( 1 ).real(), tW( 2 ).real() };
    std::sort( tRe, tRe + 3 );

    real tErr = 0.0;
    for ( uint k = 0; k < 3; ++k )
    {
        tErr = std::max( tErr, ( real ) std::abs( tRe[ k ] - Treal( k + 1 ) ) );
        tErr = std::max( tErr, ( real ) std::abs( tW( k ).imag() ) );
    }
    EXPECT_LT( tErr, S::tol() );
}

TYPED_TEST( LapackTest, GeevComplexPairsResidual )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;
    typedef std::complex< Treal > Tcplx ;

    // rotation block + decoupled row: eigenvalues { +i, -i, 2 }. For the
    // real flavors this exercises the LAPACK conjugate-pair packing of the
    // eigenvectors across two consecutive columns
    belfem::Matrix< TypeParam > tA0( 3, 3 );
    tA0( 0, 0 ) = S::mk( 0.0, 0.0 ); tA0( 0, 1 ) = S::mk( -1.0, 0.0 ); tA0( 0, 2 ) = S::mk( 0.0, 0.0 );
    tA0( 1, 0 ) = S::mk( 1.0, 0.0 ); tA0( 1, 1 ) = S::mk(  0.0, 0.0 ); tA0( 1, 2 ) = S::mk( 0.0, 0.0 );
    tA0( 2, 0 ) = S::mk( 0.0, 0.0 ); tA0( 2, 1 ) = S::mk(  0.0, 0.0 ); tA0( 2, 2 ) = S::mk( 2.0, 0.0 );

    belfem::Matrix< TypeParam > tA( tA0 );
    belfem::Vector< Tcplx > tW;
    belfem::Matrix< TypeParam > tVL, tVR;
    belfem::Vector< Treal > tWork;

    belfem::int_t tInfo = belfem::geev( tA, tW, tVL, tVR, tWork );
    EXPECT_EQ( tInfo, 0 );
    EXPECT_EQ( tVR.n_rows(), 3u );
    EXPECT_EQ( tVR.n_cols(), 3u );

    // residual || A * v_j - w_j * v_j || for every eigenpair
    real tErr = 0.0;
    for ( uint j = 0; j < 3; ++j )
    {
        Tcplx tV[ 3 ];

        if constexpr ( std::is_same< TypeParam, Treal >::value )
        {
            // real flavor: unpack the LAPACK pair convention
            if ( tW( j ).imag() > Treal( 0 ) )
            {
                for ( uint i = 0; i < 3; ++i )
                    tV[ i ] = Tcplx( tVR( i, j ), tVR( i, j + 1 ) );
            }
            else if ( tW( j ).imag() < Treal( 0 ) )
            {
                for ( uint i = 0; i < 3; ++i )
                    tV[ i ] = Tcplx( tVR( i, j - 1 ), -tVR( i, j ) );
            }
            else
            {
                for ( uint i = 0; i < 3; ++i )
                    tV[ i ] = Tcplx( tVR( i, j ), Treal( 0 ) );
            }
        }
        else
        {
            for ( uint i = 0; i < 3; ++i )
            {
                tV[ i ] = tVR( i, j );
            }
        }

        for ( uint i = 0; i < 3; ++i )
        {
            Tcplx tSum( Treal( 0 ), Treal( 0 ) );
            for ( uint k = 0; k < 3; ++k )
            {
                tSum += Tcplx( tA0( i, k ) ) * tV[ k ];
            }
            tErr = std::max( tErr, ( real ) std::abs( tSum - tW( j ) * tV[ i ] ) );
        }
    }
    EXPECT_LT( tErr, S::tol() );

    // the spectrum itself: sorted by imaginary part it must be -i, 0+2, +i
    Treal tIm[ 3 ] = { tW( 0 ).imag(), tW( 1 ).imag(), tW( 2 ).imag() };
    std::sort( tIm, tIm + 3 );
    EXPECT_NEAR( ( real ) tIm[ 0 ], -1.0, S::tol() );
    EXPECT_NEAR( ( real ) tIm[ 1 ],  0.0, S::tol() );
    EXPECT_NEAR( ( real ) tIm[ 2 ],  1.0, S::tol() );
}

// =============================================================================
// gees ( Schur decomposition )
// =============================================================================

TYPED_TEST( LapackTest, GeesSchurReconstruct )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;
    typedef std::complex< Treal > Tcplx ;

    // rotation block + decoupled row: eigenvalues { +i, -i, 2 }
    belfem::Matrix< TypeParam > tA0( 3, 3 );
    tA0( 0, 0 ) = S::mk( 0.0, 0.0 ); tA0( 0, 1 ) = S::mk( -1.0, 0.0 ); tA0( 0, 2 ) = S::mk( 0.5, 0.0 );
    tA0( 1, 0 ) = S::mk( 1.0, 0.0 ); tA0( 1, 1 ) = S::mk(  0.0, 0.0 ); tA0( 1, 2 ) = S::mk( -0.25, 0.0 );
    tA0( 2, 0 ) = S::mk( 0.0, 0.0 ); tA0( 2, 1 ) = S::mk(  0.0, 0.0 ); tA0( 2, 2 ) = S::mk( 2.0, 0.0 );

    belfem::Matrix< TypeParam > tA( tA0 );   // overwritten with T
    belfem::Vector< Tcplx > tW;
    belfem::Matrix< TypeParam > tVS;
    belfem::Vector< Treal > tWork;
    belfem::Vector< belfem::int_t > tBWork;

    belfem::int_t tInfo = belfem::gees( tA, tW, tVS, tWork, tBWork );
    EXPECT_EQ( tInfo, 0 );
    EXPECT_EQ( tVS.n_rows(), 3u );
    EXPECT_EQ( tVS.n_cols(), 3u );

    // spectrum: sorted by imaginary part it must be -i, 2, +i
    Treal tIm[ 3 ] = { tW( 0 ).imag(), tW( 1 ).imag(), tW( 2 ).imag() };
    std::sort( tIm, tIm + 3 );
    EXPECT_NEAR( ( real ) tIm[ 0 ], -1.0, S::tol() );
    EXPECT_NEAR( ( real ) tIm[ 1 ],  0.0, S::tol() );
    EXPECT_NEAR( ( real ) tIm[ 2 ],  1.0, S::tol() );

    real tErr = 0.0;

    // Schur vectors are unitary: Z^H * Z = I
    for ( uint i = 0; i < 3; ++i )
        for ( uint j = 0; j < 3; ++j )
        {
            TypeParam tSum = S::mk( 0.0, 0.0 );
            for ( uint k = 0; k < 3; ++k )
            {
                tSum += S::conj( tVS( k, i ) ) * tVS( k, j );
            }
            tErr = std::max( tErr, ( real ) std::abs(
                tSum - S::mk( i == j ? 1.0 : 0.0, 0.0 ) ) );
        }

    // reconstruction: A0 * Z == Z * T ( T sits in tA after the call )
    for ( uint i = 0; i < 3; ++i )
        for ( uint j = 0; j < 3; ++j )
        {
            TypeParam tLhs = S::mk( 0.0, 0.0 );
            TypeParam tRhs = S::mk( 0.0, 0.0 );
            for ( uint k = 0; k < 3; ++k )
            {
                tLhs += tA0( i, k ) * tVS( k, j );
                tRhs += tVS( i, k ) * tA( k, j );
            }
            tErr = std::max( tErr, ( real ) std::abs( tLhs - tRhs ) );
        }

    EXPECT_LT( tErr, S::tol() );
}

TYPED_TEST( LapackTest, GeesSortedSelect )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;
    typedef std::complex< Treal > Tcplx ;

    belfem::Matrix< TypeParam > tA0( 3, 3 );
    tA0( 0, 0 ) = S::mk( 0.0, 0.0 ); tA0( 0, 1 ) = S::mk( -1.0, 0.0 ); tA0( 0, 2 ) = S::mk( 0.5, 0.0 );
    tA0( 1, 0 ) = S::mk( 1.0, 0.0 ); tA0( 1, 1 ) = S::mk(  0.0, 0.0 ); tA0( 1, 2 ) = S::mk( -0.25, 0.0 );
    tA0( 2, 0 ) = S::mk( 0.0, 0.0 ); tA0( 2, 1 ) = S::mk(  0.0, 0.0 ); tA0( 2, 2 ) = S::mk( 2.0, 0.0 );

    // select eigenvalues with real part > 1: only lambda = 2 qualifies,
    // so it must be moved to the top left of the Schur form
    belfem::lapack::gees_select_t< TypeParam > tSelect ;
    if constexpr ( std::is_same< TypeParam, Treal >::value )
    {
        tSelect = []( const TypeParam * aWr, const TypeParam * ) -> belfem::int_t
        {
            return *aWr > TypeParam( 1 ) ? 1 : 0 ;
        };
    }
    else
    {
        tSelect = []( const TypeParam * aW ) -> belfem::int_t
        {
            return aW->real() > Treal( 1 ) ? 1 : 0 ;
        };
    }

    belfem::Matrix< TypeParam > tA( tA0 );
    belfem::Vector< Tcplx > tW;
    belfem::Matrix< TypeParam > tVS;
    belfem::Vector< Treal > tWork;
    belfem::Vector< belfem::int_t > tBWork;
    belfem::int_t tSdim = 0;

    belfem::int_t tInfo = belfem::gees(
        tA, tW, tVS, tWork, tBWork, tSelect, 'V', 'S', &tSdim );

    EXPECT_EQ( tInfo, 0 );
    EXPECT_EQ( tSdim, 1 );

    // the selected eigenvalue leads
    EXPECT_NEAR( ( real ) tW( 0 ).real(), 2.0, S::tol() );
    EXPECT_NEAR( ( real ) tW( 0 ).imag(), 0.0, S::tol() );

    // reconstruction must still hold after reordering
    real tErr = 0.0;
    for ( uint i = 0; i < 3; ++i )
        for ( uint j = 0; j < 3; ++j )
        {
            TypeParam tLhs = S::mk( 0.0, 0.0 );
            TypeParam tRhs = S::mk( 0.0, 0.0 );
            for ( uint k = 0; k < 3; ++k )
            {
                tLhs += tA0( i, k ) * tVS( k, j );
                tRhs += tVS( i, k ) * tA( k, j );
            }
            tErr = std::max( tErr, ( real ) std::abs( tLhs - tRhs ) );
        }
    EXPECT_LT( tErr, S::tol() );
}

TYPED_TEST( LapackTest, GeesSortedConjugatePair )
{
    using S = Scalar< TypeParam >;
    typedef belfem::lapack::real_t< TypeParam > Treal ;
    typedef std::complex< Treal > Tcplx ;

    belfem::Matrix< TypeParam > tA0( 3, 3 );
    tA0( 0, 0 ) = S::mk( 0.0, 0.0 ); tA0( 0, 1 ) = S::mk( -1.0, 0.0 ); tA0( 0, 2 ) = S::mk( 0.5, 0.0 );
    tA0( 1, 0 ) = S::mk( 1.0, 0.0 ); tA0( 1, 1 ) = S::mk(  0.0, 0.0 ); tA0( 1, 2 ) = S::mk( -0.25, 0.0 );
    tA0( 2, 0 ) = S::mk( 0.0, 0.0 ); tA0( 2, 1 ) = S::mk(  0.0, 0.0 ); tA0( 2, 2 ) = S::mk( 2.0, 0.0 );

    // select imag > 0, i.e. only +i. The real flavors select the whole
    // conjugate pair if EITHER member matches, so sdim == 2 there; the
    // complex flavors count exactly the matching eigenvalue, sdim == 1
    belfem::lapack::gees_select_t< TypeParam > tSelect ;
    if constexpr ( std::is_same< TypeParam, Treal >::value )
    {
        tSelect = []( const TypeParam *, const TypeParam * aWi ) -> belfem::int_t
        {
            return *aWi > TypeParam( 0 ) ? 1 : 0 ;
        };
    }
    else
    {
        tSelect = []( const TypeParam * aW ) -> belfem::int_t
        {
            return aW->imag() > Treal( 0 ) ? 1 : 0 ;
        };
    }

    belfem::Matrix< TypeParam > tA( tA0 );
    belfem::Vector< Tcplx > tW;
    belfem::Matrix< TypeParam > tVS;
    belfem::Vector< Treal > tWork;
    belfem::Vector< belfem::int_t > tBWork;
    belfem::int_t tSdim = 0;

    belfem::int_t tInfo = belfem::gees(
        tA, tW, tVS, tWork, tBWork, tSelect, 'V', 'S', &tSdim );

    EXPECT_EQ( tInfo, 0 );
    if constexpr ( std::is_same< TypeParam, Treal >::value )
    {
        EXPECT_EQ( tSdim, 2 );
    }
    else
    {
        EXPECT_EQ( tSdim, 1 );
        EXPECT_NEAR( ( real ) tW( 0 ).imag(), 1.0, S::tol() );
    }

    // reconstruction must hold after reordering
    real tErr = 0.0;
    for ( uint i = 0; i < 3; ++i )
        for ( uint j = 0; j < 3; ++j )
        {
            TypeParam tLhs = S::mk( 0.0, 0.0 );
            TypeParam tRhs = S::mk( 0.0, 0.0 );
            for ( uint k = 0; k < 3; ++k )
            {
                tLhs += tA0( i, k ) * tVS( k, j );
                tRhs += tVS( i, k ) * tA( k, j );
            }
            tErr = std::max( tErr, ( real ) std::abs( tLhs - tRhs ) );
        }
    EXPECT_LT( tErr, S::tol() );
}
