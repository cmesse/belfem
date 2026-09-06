/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for tensor contraction, rotation, Kelvin-Christoffel, and
 * invert_symmetric kernels — validated against naive reference loops.
 * See: tests_04_tensor.md §3–§6
 *
 * The hand-unrolled pointer-arithmetic kernels are the primary risk area.
 * Each kernel is tested against an obviously-correct naive loop implementation.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstring>

#include "typedefs.hpp"
#include "cl_Tensor.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_ddot.hpp"
#include "fn_identity.hpp"
#include "fn_rotate.hpp"
#include "fn_kelvin_christoffel.hpp"
#include "fn_invert_symmetric.hpp"
#include "fn_norm.hpp"
#include "fn_inv.hpp"
#include "fn_compliance_matrix.hpp"

namespace
{
    const belfem::real tEps = 1e-12;  // exact-in-theory
    const belfem::real tTol = 1e-9;   // through inversion / many FP ops

    // =========================================================================
    // Reference loop helpers — slow but obviously correct
    // =========================================================================

    // A_ijkl * B_kl = C_ij
    void ref_contract42( const belfem::real * A,
                         const belfem::real * B,
                               belfem::real * C )
    {
        for( int i = 0; i < 3; ++i )
        for( int j = 0; j < 3; ++j )
        {
            belfem::real tSum = 0.0;
            for( int k = 0; k < 3; ++k )
            for( int l = 0; l < 3; ++l )
                tSum += A[ l*27 + k*9 + j*3 + i ] * B[ l*3 + k ];
            C[ j*3 + i ] = tSum;
        }
    }

    // A_ijmn * B_mnkl = C_ijkl
    void ref_contract44( const belfem::real * A,
                         const belfem::real * B,
                               belfem::real * C )
    {
        for( int i = 0; i < 3; ++i )
        for( int j = 0; j < 3; ++j )
        for( int k = 0; k < 3; ++k )
        for( int l = 0; l < 3; ++l )
        {
            belfem::real tSum = 0.0;
            for( int m = 0; m < 3; ++m )
            for( int n = 0; n < 3; ++n )
                tSum += A[ n*27 + m*9 + j*3 + i ] * B[ l*27 + k*9 + n*3 + m ];
            C[ l*27 + k*9 + j*3 + i ] = tSum;
        }
    }

    // A_mnop = B_ijkl * R_im * R_jn * R_ko * R_lp
    void ref_rotate42( const belfem::real * B,
                       const belfem::real * R,
                             belfem::real * A )
    {
        std::fill( A, A + 81, 0.0 );
        for( int m = 0; m < 3; ++m )
        for( int n = 0; n < 3; ++n )
        for( int o = 0; o < 3; ++o )
        for( int p = 0; p < 3; ++p )
        for( int i = 0; i < 3; ++i )
        for( int j = 0; j < 3; ++j )
        for( int k = 0; k < 3; ++k )
        for( int l = 0; l < 3; ++l )
            A[ p*27 + o*9 + n*3 + m ] += B[ l*27 + k*9 + j*3 + i ]
                * R[ m*3 + i ] * R[ n*3 + j ] * R[ o*3 + k ] * R[ p*3 + l ];
    }

    // C_ik = A_ijkl * n_j * n_l
    void ref_kelvin_christoffel( const belfem::real * A,
                                 const belfem::real * n,
                                       belfem::real * C )
    {
        for( int i = 0; i < 3; ++i )
        for( int k = 0; k < 3; ++k )
        {
            belfem::real tSum = 0.0;
            for( int j = 0; j < 3; ++j )
            for( int l = 0; l < 3; ++l )
                tSum += A[ l*27 + k*9 + j*3 + i ] * n[ j ] * n[ l ];
            C[ k*3 + i ] = tSum;
        }
    }

    // Dense column-major staging of a 3x3 matrix for the reference loops.
    // Matrix::data() must NOT be indexed linearly: under Blaze the columns
    // are padded, so the inter-column stride exceeds 3 (see CLAUDE.md,
    // matrix memory layout) — only the ( i, j ) accessor is backend-safe.
    void dense33( const belfem::Matrix< belfem::real > & aM, belfem::real * aOut )
    {
        for( int j = 0; j < 3; ++j )
        {
            for( int i = 0; i < 3; ++i )
            {
                aOut[ j * 3 + i ] = aM( i, j );
            }
        }
    }

    // Fill a deterministic non-symmetric tensor for testing
    void fill_test_tensor( belfem::Tensor< belfem::real > & aTen )
    {
        for( belfem::index_t n = 0; n < 81; ++n )
        {
            aTen.data()[ n ] = std::sin( 0.3 * n + 0.7 );
        }
    }

    // Create a rotation matrix about z-axis by angle θ
    void make_rotation_z( belfem::real aTheta,
                          belfem::Matrix< belfem::real > & aR )
    {
        aR.set_size( 3, 3, 0.0 );
        belfem::real tC = std::cos( aTheta );
        belfem::real tS = std::sin( aTheta );
        aR( 0, 0 ) =  tC;
        aR( 0, 1 ) = -tS;
        aR( 1, 0 ) =  tS;
        aR( 1, 1 ) =  tC;
        aR( 2, 2 ) =  1.0;
    }

    // Create a rotation matrix about x-axis by angle θ
    void make_rotation_x( belfem::real aTheta,
                          belfem::Matrix< belfem::real > & aR )
    {
        aR.set_size( 3, 3, 0.0 );
        belfem::real tC = std::cos( aTheta );
        belfem::real tS = std::sin( aTheta );
        aR( 0, 0 ) =  1.0;
        aR( 1, 1 ) =  tC;
        aR( 1, 2 ) = -tS;
        aR( 2, 1 ) =  tS;
        aR( 2, 2 ) =  tC;
    }
}

// =============================================================================
// §3.1  contract42 — A_ijkl B_kl = C_ij  [semantic]
// =============================================================================

TEST( TensorKernel, Contract42VsReferenceLoop )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    fill_test_tensor( tA );

    belfem::Matrix< belfem::real > tB( 3, 3 );
    for( int j = 0; j < 3; ++j )
        for( int i = 0; i < 3; ++i )
            tB( i, j ) = std::cos( 0.5 * ( j * 3 + i ) + 0.3 );

    belfem::Matrix< belfem::real > tC( 3, 3, 0.0 );
    tC = tA % tB;

    belfem::real tBd[9];
    dense33( tB, tBd );

    belfem::real tRef[9];
    ref_contract42( tA.data(), tBd, tRef );

    for( int j = 0; j < 3; ++j )
    {
        for( int i = 0; i < 3; ++i )
        {
            EXPECT_NEAR( tC( i, j ), tRef[ j * 3 + i ], tEps );
        }
    }
}

TEST( TensorKernel, Contract42IdentityTensor )
{
    belfem::Tensor< belfem::real > tI( 3, 3, 3, 3 );
    belfem::tensor::identity( tI );

    // symmetric matrix
    belfem::Matrix< belfem::real > tS = { { 1.0, 2.0, 3.0 },
                                            { 2.0, 5.0, 4.0 },
                                            { 3.0, 4.0, 6.0 } };

    belfem::Matrix< belfem::real > tC = tI % tS;

    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    {
        EXPECT_NEAR( tC( i, j ), tS( i, j ), tEps );
    }
}

TEST( TensorKernel, Contract42ZeroTensor )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 0.0 );
    belfem::Matrix< belfem::real > tB( 3, 3, 1.0 );
    belfem::Matrix< belfem::real > tC = tA % tB;

    for( int i = 0; i < 9; ++i )
    {
        EXPECT_NEAR( tC.data()[ i ], 0.0, tEps );
    }
}

TEST( TensorKernel, Contract42ZeroMatrix )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 1.0 );
    belfem::Matrix< belfem::real > tB( 3, 3, 0.0 );
    belfem::Matrix< belfem::real > tC = tA % tB;

    for( int i = 0; i < 9; ++i )
    {
        EXPECT_NEAR( tC.data()[ i ], 0.0, tEps );
    }
}

TEST( TensorKernel, Contract42MemberDdotMatchesFreeFunction )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    fill_test_tensor( tA );

    belfem::Matrix< belfem::real > tB( 3, 3, 2.0 );
    belfem::Matrix< belfem::real > tC1( 3, 3 );
    belfem::Matrix< belfem::real > tC2( 3, 3 );

    tA.ddot( tB, tC1 );
    belfem::ddot( tA, tB, tC2 );

    for( int i = 0; i < 9; ++i )
    {
        EXPECT_NEAR( tC1.data()[ i ], tC2.data()[ i ], tEps );
    }
}

TEST( TensorKernel, Contract42OperatorPercentMatchesDdot )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    fill_test_tensor( tA );

    belfem::Matrix< belfem::real > tB( 3, 3, 3.0 );
    belfem::Matrix< belfem::real > tC1( 3, 3 );
    tA.ddot( tB, tC1 );

    belfem::Matrix< belfem::real > tC2 = tA % tB;

    for( int i = 0; i < 9; ++i )
    {
        EXPECT_NEAR( tC1.data()[ i ], tC2.data()[ i ], tEps );
    }
}

// =============================================================================
// §3.2  contract44 — A_ijmn B_mnkl = C_ijkl  [semantic]
// =============================================================================

TEST( TensorKernel, Contract44VsReferenceLoop )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    fill_test_tensor( tA );

    // second tensor with different pattern
    for( belfem::index_t n = 0; n < 81; ++n )
        tB.data()[ n ] = std::cos( 0.2 * n + 1.3 );

    belfem::Tensor< belfem::real > tC = tA % tB;

    belfem::real tRef[81];
    ref_contract44( tA.data(), tB.data(), tRef );

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tC.data()[ n ], tRef[ n ], tTol );
    }
}

TEST( TensorKernel, Contract44ZeroLeft )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 0.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 1.0 );
    belfem::Tensor< belfem::real > tC = tA % tB;

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tC.data()[ n ], 0.0, tEps );
    }
}

TEST( TensorKernel, Contract44MemberDdotMatchesFreeFunction )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    fill_test_tensor( tA );
    for( belfem::index_t n = 0; n < 81; ++n )
        tB.data()[ n ] = std::cos( 0.4 * n );

    belfem::Tensor< belfem::real > tC1( 3, 3, 3, 3 );
    belfem::Tensor< belfem::real > tC2( 3, 3, 3, 3 );

    tA.ddot( tB, tC1 );
    belfem::ddot( tA, tB, tC2 );

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tC1.data()[ n ], tC2.data()[ n ], tEps );
    }
}

TEST( TensorKernel, Contract44OperatorPercentMatchesDdot )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 2.0 );
    fill_test_tensor( tA );

    belfem::Tensor< belfem::real > tC1( 3, 3, 3, 3 );
    tA.ddot( tB, tC1 );

    belfem::Tensor< belfem::real > tC2 = tA % tB;

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tC1.data()[ n ], tC2.data()[ n ], tEps );
    }
}

// =============================================================================
// §3.3  contract42 / contract44  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( TensorKernelDebug, Contract42Non3333Throws )
{
    belfem::Tensor< belfem::real > tA( 2, 2, 2, 2, 0.0 );
    belfem::Matrix< belfem::real > tB( 3, 3, 0.0 );
    belfem::Matrix< belfem::real > tC( 3, 3 );
    EXPECT_THROW( tA.ddot( tB, tC ), std::runtime_error );
}

TEST( TensorKernelDebug, Contract42WrongMatrixSizeThrows )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 0.0 );
    belfem::Matrix< belfem::real > tB( 4, 3, 0.0 );
    belfem::Matrix< belfem::real > tC( 3, 3 );
    EXPECT_THROW( tA.ddot( tB, tC ), std::runtime_error );
}

TEST( TensorKernelDebug, Contract44Non3333Throws )
{
    belfem::Tensor< belfem::real > tA( 2, 2, 2, 2, 0.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 0.0 );
    belfem::Tensor< belfem::real > tC( 3, 3, 3, 3 );
    EXPECT_THROW( belfem::ddot( tA, tB, tC ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §4.1  rotate42  [semantic]
// =============================================================================

TEST( TensorKernel, Rotate42VsReferenceLoop )
{
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    fill_test_tensor( tB );

    belfem::Matrix< belfem::real > tR( 3, 3 );
    make_rotation_z( 0.7, tR );

    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    belfem::rotate( tB, tR, tA );

    belfem::real tRd[9];
    dense33( tR, tRd );

    belfem::real tRef[81];
    ref_rotate42( tB.data(), tRd, tRef );

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tA.data()[ n ], tRef[ n ], tTol );
    }
}

TEST( TensorKernel, Rotate42IdentityRotation )
{
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    fill_test_tensor( tB );

    // identity rotation matrix
    belfem::Matrix< belfem::real > tR( 3, 3, 0.0 );
    tR( 0, 0 ) = 1.0;
    tR( 1, 1 ) = 1.0;
    tR( 2, 2 ) = 1.0;

    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    belfem::rotate( tB, tR, tA );

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tA.data()[ n ], tB.data()[ n ], tTol );
    }
}

TEST( TensorKernel, Rotate42IsotropicInvariant )
{
    // isotropic tensor should be unchanged by any rotation
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    tB.fill_isotropic_elasticity( 200.0e9, 0.3 );

    belfem::Matrix< belfem::real > tR( 3, 3 );
    make_rotation_z( 0.42 * M_PI, tR );

    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    belfem::rotate( tB, tR, tA );

    // tolerance relative to tensor magnitude — values are O(1e11),
    // so absolute noise in zero entries can be O(1e-5)
    belfem::real tMaxVal = 0.0;
    for( int n = 0; n < 81; ++n )
    {
        belfem::real tAbs = std::abs( tB.data()[ n ] );
        if( tAbs > tMaxVal ) tMaxVal = tAbs;
    }
    belfem::real tAbsTol = tMaxVal * 1e-9;

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tA.data()[ n ], tB.data()[ n ], tAbsTol );
    }
}

// =============================================================================
// §5.1  Kelvin-Christoffel  [semantic]
// =============================================================================

TEST( TensorKernel, KelvinChristoffelVsReferenceLoop )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    fill_test_tensor( tA );

    belfem::Vector< belfem::real > tN = { 0.6, 0.8, 0.0 };

    belfem::Matrix< belfem::real > tC( 3, 3 );
    belfem::kelvin_christoffel( tA, tN, tC );

    belfem::real tRef[9];
    ref_kelvin_christoffel( tA.data(), tN.data(), tRef );

    for( int j = 0; j < 3; ++j )
    {
        for( int i = 0; i < 3; ++i )
        {
            EXPECT_NEAR( tC( i, j ), tRef[ j * 3 + i ], tEps );
        }
    }
}

TEST( TensorKernel, KelvinChristoffelIsotropicXDirection )
{
    // isotropic C with n = {1,0,0}
    // Γ should be diagonal with (λ+2μ, μ, μ)
    belfem::real tE  = 200.0;
    belfem::real tNu = 0.25;
    belfem::real tLambda = tE * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) );
    belfem::real tMu     = tE / ( 2.0 * ( 1.0 + tNu ) );

    belfem::Tensor< belfem::real > tC( 3, 3, 3, 3 );
    tC.fill_isotropic_elasticity( tE, tNu );

    belfem::Vector< belfem::real > tN = { 1.0, 0.0, 0.0 };
    belfem::Matrix< belfem::real > tGamma( 3, 3 );
    belfem::kelvin_christoffel( tC, tN, tGamma );

    EXPECT_NEAR( tGamma( 0, 0 ), tLambda + 2.0 * tMu, tTol );
    EXPECT_NEAR( tGamma( 1, 1 ), tMu, tTol );
    EXPECT_NEAR( tGamma( 2, 2 ), tMu, tTol );
    EXPECT_NEAR( tGamma( 0, 1 ), 0.0, tTol );
    EXPECT_NEAR( tGamma( 0, 2 ), 0.0, tTol );
    EXPECT_NEAR( tGamma( 1, 2 ), 0.0, tTol );
}

TEST( TensorKernel, KelvinChristoffelScalingProperty )
{
    // scaling n by s → Γ scales by s²
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    fill_test_tensor( tA );

    belfem::Vector< belfem::real > tN = { 0.6, 0.8, 0.0 };
    belfem::Matrix< belfem::real > tG1( 3, 3 );
    belfem::kelvin_christoffel( tA, tN, tG1 );

    belfem::real tS = 2.5;
    belfem::Vector< belfem::real > tNs = { tN( 0 ) * tS, tN( 1 ) * tS, tN( 2 ) * tS };
    belfem::Matrix< belfem::real > tG2( 3, 3 );
    belfem::kelvin_christoffel( tA, tNs, tG2 );

    for( int i = 0; i < 9; ++i )
    {
        EXPECT_NEAR( tG2.data()[ i ], tG1.data()[ i ] * tS * tS, tTol );
    }
}

// =============================================================================
// §5.2  Kelvin-Christoffel  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( TensorKernelDebug, KelvinChristoffelNon3333Throws )
{
    belfem::Tensor< belfem::real > tA( 2, 2, 2, 2, 0.0 );
    belfem::Vector< belfem::real > tN = { 1.0, 0.0, 0.0 };
    belfem::Matrix< belfem::real > tC( 3, 3 );
    EXPECT_THROW( belfem::kelvin_christoffel( tA, tN, tC ), std::runtime_error );
}

TEST( TensorKernelDebug, KelvinChristoffelWrongVectorLengthThrows )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 0.0 );
    belfem::Vector< belfem::real > tN( 2, 1.0 );
    belfem::Matrix< belfem::real > tC( 3, 3 );
    EXPECT_THROW( belfem::kelvin_christoffel( tA, tN, tC ), std::runtime_error );
}

TEST( TensorKernelDebug, KelvinChristoffelWrongMatrixSizeThrows )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 0.0 );
    belfem::Vector< belfem::real > tN = { 1.0, 0.0, 0.0 };
    belfem::Matrix< belfem::real > tC( 4, 3 );
    EXPECT_THROW( belfem::kelvin_christoffel( tA, tN, tC ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §6.1  Invert Symmetric  [semantic]
// =============================================================================

TEST( TensorKernel, InvertSymmetricRoundTrip )
{
    // C/C-SiC-like orthotropic material
    belfem::real tE1 = 60.0, tE2 = 65.0, tE3 = 20.0;
    belfem::real tG23 = 8.0, tG13 = 8.5, tG12 = 9.0;
    belfem::real tNu23 = 0.2, tNu13 = 0.25, tNu12 = 0.1;

    belfem::Matrix< belfem::real > tS( 6, 6, 0.0 );
    belfem::compliance_matrix( tE1, tE2, tE3,
                                tNu23, tNu13, tNu12,
                                tG23, tG13, tG12, tS );

    belfem::Matrix< belfem::real > tC = belfem::inv( tS );
    belfem::Tensor< belfem::real > tT( tC );

    // make a copy and invert
    belfem::Tensor< belfem::real > tM( 3, 3, 3, 3 );
    tM = tT;

    belfem::Vector< belfem::real > tWork( 72 );
    belfem::Vector< belfem::int_t >          tPivot( 36 );
    belfem::tensor::invert_symmetric( tM, tWork, tPivot );

    // C_inv % C should ≈ I_sym
    belfem::Tensor< belfem::real > tProduct = tT % tM;

    belfem::Tensor< belfem::real > tI( 3, 3, 3, 3 );
    belfem::tensor::identity( tI );

    belfem::Vector< belfem::real > tError( 81 );
    for( belfem::index_t n = 0; n < 81; ++n )
    {
        tError( n ) = tI.data()[ n ] - tProduct.data()[ n ];
    }
    EXPECT_NEAR( belfem::norm( tError ), 0.0, tTol );
}

TEST( TensorKernel, InvertSymmetricPreservesSymmetry )
{
    // isotropic material — result should have minor symmetry
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    tTen.fill_isotropic_elasticity( 200.0, 0.3 );

    belfem::Vector< belfem::real > tWork( 72 );
    belfem::Vector< belfem::int_t >          tPivot( 36 );
    belfem::tensor::invert_symmetric( tTen, tWork, tPivot );

    // minor symmetry: A(i,j,k,l) = A(j,i,k,l) = A(i,j,l,k)
    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    for( int k = 0; k < 3; ++k )
    for( int l = 0; l < 3; ++l )
    {
        EXPECT_NEAR( tTen( i, j, k, l ), tTen( j, i, k, l ), tTol );
        EXPECT_NEAR( tTen( i, j, k, l ), tTen( i, j, l, k ), tTol );
    }
}

TEST( TensorKernel, InvertSymmetricPreservesMajorSymmetry )
{
    // Codex finding #4: also check major symmetry A(i,j,k,l) = A(k,l,i,j)
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    tTen.fill_isotropic_elasticity( 200.0, 0.3 );

    belfem::Vector< belfem::real > tWork( 72 );
    belfem::Vector< belfem::int_t >          tPivot( 36 );
    belfem::tensor::invert_symmetric( tTen, tWork, tPivot );

    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    for( int k = 0; k < 3; ++k )
    for( int l = 0; l < 3; ++l )
    {
        EXPECT_NEAR( tTen( i, j, k, l ), tTen( k, l, i, j ), tTol );
    }
}

TEST( TensorKernel, Contract44IdentityRight )
{
    // A : I_sym = A for symmetric A
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    tA.fill_isotropic_elasticity( 200.0, 0.3 );

    belfem::Tensor< belfem::real > tI( 3, 3, 3, 3 );
    belfem::tensor::identity( tI );

    belfem::Tensor< belfem::real > tResult( 3, 3, 3, 3 );
    belfem::ddot( tA, tI, tResult );

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tResult.data()[ n ], tA.data()[ n ], tEps );
    }
}

TEST( TensorKernel, Contract44IdentityLeft )
{
    // I_sym : A = A for symmetric A
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    tA.fill_isotropic_elasticity( 200.0, 0.3 );

    belfem::Tensor< belfem::real > tI( 3, 3, 3, 3 );
    belfem::tensor::identity( tI );

    belfem::Tensor< belfem::real > tResult( 3, 3, 3, 3 );
    belfem::ddot( tI, tA, tResult );

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tResult.data()[ n ], tA.data()[ n ], tEps );
    }
}

TEST( TensorKernel, Rotate42TwoRotationsMatchComposition )
{
    // R2 * (R1 * A * R1^T) * R2^T == (R2*R1) * A * (R2*R1)^T
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    fill_test_tensor( tA );

    belfem::Matrix< belfem::real > tR1( 3, 3 );
    belfem::Matrix< belfem::real > tR2( 3, 3 );
    // Use DIFFERENT axes so rotations don't commute (Codex finding)
    make_rotation_z( 0.3, tR1 );
    make_rotation_x( 0.7, tR2 );

    // sequential: first R1, then R2
    belfem::Tensor< belfem::real > tTemp( 3, 3, 3, 3 );
    belfem::Tensor< belfem::real > tSeq( 3, 3, 3, 3 );
    belfem::rotate( tA, tR1, tTemp );
    belfem::rotate( tTemp, tR2, tSeq );

    // composed: for tensor rotation, sequential rotate(A,R1) then rotate(T,R2)
    // is equivalent to rotate(A, R1*R2) because the contraction convention is
    // A_mnop = R_im R_jn R_ko R_lp B_ijkl (sum over inner indices)
    belfem::Matrix< belfem::real > tRc( 3, 3, 0.0 );
    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    for( int k = 0; k < 3; ++k )
        tRc( i, j ) += tR1( i, k ) * tR2( k, j );

    belfem::Tensor< belfem::real > tComp( 3, 3, 3, 3 );
    belfem::rotate( tA, tRc, tComp );

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tSeq.data()[ n ], tComp.data()[ n ], tEps );
    }
}

TEST( TensorKernel, Contract42IsotropicHydrostatic )
{
    // C:ε = σ. For hydrostatic strain ε = δ, isotropic C gives σ = (3λ+2μ)δ
    belfem::real tE = 200.0, tNu = 0.3;
    belfem::real tLambda = tE * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) );
    belfem::real tMu     = tE / ( 2.0 * ( 1.0 + tNu ) );

    belfem::Tensor< belfem::real > tC( 3, 3, 3, 3 );
    tC.fill_isotropic_elasticity( tE, tNu );

    // identity strain ε_kl = δ_kl
    belfem::Matrix< belfem::real > tEps( 3, 3, 0.0 );
    tEps( 0, 0 ) = 1.0;
    tEps( 1, 1 ) = 1.0;
    tEps( 2, 2 ) = 1.0;

    belfem::Matrix< belfem::real > tSigma( 3, 3, 0.0 );
    belfem::ddot( tC, tEps, tSigma );

    belfem::real tExpected = 3.0 * tLambda + 2.0 * tMu;
    EXPECT_NEAR( tSigma( 0, 0 ), tExpected, tTol );
    EXPECT_NEAR( tSigma( 1, 1 ), tExpected, tTol );
    EXPECT_NEAR( tSigma( 2, 2 ), tExpected, tTol );
    EXPECT_NEAR( tSigma( 0, 1 ), 0.0, tTol );
}

TEST( TensorKernel, Contract42IsotropicPureShear )
{
    // pure shear ε_12 = ε_21 = 0.5 → σ_12 = 2μ * 0.5 = μ
    belfem::real tE = 200.0, tNu = 0.3;
    belfem::real tMu = tE / ( 2.0 * ( 1.0 + tNu ) );

    belfem::Tensor< belfem::real > tC( 3, 3, 3, 3 );
    tC.fill_isotropic_elasticity( tE, tNu );

    belfem::Matrix< belfem::real > tEps( 3, 3, 0.0 );
    tEps( 0, 1 ) = 0.5;
    tEps( 1, 0 ) = 0.5;

    belfem::Matrix< belfem::real > tSigma( 3, 3, 0.0 );
    belfem::ddot( tC, tEps, tSigma );

    EXPECT_NEAR( tSigma( 0, 1 ), tMu, tTol );
    EXPECT_NEAR( tSigma( 1, 0 ), tMu, tTol );
    EXPECT_NEAR( tSigma( 0, 0 ), 0.0, tTol );
}

TEST( TensorKernel, KelvinChristoffelSymmetricForSymmetricTensor )
{
    // Γ_ik = C_ijkl n_j n_l should be symmetric for symmetric C
    belfem::Tensor< belfem::real > tC( 3, 3, 3, 3 );
    tC.fill_isotropic_elasticity( 200.0, 0.3 );

    belfem::Vector< belfem::real > tN = { 0.5, 0.6, std::sqrt( 1.0 - 0.25 - 0.36 ) };
    belfem::Matrix< belfem::real > tGamma( 3, 3, 0.0 );
    belfem::kelvin_christoffel( tC, tN, tGamma );

    EXPECT_NEAR( tGamma( 0, 1 ), tGamma( 1, 0 ), tTol );
    EXPECT_NEAR( tGamma( 0, 2 ), tGamma( 2, 0 ), tTol );
    EXPECT_NEAR( tGamma( 1, 2 ), tGamma( 2, 1 ), tTol );
}

TEST( TensorKernel, Rotate42HighLevelMatchesLowLevel )
{
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    fill_test_tensor( tB );

    belfem::Matrix< belfem::real > tR( 3, 3 );
    make_rotation_z( 0.42 * M_PI, tR );

    // high-level API
    belfem::Tensor< belfem::real > tA1( 3, 3, 3, 3 );
    belfem::rotate( tB, tR, tA1 );

    // low-level API
    belfem::Tensor< belfem::real > tA2( 3, 3, 3, 3 );
    belfem::tensor::rotate42( tB.data(), tR, tA2.data() );

    for( int n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tA1.data()[ n ], tA2.data()[ n ], tEps );
    }
}

// =============================================================================
// §6.2  Invert Symmetric  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( TensorKernelDebug, InvertSymmetricNon3333Throws )
{
    belfem::Tensor< belfem::real > tTen( 2, 2, 2, 2, 0.0 );
    belfem::Vector< belfem::real > tWork( 72 );
    belfem::Vector< belfem::int_t >          tPivot( 36 );
    EXPECT_THROW( belfem::tensor::invert_symmetric( tTen, tWork, tPivot ),
                  std::runtime_error );
}

TEST( TensorKernelDebug, InvertSymmetricWorkVectorTooShortThrows )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 1.0 );
    belfem::Vector< belfem::real > tWork( 50 );   // needs 72
    belfem::Vector< belfem::int_t >          tPivot( 36 );
    EXPECT_THROW( belfem::tensor::invert_symmetric( tTen, tWork, tPivot ),
                  std::runtime_error );
}

TEST( TensorKernelDebug, InvertSymmetricPivotVectorTooShortThrows )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 1.0 );
    belfem::Vector< belfem::real > tWork( 72 );
    belfem::Vector< belfem::int_t >          tPivot( 20 );   // needs 36
    EXPECT_THROW( belfem::tensor::invert_symmetric( tTen, tWork, tPivot ),
                  std::runtime_error );
}

#endif // NDEBUG
