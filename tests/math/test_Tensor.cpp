/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Tensor<T> class semantics, constitutive fill, and
 * mat↔ten conversion.
 * See: tests_04_tensor.md §1–§2
 *
 * Tensor uses manual malloc/free — all tests are Valgrind candidates.
 * Flat storage: data[L*27 + K*9 + J*3 + I] == operator()(I,J,K,L) for 3×3×3×3.
 * operator== is EXACT elementwise. Use EXPECT_NEAR for anything with arithmetic.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Tensor.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_identity.hpp"
#include "fn_inv.hpp"
#include "fn_compliance_matrix.hpp"

namespace
{
    const belfem::real tEps = 1e-12;  // exact-in-theory
    const belfem::real tTol = 1e-9;   // through inversion / many FP ops
}

// =============================================================================
// §1.1  Construction & Destruction  [semantic]
// =============================================================================

TEST( Tensor, SizedConstructor )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    EXPECT_EQ( tTen.capacity(), 81u );
    EXPECT_TRUE( tTen.is_3333() );
}

TEST( Tensor, SizedConstructorNon3333 )
{
    belfem::Tensor< belfem::real > tTen( 2, 2, 2, 2 );
    EXPECT_EQ( tTen.capacity(), 16u );
    EXPECT_FALSE( tTen.is_3333() );
}

TEST( Tensor, FillValueConstructor )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 7.0 );
    for( belfem::index_t n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tTen.data()[ n ], 7.0, tEps );
    }
}

TEST( Tensor, FromElasticityMatrix )
{
    // create a known 6×6 matrix
    belfem::Matrix< belfem::real > tC( 6, 6, 0.0 );
    tC( 0, 0 ) = 100.0;
    tC( 1, 1 ) = 100.0;
    tC( 2, 2 ) = 100.0;
    tC( 3, 3 ) = 50.0;
    tC( 4, 4 ) = 50.0;
    tC( 5, 5 ) = 50.0;

    belfem::Tensor< belfem::real > tTen( tC );
    EXPECT_TRUE( tTen.is_3333() );

    // verify a diagonal entry: C(0,0) maps to A(0,0,0,0)
    EXPECT_NEAR( tTen( 0, 0, 0, 0 ), 100.0, tEps );
}

TEST( Tensor, CopyConstructorDeepCopies )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 1.0 );
    belfem::Tensor< belfem::real > tB( tA );
    tB( 0, 0, 0, 0 ) = 99.0;
    EXPECT_NEAR( tA( 0, 0, 0, 0 ), 1.0, tEps );
}

TEST( Tensor, MoveConstructorTransfers )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 2.0 );
    belfem::Tensor< belfem::real > tB( std::move( tA ) );
    EXPECT_NEAR( tB( 0, 0, 0, 0 ), 2.0, tEps );
    EXPECT_EQ( tA.data(), nullptr );
}

TEST( Tensor, CopyAssignment )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 3.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    tB = tA;
    EXPECT_TRUE( tA == tB );
    tB( 0, 0, 0, 0 ) = 99.0;
    EXPECT_NEAR( tA( 0, 0, 0, 0 ), 3.0, tEps );
}

TEST( Tensor, MoveAssignment )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 4.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    tB = std::move( tA );
    EXPECT_NEAR( tB( 0, 0, 0, 0 ), 4.0, tEps );
    EXPECT_EQ( tA.data(), nullptr );
}

TEST( Tensor, SelfCopyAssignment )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 7.0 );
    tA = tA;
    EXPECT_NEAR( tA( 0, 0, 0, 0 ), 7.0, tEps );
}

TEST( Tensor, SelfMoveAssignment )
{
    // BUG-T1 regression: self-move used to free(mData) before reading it
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 7.0 );
#pragma GCC diagnostic push
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif
    tA = std::move( tA );
#pragma GCC diagnostic pop
    EXPECT_NE( tA.data(), nullptr );
    EXPECT_NEAR( tA( 0, 0, 0, 0 ), 7.0, tEps );
}

TEST( Tensor, ScalarAssignment )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    tTen = 5.0;
    for( belfem::index_t n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tTen.data()[ n ], 5.0, tEps );
    }
}

TEST( Tensor, MatrixAssignment )
{
    belfem::Matrix< belfem::real > tC( 6, 6, 0.0 );
    tC( 0, 0 ) = 200.0;
    tC( 3, 3 ) = 75.0;

    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    tTen = tC;
    EXPECT_NEAR( tTen( 0, 0, 0, 0 ), 200.0, tEps );
}

// =============================================================================
// §1.2  Construction  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( TensorDebug, FromMatrixWrongSizeThrows )
{
    belfem::Matrix< belfem::real > tC( 5, 6, 0.0 );
    EXPECT_THROW( belfem::Tensor< belfem::real > tTen( tC ), std::runtime_error );
}

TEST( TensorDebug, CopyAssignmentSizeMismatchThrows )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 1.0 );
    belfem::Tensor< belfem::real > tB( 2, 2, 2, 2, 1.0 );
    EXPECT_THROW( tB = tA, std::runtime_error );
}

TEST( TensorDebug, MoveAssignmentSizeMismatchThrows )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 1.0 );
    belfem::Tensor< belfem::real > tB( 2, 2, 2, 2, 1.0 );
    EXPECT_THROW( tB = std::move( tA ), std::runtime_error );
}

TEST( TensorDebug, MatrixAssignmentNon3333Throws )
{
    belfem::Matrix< belfem::real > tC( 6, 6, 0.0 );
    belfem::Tensor< belfem::real > tTen( 2, 2, 2, 2 );
    EXPECT_THROW( tTen = tC, std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §1.3  Access  [semantic]
// =============================================================================

TEST( Tensor, ParenthesisReadWrite )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 0.0 );
    tTen( 1, 2, 0, 1 ) = 42.0;
    EXPECT_NEAR( tTen( 1, 2, 0, 1 ), 42.0, tEps );
}

TEST( Tensor, FlatLayoutMatchesOperator )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    belfem::real tVal = 1.0;
    for( int l = 0; l < 3; ++l )
    for( int k = 0; k < 3; ++k )
    for( int j = 0; j < 3; ++j )
    for( int i = 0; i < 3; ++i )
    {
        tTen( i, j, k, l ) = tVal;
        tVal += 1.0;
    }

    for( int l = 0; l < 3; ++l )
    for( int k = 0; k < 3; ++k )
    for( int j = 0; j < 3; ++j )
    for( int i = 0; i < 3; ++i )
    {
        EXPECT_NEAR( tTen.data()[ l * 27 + k * 9 + j * 3 + i ],
                     tTen( i, j, k, l ), tEps );
    }
}

TEST( Tensor, DataPointerNonNull )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    EXPECT_NE( tTen.data(), nullptr );
}

TEST( Tensor, CapacityMatchesProduct )
{
    belfem::Tensor< belfem::real > tTen( 2, 3, 4, 5 );
    EXPECT_EQ( tTen.capacity(), 120u );
}

TEST( Tensor, Is3333 )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    EXPECT_TRUE( tA.is_3333() );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 2 );
    EXPECT_FALSE( tB.is_3333() );
}

// =============================================================================
// §1.4  Access  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( TensorDebug, IndexIOutOfBoundsThrows )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    EXPECT_THROW( tTen( 3, 0, 0, 0 ), std::runtime_error );
}

TEST( TensorDebug, IndexJOutOfBoundsThrows )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    EXPECT_THROW( tTen( 0, 3, 0, 0 ), std::runtime_error );
}

TEST( TensorDebug, IndexKOutOfBoundsThrows )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    EXPECT_THROW( tTen( 0, 0, 3, 0 ), std::runtime_error );
}

TEST( TensorDebug, IndexLOutOfBoundsThrows )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    EXPECT_THROW( tTen( 0, 0, 0, 3 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §1.5  Operators  [semantic]
// =============================================================================

TEST( Tensor, PlusEqualsScalar )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 1.0 );
    tTen += 2.0;
    EXPECT_NEAR( tTen( 0, 0, 0, 0 ), 3.0, tEps );
    EXPECT_NEAR( tTen( 2, 2, 2, 2 ), 3.0, tEps );
}

TEST( Tensor, MinusEqualsScalar )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 5.0 );
    tTen -= 2.0;
    EXPECT_NEAR( tTen( 1, 1, 1, 1 ), 3.0, tEps );
}

TEST( Tensor, TimesEqualsScalar )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 3.0 );
    tTen *= 2.0;
    EXPECT_NEAR( tTen( 0, 0, 0, 0 ), 6.0, tEps );
}

TEST( Tensor, DivideEqualsScalar )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 6.0 );
    tTen /= 2.0;
    EXPECT_NEAR( tTen( 0, 0, 0, 0 ), 3.0, tEps );
}

TEST( Tensor, PlusEqualsTensor )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 1.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 2.0 );
    tA += tB;
    EXPECT_NEAR( tA( 0, 0, 0, 0 ), 3.0, tEps );
}

TEST( Tensor, MinusEqualsTensor )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 5.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 2.0 );
    tA -= tB;
    EXPECT_NEAR( tA( 0, 0, 0, 0 ), 3.0, tEps );
}

TEST( Tensor, BinaryPlus )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 1.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 2.0 );
    belfem::Tensor< belfem::real > tC = tA + tB;
    EXPECT_NEAR( tC( 0, 0, 0, 0 ), 3.0, tEps );
}

TEST( Tensor, BinaryMinus )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 5.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 2.0 );
    belfem::Tensor< belfem::real > tC = tA - tB;
    EXPECT_NEAR( tC( 0, 0, 0, 0 ), 3.0, tEps );
}

TEST( Tensor, EqualityExact )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 7.0 );
    belfem::Tensor< belfem::real > tB( tA );
    EXPECT_TRUE( tA == tB );
}

TEST( Tensor, EqualityDifferent )
{
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 7.0 );
    belfem::Tensor< belfem::real > tB( tA );
    tB( 1, 1, 1, 1 ) = 8.0;
    EXPECT_FALSE( tA == tB );
}

// =============================================================================
// §1.6  Operators  [debug]
// =============================================================================

#ifndef NDEBUG

// The addition operators were generalized from the old 3x3x3x3-only rule to
// any MATCHING shapes ( cl_Tensor.hpp: per-dimension asserts ); equal-shape
// addition is legal, mismatched shapes must throw.

TEST( TensorDebug, PlusEqualsShapeMismatchThrows )
{
    belfem::Tensor< belfem::real > tA( 2, 2, 2, 2, 1.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 1.0 );
    EXPECT_THROW( tA += tB, std::runtime_error );
}

TEST( TensorDebug, PlusEqualsMatchingShapeIsLegal )
{
    belfem::Tensor< belfem::real > tA( 2, 2, 2, 2, 1.0 );
    belfem::Tensor< belfem::real > tB( 2, 2, 2, 2, 1.0 );
    EXPECT_NO_THROW( tA += tB );
    EXPECT_NEAR( tA( 0, 0, 0, 0 ), 2.0, 1e-12 );
}

TEST( TensorDebug, BinaryPlusShapeMismatchThrows )
{
    belfem::Tensor< belfem::real > tA( 2, 2, 2, 2, 1.0 );
    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3, 1.0 );
    EXPECT_THROW( auto tC = tA + tB, std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §2.1  Isotropic Fill  [semantic]
// =============================================================================

TEST( TensorFill, FillAbProducesCorrectStructure )
{
    belfem::real tA = 10.0;
    belfem::real tB = 6.0;
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    tTen.fill( tA, tB );

    belfem::real tDiag     = tA + tB * 4.0 / 3.0;   // a + 4b/3
    belfem::real tOffDiag  = tA - tB * 2.0 / 3.0;   // a - 2b/3

    // diagonal: (0,0,0,0), (1,1,1,1), (2,2,2,2)
    EXPECT_NEAR( tTen( 0, 0, 0, 0 ), tDiag, tEps );
    EXPECT_NEAR( tTen( 1, 1, 1, 1 ), tDiag, tEps );
    EXPECT_NEAR( tTen( 2, 2, 2, 2 ), tDiag, tEps );

    // off-diagonal coupling: (0,0,1,1), (0,0,2,2)
    EXPECT_NEAR( tTen( 0, 0, 1, 1 ), tOffDiag, tEps );
    EXPECT_NEAR( tTen( 0, 0, 2, 2 ), tOffDiag, tEps );

    // shear: (0,1,0,1), (0,1,1,0)
    EXPECT_NEAR( tTen( 0, 1, 0, 1 ), tB, tEps );
    EXPECT_NEAR( tTen( 0, 1, 1, 0 ), tB, tEps );

    // zero entry
    EXPECT_NEAR( tTen( 0, 0, 0, 1 ), 0.0, tEps );
}

TEST( TensorFill, FillIsotropicElasticity )
{
    belfem::real tE  = 200.0e9;
    belfem::real tNu = 0.3;

    belfem::real tK = tE / ( 3.0 * ( 1.0 - 2.0 * tNu ) );
    belfem::real tG = tE / ( 2.0 * ( 1.0 + tNu ) );

    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3 );
    tA.fill( tK, tG );

    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    tB.fill_isotropic_elasticity( tE, tNu );

    // both should be identical
    for( belfem::index_t n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tA.data()[ n ], tB.data()[ n ], tEps );
    }
}

TEST( TensorFill, FillIsotropicKnownSteel )
{
    // E = 200 GPa, ν = 0.3
    belfem::real tE  = 200.0e9;
    belfem::real tNu = 0.3;

    // Lamé parameters
    belfem::real tLambda = tE * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) );
    belfem::real tMu     = tE / ( 2.0 * ( 1.0 + tNu ) );

    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    tTen.fill_isotropic_elasticity( tE, tNu );

    // C11 = λ + 2μ
    EXPECT_NEAR( tTen( 0, 0, 0, 0 ), tLambda + 2.0 * tMu, tTol );
    // C12 = λ
    EXPECT_NEAR( tTen( 0, 0, 1, 1 ), tLambda, tTol );
    // C44 = μ  (shear)
    EXPECT_NEAR( tTen( 0, 1, 0, 1 ), tMu, tTol );
}

// =============================================================================
// §2.2  Isotropic Fill  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( TensorDebug, FillAbNon3333Throws )
{
    belfem::Tensor< belfem::real > tTen( 2, 2, 2, 2 );
    EXPECT_THROW( tTen.fill( 1.0, 1.0 ), std::runtime_error );
}

TEST( TensorDebug, FillIsotropicNon3333Throws )
{
    belfem::Tensor< belfem::real > tTen( 2, 2, 2, 2 );
    EXPECT_THROW( tTen.fill_isotropic_elasticity( 200e9, 0.3 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §2.3  Orthotropic Fill  [semantic]
// =============================================================================

TEST( TensorFill, FillOrthotropicProducesSymmetricVoigtMatrix )
{
    // (idea from ChatGPT)
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3 );
    tTen.fill_orthotropic_elasticity( 60.0, 65.0, 20.0,
                                       0.2, 0.25, 0.1,
                                       8.0, 8.5, 9.0 );

    belfem::Matrix< belfem::real > tC( 6, 6 );
    tTen.to_matrix( tC );

    for( size_t i = 0; i < 6; ++i )
    for( size_t j = 0; j < 6; ++j )
    {
        EXPECT_NEAR( tC( i, j ), tC( j, i ), tTol );
    }
}

TEST( TensorFill, FillOrthotropicIsotropicLimit )
{
    // use equal properties → should match isotropic
    belfem::real tE  = 100.0;
    belfem::real tNu = 0.25;
    belfem::real tG  = tE / ( 2.0 * ( 1.0 + tNu ) );

    belfem::Tensor< belfem::real > tIso( 3, 3, 3, 3 );
    tIso.fill_isotropic_elasticity( tE, tNu );

    belfem::Tensor< belfem::real > tOrtho( 3, 3, 3, 3 );
    tOrtho.fill_orthotropic_elasticity( tE, tE, tE, tNu, tNu, tNu, tG, tG, tG );

    for( belfem::index_t n = 0; n < 81; ++n )
    {
        EXPECT_NEAR( tIso.data()[ n ], tOrtho.data()[ n ], tTol );
    }
}

// =============================================================================
// §2.4  Identity Tensor  [semantic]
// =============================================================================

TEST( TensorFill, IdentityTensorValues )
{
    belfem::Tensor< belfem::real > tI( 3, 3, 3, 3 );
    belfem::tensor::identity( tI );

    // I_ijkl = 0.5 * (δ_ik δ_jl + δ_il δ_jk)
    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    for( int k = 0; k < 3; ++k )
    for( int l = 0; l < 3; ++l )
    {
        belfem::real tExpected = 0.5 * (
            ( i == k ? 1.0 : 0.0 ) * ( j == l ? 1.0 : 0.0 ) +
            ( i == l ? 1.0 : 0.0 ) * ( j == k ? 1.0 : 0.0 ) );
        EXPECT_NEAR( tI( i, j, k, l ), tExpected, tEps );
    }
}

TEST( TensorFill, IdentityContractionWithMatrix )
{
    // I % S ≈ S for a symmetric matrix
    belfem::Tensor< belfem::real > tI( 3, 3, 3, 3 );
    belfem::tensor::identity( tI );

    belfem::Matrix< belfem::real > tS = { { 1.0, 2.0, 3.0 },
                                            { 2.0, 5.0, 4.0 },
                                            { 3.0, 4.0, 6.0 } };

    belfem::Matrix< belfem::real > tResult = tI % tS;

    for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
    {
        EXPECT_NEAR( tResult( i, j ), tS( i, j ), tEps );
    }
}

// =============================================================================
// §2.5  Compliance Matrix  [semantic]
// =============================================================================

TEST( TensorConversion, ComplianceMatrixIsotropic )
{
    // (idea from ChatGPT) — verify specific known entries
    belfem::real tE  = 200.0e9;
    belfem::real tNu = 0.3;
    belfem::real tG  = tE / ( 2.0 * ( 1.0 + tNu ) );

    belfem::Matrix< belfem::real > tS( 6, 6, 0.0 );
    belfem::compliance_matrix( tE, tE, tE, tNu, tNu, tNu, tG, tG, tG, tS );

    EXPECT_NEAR( tS( 0, 0 ),  1.0 / tE,   tTol );
    EXPECT_NEAR( tS( 1, 1 ),  1.0 / tE,   tTol );
    EXPECT_NEAR( tS( 0, 1 ), -tNu / tE,   tTol );
    EXPECT_NEAR( tS( 3, 3 ),  1.0 / tG,   tTol );
}

TEST( TensorConversion, ComplianceMatrixSymmetric )
{
    // (idea from ChatGPT)
    belfem::Matrix< belfem::real > tS( 6, 6, 0.0 );
    belfem::compliance_matrix( 60.0, 65.0, 20.0,
                                0.2, 0.25, 0.1,
                                8.0, 8.5, 9.0, tS );

    for( size_t i = 0; i < 6; ++i )
    for( size_t j = 0; j < 6; ++j )
    {
        EXPECT_NEAR( tS( i, j ), tS( j, i ), tTol );
    }
}

TEST( TensorConversion, StiffnessTimesComplianceIsIdentity )
{
    // (idea from ChatGPT) — inv(S) * S ≈ I(6×6)
    belfem::Matrix< belfem::real > tS( 6, 6, 0.0 );
    belfem::compliance_matrix( 60.0, 65.0, 20.0,
                                0.2, 0.25, 0.1,
                                8.0, 8.5, 9.0, tS );

    belfem::Matrix< belfem::real > tC = belfem::inv( tS );
    belfem::Matrix< belfem::real > tI( tC * tS );

    for( size_t i = 0; i < 6; ++i )
    for( size_t j = 0; j < 6; ++j )
    {
        belfem::real tExpected = ( i == j ) ? 1.0 : 0.0;
        EXPECT_NEAR( tI( i, j ), tExpected, tTol );
    }
}

// =============================================================================
// §2.6  Mat ↔ Ten Conversion  [semantic]
// =============================================================================

TEST( TensorConversion, MatToTenToMatRoundTrip )
{
    // non-trivial asymmetric-looking 6×6
    belfem::Matrix< belfem::real > tC( 6, 6, 0.0 );
    belfem::real tVal = 1.0;
    for( size_t i = 0; i < 6; ++i )
    for( size_t j = 0; j < 6; ++j )
    {
        tC( i, j ) = tVal;
        tVal += 1.0;
    }

    belfem::Tensor< belfem::real > tTen( tC );
    belfem::Matrix< belfem::real > tC2( 6, 6 );
    tTen.to_matrix( tC2 );

    for( size_t i = 0; i < 6; ++i )
    for( size_t j = 0; j < 6; ++j )
    {
        EXPECT_NEAR( tC2( i, j ), tC( i, j ), tEps );
    }
}

TEST( TensorConversion, ConstructorFromMatrixMatchesMatToTen )
{
    belfem::Matrix< belfem::real > tC( 6, 6, 0.0 );
    tC( 0, 0 ) = 10.0;
    tC( 3, 3 ) = 5.0;

    belfem::Tensor< belfem::real > tA( tC );

    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    belfem::tensor::mat_to_ten( tC, tB.data() );

    EXPECT_TRUE( tA == tB );
}

TEST( TensorConversion, ConstructorFromMatrixMatchesMatrixAssignment )
{
    // Tensor(C) should produce same result as tTen = C (idea from ChatGPT)
    belfem::Matrix< belfem::real > tC( 6, 6, 0.0 );
    belfem::real tVal = 1.0;
    for( size_t i = 0; i < 6; ++i )
    for( size_t j = 0; j < 6; ++j )
    {
        tC( i, j ) = tVal;
        tVal += 1.0;
    }

    belfem::Tensor< belfem::real > tCtor( tC );
    belfem::Tensor< belfem::real > tAssign( 3, 3, 3, 3 );
    tAssign = tC;

    EXPECT_TRUE( tCtor == tAssign );
}

TEST( TensorConversion, SymmetricTensorMappingRoundTrip )
{
    // Voigt mapping only preserves minor-symmetric tensors (A_ijkl = A_jikl = A_ijlk).
    // Use isotropic fill which guarantees full symmetry.
    belfem::Tensor< belfem::real > tA( 3, 3, 3, 3, 0.0 );
    tA.fill_isotropic_elasticity( 200.0e9, 0.3 );

    belfem::Matrix< belfem::real > tC( 6, 6 );
    tA.to_matrix( tC );

    belfem::Tensor< belfem::real > tB( 3, 3, 3, 3 );
    tB = tC;

    EXPECT_TRUE( tA == tB );
}

// =============================================================================
// §2.7  Conversion  [debug]
// =============================================================================

// =============================================================================
// Voigt convention pinning (Codex finding #2)
// =============================================================================

TEST( TensorConversion, IsotropicTensorToVoigtMatrix )
{
    // Pin the Voigt convention: for isotropic steel, the 6×6 stiffness matrix
    // should have known analytical entries. This is NOT a round-trip — it directly
    // checks the ten_to_mat mapping against independently computed values.
    belfem::real tE  = 200.0e9;
    belfem::real tNu = 0.3;
    belfem::real tLambda = tE * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) );
    belfem::real tMu     = tE / ( 2.0 * ( 1.0 + tNu ) );

    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 0.0 );
    tTen.fill_isotropic_elasticity( tE, tNu );

    belfem::Matrix< belfem::real > tC( 6, 6, 0.0 );
    tTen.to_matrix( tC );

    // diagonal: C_11 = C_22 = C_33 = λ + 2μ
    belfem::real tDiag = tLambda + 2.0 * tMu;
    EXPECT_NEAR( tC( 0, 0 ), tDiag,   tTol * tDiag );
    EXPECT_NEAR( tC( 1, 1 ), tDiag,   tTol * tDiag );
    EXPECT_NEAR( tC( 2, 2 ), tDiag,   tTol * tDiag );

    // off-diagonal: C_12 = C_13 = C_23 = λ
    EXPECT_NEAR( tC( 0, 1 ), tLambda, tTol * tDiag );
    EXPECT_NEAR( tC( 0, 2 ), tLambda, tTol * tDiag );
    EXPECT_NEAR( tC( 1, 2 ), tLambda, tTol * tDiag );

    // shear: C_44 = C_55 = C_66 = μ
    EXPECT_NEAR( tC( 3, 3 ), tMu,     tTol * tDiag );
    EXPECT_NEAR( tC( 4, 4 ), tMu,     tTol * tDiag );
    EXPECT_NEAR( tC( 5, 5 ), tMu,     tTol * tDiag );

    // off-diagonal shear should be zero
    EXPECT_NEAR( tC( 0, 3 ), 0.0,     tTol * tDiag );
    EXPECT_NEAR( tC( 3, 4 ), 0.0,     tTol * tDiag );
}

#ifndef NDEBUG

TEST( TensorDebug, ToMatrixWrongSizeThrows )
{
    belfem::Tensor< belfem::real > tTen( 3, 3, 3, 3, 0.0 );
    belfem::Matrix< belfem::real > tM( 5, 6 );
    EXPECT_THROW( tTen.to_matrix( tM ), std::runtime_error );
}

TEST( TensorDebug, ToMatrixNon3333Throws )
{
    belfem::Tensor< belfem::real > tTen( 2, 2, 2, 2, 0.0 );
    belfem::Matrix< belfem::real > tM( 6, 6 );
    EXPECT_THROW( tTen.to_matrix( tM ), std::runtime_error );
}

#endif // NDEBUG
