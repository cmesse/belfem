/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for SpMatrix: construction from dense, element access,
 * CSR/CSC structure, multiply, transpose, COO indices, indexing base.
 * See: tests_08_sparse.md §1–§4
 */

#include <gtest/gtest.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <set>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"

namespace
{
    const belfem::real tEps = 1e-12;
    const belfem::real tTol = 1e-9;

    // Standard 4×4 tridiagonal test matrix:
    //  2 -1  0  0
    // -1  2 -1  0
    //  0 -1  2 -1
    //  0  0 -1  2
    belfem::Matrix< belfem::real > make_tridiag4()
    {
        belfem::uint tN = 4;
        belfem::Matrix< belfem::real > tK( tN, tN, 0.0 );
        for( belfem::uint k = 0; k < tN; ++k )
        {
            tK( k, k ) = 2.0;
            if( k > 0 )     tK( k, k - 1 ) = -1.0;
            if( k < tN - 1 ) tK( k, k + 1 ) = -1.0;
        }
        return tK;
    }

    // 3×3 fully dense test matrix (no zeros)
    belfem::Matrix< belfem::real > make_full3()
    {
        belfem::Matrix< belfem::real > tM( 3, 3 );
        tM( 0, 0 ) = 1.0; tM( 0, 1 ) = 2.0; tM( 0, 2 ) = 3.0;
        tM( 1, 0 ) = 4.0; tM( 1, 1 ) = 5.0; tM( 1, 2 ) = 6.0;
        tM( 2, 0 ) = 7.0; tM( 2, 1 ) = 8.0; tM( 2, 2 ) = 9.0;
        return tM;
    }

    // 3×5 rectangular test matrix
    belfem::Matrix< belfem::real > make_rect35()
    {
        belfem::Matrix< belfem::real > tM( 3, 5, 0.0 );
        tM( 0, 0 ) = 1.0; tM( 0, 2 ) = 2.0; tM( 0, 4 ) = 3.0;
        tM( 1, 1 ) = 4.0; tM( 1, 3 ) = 5.0;
        tM( 2, 0 ) = 6.0; tM( 2, 4 ) = 7.0;
        return tM;
    }

    // Same sparsity as make_tridiag4(), values scaled — the payload used by
    // the HDF5 tests, so a successful load is visible in the values.
    belfem::Matrix< belfem::real > make_tridiag4_scaled( const belfem::real aFactor )
    {
        belfem::Matrix< belfem::real > tK = make_tridiag4();
        for( belfem::uint j = 0; j < 4; ++j )
        {
            for( belfem::uint i = 0; i < 4; ++i )
            {
                tK( i, j ) *= aFactor;
            }
        }
        return tK;
    }

    // 4x4, nnz = 10, row counts 2,3,3,2 — the SAME pointer array as
    // make_tridiag4(), but different column indices, so a CSR load of this
    // pattern fails the index compare and nothing else.
    belfem::Matrix< belfem::real > make_pattern4_indices()
    {
        belfem::Matrix< belfem::real > tM( 4, 4, 0.0 );
        tM( 0, 0 ) = 1.0; tM( 0, 2 ) = 1.0;
        tM( 1, 0 ) = 1.0; tM( 1, 1 ) = 1.0; tM( 1, 3 ) = 1.0;
        tM( 2, 0 ) = 1.0; tM( 2, 2 ) = 1.0; tM( 2, 3 ) = 1.0;
        tM( 3, 1 ) = 1.0; tM( 3, 3 ) = 1.0;
        return tM;
    }

    // 4x4, nnz = 10, row counts 3,2,3,2 — different pointer array, so a CSR
    // load of this pattern fails the pointer compare before the indices.
    belfem::Matrix< belfem::real > make_pattern4_pointers()
    {
        belfem::Matrix< belfem::real > tM( 4, 4, 0.0 );
        tM( 0, 0 ) = 1.0; tM( 0, 1 ) = 1.0; tM( 0, 2 ) = 1.0;
        tM( 1, 0 ) = 1.0; tM( 1, 1 ) = 1.0;
        tM( 2, 1 ) = 1.0; tM( 2, 2 ) = 1.0; tM( 2, 3 ) = 1.0;
        tM( 3, 2 ) = 1.0; tM( 3, 3 ) = 1.0;
        return tM;
    }
}

// =============================================================================
// §1.1  Dense Constructor  [semantic]
// =============================================================================

TEST( SpMatrixConstruction, ConstructCSRFromDense )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    EXPECT_EQ( tM.type(), belfem::SpMatrixType::CSR );
    EXPECT_EQ( tM.n_rows(), 4u );
    EXPECT_EQ( tM.n_cols(), 4u );
}

TEST( SpMatrixConstruction, ConstructCSCFromDense )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSC );

    EXPECT_EQ( tM.type(), belfem::SpMatrixType::CSC );
    EXPECT_EQ( tM.n_rows(), 4u );
    EXPECT_EQ( tM.n_cols(), 4u );
}

TEST( SpMatrixConstruction, NonzeroCountCorrect )
{
    // tridiag4 has 4 diagonal + 3+3 off-diagonal = 10 nonzeros
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    EXPECT_EQ( tM.number_of_nonzeros(), 10u );
}

TEST( SpMatrixConstruction, DiagonalMatrixCSR )
{
    belfem::uint tN = 5;
    belfem::Matrix< belfem::real > tDense( tN, tN, 0.0 );
    for( belfem::uint i = 0; i < tN; ++i )
    {
        tDense( i, i ) = static_cast< belfem::real >( i + 1 );
    }
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    EXPECT_EQ( tM.number_of_nonzeros(), tN );
}

TEST( SpMatrixConstruction, FullDenseMatrixCSR )
{
    belfem::Matrix< belfem::real > tDense = make_full3();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    EXPECT_EQ( tM.number_of_nonzeros(), 9u );
}

TEST( SpMatrixConstruction, RectangularCSR )
{
    belfem::Matrix< belfem::real > tDense = make_rect35();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    EXPECT_EQ( tM.n_rows(), 3u );
    EXPECT_EQ( tM.n_cols(), 5u );
    EXPECT_EQ( tM.number_of_nonzeros(), 7u );
}

TEST( SpMatrixConstruction, RectangularCSC )
{
    belfem::Matrix< belfem::real > tDense = make_rect35();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSC );

    EXPECT_EQ( tM.n_rows(), 3u );
    EXPECT_EQ( tM.n_cols(), 5u );
    EXPECT_EQ( tM.number_of_nonzeros(), 7u );
}

// =============================================================================
// §1.2  Dense Round-Trip  [semantic]
// =============================================================================

TEST( SpMatrixRoundTrip, CSRReadBackMatchesDense )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    // need indexing base set for operator() to work
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    const belfem::SpMatrix & tConst = tM;
    for( belfem::index_t i = 0; i < 4; ++i )
    for( belfem::index_t j = 0; j < 4; ++j )
    {
        EXPECT_NEAR( tConst( i, j ), tDense( i, j ), tEps );
    }
}

TEST( SpMatrixRoundTrip, CSCReadBackMatchesDense )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSC );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    const belfem::SpMatrix & tConst = tM;
    for( belfem::index_t i = 0; i < 4; ++i )
    for( belfem::index_t j = 0; j < 4; ++j )
    {
        EXPECT_NEAR( tConst( i, j ), tDense( i, j ), tEps );
    }
}

TEST( SpMatrixRoundTrip, CSRAndCSCAgree )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();

    belfem::SpMatrix tCSR( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tCSC( tDense, belfem::SpMatrixType::CSC );

    tCSR.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    tCSC.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    const belfem::SpMatrix & tConstCSR = tCSR;
    const belfem::SpMatrix & tConstCSC = tCSC;

    for( belfem::index_t i = 0; i < 4; ++i )
    for( belfem::index_t j = 0; j < 4; ++j )
    {
        EXPECT_NEAR( tConstCSR( i, j ), tConstCSC( i, j ), tEps );
    }
}

// =============================================================================
// §1.3  Pointer/Index Structure  [semantic]
// =============================================================================

TEST( SpMatrixStructure, CSRPointersMonotone )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    for( belfem::index_t i = 0; i < tM.n_rows(); ++i )
    {
        EXPECT_LE( tM.pointers()[ i ], tM.pointers()[ i + 1 ] );
    }
}

TEST( SpMatrixStructure, CSRPointersLastEqualsNNZ )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    EXPECT_EQ( tM.pointers()[ tM.n_rows() ],
               static_cast< belfem::int_t >( tM.number_of_nonzeros() ) );
}

TEST( SpMatrixStructure, CSCPointersLastEqualsNNZ )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSC );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    EXPECT_EQ( tM.pointers()[ tM.n_cols() ],
               static_cast< belfem::int_t >( tM.number_of_nonzeros() ) );
}

TEST( SpMatrixStructure, CSRIndicesSorted )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    const belfem::int_t * tPtr = tM.pointers();
    const belfem::int_t * tIdx = tM.indices();

    for( belfem::index_t i = 0; i < tM.n_rows(); ++i )
    {
        for( belfem::int_t k = tPtr[ i ] + 1; k < tPtr[ i + 1 ]; ++k )
        {
            EXPECT_LT( tIdx[ k - 1 ], tIdx[ k ] );
        }
    }
}

// =============================================================================
// §2.1  Element Access  [semantic]
// =============================================================================

TEST( SpMatrixAccess, ConstAccessStructuralZero )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    const belfem::SpMatrix & tConst = tM;
    // (0,2) is a structural zero in tridiag4
    EXPECT_NEAR( tConst( 0, 2 ), 0.0, tEps );
}

TEST( SpMatrixAccess, FillSetsAllValues )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    tM.fill( 3.14 );

    // all nonzero positions should be 3.14
    EXPECT_NEAR( tM( 0, 0 ), 3.14, tEps );
    EXPECT_NEAR( tM( 0, 1 ), 3.14, tEps );
    EXPECT_NEAR( tM( 1, 0 ), 3.14, tEps );

    // structural zeros still return 0.0
    const belfem::SpMatrix & tConst = tM;
    EXPECT_NEAR( tConst( 0, 2 ), 0.0, tEps );
}

TEST( SpMatrixAccess, IndexReturnsNNZForZeroPosition )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    // (0,3) is structural zero → index should return nnz
    int tIdx = tM.index( 0, 3 );
    EXPECT_EQ( tIdx, static_cast< int >( tM.number_of_nonzeros() ) );
}

TEST( SpMatrixAccess, IndexReturnValidForNonzero )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    // (0,0) is nonzero → index should be < nnz
    int tIdx = tM.index( 0, 0 );
    EXPECT_LT( tIdx, static_cast< int >( tM.number_of_nonzeros() ) );
    EXPECT_GE( tIdx, 0 );
}

// --- §2.2 Element Access [debug] ---

#ifndef NDEBUG

TEST( SpMatrixDebug, WritableAccessStructuralZeroThrows )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    // (0,3) is structural zero — writable access should assert
    EXPECT_THROW( tM( 0, 3 ) = 1.0, std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §2.3  Indexing Base Conversion  [semantic]
// =============================================================================

TEST( SpMatrixIndexing, CppToFortranRoundTrip )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_EQ( tM.indexing_base(), 0 );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
    EXPECT_EQ( tM.indexing_base(), 1 );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_EQ( tM.indexing_base(), 0 );
}

TEST( SpMatrixIndexing, FortranBasePointersShifted )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_EQ( tM.pointers()[ 0 ], 0 );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
    EXPECT_EQ( tM.pointers()[ 0 ], 1 );
}

TEST( SpMatrixIndexing, ElementAccessAfterFortranRoundTrip )
{
    // Convert to Fortran, then back to Cpp — verify data survives the round-trip
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
    EXPECT_EQ( tM.indexing_base(), 1 );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_EQ( tM.indexing_base(), 0 );

    const belfem::SpMatrix & tConst = tM;
    EXPECT_NEAR( tConst( 0, 0 ), 2.0, tEps );
    EXPECT_NEAR( tConst( 0, 1 ), -1.0, tEps );
    EXPECT_NEAR( tConst( 0, 2 ), 0.0, tEps );
}

TEST( SpMatrixIndexing, DoubleConversionIsIdempotent )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    EXPECT_EQ( tM.pointers()[ 0 ], 0 );

    const belfem::SpMatrix & tConst = tM;
    EXPECT_NEAR( tConst( 0, 0 ), 2.0, tEps );
}

// =============================================================================
// §3.1  Multiply  [semantic]
// =============================================================================

TEST( SpMatrixMultiply, MultiplyCSRIdentity )
{
    belfem::uint tN = 4;
    belfem::Matrix< belfem::real > tI( tN, tN, 0.0 );
    for( belfem::uint i = 0; i < tN; ++i ) tI( i, i ) = 1.0;

    belfem::SpMatrix tM( tI, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tY( tN, 0.0 );

    tM.multiply( tX, tY );

    for( belfem::uint i = 0; i < tN; ++i )
    {
        EXPECT_NEAR( tY( i ), tX( i ), tEps );
    }
}

TEST( SpMatrixMultiply, MultiplyCSRTridiagonal )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tY( 4, 0.0 );

    // compute expected: b = K * x via dense
    belfem::Vector< belfem::real > tExpected( tDense * tX );

    tM.multiply( tX, tY );

    for( belfem::uint i = 0; i < 4; ++i )
    {
        EXPECT_NEAR( tY( i ), tExpected( i ), tEps );
    }
}

TEST( SpMatrixMultiply, MultiplyCSCTridiagonal )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSC );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tY( 4, 0.0 );
    belfem::Vector< belfem::real > tExpected( tDense * tX );

    tM.multiply( tX, tY );

    for( belfem::uint i = 0; i < 4; ++i )
    {
        EXPECT_NEAR( tY( i ), tExpected( i ), tEps );
    }
}

TEST( SpMatrixMultiply, MultiplyRectangular )
{
    belfem::Matrix< belfem::real > tDense = make_rect35();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    belfem::Vector< belfem::real > tY( 3, 0.0 );
    belfem::Vector< belfem::real > tExpected( tDense * tX );

    tM.multiply( tX, tY );

    for( belfem::uint i = 0; i < 3; ++i )
    {
        EXPECT_NEAR( tY( i ), tExpected( i ), tEps );
    }
}

TEST( SpMatrixMultiply, MultiplyTransposed )
{
    // (idea from ChatGPT) — test aTransposedFlag=true path
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tY( 4, 0.0 );

    // compute expected: y = Aᵀ * x (tridiag is symmetric so Aᵀ*x == A*x)
    belfem::Vector< belfem::real > tExpected( tDense * tX );

    tM.multiply( tX, tY, 1.0, 0.0, true );

    for( belfem::uint i = 0; i < 4; ++i )
    {
        EXPECT_NEAR( tY( i ), tExpected( i ), tEps );
    }
}

TEST( SpMatrixMultiply, MultiplyTransposedRectangularCSR )
{
    // non-symmetric rectangular case: y = Aᵀ * x with A 3x5, so a silently
    // non-transposed product cannot pass by accident
    belfem::Matrix< belfem::real > tDense = make_rect35();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, -2.0, 3.0 };
    belfem::Vector< belfem::real > tY( 5, 0.0 );

    // expected: y_j = sum_i A(i,j) * x_i
    belfem::Vector< belfem::real > tExpected( 5, 0.0 );
    for( belfem::uint j = 0; j < 5; ++j )
    {
        for( belfem::uint i = 0; i < 3; ++i )
        {
            tExpected( j ) += tDense( i, j ) * tX( i );
        }
    }

    tM.multiply( tX, tY, 1.0, 0.0, true );

    for( belfem::uint j = 0; j < 5; ++j )
    {
        EXPECT_NEAR( tY( j ), tExpected( j ), tEps );
    }
}

TEST( SpMatrixMultiply, MultiplyTransposedRectangularCSC )
{
    belfem::Matrix< belfem::real > tDense = make_rect35();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSC );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, -2.0, 3.0 };
    belfem::Vector< belfem::real > tY( 5, 0.0 );

    belfem::Vector< belfem::real > tExpected( 5, 0.0 );
    for( belfem::uint j = 0; j < 5; ++j )
    {
        for( belfem::uint i = 0; i < 3; ++i )
        {
            tExpected( j ) += tDense( i, j ) * tX( i );
        }
    }

    tM.multiply( tX, tY, 1.0, 0.0, true );

    for( belfem::uint j = 0; j < 5; ++j )
    {
        EXPECT_NEAR( tY( j ), tExpected( j ), tEps );
    }
}

TEST( SpMatrixMultiply, MultiplyAlphaBeta )
{
    // y := alpha * A * x + beta * y with y preset, alpha = 2, beta = 3
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tY = { 1.0, -1.0, 0.5, 2.0 };

    belfem::Vector< belfem::real > tAx( tDense * tX );
    belfem::Vector< belfem::real > tExpected( 4 );
    for( belfem::uint i = 0; i < 4; ++i )
    {
        tExpected( i ) = 2.0 * tAx( i ) + 3.0 * tY( i );
    }

    tM.multiply( tX, tY, 2.0, 3.0 );

    for( belfem::uint i = 0; i < 4; ++i )
    {
        EXPECT_NEAR( tY( i ), tExpected( i ), tEps );
    }
}

TEST( SpMatrixMultiply, MultiplyRestoresIndexingBase )
{
    // Contract: multiply() must restore the indexing base to what it was
    // before the call. Internal Fortran conversion is transparent to the caller.
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    belfem::Vector< belfem::real > tX = { 1.0, 2.0, 3.0, 4.0 };
    belfem::Vector< belfem::real > tY( 4, 0.0 );

    tM.multiply( tX, tY );

    // indexing base must be restored to Cpp (0)
    EXPECT_EQ( tM.indexing_base(), 0 );

    // element access should work directly — no manual restoration needed
    const belfem::SpMatrix & tConst = tM;
    EXPECT_NEAR( tConst( 0, 0 ), 2.0, tEps );
    EXPECT_NEAR( tConst( 0, 1 ), -1.0, tEps );
    EXPECT_NEAR( tConst( 1, 0 ), -1.0, tEps );
    EXPECT_NEAR( tConst( 0, 2 ), 0.0, tEps );   // structural zero
}

// =============================================================================
// §3.2  Transpose  [semantic]
// =============================================================================

TEST( SpMatrixTranspose, TransposeSwapsType )
{
    belfem::Matrix< belfem::real > tDense = make_rect35();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.transpose();

    EXPECT_EQ( tM.type(), belfem::SpMatrixType::CSC );
}

TEST( SpMatrixTranspose, TransposeSwapsDimensions )
{
    belfem::Matrix< belfem::real > tDense = make_rect35();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.transpose();

    EXPECT_EQ( tM.n_rows(), 5u );
    EXPECT_EQ( tM.n_cols(), 3u );
}

TEST( SpMatrixTranspose, TransposePreservesValues )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    tM.transpose();
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    const belfem::SpMatrix & tConst = tM;
    // tridiag is symmetric → A(i,j) == Aᵀ(i,j) == A(j,i)
    for( belfem::index_t i = 0; i < 4; ++i )
    for( belfem::index_t j = 0; j < 4; ++j )
    {
        EXPECT_NEAR( tConst( i, j ), tDense( j, i ), tEps );
    }
}

TEST( SpMatrixTranspose, TransposeNonSymmetricRectangular )
{
    // Codex finding: previous transpose tests used symmetric fixture.
    // Use non-symmetric rectangular matrix to truly verify transpose.
    belfem::Matrix< belfem::real > tDense = make_rect35();  // 3×5
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    tM.transpose();
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    // after transpose: 5×3
    EXPECT_EQ( tM.n_rows(), 5u );
    EXPECT_EQ( tM.n_cols(), 3u );

    const belfem::SpMatrix & tConst = tM;
    for( belfem::index_t i = 0; i < 5; ++i )
    for( belfem::index_t j = 0; j < 3; ++j )
    {
        EXPECT_NEAR( tConst( i, j ), tDense( j, i ), tEps );
    }
}

TEST( SpMatrixTranspose, DoubleTransposeRoundTrip )
{
    belfem::Matrix< belfem::real > tDense = make_rect35();  // 3×5
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    tM.transpose();
    tM.transpose();
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );

    // back to 3×5
    EXPECT_EQ( tM.n_rows(), 3u );
    EXPECT_EQ( tM.n_cols(), 5u );

    const belfem::SpMatrix & tConst = tM;
    for( belfem::index_t i = 0; i < 3; ++i )
    for( belfem::index_t j = 0; j < 5; ++j )
    {
        EXPECT_NEAR( tConst( i, j ), tDense( i, j ), tEps );
    }
}

// =============================================================================
// §3.3  COO Indices  [semantic]
// =============================================================================

TEST( SpMatrixCOO, CreateCOOForCSR )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    EXPECT_FALSE( tM.have_coo_indices() );

    tM.create_coo_indices();

    EXPECT_TRUE( tM.have_coo_indices() );
    EXPECT_NE( tM.rows(), nullptr );

    // (idea from ChatGPT) — verify COO entries map to correct values
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    const belfem::SpMatrix & tConst = tM;
    for( belfem::index_t k = 0; k < tM.number_of_nonzeros(); ++k )
    {
        EXPECT_NEAR( tConst( tM.rows()[ k ], tM.cols()[ k ] ),
                     tM.data( k ), tEps );
    }
}

TEST( SpMatrixCOO, CreateCOOForCSC )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSC );

    tM.create_coo_indices();

    EXPECT_TRUE( tM.have_coo_indices() );
    EXPECT_NE( tM.cols(), nullptr );
}

TEST( SpMatrixCOO, FreeCOOIndices )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.create_coo_indices();
    ASSERT_TRUE( tM.have_coo_indices() );

    tM.free_coo_indices();

    EXPECT_FALSE( tM.have_coo_indices() );
}

// =============================================================================
// §5  Parent/Child Structure Sharing  [semantic]
// =============================================================================

TEST( SpMatrixParentChild, ChildSharesStructureAndOwnsValues )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    // same shape and pattern
    EXPECT_EQ( tM.n_rows(), tK.n_rows() );
    EXPECT_EQ( tM.n_cols(), tK.n_cols() );
    EXPECT_EQ( tM.number_of_nonzeros(), tK.number_of_nonzeros() );
    EXPECT_EQ( tM.type(), tK.type() );

    // structure arrays are aliased, values are not
    EXPECT_EQ( tM.pointers(), tK.pointers() );
    EXPECT_EQ( tM.cols(),     tK.cols() );
    EXPECT_NE( tM.data(),     tK.data() );

    // child values are zero-initialized, parent values untouched
    for( belfem::index_t k = 0; k < tM.number_of_nonzeros(); ++k )
    {
        EXPECT_EQ( tM.data( k ), 0.0 );
    }
    EXPECT_NEAR( tK( 0, 0 ), 2.0, tEps );
}

TEST( SpMatrixParentChild, ValuesAreIndependent )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    tM.fill( 5.0 );

    // parent values must be unaffected
    EXPECT_NEAR( tK( 0, 0 ),  2.0, tEps );
    EXPECT_NEAR( tK( 0, 1 ), -1.0, tEps );
    EXPECT_NEAR( tM( 0, 0 ),  5.0, tEps );
}

TEST( SpMatrixParentChild, ParentInitiatedBaseFlipReachesChild )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );
    tM.fill( 5.0 );

    // flip on the parent: the shared arrays convert once, the child's index
    // function must follow
    tK.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
    EXPECT_EQ( tM.indexing_base(), 1 );

    tK.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_EQ( tM.indexing_base(), 0 );

    // child element access must still be consistent after the round trip
    const belfem::SpMatrix & tConst = tM;
    EXPECT_NEAR( tConst( 1, 1 ), 5.0, tEps );
    EXPECT_NEAR( tConst( 0, 3 ), 0.0, tEps );  // structural zero
}

TEST( SpMatrixParentChild, ChildInitiatedBaseFlipReachesParent )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
    EXPECT_EQ( tK.indexing_base(), 1 );

    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_EQ( tK.indexing_base(), 0 );

    const belfem::SpMatrix & tConst = tK;
    EXPECT_NEAR( tConst( 0, 0 ), 2.0, tEps );
}

TEST( SpMatrixParentChild, CooParentFirst )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    // parent creates first; child must end up aliasing, not allocating
    tK.create_coo_indices();
    EXPECT_TRUE( tM.have_coo_indices() );
    EXPECT_EQ( tM.rows(), tK.rows() );

    // freeing from either side clears both sides
    tM.free_coo_indices();
    EXPECT_FALSE( tK.have_coo_indices() );
    EXPECT_FALSE( tM.have_coo_indices() );
    EXPECT_EQ( tK.rows(), nullptr );
    EXPECT_EQ( tM.rows(), nullptr );
}

TEST( SpMatrixParentChild, CooChildFirst )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    // child creates first; parent computes, child aliases
    tM.create_coo_indices();
    EXPECT_TRUE( tK.have_coo_indices() );
    EXPECT_EQ( tM.rows(), tK.rows() );
    EXPECT_NE( tM.rows(), nullptr );

    // free from the parent side; the child must drop its alias too
    tK.free_coo_indices();
    EXPECT_FALSE( tM.have_coo_indices() );
    EXPECT_EQ( tM.rows(), nullptr );

    // re-create after a free cycle must re-alias cleanly
    tK.create_coo_indices();
    EXPECT_TRUE( tM.have_coo_indices() );
    EXPECT_EQ( tM.rows(), tK.rows() );
}

TEST( SpMatrixParentChild, MemoryAccounting )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    // child owns only its values
    EXPECT_EQ( tM.memory(),
               tM.number_of_nonzeros() * sizeof( belfem::real ) );

    // parent owns values, pointer array and structural indices
    EXPECT_EQ( tK.memory(),
               tK.number_of_nonzeros() * ( sizeof( belfem::real ) + sizeof( belfem::int_t ) )
             + tK.n_pointers() * sizeof( belfem::int_t ) );
}

TEST( SpMatrixParentChild, ChildDestroyedFirst )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix * tK = new belfem::SpMatrix( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix * tM = new belfem::SpMatrix( tK );

    delete tM;

    // the parent must be fully usable afterwards, including base flips
    // ( which would touch a stale child link if unlinking failed )
    tK->set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
    tK->set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    const belfem::SpMatrix & tConst = *tK;
    EXPECT_NEAR( tConst( 2, 2 ), 2.0, tEps );

    // and it must accept a new child again
    belfem::SpMatrix * tM2 = new belfem::SpMatrix( tK );
    EXPECT_EQ( tM2->pointers(), tK->pointers() );
    delete tM2;
    delete tK;
}

TEST( SpMatrixParentChild, ParentDestroyedFirst )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix * tK = new belfem::SpMatrix( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix * tM = new belfem::SpMatrix( tK );
    tM->fill( 3.0 );

    // ownership of the structure transfers to the child
    delete tK;

    EXPECT_EQ( tM->n_rows(), 4u );
    EXPECT_EQ( tM->number_of_nonzeros(), 10u );
    const belfem::SpMatrix & tConst = *tM;
    EXPECT_NEAR( tConst( 1, 0 ), 3.0, tEps );

    // multiply exercises pointers, indices and values together
    belfem::Vector< belfem::real > tX( 4, 1.0 );
    belfem::Vector< belfem::real > tY( 4, 0.0 );
    tM->multiply( tX, tY );
    EXPECT_NEAR( tY( 0 ), 6.0, tEps );   // two entries in row 0
    EXPECT_NEAR( tY( 1 ), 9.0, tEps );   // three entries in row 1

    delete tM;
}

TEST( SpMatrixParentChild, CopyFromChildIsUnlinkedDeepCopy )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );
    tM.fill( 4.0 );

    belfem::SpMatrix tCopy;
    tCopy = tM;

    // deep copy: own structure, own values
    EXPECT_NE( tCopy.pointers(), tK.pointers() );
    EXPECT_NE( tCopy.data(), tM.data() );
    EXPECT_EQ( tCopy.number_of_nonzeros(), tM.number_of_nonzeros() );
    const belfem::SpMatrix & tConst = tCopy;
    EXPECT_NEAR( tConst( 0, 1 ), 4.0, tEps );

    // mutating the copy must not touch parent or child
    tCopy.fill( 9.0 );
    EXPECT_NEAR( tM( 0, 1 ), 4.0, tEps );
    EXPECT_NEAR( tK( 0, 1 ), -1.0, tEps );
}

TEST( SpMatrixParentChild, CscChildSharesCooColumns )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSC );
    belfem::SpMatrix tM( &tK );

    // for CSC the structural array is mRows and the COO array is mColumns
    EXPECT_EQ( tM.rows(), tK.rows() );

    tK.create_coo_indices();
    EXPECT_TRUE( tM.have_coo_indices() );
    EXPECT_EQ( tM.cols(), tK.cols() );
    EXPECT_NE( tM.cols(), nullptr );

    tM.free_coo_indices();
    EXPECT_FALSE( tK.have_coo_indices() );
    EXPECT_EQ( tK.cols(), nullptr );
    EXPECT_EQ( tM.cols(), nullptr );
}

TEST( SpMatrixParentChild, CooExistsBeforeChildConstruction )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    tK.create_coo_indices();

    // the child must inherit the flag and the alias at construction
    belfem::SpMatrix tM( &tK );
    EXPECT_TRUE( tM.have_coo_indices() );
    EXPECT_EQ( tM.rows(), tK.rows() );
    EXPECT_NE( tM.rows(), nullptr );
}

TEST( SpMatrixParentChild, ParentInFortranBaseAtChildConstruction )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    tK.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );

    // the child must come up in the parent's base
    belfem::SpMatrix tM( &tK );
    EXPECT_EQ( tM.indexing_base(), 1 );

    // and a flip back through the child must reach both sides
    tM.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_EQ( tK.indexing_base(), 0 );
    const belfem::SpMatrix & tConst = tK;
    EXPECT_NEAR( tConst( 0, 0 ), 2.0, tEps );
}

TEST( SpMatrixParentChild, MemoryAccountingWithCoo )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    tK.create_coo_indices();

    // the coo array is charged to the parent only
    EXPECT_EQ( tM.memory(),
               tM.number_of_nonzeros() * sizeof( belfem::real ) );
    EXPECT_EQ( tK.memory(),
               tK.number_of_nonzeros() * ( sizeof( belfem::real ) + 2 * sizeof( belfem::int_t ) )
             + tK.n_pointers() * sizeof( belfem::int_t ) );
}

// =============================================================================
// §6  Parent/Child HDF5 Round Trip  [semantic]
//
// Covers SpMatrix::load() on LINKED matrices — the verify-only path used by
// SolverData::load_system, where the structure is never replaced and only the
// values are read. The error paths are BELFEM_ERROR and therefore catchable:
// tests/sparse/test_sparse_main.cpp arms set_throw_on_error( true ).
// =============================================================================

#ifdef BELFEM_HDF5

#include <cstdio>
#include "filetools.hpp"

namespace
{
    void
    remove_file_if_exists( const std::string & aPath )
    {
        if ( belfem::file_exists( aPath ) )
        {
            std::remove( aPath.c_str() );
        }
    }

    // elementwise structure comparison, for matrices that own separate arrays
    void
    expect_same_structure( const belfem::SpMatrix & aA, const belfem::SpMatrix & aB )
    {
        ASSERT_EQ( aA.n_pointers(), aB.n_pointers() );
        ASSERT_EQ( aA.number_of_nonzeros(), aB.number_of_nonzeros() );

        for( belfem::index_t k = 0; k < aA.n_pointers(); ++k )
        {
            EXPECT_EQ( aA.pointers()[ k ], aB.pointers()[ k ] );
        }
        for( belfem::index_t k = 0; k < aA.number_of_nonzeros(); ++k )
        {
            EXPECT_EQ( aA.indices()[ k ], aB.indices()[ k ] );
        }
    }
}

class SpMatrixHDF5Test : public ::testing::Test
{
protected:

    std::string mPath;

    // unique per test, no path separators: the file is created in the cwd
    void
    SetUp() override
    {
        const ::testing::TestInfo * tInfo =
                ::testing::UnitTest::GetInstance()->current_test_info();

        mPath = std::string( "test_spmatrix_" )
              + tInfo->test_suite_name()
              + "_"
              + tInfo->name()
              + ".hdf5";

        remove_file_if_exists( mPath );
    }

    void
    TearDown() override
    {
        remove_file_if_exists( mPath );
    }

    // helper: write a dense pattern to mPath as a sparse matrix of aType
    void
    save_matrix( const belfem::Matrix< belfem::real > & aDense,
                 const belfem::SpMatrixType aType = belfem::SpMatrixType::CSR )
    {
        belfem::SpMatrix tSource( aDense, aType );
        tSource.save( mPath, "Matrix", belfem::FileMode::NEW );
    }
};

// -----------------------------------------------------------------------------
// regression anchor: the UNLINKED branch must still rebuild from file
// -----------------------------------------------------------------------------

TEST_F( SpMatrixHDF5Test, UnlinkedRoundTrip )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tSource( tDense, belfem::SpMatrixType::CSR );
    tSource.save( mPath, "Matrix", belfem::FileMode::NEW );

    belfem::SpMatrix tLoaded;
    tLoaded.load( mPath, "Matrix" );

    EXPECT_EQ( tLoaded.type(), belfem::SpMatrixType::CSR );
    EXPECT_EQ( tLoaded.n_rows(), tSource.n_rows() );
    EXPECT_EQ( tLoaded.n_cols(), tSource.n_cols() );
    expect_same_structure( tLoaded, tSource );

    for( belfem::index_t k = 0; k < tLoaded.number_of_nonzeros(); ++k )
    {
        EXPECT_NEAR( tLoaded.data( k ), tSource.data( k ), tEps );
    }
}

// -----------------------------------------------------------------------------
// linked loads
// -----------------------------------------------------------------------------

TEST_F( SpMatrixHDF5Test, ChildLoadsMatchingFile )
{
    // file holds the same pattern with x10 values
    this->save_matrix( make_tridiag4_scaled( 10.0 ) );

    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    tM.load( mPath, "Matrix" );

    // the child received the file values
    EXPECT_NEAR( tM( 0, 0 ), 20.0, tEps );
    EXPECT_NEAR( tM( 0, 1 ), -10.0, tEps );

    // the parent is untouched
    EXPECT_NEAR( tK( 0, 0 ), 2.0, tEps );
    EXPECT_NEAR( tK( 0, 1 ), -1.0, tEps );

    // the structure is still shared, the values are not
    EXPECT_EQ( tM.pointers(), tK.pointers() );
    EXPECT_EQ( tM.indices(),  tK.indices() );
    EXPECT_NE( tM.data(),     tK.data() );

    // and the pair is still linked: a flip on the parent reaches the child
    tK.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
    EXPECT_EQ( tM.indexing_base(), 1 );
    tK.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_NEAR( tM( 1, 1 ), 20.0, tEps );
}

TEST_F( SpMatrixHDF5Test, LinkedParentLoadsValuesOnly )
{
    // this is the mechanism behind SolverData::load_system, which loads into
    // mSystemMatrix while mJacobianMatrix is attached as its child
    this->save_matrix( make_tridiag4_scaled( 10.0 ) );

    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );
    tM.fill( 7.0 );

    belfem::int_t * tPointersBefore = tK.pointers();
    belfem::int_t * tIndicesBefore  = tK.indices();

    tK.load( mPath, "Matrix" );

    // the parent took the file values
    EXPECT_NEAR( tK( 0, 0 ), 20.0, tEps );

    // the child kept its own values
    EXPECT_NEAR( tM( 0, 0 ), 7.0, tEps );

    // no structure was reallocated, and the link still holds
    EXPECT_EQ( tK.pointers(), tPointersBefore );
    EXPECT_EQ( tK.indices(),  tIndicesBefore );
    EXPECT_EQ( tM.pointers(), tK.pointers() );
    EXPECT_NE( tM.data(),     tK.data() );

    tK.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
    EXPECT_EQ( tM.indexing_base(), 1 );
    tK.set_indexing_base( belfem::SpMatrixIndexingBase::Cpp );
    EXPECT_NEAR( tM( 2, 2 ), 7.0, tEps );
}

TEST_F( SpMatrixHDF5Test, ChildLoadsOtherBaseFile )
{
    // the file is written by a SEPARATE, unlinked matrix in Fortran base;
    // flipping the pair itself would convert the shared arrays and leave
    // nothing for the shift to absorb
    {
        belfem::SpMatrix tSource( make_tridiag4_scaled( 10.0 ), belfem::SpMatrixType::CSR );
        tSource.set_indexing_base( belfem::SpMatrixIndexingBase::Fortran );
        tSource.save( mPath, "Matrix", belfem::FileMode::NEW );
    }

    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );
    ASSERT_EQ( tK.indexing_base(), 0 );

    // the compare tolerates the base offset, so this must succeed
    tM.load( mPath, "Matrix" );

    // the shared arrays stay in c++ base, the values arrived
    EXPECT_EQ( tM.indexing_base(), 0 );
    EXPECT_EQ( tK.indexing_base(), 0 );
    EXPECT_NEAR( tM( 0, 0 ), 20.0, tEps );
    EXPECT_NEAR( tK( 0, 0 ),  2.0, tEps );
}

// -----------------------------------------------------------------------------
// rejected loads — each pattern isolates exactly one check
// -----------------------------------------------------------------------------

TEST_F( SpMatrixHDF5Test, ChildRejectsIndexMismatch )
{
    // same dimensions, same nnz, same pointer array, different columns
    this->save_matrix( make_pattern4_indices() );

    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    EXPECT_THROW( tM.load( mPath, "Matrix" ), std::runtime_error );

    // the checks run before the values are read, so nothing was written,
    // and the failed load must not have touched the shared structure
    for( belfem::index_t k = 0; k < tM.number_of_nonzeros(); ++k )
    {
        EXPECT_EQ( tM.data( k ), 0.0 );
    }
    EXPECT_NEAR( tK( 0, 0 ), 2.0, tEps );
    EXPECT_EQ( tM.pointers(), tK.pointers() );
    EXPECT_EQ( tM.indices(),  tK.indices() );
}

TEST_F( SpMatrixHDF5Test, ChildRejectsPointerMismatch )
{
    // same dimensions and nnz, different row distribution
    this->save_matrix( make_pattern4_pointers() );

    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    EXPECT_THROW( tM.load( mPath, "Matrix" ), std::runtime_error );

    for( belfem::index_t k = 0; k < tM.number_of_nonzeros(); ++k )
    {
        EXPECT_EQ( tM.data( k ), 0.0 );
    }
}

TEST_F( SpMatrixHDF5Test, ChildRejectsDimensionMismatch )
{
    this->save_matrix( make_full3() );

    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    EXPECT_THROW( tM.load( mPath, "Matrix" ), std::runtime_error );

    EXPECT_EQ( tM.n_rows(), 4u );
    EXPECT_EQ( tM.number_of_nonzeros(), 10u );
}

TEST_F( SpMatrixHDF5Test, ChildRejectsFormatMismatch )
{
    // the tridiagonal matrix is symmetric, so its CSC dump carries the very
    // same pointer and index arrays as the CSR one — the format string in the
    // file is the only thing that distinguishes them, and the type check is
    // the only guard that reads it
    this->save_matrix( make_tridiag4_scaled( 10.0 ), belfem::SpMatrixType::CSC );

    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );

    EXPECT_THROW( tM.load( mPath, "Matrix" ), std::runtime_error );

    EXPECT_EQ( tM.type(), belfem::SpMatrixType::CSR );
    for( belfem::index_t k = 0; k < tM.number_of_nonzeros(); ++k )
    {
        EXPECT_EQ( tM.data( k ), 0.0 );
    }
}

// -----------------------------------------------------------------------------
// saving from a child, and the orphaned owner
// -----------------------------------------------------------------------------

TEST_F( SpMatrixHDF5Test, SaveFromChildRoundTrips )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );
    tM.fill( 3.5 );

    // a child must write its inherited structure and its own values
    tM.save( mPath, "Matrix", belfem::FileMode::NEW );

    belfem::SpMatrix tLoaded;
    tLoaded.load( mPath, "Matrix" );

    EXPECT_EQ( tLoaded.n_rows(), 4u );
    EXPECT_EQ( tLoaded.number_of_nonzeros(), 10u );
    expect_same_structure( tLoaded, tK );

    for( belfem::index_t k = 0; k < tLoaded.number_of_nonzeros(); ++k )
    {
        EXPECT_NEAR( tLoaded.data( k ), 3.5, tEps );
    }
}

TEST_F( SpMatrixHDF5Test, OrphanedOwnerLoadRebuilds )
{
    // the ownership transfer happens when the PARENT dies first, which needs
    // heap allocation: on the stack the child would be destroyed first
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix * tK = new belfem::SpMatrix( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix * tM = new belfem::SpMatrix( tK );

    delete tK;

    // the orphan is an ordinary unlinked matrix again, so load() rebuilds it
    // from file — even with a different sparsity pattern
    this->save_matrix( make_pattern4_pointers() );
    tM->load( mPath, "Matrix" );

    EXPECT_EQ( tM->n_rows(), 4u );
    EXPECT_EQ( tM->number_of_nonzeros(), 10u );

    // row 0 of the new pattern has three entries, the tridiagonal had two
    EXPECT_EQ( tM->pointers()[ 1 ], 3 );
    const belfem::SpMatrix & tConst = *tM;
    EXPECT_NEAR( tConst( 0, 2 ), 1.0, tEps );

    delete tM;
}

#endif

// =============================================================================
// §7  Transpose Rebinding and Linked Multiply  [semantic]
// =============================================================================

TEST( SpMatrixTranspose, AccessAfterTransposeWithoutRebase )
{
    // the lock for the stale-index-function defect: transpose() must leave a
    // usable matrix WITHOUT the caller re-asserting the indexing base — the
    // older transpose tests all did that, which masked the bug
    belfem::Matrix< belfem::real > tDense = make_rect35();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.transpose();

    EXPECT_EQ( tM.type(), belfem::SpMatrixType::CSC );
    EXPECT_EQ( tM.n_rows(), 5u );
    EXPECT_EQ( tM.n_cols(), 3u );

    const belfem::SpMatrix & tConst = tM;
    for( belfem::index_t i = 0; i < 5; ++i )
    for( belfem::index_t j = 0; j < 3; ++j )
    {
        EXPECT_NEAR( tConst( i, j ), tDense( j, i ), tEps );
    }
}

TEST( SpMatrixTranspose, TransposeDropsCooIndices )
{
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tM( tDense, belfem::SpMatrixType::CSR );

    tM.create_coo_indices();
    ASSERT_TRUE( tM.have_coo_indices() );
    ASSERT_NE( tM.rows(), nullptr );

    // the transpose consumes the coo array ( it becomes the structural one of
    // the new type ), so the flag must not survive
    tM.transpose();

    EXPECT_FALSE( tM.have_coo_indices() );
    EXPECT_EQ( tM.cols(), nullptr );      // coo of the new CSC type is absent
    EXPECT_NE( tM.rows(), nullptr );      // structural array of CSC is present

    // and a fresh creation on the transposed matrix must work
    tM.create_coo_indices();
    EXPECT_TRUE( tM.have_coo_indices() );
    EXPECT_NE( tM.cols(), nullptr );
}

TEST( SpMatrixParentChild, LinkedMultiply )
{
    // the production pattern: SolverData multiplies while the pair is linked.
    // multiply() flips the SHARED arrays to fortran base and back, so both
    // index functions must survive the round trip on both call sides.
    belfem::Matrix< belfem::real > tDense = make_tridiag4();
    belfem::SpMatrix tK( tDense, belfem::SpMatrixType::CSR );
    belfem::SpMatrix tM( &tK );
    tM.fill( 3.0 );

    belfem::Vector< belfem::real > tX( 4, 1.0 );
    belfem::Vector< belfem::real > tY( 4, 0.0 );

    // multiply on the parent: tridiagonal row sums are 1, 0, 0, 1
    tK.multiply( tX, tY );
    EXPECT_NEAR( tY( 0 ), 1.0, tEps );
    EXPECT_NEAR( tY( 1 ), 0.0, tEps );
    EXPECT_NEAR( tY( 2 ), 0.0, tEps );
    EXPECT_NEAR( tY( 3 ), 1.0, tEps );

    // base restored on both sides
    EXPECT_EQ( tK.indexing_base(), 0 );
    EXPECT_EQ( tM.indexing_base(), 0 );

    // multiply on the child: rows carry 2, 3, 3, 2 entries of 3.0 each
    tM.multiply( tX, tY );
    EXPECT_NEAR( tY( 0 ), 6.0, tEps );
    EXPECT_NEAR( tY( 1 ), 9.0, tEps );
    EXPECT_NEAR( tY( 2 ), 9.0, tEps );
    EXPECT_NEAR( tY( 3 ), 6.0, tEps );

    EXPECT_EQ( tK.indexing_base(), 0 );
    EXPECT_EQ( tM.indexing_base(), 0 );

    // element access still consistent on both sides, structure still shared
    const belfem::SpMatrix & tConstK = tK;
    const belfem::SpMatrix & tConstM = tM;
    EXPECT_NEAR( tConstK( 1, 1 ), 2.0, tEps );
    EXPECT_NEAR( tConstM( 1, 1 ), 3.0, tEps );
    EXPECT_EQ( tM.pointers(), tK.pointers() );
    EXPECT_EQ( tM.indices(),  tK.indices() );
}

// =============================================================================
// §5  Accessor and base handling  [semantic]
// =============================================================================
//
// position() replaced a member-function-pointer dispatch over four
// ( type x base ) search variants. It reads the base from pointers()[ 0 ] on
// every call, so it must agree with a dense reference for CSR and CSC in both
// bases, and multiply() must no longer rewrite the index arrays at all.

namespace
{
    //! Reference position of ( i, j ) derived from the raw CSR/CSC arrays,
    //! independent of position() itself.
    belfem::int_t
    reference_position( const belfem::SpMatrix & aMatrix,
                        const belfem::index_t    aRow,
                        const belfem::index_t    aCol )
    {
        const belfem::int_t tBase = aMatrix.pointers()[ 0 ];
        const bool tIsCsr = ( aMatrix.type() == belfem::SpMatrixType::CSR );

        const belfem::index_t tSlice  = tIsCsr ? aRow : aCol ;
        const belfem::int_t   tTarget =
            ( belfem::int_t )( tIsCsr ? aCol : aRow ) + tBase ;

        const belfem::int_t * tIndices = tIsCsr ? aMatrix.cols() : aMatrix.rows() ;

        for ( belfem::int_t k = aMatrix.pointers()[ tSlice ]     - tBase ;
                            k < aMatrix.pointers()[ tSlice + 1 ] - tBase ; ++k )
        {
            if ( tIndices[ k ] == tTarget ) return k ;
        }
        return aMatrix.number_of_nonzeros() ;
    }

    //! A pattern with empty rows, a full row, and isolated entries, so the
    //! search hits first/last/absent columns and zero-length slices.
    belfem::Matrix< belfem::real >
    make_awkward_pattern( const belfem::uint aN )
    {
        belfem::Matrix< belfem::real > aDense( aN, aN, 0.0 );

        for ( belfem::uint j = 0; j < aN; ++j )
        {
            aDense( 0, j ) = 1.0 + j ;             // full first row
        }
        for ( belfem::uint i = 2; i < aN; ++i )
        {
            aDense( i, i ) = 2.0 + i ;             // diagonal
            if ( i + 2 < aN ) aDense( i, i + 2 ) = -0.5 ;
        }
        aDense( aN - 1, 0 ) = 7.0 ;                // last row, first column
        // row 1 stays empty on purpose
        return aDense ;
    }
}

TEST( SpMatrix, PositionMatchesReferenceInBothBases )
{
    const belfem::uint tN = 9 ;
    belfem::Matrix< belfem::real > tDense = make_awkward_pattern( tN );

    for ( int tTypeIdx = 0; tTypeIdx < 2; ++tTypeIdx )
    {
        const belfem::SpMatrixType tType = ( tTypeIdx == 0 )
            ? belfem::SpMatrixType::CSR : belfem::SpMatrixType::CSC ;

        belfem::SpMatrix tA( tDense, tType );

        for ( int tBaseIdx = 0; tBaseIdx < 2; ++tBaseIdx )
        {
            tA.set_indexing_base( tBaseIdx == 0
                ? belfem::SpMatrixIndexingBase::Cpp
                : belfem::SpMatrixIndexingBase::Fortran );

            for ( belfem::index_t i = 0; i < tN; ++i )
            {
                for ( belfem::index_t j = 0; j < tN; ++j )
                {
                    const belfem::int_t tExpected = reference_position( tA, i, j );

                    EXPECT_EQ( tA.position( i, j ), tExpected )
                        << "type " << tTypeIdx << " base " << tBaseIdx
                        << " entry ( " << i << ", " << j << " )" ;

                    // index() is a thin wrapper and must agree exactly
                    EXPECT_EQ( tA.index( i, j ), tA.position( i, j ) );

                    // a structural zero must report the miss sentinel
                    if ( tDense( i, j ) == 0.0 )
                    {
                        EXPECT_EQ( tA.position( i, j ),
                                   ( belfem::int_t ) tA.number_of_nonzeros() );
                    }
                }
            }
        }
    }
}

TEST( SpMatrix, PositionHandlesLongAndShortSlices )
{
    // a dense row next to a one-entry row, so the search is exercised at
    // both ends of the slice-length range
    const belfem::uint tN = 96 ;

    belfem::Matrix< belfem::real > tDense( tN, tN, 0.0 );
    for ( belfem::uint j = 0; j < tN; ++j )
    {
        tDense( 0, j ) = 1.0 + j ;
    }
    tDense( 1, 1 ) = 5.0 ;      // short slice, linear branch

    belfem::SpMatrix tA( tDense, belfem::SpMatrixType::CSR );

    for ( belfem::index_t j = 0; j < tN; ++j )
    {
        EXPECT_EQ( tA.position( 0, j ), reference_position( tA, 0, j ) );
    }
    EXPECT_EQ( tA.position( 1, 1 ), reference_position( tA, 1, 1 ) );
    EXPECT_EQ( tA.position( 1, 0 ), ( belfem::int_t ) tA.number_of_nonzeros() );
}

TEST( SpMatrix, PositionsInSliceMatchesPerEntryLookup )
{
    const belfem::uint tN = 9 ;
    belfem::Matrix< belfem::real > tDense = make_awkward_pattern( tN );

    for ( int tTypeIdx = 0; tTypeIdx < 2; ++tTypeIdx )
    {
        const belfem::SpMatrixType tType = ( tTypeIdx == 0 )
            ? belfem::SpMatrixType::CSR : belfem::SpMatrixType::CSC ;

        belfem::SpMatrix tA( tDense, tType );
        const bool tIsCsr = ( tType == belfem::SpMatrixType::CSR );

        for ( int tBaseIdx = 0; tBaseIdx < 2; ++tBaseIdx )
        {
            tA.set_indexing_base( tBaseIdx == 0
                ? belfem::SpMatrixIndexingBase::Cpp
                : belfem::SpMatrixIndexingBase::Fortran );

            for ( belfem::index_t tSlice = 0; tSlice < tN; ++tSlice )
            {
                // every other index, so the list mixes present and absent
                // entries while staying strictly ascending
                belfem::Cell< belfem::int_t > tCols ;
                for ( belfem::int_t k = 0; k < ( belfem::int_t ) tN; k += 2 )
                {
                    tCols.push( k );
                }

                belfem::Cell< belfem::int_t > tPos ;
                tPos.set_size( tCols.size(), -1 );

                tA.positions_in_slice( tSlice, tCols.data(),
                                       ( belfem::uint ) tCols.size(),
                                       tPos.data() );

                for ( belfem::index_t k = 0; k < tCols.size(); ++k )
                {
                    const belfem::index_t tRow = tIsCsr ? tSlice
                                               : ( belfem::index_t ) tCols( k );
                    const belfem::index_t tCol = tIsCsr
                                               ? ( belfem::index_t ) tCols( k )
                                               : tSlice ;

                    EXPECT_EQ( tPos( k ), tA.position( tRow, tCol ) )
                        << "type " << tTypeIdx << " base " << tBaseIdx
                        << " slice " << tSlice << " entry " << k ;
                }
            }
        }
    }
}

TEST( SpMatrix, MultiplyIsBaseNeutral )
{
    belfem::Matrix< belfem::real > tDense = make_awkward_pattern( 9 );

    for ( int tTypeIdx = 0; tTypeIdx < 2; ++tTypeIdx )
    {
        const belfem::SpMatrixType tType = ( tTypeIdx == 0 )
            ? belfem::SpMatrixType::CSR : belfem::SpMatrixType::CSC ;

        belfem::SpMatrix tA( tDense, tType );

        belfem::Vector< belfem::real > tX( 9 );
        for ( belfem::uint k = 0; k < 9; ++k )
        {
            tX( k ) = 1.0 + 0.25 * k ;
        }

        // dense reference
        belfem::Vector< belfem::real > tRef( 9, 0.0 );
        for ( belfem::uint i = 0; i < 9; ++i )
        {
            for ( belfem::uint j = 0; j < 9; ++j )
            {
                tRef( i ) += tDense( i, j ) * tX( j );
            }
        }

        for ( int tBaseIdx = 0; tBaseIdx < 2; ++tBaseIdx )
        {
            const belfem::SpMatrixIndexingBase tBase = ( tBaseIdx == 0 )
                ? belfem::SpMatrixIndexingBase::Cpp
                : belfem::SpMatrixIndexingBase::Fortran ;

            tA.set_indexing_base( tBase );
            const belfem::int_t tBaseBefore = tA.indexing_base();

            belfem::Vector< belfem::real > tY( 9, 0.0 );
            tA.multiply( tX, tY );

            for ( belfem::uint i = 0; i < 9; ++i )
            {
                EXPECT_NEAR( tY( i ), tRef( i ), tTol )
                    << "type " << tTypeIdx << " base " << tBaseIdx
                    << " row " << i ;
            }

            // multiply() must not touch the base any more
            EXPECT_EQ( tA.indexing_base(), tBaseBefore );
        }
    }
}

TEST( SpMatrix, MultiplyAlphaBetaIsBaseNeutral )
{
    belfem::Matrix< belfem::real > tDense = make_awkward_pattern( 9 );
    belfem::SpMatrix tA( tDense, belfem::SpMatrixType::CSR );

    belfem::Vector< belfem::real > tX( 9 );
    for ( belfem::uint k = 0; k < 9; ++k )
    {
        tX( k ) = 0.5 + 0.1 * k ;
    }

    const belfem::real tAlpha = 2.5 ;
    const belfem::real tBeta  = -0.75 ;

    for ( int tBaseIdx = 0; tBaseIdx < 2; ++tBaseIdx )
    {
        tA.set_indexing_base( tBaseIdx == 0
            ? belfem::SpMatrixIndexingBase::Cpp
            : belfem::SpMatrixIndexingBase::Fortran );

        const belfem::int_t tBaseBefore = tA.indexing_base();

        belfem::Vector< belfem::real > tY( 9 );
        belfem::Vector< belfem::real > tRef( 9 );
        for ( belfem::uint k = 0; k < 9; ++k )
        {
            tY( k )   = 3.0 - 0.2 * k ;
            tRef( k ) = tBeta * tY( k );
        }
        for ( belfem::uint i = 0; i < 9; ++i )
        {
            for ( belfem::uint j = 0; j < 9; ++j )
            {
                tRef( i ) += tAlpha * tDense( i, j ) * tX( j );
            }
        }

        tA.multiply( tX, tY, tAlpha, tBeta, false );

        for ( belfem::uint i = 0; i < 9; ++i )
        {
            EXPECT_NEAR( tY( i ), tRef( i ), tTol ) << "base " << tBaseIdx
                                                    << " row " << i ;
        }
        EXPECT_EQ( tA.indexing_base(), tBaseBefore );
    }
}

TEST( SpMatrix, WriteToAbsentEntryHitsTheSentinelSlot )
{
    // allocate_values() reserves one slot past the nonzeros so that a release
    // build assembling an entry outside the pattern lands in a dump slot
    // instead of one element past the allocation. Run under ASan to mean
    // anything - the assertion here only checks that nothing else moved.
    belfem::Matrix< belfem::real > tDense = make_awkward_pattern( 9 );
    belfem::SpMatrix tA( tDense, belfem::SpMatrixType::CSR );

    belfem::Vector< belfem::real > tBefore( tA.number_of_nonzeros() );
    for ( belfem::index_t k = 0; k < ( belfem::index_t ) tA.number_of_nonzeros(); ++k )
    {
        tA.data()[ k ] = 1.0 + k ;
        tBefore( k )   = tA.data()[ k ];
    }

    // ( 1, 3 ) is a structural zero - row 1 is empty
    const belfem::int_t tPos = tA.position( 1, 3 );
    ASSERT_EQ( tPos, ( belfem::int_t ) tA.number_of_nonzeros() );

    tA.data()[ tPos ] += 12345.0 ;

    for ( belfem::index_t k = 0; k < ( belfem::index_t ) tA.number_of_nonzeros(); ++k )
    {
        EXPECT_EQ( tA.data()[ k ], tBefore( k ) ) << "entry " << k << " moved";
    }
}

TEST( SpMatrix, CscMultiplyMatchesReferenceAcrossThreadCounts )
{
    // matvec_csc scatters into y; the per-nonzero atomic was replaced by an
    // array reduction, so the result must be thread-count independent
    const belfem::uint tN = 64 ;

    belfem::Matrix< belfem::real > tDense( tN, tN, 0.0 );
    for ( belfem::uint i = 0; i < tN; ++i )
    {
        tDense( i, i ) = 4.0 + i ;
        if ( i > 0 )      tDense( i, i - 1 ) = -1.0 ;
        if ( i + 1 < tN ) tDense( i, i + 1 ) = -2.0 ;
        tDense( i, ( i * 7 ) % tN ) += 0.5 ;      // scattered off-diagonal
    }

    belfem::SpMatrix tA( tDense, belfem::SpMatrixType::CSC );

    belfem::Vector< belfem::real > tX( tN );
    for ( belfem::uint k = 0; k < tN; ++k )
    {
        tX( k ) = 1.0 - 0.01 * k ;
    }

    belfem::Vector< belfem::real > tRef( tN, 0.0 );
    for ( belfem::uint i = 0; i < tN; ++i )
    {
        for ( belfem::uint j = 0; j < tN; ++j )
        {
            tRef( i ) += tDense( i, j ) * tX( j );
        }
    }

    const int tThreadCounts[] = { 1, 2, 8 };

    for ( int tThreads : tThreadCounts )
    {
#ifdef _OPENMP
        omp_set_num_threads( tThreads );
#else
        ( void ) tThreads ;
#endif
        belfem::Vector< belfem::real > tY( tN, 0.0 );
        tA.multiply( tX, tY );

        for ( belfem::uint i = 0; i < tN; ++i )
        {
            EXPECT_NEAR( tY( i ), tRef( i ), tTol )
                << "threads " << tThreads << " row " << i ;
        }
    }
}
