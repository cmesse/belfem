/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Matrix<T> wrapper and matrix-centric free functions/operators
 * See: tests_00_strategy.md, tests_02_linalg.md §2, §3.3–3.4, §4.4–4.9
 *
 * Matrix<T> wraps arma::Mat<T> or blaze::DynamicMatrix<T>.
 * Column-major storage: data()[i + j * n_rows()] == operator()(i, j).
 *
 * Includes BUG-L1 regression test for submat() bounds checking.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <string>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_inv2.hpp"
#include "fn_inv3.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"
#include "op_MatrixPlus.hpp"
#include "op_MatrixMinus.hpp"
#include "op_MatrixTimes.hpp"
// operator== for Matrix included transitively via cl_Matrix.hpp

namespace
{
    const belfem::real tEps = 1e-12;  // for exact-in-theory results
}

// =============================================================================
// §2.1  Construction & Destruction  [semantic]
// =============================================================================

TEST( Matrix, DefaultConstructorEmpty )
{
    belfem::Matrix< belfem::real > tMat;
    EXPECT_EQ( tMat.n_rows(), 0u );
    EXPECT_EQ( tMat.n_cols(), 0u );
}

TEST( Matrix, SizedConstructor )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    EXPECT_EQ( tMat.n_rows(), 3u );
    EXPECT_EQ( tMat.n_cols(), 4u );
}

TEST( Matrix, SizedWithFillValue )
{
    belfem::Matrix< belfem::real > tMat( 3, 4, 2.0 );
    EXPECT_EQ( tMat.n_rows(), 3u );
    EXPECT_EQ( tMat.n_cols(), 4u );
    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 4; ++j )
        {
            EXPECT_NEAR( tMat( i, j ), 2.0, tEps );
        }
    }
}

TEST( Matrix, InitializerListConstructor )
{
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0, 3.0 },
                                             { 4.0, 5.0, 6.0 } };
    EXPECT_EQ( tMat.n_rows(), 2u );
    EXPECT_EQ( tMat.n_cols(), 3u );
    EXPECT_NEAR( tMat( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tMat( 1, 0 ), 4.0, tEps );
    EXPECT_NEAR( tMat( 0, 2 ), 3.0, tEps );
    EXPECT_NEAR( tMat( 1, 2 ), 6.0, tEps );
}

TEST( Matrix, CopyConstructorDeepCopies )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    belfem::Matrix< belfem::real > tB( tA );
    tB( 0, 0 ) = 99.0;
    EXPECT_NEAR( tA( 0, 0 ), 1.0, tEps );
}

TEST( Matrix, MoveConstructor )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    belfem::Matrix< belfem::real > tB( std::move( tA ) );
    EXPECT_EQ( tB.n_rows(), 2u );
    EXPECT_EQ( tB.n_cols(), 2u );
    EXPECT_NEAR( tB( 0, 0 ), 1.0, tEps );
    EXPECT_NO_THROW( ( void ) tA.n_rows() );
}

TEST( Matrix, CopyAssignment )
{
    belfem::Matrix< belfem::real > tA = { { 1.0 }, { 2.0 } };
    belfem::Matrix< belfem::real > tB;
    tB = tA;
    tB( 0, 0 ) = 99.0;
    EXPECT_NEAR( tA( 0, 0 ), 1.0, tEps );
}

TEST( Matrix, MoveAssignment )
{
    belfem::Matrix< belfem::real > tA = { { 1.0 }, { 2.0 } };
    belfem::Matrix< belfem::real > tB;
    tB = std::move( tA );
    EXPECT_EQ( tB.n_rows(), 2u );
    EXPECT_NEAR( tB( 0, 0 ), 1.0, tEps );
}

TEST( Matrix, SelfCopyAssignment )
{
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    tMat = tMat;
    EXPECT_NEAR( tMat( 1, 1 ), 4.0, tEps );
}

TEST( Matrix, SelfMoveAssignment )
{
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0 }, { 3.0, 4.0 } };
#pragma GCC diagnostic push
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif
    tMat = std::move( tMat );
#pragma GCC diagnostic pop
    EXPECT_NEAR( tMat( 1, 1 ), 4.0, tEps );
}

TEST( Matrix, AssignFromScalar )
{
    belfem::Matrix< belfem::real > tMat( 2, 3 );
    tMat = 5.0;
    for( size_t i = 0; i < 2; ++i )
    {
        for( size_t j = 0; j < 3; ++j )
        {
            EXPECT_NEAR( tMat( i, j ), 5.0, tEps );
        }
    }
}

// =============================================================================
// §2.2  Memory & Access  [semantic]
// =============================================================================

TEST( Matrix, DataPointerNonNull )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    EXPECT_NE( tMat.data(), nullptr );
}

TEST( Matrix, ColumnMajorLayout )
{
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0, 3.0 },
                                             { 4.0, 5.0, 6.0 } };
    // column-major: column j starts at data[ j * spacing() ] and is
    // contiguous. spacing() is the storage stride between columns ( the
    // BLAS leading dimension ): n_rows() under Armadillo, but >= n_rows()
    // under Blaze, which pads every column to the SIMD width. Indexing
    // with n_rows() is only correct when the two coincide.
    EXPECT_GE( tMat.spacing(), tMat.n_rows() );
    for( size_t j = 0; j < tMat.n_cols(); ++j )
    {
        for( size_t i = 0; i < tMat.n_rows(); ++i )
        {
            EXPECT_NEAR( tMat.data()[ i + j * tMat.spacing() ],
                         tMat( i, j ), tEps );
        }
    }
}

TEST( Matrix, ConstDataPointer )
{
    belfem::Matrix< belfem::real > tMat( 2, 2, 1.0 );
    const belfem::Matrix< belfem::real > & tRef = tMat;
    EXPECT_EQ( tRef.data(), tMat.data() );
}

TEST( Matrix, MatrixDataExposesBackend )
{
    // modify through matrix_data(), verify wrapper sees changes (idea from ChatGPT)
    belfem::Matrix< belfem::real > tMat( 2, 2, 0.0 );

    auto & tRaw = tMat.matrix_data();
    tRaw( 0, 0 ) = 1.0;
    tRaw( 0, 1 ) = 2.0;
    tRaw( 1, 0 ) = 3.0;
    tRaw( 1, 1 ) = 4.0;

    EXPECT_NEAR( tMat( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tMat( 1, 1 ), 4.0, tEps );
}

TEST( Matrix, ParenthesisReadWrite )
{
    belfem::Matrix< belfem::real > tMat( 2, 2 );
    tMat( 0, 0 ) = 1.0;
    tMat( 0, 1 ) = 2.0;
    tMat( 1, 0 ) = 3.0;
    tMat( 1, 1 ) = 4.0;
    EXPECT_NEAR( tMat( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tMat( 1, 1 ), 4.0, tEps );
}

TEST( Matrix, ConstParenthesisAccess )
{
    belfem::Matrix< belfem::real > tMat = { { 3.14 } };
    const belfem::Matrix< belfem::real > & tRef = tMat;
    EXPECT_NEAR( tRef( 0, 0 ), 3.14, tEps );
}

TEST( Matrix, CapacityMatchesProduct )
{
    // capacity() is the physical allocation: exactly rows*cols under
    // Armadillo, but >= rows*cols under Blaze, whose columns are padded
    // for SIMD alignment ( see lapack::leading_dimension() )
    belfem::Matrix< belfem::real > tMat( 3, 5 );
    EXPECT_GE( tMat.capacity(), 15u );
#ifdef BELFEM_ARMADILLO
    EXPECT_EQ( tMat.capacity(), 15u );
#endif
}

// =============================================================================
// §2.3  Sizing  [semantic]
// =============================================================================

TEST( Matrix, FillSetsAllElements )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    tMat.fill( 3.0 );
    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 4; ++j )
        {
            EXPECT_NEAR( tMat( i, j ), 3.0, tEps );
        }
    }
}

TEST( Matrix, SetSizeChanges )
{
    belfem::Matrix< belfem::real > tMat;
    tMat.set_size( 5, 7 );
    EXPECT_EQ( tMat.n_rows(), 5u );
    EXPECT_EQ( tMat.n_cols(), 7u );
}

TEST( Matrix, SetSizeWithValue )
{
    belfem::Matrix< belfem::real > tMat;
    tMat.set_size( 5, 7, 1.0 );
    EXPECT_EQ( tMat.n_rows(), 5u );
    EXPECT_EQ( tMat.n_cols(), 7u );
    EXPECT_NEAR( tMat( 4, 6 ), 1.0, tEps );
}

// =============================================================================
// §2.4  Row / Column / Submatrix Views  [semantic]
// =============================================================================

TEST( Matrix, RowViewReadsCorrectly )
{
    belfem::Matrix< belfem::real > tMat = { { 16.0, 2.0, 3.0, 13.0 },
                                             {  5.0, 11.0, 10.0, 8.0 },
                                             {  9.0, 7.0, 6.0, 12.0 } };
    belfem::Vector< belfem::real > tRow( tMat.row( 1 ) );
    EXPECT_EQ( tRow.length(), 4u );
    EXPECT_NEAR( tRow( 0 ), 5.0,  tEps );
    EXPECT_NEAR( tRow( 1 ), 11.0, tEps );
    EXPECT_NEAR( tRow( 2 ), 10.0, tEps );
    EXPECT_NEAR( tRow( 3 ), 8.0,  tEps );
}

TEST( Matrix, RowViewWriteModifiesParent )
{
    // view is aliased — modifying it changes the parent (idea from ChatGPT)
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0, 3.0 },
                                             { 4.0, 5.0, 6.0 } };
    auto tRow = tMat.row( 0 );
    tRow[ 1 ] = 99.0;   // operator[] works on both arma and blaze views
    EXPECT_NEAR( tMat( 0, 1 ), 99.0, tEps );
}

TEST( Matrix, ColViewReadsCorrectly )
{
    belfem::Matrix< belfem::real > tMat = { { 16.0, 2.0, 3.0, 13.0 },
                                             {  5.0, 11.0, 10.0, 8.0 },
                                             {  9.0, 7.0, 6.0, 12.0 } };
    belfem::Vector< belfem::real > tCol( tMat.col( 2 ) );
    EXPECT_EQ( tCol.length(), 3u );
    EXPECT_NEAR( tCol( 0 ), 3.0,  tEps );
    EXPECT_NEAR( tCol( 1 ), 10.0, tEps );
    EXPECT_NEAR( tCol( 2 ), 6.0,  tEps );
}

TEST( Matrix, ColViewWriteModifiesParent )
{
    // view is aliased — modifying it changes the parent (idea from ChatGPT)
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0, 3.0 },
                                             { 4.0, 5.0, 6.0 } };
    auto tCol = tMat.col( 1 );
    tCol[ 0 ] = 77.0;   // operator[] works on both arma and blaze views
    EXPECT_NEAR( tMat( 0, 1 ), 77.0, tEps );
}

TEST( Matrix, SetRowFromVector )
{
    belfem::Matrix< belfem::real > tMat( 3, 4, 0.0 );
    belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0, 4.0 };
    tMat.set_row( 1, tVec );
    EXPECT_NEAR( tMat( 1, 0 ), 1.0, tEps );
    EXPECT_NEAR( tMat( 1, 1 ), 2.0, tEps );
    EXPECT_NEAR( tMat( 1, 2 ), 3.0, tEps );
    EXPECT_NEAR( tMat( 1, 3 ), 4.0, tEps );
}

TEST( Matrix, SetColFromVector )
{
    belfem::Matrix< belfem::real > tMat( 3, 4, 0.0 );
    belfem::Vector< belfem::real > tVec = { 10.0, 20.0, 30.0 };
    tMat.set_col( 2, tVec );
    EXPECT_NEAR( tMat( 0, 2 ), 10.0, tEps );
    EXPECT_NEAR( tMat( 1, 2 ), 20.0, tEps );
    EXPECT_NEAR( tMat( 2, 2 ), 30.0, tEps );
}

TEST( Matrix, SubmatReadsCorrectly )
{
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0, 3.0, 4.0 },
                                             { 5.0, 6.0, 7.0, 8.0 },
                                             { 9.0, 10.0, 11.0, 12.0 },
                                             { 13.0, 14.0, 15.0, 16.0 } };
    // extract 2x2 submatrix from rows 1-2, cols 1-2
    belfem::Matrix< belfem::real > tSub( tMat.submat( 1, 1, 2, 2 ) );
    EXPECT_NEAR( tSub( 0, 0 ), 6.0,  tEps );
    EXPECT_NEAR( tSub( 0, 1 ), 7.0,  tEps );
    EXPECT_NEAR( tSub( 1, 0 ), 10.0, tEps );
    EXPECT_NEAR( tSub( 1, 1 ), 11.0, tEps );
}

// --- BUG-L1 regression test ---
// submat() bounds check compares row indices against n_cols() instead of
// n_rows(), and reuses aLastRow where aLastCol is intended.
// A non-square matrix (3×5) where row index 2 < 5 (n_cols) but 2 < 3 (n_rows)
// exposes the bug: the buggy assert passes but the correct assert also passes.
// If the bug were fixed to compare against n_rows(), a row index of 3 on a
// 3-row matrix would correctly fail.

TEST( Matrix, SubmatOnNonSquareMatrix )
{
    belfem::Matrix< belfem::real > tMat( 3, 5, 0.0 );
    // fill with identifiable values
    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 5; ++j )
        {
            tMat( i, j ) = static_cast< belfem::real >( i * 10 + j );
        }
    }

    // extract submat(0, 0, 2, 3) — rows 0-2, cols 0-3
    // should succeed: all row indices < 3 (n_rows), all col indices < 5 (n_cols)
    belfem::Matrix< belfem::real > tSub( tMat.submat( 0, 0, 2, 3 ) );
    EXPECT_EQ( tSub.n_rows(), 3u );
    EXPECT_EQ( tSub.n_cols(), 4u );
    EXPECT_NEAR( tSub( 2, 3 ), 23.0, tEps );
}

// =============================================================================
// §2.5  Row / Column / Submatrix  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( MatrixDebug, RowOutOfBoundsThrows )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    EXPECT_THROW( tMat.row( 3 ), std::runtime_error );
}

TEST( MatrixDebug, ColOutOfBoundsThrows )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    EXPECT_THROW( tMat.col( 4 ), std::runtime_error );
}

TEST( MatrixDebug, ParenthesisRowOutOfBoundsThrows )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    EXPECT_THROW( tMat( 3, 0 ), std::runtime_error );
}

TEST( MatrixDebug, ParenthesisColOutOfBoundsThrows )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    EXPECT_THROW( tMat( 0, 4 ), std::runtime_error );
}

TEST( MatrixDebug, SetRowLengthMismatchThrows )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    belfem::Vector< belfem::real > tVec( 3 );   // wrong: need 4
    EXPECT_THROW( tMat.set_row( 0, tVec ), std::runtime_error );
}

TEST( MatrixDebug, SetColLengthMismatchThrows )
{
    belfem::Matrix< belfem::real > tMat( 3, 4 );
    belfem::Vector< belfem::real > tVec( 4 );   // wrong: need 3
    EXPECT_THROW( tMat.set_col( 0, tVec ), std::runtime_error );
}

// --- BUG-L1 regression: submat bounds on non-square matrix ---
// Before fix: all 4 asserts compared against n_cols(), and used aLastRow
// instead of aLastCol. On a 3×5 matrix, row index 3 would pass the buggy
// check (3 < 5) but must fail the correct check (3 < 3).

TEST( MatrixDebug, SubmatRowOutOfBoundsThrows )
{
    belfem::Matrix< belfem::real > tMat( 3, 5, 0.0 );
    // aFirstRow=3 on 3-row matrix → out of bounds
    EXPECT_THROW( tMat.submat( 3, 0, 3, 4 ), std::runtime_error );
}

TEST( MatrixDebug, SubmatColOutOfBoundsThrows )
{
    belfem::Matrix< belfem::real > tMat( 5, 3, 0.0 );
    // aFirstCol=3 on 3-col matrix → out of bounds
    // Before fix: the 3rd assert checked aFirstRow (0) < n_cols() (3) — passed.
    // After fix: checks aFirstCol (3) < n_cols() (3) — correctly fails.
    EXPECT_THROW( tMat.submat( 0, 3, 4, 3 ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §2.6  Compound Operators  [semantic]
// =============================================================================

TEST( Matrix, PlusEqualsScalar )
{
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    tMat += 2.0;
    EXPECT_NEAR( tMat( 0, 0 ), 3.0, tEps );
    EXPECT_NEAR( tMat( 1, 1 ), 6.0, tEps );
}

TEST( Matrix, PlusEqualsMatrix )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    belfem::Matrix< belfem::real > tB = { { 10.0, 20.0 }, { 30.0, 40.0 } };
    tA += tB;
    EXPECT_NEAR( tA( 0, 0 ), 11.0, tEps );
    EXPECT_NEAR( tA( 1, 1 ), 44.0, tEps );
}

TEST( Matrix, MinusEqualsScalar )
{
    belfem::Matrix< belfem::real > tMat = { { 5.0, 7.0 }, { 9.0, 11.0 } };
    tMat -= 2.0;
    EXPECT_NEAR( tMat( 0, 0 ), 3.0, tEps );
    EXPECT_NEAR( tMat( 1, 1 ), 9.0, tEps );
}

TEST( Matrix, MinusEqualsMatrix )
{
    belfem::Matrix< belfem::real > tA = { { 10.0, 20.0 }, { 30.0, 40.0 } };
    belfem::Matrix< belfem::real > tB = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    tA -= tB;
    EXPECT_NEAR( tA( 0, 0 ), 9.0,  tEps );
    EXPECT_NEAR( tA( 1, 1 ), 36.0, tEps );
}

TEST( Matrix, TimesEqualsScalar )
{
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    tMat *= 3.0;
    EXPECT_NEAR( tMat( 0, 0 ), 3.0,  tEps );
    EXPECT_NEAR( tMat( 1, 1 ), 12.0, tEps );
}

TEST( Matrix, TimesEqualsMatrix )
{
    // in-place matrix multiplication (idea from ChatGPT)
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    belfem::Matrix< belfem::real > tB = { { 2.0, 0.0 }, { 1.0, 2.0 } };
    tA *= tB;
    EXPECT_NEAR( tA( 0, 0 ), 4.0,  tEps );
    EXPECT_NEAR( tA( 0, 1 ), 4.0,  tEps );
    EXPECT_NEAR( tA( 1, 0 ), 10.0, tEps );
    EXPECT_NEAR( tA( 1, 1 ), 8.0,  tEps );
}

TEST( Matrix, DivideEqualsScalar )
{
    belfem::Matrix< belfem::real > tMat = { { 2.0, 4.0 }, { 6.0, 8.0 } };
    tMat /= 2.0;
    EXPECT_NEAR( tMat( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tMat( 1, 1 ), 4.0, tEps );
}

// =============================================================================
// §2.7  Print  [semantic] (smoke test only)
// =============================================================================

TEST( Matrix, PrintProducesOutput )
{
    // capture stdout, verify non-empty (idea from ChatGPT)
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    testing::internal::CaptureStdout();
    tMat.print( "tMat" );
    std::string tOut = testing::internal::GetCapturedStdout();
    EXPECT_FALSE( tOut.empty() );
}

// =============================================================================
// §3.3  Matrix Binary Operators  [semantic]
// =============================================================================

TEST( Matrix, MatrixPlusMatrix )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    belfem::Matrix< belfem::real > tB = { { 5.0, 6.0 }, { 7.0, 8.0 } };
    belfem::Matrix< belfem::real > tC( tA + tB );
    EXPECT_NEAR( tC( 0, 0 ), 6.0,  tEps );
    EXPECT_NEAR( tC( 0, 1 ), 8.0,  tEps );
    EXPECT_NEAR( tC( 1, 0 ), 10.0, tEps );
    EXPECT_NEAR( tC( 1, 1 ), 12.0, tEps );
}

TEST( Matrix, MatrixMinusMatrix )
{
    belfem::Matrix< belfem::real > tA = { { 6.0, 8.0 }, { 10.0, 12.0 } };
    belfem::Matrix< belfem::real > tB = { { 5.0, 6.0 }, { 7.0, 8.0 } };
    belfem::Matrix< belfem::real > tC( tA - tB );
    EXPECT_NEAR( tC( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tC( 1, 1 ), 4.0, tEps );
}

TEST( Matrix, MatrixTimesMatrix )
{
    // magic square property: A*A has known entries
    belfem::Matrix< belfem::real > tA = { { 8.0, 1.0, 6.0 },
                                           { 3.0, 5.0, 7.0 },
                                           { 4.0, 9.0, 2.0 } };
    belfem::Matrix< belfem::real > tC( tA * tA );
    EXPECT_NEAR( tC( 0, 0 ), 91.0, tEps );
    EXPECT_NEAR( tC( 1, 1 ), 91.0, tEps );
    EXPECT_NEAR( tC( 2, 2 ), 91.0, tEps );
    EXPECT_NEAR( tC( 0, 1 ), 67.0, tEps );
}

TEST( Matrix, MatrixTimesVector )
{
    belfem::Matrix< belfem::real > tA = { { 8.0, 1.0, 6.0 },
                                           { 3.0, 5.0, 7.0 },
                                           { 4.0, 9.0, 2.0 } };
    belfem::Vector< belfem::real > tX = { 7.0, 3.0, 5.0 };
    belfem::Vector< belfem::real > tY( tA * tX );
    EXPECT_NEAR( tY( 0 ), 89.0, tEps );
    EXPECT_NEAR( tY( 1 ), 71.0, tEps );
    EXPECT_NEAR( tY( 2 ), 65.0, tEps );
}

TEST( Matrix, MatrixTimesScalar )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    belfem::Matrix< belfem::real > tC( tA * 2.0 );
    EXPECT_NEAR( tC( 0, 0 ), 2.0, tEps );
    EXPECT_NEAR( tC( 1, 1 ), 8.0, tEps );
}

TEST( Matrix, MatrixEqualityTrue )
{
    belfem::Matrix< belfem::real > tA = { { 1.0 }, { 2.0 } };
    belfem::Matrix< belfem::real > tB = { { 1.0 }, { 2.0 } };
    EXPECT_TRUE( tA == tB );
}

TEST( Matrix, MatrixEqualityFalse )
{
    belfem::Matrix< belfem::real > tA = { { 1.0 }, { 2.0 } };
    belfem::Matrix< belfem::real > tB = { { 1.0 }, { 3.0 } };
    EXPECT_FALSE( tA == tB );
}

TEST( Matrix, MatrixEqualityScalar )
{
    belfem::Matrix< belfem::real > tA = { { 4.0 }, { 4.0 }, { 4.0 } };
    EXPECT_TRUE( tA == 4.0 );
    belfem::Matrix< belfem::real > tB = { { 4.0 }, { 3.0 } };
    EXPECT_FALSE( tB == 4.0 );
}

// =============================================================================
// §3.4  Mixed Expression Smoke Test  [semantic]
// =============================================================================

TEST( Matrix, ChainedMatrixVectorProduct )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    belfem::Matrix< belfem::real > tB = { { 5.0, 6.0 }, { 7.0, 8.0 } };
    belfem::Vector< belfem::real > tX = { 1.0, 1.0 };

    // A * B * x should match A * (B * x)
    belfem::Vector< belfem::real > tBx( tB * tX );
    belfem::Vector< belfem::real > tABx( tA * tBx );

    belfem::Matrix< belfem::real > tAB( tA * tB );
    belfem::Vector< belfem::real > tABx2( tAB * tX );

    for( size_t i = 0; i < 2; ++i )
    {
        EXPECT_NEAR( tABx( i ), tABx2( i ), tEps );
    }
}

// =============================================================================
// §4.5  Determinant  [semantic]
// =============================================================================

TEST( Matrix, DetIdentity2x2 )
{
    belfem::Matrix< belfem::real > tI = { { 1.0, 0.0 }, { 0.0, 1.0 } };
    EXPECT_NEAR( belfem::det( tI ), 1.0, tEps );
}

TEST( Matrix, DetIdentity3x3 )
{
    belfem::Matrix< belfem::real > tI = { { 1.0, 0.0, 0.0 },
                                           { 0.0, 1.0, 0.0 },
                                           { 0.0, 0.0, 1.0 } };
    EXPECT_NEAR( belfem::det( tI ), 1.0, tEps );
}

TEST( Matrix, DetKnownMatrix )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    EXPECT_NEAR( belfem::det( tA ), -2.0, tEps );
}

TEST( Matrix, DetSingular )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 2.0, 4.0 } };
    EXPECT_NEAR( belfem::det( tA ), 0.0, tEps );
}

// =============================================================================
// §4.6  Inverse  [semantic]
// =============================================================================

TEST( Matrix, InvIdentity )
{
    belfem::Matrix< belfem::real > tI = { { 1.0, 0.0 }, { 0.0, 1.0 } };
    belfem::Matrix< belfem::real > tInv( belfem::inv( tI ) );
    EXPECT_NEAR( tInv( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tInv( 0, 1 ), 0.0, tEps );
    EXPECT_NEAR( tInv( 1, 0 ), 0.0, tEps );
    EXPECT_NEAR( tInv( 1, 1 ), 1.0, tEps );
}

TEST( Matrix, InvRoundTrip )
{
    belfem::Matrix< belfem::real > tA = { { 2.0, 1.0 }, { 1.0, 3.0 } };
    belfem::Matrix< belfem::real > tInv( belfem::inv( tA ) );
    belfem::Matrix< belfem::real > tProduct( tA * tInv );
    EXPECT_NEAR( tProduct( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tProduct( 0, 1 ), 0.0, tEps );
    EXPECT_NEAR( tProduct( 1, 0 ), 0.0, tEps );
    EXPECT_NEAR( tProduct( 1, 1 ), 1.0, tEps );
}

// =============================================================================
// §4.7  Inv2 and Inv3 (BELFEM-authored)  [semantic]
// =============================================================================

TEST( Matrix, Inv2ReturnsCorrectInverse )
{
    belfem::Matrix< belfem::real > tA = { { 4.0, 7.0 }, { 2.0, 6.0 } };
    belfem::Matrix< belfem::real > tB( 2, 2 );

    belfem::real tDetA = belfem::inv2( tA, tB );

    // verify A * B ≈ I
    belfem::Matrix< belfem::real > tProduct( tA * tB );
    EXPECT_NEAR( tProduct( 0, 0 ), 1.0, tEps );
    EXPECT_NEAR( tProduct( 0, 1 ), 0.0, tEps );
    EXPECT_NEAR( tProduct( 1, 0 ), 0.0, tEps );
    EXPECT_NEAR( tProduct( 1, 1 ), 1.0, tEps );

    // verify returned determinant
    EXPECT_NEAR( tDetA, belfem::det( tA ), tEps );
}

TEST( Matrix, Inv3ReturnsCorrectInverse )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0, 3.0 },
                                           { 0.0, 1.0, 4.0 },
                                           { 5.0, 6.0, 0.0 } };
    belfem::Matrix< belfem::real > tB( 3, 3 );

    belfem::real tDetA = belfem::inv3( tA, tB );

    // verify A * B ≈ I
    belfem::Matrix< belfem::real > tProduct( tA * tB );
    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 3; ++j )
        {
            belfem::real tExpected = ( i == j ) ? 1.0 : 0.0;
            EXPECT_NEAR( tProduct( i, j ), tExpected, tEps );
        }
    }

    // verify returned determinant
    EXPECT_NEAR( tDetA, belfem::det( tA ), tEps );
}

// =============================================================================
// §4.7  Inv2/Inv3  [debug]
// =============================================================================

#ifndef NDEBUG

TEST( MatrixDebug, Inv2SingularThrows )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0 }, { 2.0, 4.0 } };
    belfem::Matrix< belfem::real > tB( 2, 2 );
    EXPECT_THROW( belfem::inv2( tA, tB ), std::runtime_error );
}

TEST( MatrixDebug, Inv3SingularThrows )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0, 3.0 },
                                           { 4.0, 5.0, 6.0 },
                                           { 7.0, 8.0, 9.0 } };
    belfem::Matrix< belfem::real > tB( 3, 3 );
    EXPECT_THROW( belfem::inv3( tA, tB ), std::runtime_error );
}

#endif // NDEBUG

// =============================================================================
// §4.9  Transpose  [semantic]
// =============================================================================

TEST( Matrix, TransposeSwapsDimensions )
{
    belfem::Matrix< belfem::real > tA( 3, 4, 1.0 );
    belfem::Matrix< belfem::real > tB( belfem::trans( tA ) );
    EXPECT_EQ( tB.n_rows(), 4u );
    EXPECT_EQ( tB.n_cols(), 3u );
}

TEST( Matrix, TransposeValues )
{
    belfem::Matrix< belfem::real > tA = { { 8.0, 1.0, 6.0 },
                                           { 3.0, 5.0, 7.0 } };
    belfem::Matrix< belfem::real > tB( belfem::trans( tA ) );
    EXPECT_NEAR( tB( 0, 0 ), 8.0, tEps );
    EXPECT_NEAR( tB( 1, 0 ), 1.0, tEps );
    EXPECT_NEAR( tB( 2, 0 ), 6.0, tEps );
    EXPECT_NEAR( tB( 0, 1 ), 3.0, tEps );
    EXPECT_NEAR( tB( 1, 1 ), 5.0, tEps );
    EXPECT_NEAR( tB( 2, 1 ), 7.0, tEps );
}

TEST( Matrix, DoubleTransposeIdentity )
{
    belfem::Matrix< belfem::real > tA = { { 1.0, 2.0, 3.0 },
                                           { 4.0, 5.0, 6.0 } };
    belfem::Matrix< belfem::real > tB( belfem::trans( belfem::trans( tA ) ) );
    for( size_t i = 0; i < 2; ++i )
    {
        for( size_t j = 0; j < 3; ++j )
        {
            EXPECT_NEAR( tB( i, j ), tA( i, j ), tEps );
        }
    }
}

// =============================================================================
// §4.11  Min, Max for Matrix  [semantic]
// =============================================================================

TEST( Matrix, MinMatrix )
{
    belfem::Matrix< belfem::real > tMat = { { 3.0, -1.0 }, { 2.0, 5.0 } };
    EXPECT_NEAR( belfem::min( tMat ), -1.0, tEps );
}

TEST( Matrix, MaxMatrix )
{
    belfem::Matrix< belfem::real > tMat = { { 3.0, -1.0 }, { 2.0, 5.0 } };
    EXPECT_NEAR( belfem::max( tMat ), 5.0, tEps );
}

TEST( Matrix, ScalarTimesMatrix )
{
    belfem::Matrix< belfem::real > tMat = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    belfem::Matrix< belfem::real > tResult( 2.0 * tMat );

    EXPECT_NEAR( tResult( 0, 0 ), 2.0, tEps );
    EXPECT_NEAR( tResult( 0, 1 ), 4.0, tEps );
    EXPECT_NEAR( tResult( 1, 0 ), 6.0, tEps );
    EXPECT_NEAR( tResult( 1, 1 ), 8.0, tEps );
}
