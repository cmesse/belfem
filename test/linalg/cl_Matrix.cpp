#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"
#include "fn_trans.hpp"

using namespace belfem;

// initialize matrix
Matrix< real > tMatrix = { {1, 2, 3}, {4, 5, 6} };

TEST( LINALG, matrix_initialization )
{
    // test rows
    EXPECT_EQ( tMatrix.n_rows(), (uint) 2 );

    // test columns
    EXPECT_EQ( tMatrix.n_cols(), (uint) 3 );

    // test content
    EXPECT_EQ( tMatrix( 0,0 ), 1 );
    EXPECT_EQ( tMatrix( 1,0 ), 4 );
    EXPECT_EQ( tMatrix( 0,1 ), 2 );
    EXPECT_EQ( tMatrix( 1,1 ), 5 );
    EXPECT_EQ( tMatrix( 0,2 ), 3 );
    EXPECT_EQ( tMatrix( 1,2 ), 6 );

    // test writing function
    tMatrix( 0, 0 ) = -1;
    EXPECT_EQ( tMatrix( 0,0 ), -1 );
}

TEST( LINALG, MATRIX_MINMAX )
{
    EXPECT_EQ( min( tMatrix ), -1 );
    EXPECT_EQ( max( tMatrix ), 6 );
}

TEST( LINALG, matrix_operators )
{
    Matrix< real > tA = { { 8, 1, 6 }, { 3, 5, 7 }, { 4, 9, 2 } };

    // test addition
    Matrix< real > tC( tA + tA );
    EXPECT_EQ( tC( 0,0 ), 16 );
    EXPECT_EQ( tC( 1,0 ),  6 );
    EXPECT_EQ( tC( 2,0 ),  8 );
    EXPECT_EQ( tC( 0,1 ),  2 );
    EXPECT_EQ( tC( 1,1 ), 10 );
    EXPECT_EQ( tC( 2,1 ), 18 );
    EXPECT_EQ( tC( 0,2 ), 12 );
    EXPECT_EQ( tC( 1,2 ), 14 );
    EXPECT_EQ( tC( 2,2 ),  4 );

    // test subtraction
    Matrix< real > tD( tC - tA );
    EXPECT_EQ( tD( 0,0 ),  8 );
    EXPECT_EQ( tD( 1,0 ),  3 );
    EXPECT_EQ( tD( 2,0 ),  4 );
    EXPECT_EQ( tD( 0,1 ),  1 );
    EXPECT_EQ( tD( 1,1 ),  5 );
    EXPECT_EQ( tD( 2,1 ),  9 );
    EXPECT_EQ( tD( 0,2 ),  6 );
    EXPECT_EQ( tD( 1,2 ),  7 );
    EXPECT_EQ( tD( 2,2 ),  2 );

    // test internal addition
    tD += tA;
    EXPECT_EQ( tD( 0,0 ), 16 );
    EXPECT_EQ( tD( 1,0 ),  6 );
    EXPECT_EQ( tD( 2,0 ),  8 );
    EXPECT_EQ( tD( 0,1 ),  2 );
    EXPECT_EQ( tD( 1,1 ), 10 );
    EXPECT_EQ( tD( 2,1 ), 18 );
    EXPECT_EQ( tD( 0,2 ), 12 );
    EXPECT_EQ( tD( 1,2 ), 14 );
    EXPECT_EQ( tD( 2,2 ),  4 );

    // test multiplication
    Matrix< real > tE( tA*tA );
    EXPECT_EQ( tE( 0,0 ), 91 );
    EXPECT_EQ( tE( 1,0 ), 67 );
    EXPECT_EQ( tE( 2,0 ), 67 );
    EXPECT_EQ( tE( 0,1 ), 67 );
    EXPECT_EQ( tE( 1,1 ), 91 );
    EXPECT_EQ( tE( 2,1 ), 67 );
    EXPECT_EQ( tE( 0,2 ), 67 );
    EXPECT_EQ( tE( 1,2 ), 67 );
    EXPECT_EQ( tE( 2,2 ), 91 );

    Vector< real > tX = { 7, 3, 5 };

    Vector< real > tY( tA* tX );
    EXPECT_EQ( tY( 0 ), 89 );
    EXPECT_EQ( tY( 1 ), 71 );
    EXPECT_EQ( tY( 2 ), 65 );
}

TEST( LINALG, matrix_equal_equal )
{
    Matrix<real> tA = { { 4 }, { 3 }, { 2 }, { 1 } };
    Matrix<real> tB = { { 4 }, { 3 }, { 2 }, { 1 } };
    Matrix<real> tC = { { 4 }, { 4 }, { 4 }, { 4 } };

    EXPECT_TRUE(  tA == tB );
    EXPECT_FALSE( tA == tC );

    EXPECT_TRUE(  tC == 4.0 );
    EXPECT_FALSE( tA == 4.0 );
}

TEST( LINALG, vector_from_matrix )
{
    Matrix< real > tA = { { 16,  2,  3, 13 },
                          {  5, 11, 10, 8 },
                          {  9,  7, 6,  12} };


    Vector< real > tCol = tA.col( 2 );

    EXPECT_EQ( tCol.length(), (uint) 3 );

    EXPECT_EQ( tCol( 0 ),  3.0 );
    EXPECT_EQ( tCol( 1 ), 10.0 );
    EXPECT_EQ( tCol( 2 ),  6.0 );

    Vector< real > tRow = tA.row( 1 );

    EXPECT_EQ( tRow.length(), (uint) 4 );

    EXPECT_EQ( tRow( 0 ),  5.0 );
    EXPECT_EQ( tRow( 1 ), 11.0 );
    EXPECT_EQ( tRow( 2 ), 10.0 );
    EXPECT_EQ( tRow( 3 ),  8.0 );

    tA.set_col( 1, tCol );
    EXPECT_EQ( tA( 0, 1 ),  3.0 );
    EXPECT_EQ( tA( 1, 1 ), 10.0 );
    EXPECT_EQ( tA( 2, 1 ),  6.0 );

    tA.set_row( 2, tRow );
    EXPECT_EQ( tA( 2, 0 ),  5.0 );
    EXPECT_EQ( tA( 2, 1 ), 11.0 );
    EXPECT_EQ( tA( 2, 2 ), 10.0 );
    EXPECT_EQ( tA( 2, 3 ),  8.0 );
}

TEST( LINALG, trans )
{
    Matrix< real > tA = { { 8, 1, 6 }, { 3, 5, 7 } };

    Matrix< real > tB = trans( tA );

    EXPECT_EQ( tB( 0, 0 ),  8.0 );
    EXPECT_EQ( tB( 1, 0 ),  1.0 );
    EXPECT_EQ( tB( 2, 0 ),  6.0 );
    EXPECT_EQ( tB( 0, 1 ),  3.0 );
    EXPECT_EQ( tB( 1, 1 ),  5.0 );
    EXPECT_EQ( tB( 2, 1 ),  7.0 );

    Matrix< real > tC = trans( trans( tA ) + tB );

    EXPECT_EQ( tC( 0, 0 ),  16.0 );
    EXPECT_EQ( tC( 1, 0 ),   6.0 );
    EXPECT_EQ( tC( 0, 1 ),   2.0 );
    EXPECT_EQ( tC( 1, 1 ),  10.0 );
    EXPECT_EQ( tC( 0, 2 ),  12.0 );
    EXPECT_EQ( tC( 1, 2 ),  14.0 );

}