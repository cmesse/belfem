/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * CSR storage and projections of the transformation matrix.
 *
 * Tmatrix keeps only the entries of a dense T with |T_ij| > BELFEM_EPSILON.
 * The dense copy it hands back, the vector projection T^T x and the matrix
 * projection T^T A T are checked against the same explicit dense oracle,
 * which carries 0.0 where the class must have dropped an entry.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_FEM_Tmatrix.hpp"

namespace
{
    // the dense input: row 0 holds two kept entries, row 2 holds none, and
    // the three epsilon-scaled entries sit below, at and above the threshold
    void
    fill_input( belfem::Matrix< belfem::real > & aT )
    {
        aT.set_size( 4, 3, 0.0 );
        aT( 0, 0 ) = -0.3 ;
        aT( 0, 1 ) =  0.25 ;
        aT( 0, 2 ) =  0.0 ;
        aT( 1, 0 ) =  0.0 ;
        aT( 1, 1 ) =  0.5 * belfem::BELFEM_EPSILON ;
        aT( 1, 2 ) =  0.7 ;
        aT( 2, 0 ) =  belfem::BELFEM_EPSILON ;
        aT( 2, 1 ) =  0.0 ;
        aT( 2, 2 ) =  0.0 ;
        aT( 3, 0 ) =  0.0 ;
        aT( 3, 1 ) =  2.0 * belfem::BELFEM_EPSILON ;
        aT( 3, 2 ) =  0.0 ;
    }

    // what the class must store: the input with the entries at and below
    // BELFEM_EPSILON replaced by 0.0
    void
    fill_expected( belfem::Matrix< belfem::real > & aT )
    {
        aT.set_size( 4, 3, 0.0 );
        aT( 0, 0 ) = -0.3 ;
        aT( 0, 1 ) =  0.25 ;
        aT( 1, 2 ) =  0.7 ;
        aT( 3, 1 ) =  2.0 * belfem::BELFEM_EPSILON ;
    }

    // a symmetric 4 x 4 operand for the matrix projection
    void
    fill_symmetric( belfem::Matrix< belfem::real > & aA )
    {
        aA.set_size( 4, 4, 0.0 );
        aA( 0, 0 ) =  4.0 ;
        aA( 1, 1 ) =  3.0 ;
        aA( 2, 2 ) =  2.5 ;
        aA( 3, 3 ) =  1.5 ;
        aA( 0, 1 ) = aA( 1, 0 ) = -1.0 ;
        aA( 0, 2 ) = aA( 2, 0 ) =  0.5 ;
        aA( 0, 3 ) = aA( 3, 0 ) =  2.0 ;
        aA( 1, 2 ) = aA( 2, 1 ) =  0.75 ;
        aA( 1, 3 ) = aA( 3, 1 ) = -0.25 ;
        aA( 2, 3 ) = aA( 3, 2 ) =  1.25 ;
    }
}

TEST( Tmatrix, GetMatrixDropsEntriesAtAndBelowEpsilon )
{
    belfem::Matrix< belfem::real > tInput ;
    fill_input( tInput );

    belfem::Matrix< belfem::real > tExpected ;
    fill_expected( tExpected );

    belfem::fem::Tmatrix tT( tInput );

    belfem::Matrix< belfem::real > tStored ;
    tT.get_matrix( tStored );

    ASSERT_EQ( tStored.n_rows(), 4u );
    ASSERT_EQ( tStored.n_cols(), 3u );

    // the kept entries are copied, not rounded, so the comparison is exact
    for ( belfem::uint j = 0; j < 3; ++j )
    {
        for ( belfem::uint i = 0; i < 4; ++i )
        {
            EXPECT_EQ( tStored( i, j ), tExpected( i, j ) ) << "i = " << i << " j = " << j ;
        }
    }
}

TEST( Tmatrix, VectorProjectIsTransposeTimesVector )
{
    belfem::Matrix< belfem::real > tInput ;
    fill_input( tInput );

    belfem::Matrix< belfem::real > tExpected ;
    fill_expected( tExpected );

    belfem::fem::Tmatrix tT( tInput );

    belfem::Vector< belfem::real > tX( 4 );
    tX( 0 ) =  1.0 ;
    tX( 1 ) = -2.0 ;
    tX( 2 ) =  0.5 ;
    tX( 3 ) =  3.0 ;

    belfem::Vector< belfem::real > tY ;
    tT.project( tX, tY );

    ASSERT_EQ( tY.length(), 3u );

    // y_j = sum_i T_ij x_i over the dense oracle
    for ( belfem::uint j = 0; j < 3; ++j )
    {
        belfem::real tYj = 0.0 ;
        for ( belfem::uint i = 0; i < 4; ++i )
        {
            tYj += tExpected( i, j ) * tX( i );
        }
        EXPECT_NEAR( tY( j ), tYj, 1e-14 ) << "j = " << j ;
    }
}

TEST( Tmatrix, MatrixProjectIsTransposeTimesMatrixTimesT )
{
    belfem::Matrix< belfem::real > tInput ;
    fill_input( tInput );

    belfem::Matrix< belfem::real > tExpected ;
    fill_expected( tExpected );

    belfem::Matrix< belfem::real > tA ;
    fill_symmetric( tA );

    belfem::fem::Tmatrix tT( tInput );

    belfem::Matrix< belfem::real > tB ;
    tT.project( tA, tB );

    ASSERT_EQ( tB.n_rows(), 3u );
    ASSERT_EQ( tB.n_cols(), 3u );

    // B_jl = sum_i sum_k T_ij A_ik T_kl over the dense oracle
    for ( belfem::uint l = 0; l < 3; ++l )
    {
        for ( belfem::uint j = 0; j < 3; ++j )
        {
            belfem::real tBjl = 0.0 ;
            for ( belfem::uint k = 0; k < 4; ++k )
            {
                for ( belfem::uint i = 0; i < 4; ++i )
                {
                    tBjl += tExpected( i, j ) * tA( i, k ) * tExpected( k, l );
                }
            }
            EXPECT_NEAR( tB( j, l ), tBjl, 1e-13 ) << "j = " << j << " l = " << l ;
        }
    }

    // a symmetric operand projects to a symmetric result
    for ( belfem::uint l = 0; l < 3; ++l )
    {
        for ( belfem::uint j = 0; j < l; ++j )
        {
            EXPECT_NEAR( tB( j, l ), tB( l, j ), 1e-13 ) << "j = " << j << " l = " << l ;
        }
    }
}

// a nonsymmetric operand separates T^T A T from T^T A^T T
TEST( Tmatrix, MatrixProjectKeepsTheOperandOrientation )
{
    belfem::Matrix< belfem::real > tInput ;
    fill_input( tInput );

    belfem::Matrix< belfem::real > tExpected ;
    fill_expected( tExpected );

    belfem::Matrix< belfem::real > tA ;
    fill_symmetric( tA );
    tA( 0, 1 ) += 0.6 ;
    tA( 2, 3 ) -= 0.4 ;
    tA( 3, 0 ) += 1.1 ;
    ASSERT_NE( tA( 0, 1 ), tA( 1, 0 ) );

    belfem::fem::Tmatrix tT( tInput );

    belfem::Matrix< belfem::real > tB ;
    tT.project( tA, tB );

    ASSERT_EQ( tB.n_rows(), 3u );
    ASSERT_EQ( tB.n_cols(), 3u );

    for ( belfem::uint l = 0; l < 3; ++l )
    {
        for ( belfem::uint j = 0; j < 3; ++j )
        {
            belfem::real tBjl = 0.0 ;
            belfem::real tBjlTransposed = 0.0 ;
            for ( belfem::uint k = 0; k < 4; ++k )
            {
                for ( belfem::uint i = 0; i < 4; ++i )
                {
                    tBjl           += tExpected( i, j ) * tA( i, k ) * tExpected( k, l );
                    tBjlTransposed += tExpected( i, j ) * tA( k, i ) * tExpected( k, l );
                }
            }
            EXPECT_NEAR( tB( j, l ), tBjl, 1e-13 ) << "j = " << j << " l = " << l ;
            if ( std::abs( tBjl - tBjlTransposed ) > 1e-6 )
            {
                EXPECT_GT( std::abs( tB( j, l ) - tBjlTransposed ), 1e-6 )
                    << "j = " << j << " l = " << l ;
            }
        }
    }
}
