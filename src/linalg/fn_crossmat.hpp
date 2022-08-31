//
// Created by Christian Messe on 21.01.22.
//

#ifndef BELFEM_FN_CROSSMAT_HPP
#define BELFEM_FN_CROSSMAT_HPP

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "assert.hpp"
namespace belfem
{
    /**
     * crosses a normal vector with a matrix, 2D version
     */
    inline void
    crossmat(
            const Vector< real > & aN,
            const Matrix< real > & aA,
                  Vector< real > & aNxA )
    {

        uint tN = aA.n_cols() ;

        BELFEM_ASSERT( aN.length() == 2,
                      "Length of normal vector must be 2 (is %u)",
                      ( unsigned int ) aN.length() );

        BELFEM_ASSERT( aNxA.length() == tN,
                      "Length of solution vector does not match (is %u but expect %u )",
                      ( unsigned int ) aNxA.length(),
                      ( unsigned int ) tN );

        for( uint k=0; k<tN; ++k )
        {
            aNxA( k ) = aN( 0 ) * aA( 1, k ) - aN( 1 ) * aA( 0, k ) ;
        }

        // remove dust ( can be removed in the future if too slow, just for beautification )
        real tInvNorm = 1.0/norm( aNxA );
        for( uint k=0; k<tN; ++k )
        {
            if( std::abs(  aNxA( k ) * tInvNorm ) < BELFEM_EPSILON )
            {
                aNxA( k ) = 0.0 ;
            }
        }

    }


    inline void
    crossmat(
            const Vector< real > & aN,
            const Matrix< real > & aA,
            const real             aScale,
            Vector< real > & aNxA )
    {

        uint tN = aA.n_cols() ;

        BELFEM_ASSERT( aN.length() == 2,
                       "Length of normal vector must be 2 (is %u)",
                       ( unsigned int ) aN.length() );

        BELFEM_ASSERT( aNxA.length() == tN,
                       "Length of solution vector does not match (is %u but expect %u )",
                       ( unsigned int ) aNxA.length(),
                       ( unsigned int ) tN );

        for( uint k=0; k<tN; ++k )
        {
            aNxA( k ) += aScale * ( aN( 0 ) * aA( 1, k ) - aN( 1 ) * aA( 0, k ) );
        }

        // remove dust ( can be removed in the future if too slow, just for beautification )
        real tInvNorm = 1.0/norm( aNxA );
        for( uint k=0; k<tN; ++k )
        {
            if( std::abs(  aNxA( k ) * tInvNorm ) < BELFEM_EPSILON )
            {
                aNxA( k ) = 0.0 ;
            }
        }
    }

    // 3D cross
    inline void
    crossmat(
            const Vector< real > & aN,
            const Matrix< real > & aA,
                  Matrix< real > & aNxA )
    {
        uint tN = aA.n_cols() ;

        BELFEM_ASSERT( aN.length() == 3,
                      "Length of normal vector must be 3 (is %u)",
                      ( unsigned int ) aN.length() );

        BELFEM_ASSERT( aNxA.n_cols() == tN,
                      "Number of columns of solution matrix does not match (is %u but expect %u )",
                      ( unsigned int ) aNxA.n_cols(),
                      ( unsigned int ) tN );

        BELFEM_ASSERT( aNxA.n_rows() == 3,
                      "Number of rows of solution matrix does not match (is %u but expect 3)",
                      ( unsigned int ) aNxA.n_rows()  );

        for( uint k=0; k<tN; ++k )
        {
            aNxA( 0, k ) = aN( 1 ) * aA( 2, k ) - aN( 2 ) * aA( 1, k ) ;
            aNxA( 1, k ) = aN( 2 ) * aA( 0, k ) - aN( 0 ) * aA( 2, k ) ;
            aNxA( 2, k ) = aN( 0 ) * aA( 1, k ) - aN( 1 ) * aA( 0, k ) ;

        }
    }

    // 3D cross
    inline void
    crossmat(
            const Vector< real > & aN,
            const Matrix< real > & aA,
            const real             aScale,
            Matrix< real > & aNxA )
    {
        uint tN = aA.n_cols() ;

        BELFEM_ASSERT( aN.length() == 3,
                      "Length of normal vector must be 3 (is %u)",
                      ( unsigned int ) aN.length() );

        BELFEM_ASSERT( aNxA.n_cols() == tN,
                      "Number of columns of solution matrix does not match (is %u but expect %u )",
                      ( unsigned int ) aNxA.n_cols(),
                      ( unsigned int ) tN );

        BELFEM_ASSERT( aNxA.n_rows() == 3,
                      "Number of rows of solution matrix does not match (is %u but expect 3)",
                      ( unsigned int ) aNxA.n_rows()  );

        for( uint k=0; k<tN; ++k )
        {
            aNxA( 0, k ) += aScale * ( aN( 1 ) * aA( 2, k ) - aN( 2 ) * aA( 1, k ) );
            aNxA( 1, k ) += aScale * ( aN( 2 ) * aA( 0, k ) - aN( 0 ) * aA( 2, k ) );
            aNxA( 2, k ) += aScale * ( aN( 0 ) * aA( 1, k ) - aN( 1 ) * aA( 0, k ) );
        }
    }
}
#endif //BELFEM_FN_CROSSMAT_HPP
