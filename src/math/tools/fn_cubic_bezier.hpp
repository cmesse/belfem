//
// Created by Christian Messe on 03.10.19.
//

#ifndef BELFEM_FN_CUBIC_BEZIER_HPP
#define BELFEM_FN_CUBIC_BEZIER_HPP

#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"


namespace belfem
{
//------------------------------------------------------------------------------

    template < typename T >
    void
    cubic_bezier(
            const Matrix< T > & aPoints,
            Vector< T >       & aWork,
            const T           & aXi,
            Vector< T >       & aPoint )
    {

        BELFEM_ASSERT( aPoint.length() == aPoints.n_rows(), "Length of Point vector does not match" );
        BELFEM_ASSERT( aPoints.n_cols() == 4, "Matrix must have 4 columns");
        BELFEM_ASSERT( aWork.length() == 4, "Work Vector must have 4 columns");

        aWork( 0 ) = std::pow( 1.0 - aXi, 3 );
        aWork( 1 ) = 3.0 * ( ( aXi *(aXi -1.0)-1.0) * aXi  + 1.0 );
        aWork( 2 ) =  3.0 * ( ( 1.0 - (aXi +1.0)*aXi  ) * aXi  + 1.0 );
        aWork( 3 ) =( ( aXi  * ( 3.0 + aXi  ) + 3.0 ) * aXi  + 1.0 );
        aWork *= 0.125;
        aPoint = aPoints * aWork;

        // cleanup point
        for( T & tX : aPoint )
        {
            if( std::abs( tX ) < BELFEM_EPSILON )
            {
                tX = 0.0 ;
            }
        }
    }

//------------------------------------------------------------------------------

    template < typename T >
    void
    cubic_bezier_derivative(
                  const Matrix< T > & aPoints,
                  Vector< T >   & aWork,
                  const T       & aXi,
                  Vector< T >   & aPoint )
    {

        BELFEM_ASSERT( aPoint.length() == aPoints.n_rows(), "Length of Point vector does not match" );
        BELFEM_ASSERT( aPoints.n_cols() == 4, "Matrix must have 4 columns");
        BELFEM_ASSERT( aWork.length() == 4, "Work Vector must have 4 columns");


        aWork( 0 ) = - ( 1.0 - aXi) * ( 1.0 - aXi );
        aWork( 1 ) =  ( ( 3.0 * aXi - 2.0 ) * aXi - 1.0 );
        aWork( 2 ) =   1.0 - aXi * ( 2.0 + 3.0 * aXi );
        aWork( 3 ) =   ( 1.0 + aXi) * ( 1.0 + aXi );
        aWork *= 0.375;
        aPoint = aPoints * aWork;

        // cleanup point
        for( T & tX : aPoint )
        {
            if( std::abs( tX ) < BELFEM_EPSILON )
            {
                tX = 0.0 ;
            }
        }
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_BEZIER_HPP
