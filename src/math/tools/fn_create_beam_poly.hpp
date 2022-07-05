//
// Created by Christian Messe on 27.08.19.
//

#ifndef BELFEM_FN_CREATE_BEAM_POLY_HPP
#define BELFEM_FN_CREATE_BEAM_POLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_gesv.hpp"

namespace belfem
{
    inline void
    create_beam_poly( const real aX1, const real aF1, const real adF1dX,
                      const real aX2, const real aF2, const real adF2dX,
                      Vector< real > & aCoefficients )
    {
        // create the vandermonde matrix
        Matrix< real > tVandermonde( 4, 4 );

        // fill the matrix
        tVandermonde( 0, 0 ) = aX1*aX1*aX1;
        tVandermonde( 1, 0 ) = 3.0 * aX1 * aX1;
        tVandermonde( 2, 0 ) = aX2*aX2*aX2;
        tVandermonde( 3, 0 ) = 3.0 * aX2 * aX2;
        tVandermonde( 0, 1 ) = aX1*aX1;
        tVandermonde( 1, 1 ) = 2.0 * aX1;
        tVandermonde( 2, 1 ) = aX2*aX2;
        tVandermonde( 3, 1 ) = 2.0 * aX2 ;
        tVandermonde( 0, 2 ) = aX1;
        tVandermonde( 1, 2 ) = 1.0;
        tVandermonde( 2, 2 ) = aX2;
        tVandermonde( 3, 2 ) = 1.0;
        tVandermonde( 0, 3 ) = 1.0;
        tVandermonde( 1, 3 ) = 0.0;
        tVandermonde( 2, 3 ) = 1.0;
        tVandermonde( 3, 3 ) = 0.0;

        // create the right hand side
        aCoefficients.set_size( 4 );

        aCoefficients( 0 ) = aF1;
        aCoefficients( 1 ) = adF1dX;
        aCoefficients( 2 ) = aF2;
        aCoefficients( 3 ) = adF2dX;

        // allocate the pivot vector
        Vector< int > tPivot( 4 );

        // solve the system and return the coefficients
        gesv( tVandermonde, aCoefficients, tPivot );
    }
}

#endif //BELFEM_FN_CREATE_BEAM_POLY_HPP
