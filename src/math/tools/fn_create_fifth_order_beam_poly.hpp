//
// Created by Christian Messe on 28.08.19.
//

#ifndef BELFEM_FN_CREATE_FIFTH_ORDER_BEAM_POLY_HPP
#define BELFEM_FN_CREATE_FIFTH_ORDER_BEAM_POLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_gesv.hpp"

namespace belfem
{
    void
    create_fifth_order_beam_poly( const real aX1, const real aF1, const real adF1dX, const real ad2F1dX,
                                  const real aX2, const real aF2, const real adF2dX, const real ad2F2dX,
                                  Vector <real> & aCoefficients )
    {
        // create the vandermonde matrix
        Matrix <real> tVandermonde( 6, 6 );

        // fill the matrix
        tVandermonde( 0, 0 ) = aX1 * aX1 * aX1 * aX1 * aX1;
        tVandermonde( 1, 0 ) = 5 * aX1 * aX1 * aX1 * aX1;
        tVandermonde( 2, 0 ) = 15 * aX1 * aX1 * aX1;
        tVandermonde( 3, 0 ) = aX2 * aX2 * aX2 * aX2 * aX2;
        tVandermonde( 4, 0 ) = 5 * aX2 * aX2 * aX2 * aX2;
        tVandermonde( 5, 0 ) = 15 * aX2 * aX2 * aX2;

        tVandermonde( 0, 1 ) = aX1 * aX1 * aX1 * aX1;
        tVandermonde( 1, 1 ) = 4 * aX1 * aX1 * aX1;
        tVandermonde( 2, 1 ) = 12 * aX1 * aX1;
        tVandermonde( 3, 1 ) = aX2 * aX2 * aX2 * aX2;
        tVandermonde( 4, 1 ) = 4 * aX2 * aX2 * aX2;
        tVandermonde( 5, 1 ) = 12 * aX2 * aX2;

        tVandermonde( 0, 2 ) = aX1 * aX1 * aX1;
        tVandermonde( 1, 2 ) = 3 * aX1 * aX1;
        tVandermonde( 2, 2 ) = 6 * aX1;
        tVandermonde( 3, 2 ) = aX2 * aX2 * aX2;
        tVandermonde( 4, 2 ) = 3 * aX2 * aX2;
        tVandermonde( 5, 2 ) = 6 * aX2;

        tVandermonde( 0, 3 ) = aX1 * aX1;
        tVandermonde( 1, 3 ) = 2 * aX1;
        tVandermonde( 2, 3 ) = 2.0;
        tVandermonde( 3, 3 ) = aX2 * aX2;
        tVandermonde( 4, 3 ) = 2 * aX2;
        tVandermonde( 5, 3 ) = 2.0;

        tVandermonde( 0, 4 ) = aX1;
        tVandermonde( 1, 4 ) = 1.0;
        tVandermonde( 2, 4 ) = 0;
        tVandermonde( 3, 4 ) = aX2;
        tVandermonde( 4, 4 ) = 1.0;
        tVandermonde( 5, 4 ) = 0.0;

        tVandermonde( 0, 5 ) = 1.0;
        tVandermonde( 1, 5 ) = 0.0;
        tVandermonde( 2, 5 ) = 0.0;
        tVandermonde( 3, 5 ) = 1.0;
        tVandermonde( 4, 5 ) = 0.0;
        tVandermonde( 5, 5 ) = 0.0;

        // create the right hand side
        aCoefficients.set_size( 6 );

        aCoefficients( 0 ) = aF1;
        aCoefficients( 1 ) = adF1dX;
        aCoefficients( 2 ) = ad2F1dX;
        aCoefficients( 3 ) = aF2;
        aCoefficients( 4 ) = adF2dX;
        aCoefficients( 5 ) = ad2F2dX;

        // allocate the pivot vector
        Vector<int> tPivot( 6 );

        // solve the system and return the coefficients
        gesv( tVandermonde, aCoefficients, tPivot );
    }
}
#endif //BELFEM_FN_CREATE_FIFTH_ORDER_BEAM_POLY_HPP
