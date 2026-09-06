/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_CREATE_FIFTH_ORDER_BEAM_POLY_HPP
#define BELFEM_FN_CREATE_FIFTH_ORDER_BEAM_POLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "../../linalg/lapack/fn_gesv.hpp"

namespace belfem
{
    inline void
    create_fifth_order_beam_poly( const real aX1, const real aF1, const real adF1dX, const real ad2F1dX,
                                  const real aX2, const real aF2, const real adF2dX, const real ad2F2dX,
                                  Vector <real> & aCoefficients, Matrix< real > & aWork, Vector< int_t > & aPivot )
    {
        BELFEM_ASSERT( aWork.n_rows() == 6 && aWork.n_cols() == 6, "Work Matrix must be 6x6" );
        BELFEM_ASSERT( aPivot.length() >= 6, "Pivot Vector must be at least of length 6" );

        // populate the matrix
        aWork( 0, 0 ) = aX1 * aX1 * aX1 * aX1 * aX1;
        aWork( 1, 0 ) = 5 * aX1 * aX1 * aX1 * aX1;
        aWork( 2, 0 ) = 20 * aX1 * aX1 * aX1;
        aWork( 3, 0 ) = aX2 * aX2 * aX2 * aX2 * aX2;
        aWork( 4, 0 ) = 5 * aX2 * aX2 * aX2 * aX2;
        aWork( 5, 0 ) = 20 * aX2 * aX2 * aX2;

        aWork( 0, 1 ) = aX1 * aX1 * aX1 * aX1;
        aWork( 1, 1 ) = 4 * aX1 * aX1 * aX1;
        aWork( 2, 1 ) = 12 * aX1 * aX1;
        aWork( 3, 1 ) = aX2 * aX2 * aX2 * aX2;
        aWork( 4, 1 ) = 4 * aX2 * aX2 * aX2;
        aWork( 5, 1 ) = 12 * aX2 * aX2;

        aWork( 0, 2 ) = aX1 * aX1 * aX1;
        aWork( 1, 2 ) = 3 * aX1 * aX1;
        aWork( 2, 2 ) = 6 * aX1;
        aWork( 3, 2 ) = aX2 * aX2 * aX2;
        aWork( 4, 2 ) = 3 * aX2 * aX2;
        aWork( 5, 2 ) = 6 * aX2;

        aWork( 0, 3 ) = aX1 * aX1;
        aWork( 1, 3 ) = 2 * aX1;
        aWork( 2, 3 ) = 2.0;
        aWork( 3, 3 ) = aX2 * aX2;
        aWork( 4, 3 ) = 2 * aX2;
        aWork( 5, 3 ) = 2.0;

        aWork( 0, 4 ) = aX1;
        aWork( 1, 4 ) = 1.0;
        aWork( 2, 4 ) = 0;
        aWork( 3, 4 ) = aX2;
        aWork( 4, 4 ) = 1.0;
        aWork( 5, 4 ) = 0.0;

        aWork( 0, 5 ) = 1.0;
        aWork( 1, 5 ) = 0.0;
        aWork( 2, 5 ) = 0.0;
        aWork( 3, 5 ) = 1.0;
        aWork( 4, 5 ) = 0.0;
        aWork( 5, 5 ) = 0.0;


        // create the right hand side
        aCoefficients.set_size( 6 );

        aCoefficients( 0 ) = aF1;
        aCoefficients( 1 ) = adF1dX;
        aCoefficients( 2 ) = ad2F1dX;
        aCoefficients( 3 ) = aF2;
        aCoefficients( 4 ) = adF2dX;
        aCoefficients( 5 ) = ad2F2dX;

        // solve the system and return the coefficients
        gesv( aWork, aCoefficients, aPivot );
    }

    inline void
    create_fifth_order_beam_poly( const real aX1, const real aF1, const real adF1dX, const real ad2F1dX,
                                  const real aX2, const real aF2, const real adF2dX, const real ad2F2dX,
                                  Vector <real> & aCoefficients )
    {
        // create the vandermonde matrix
        Matrix <real> tVandermonde( 6, 6 );

        Vector< int_t > tPivot( 6 );

        create_fifth_order_beam_poly( aX1, aF1, adF1dX, ad2F1dX, aX2, aF2, adF2dX, ad2F2dX, aCoefficients, tVandermonde, tPivot );
    }

}
#endif //BELFEM_FN_CREATE_FIFTH_ORDER_BEAM_POLY_HPP
