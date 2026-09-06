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

#ifndef BELFEM_FN_CREATE_BEAM_POLY_HPP
#define BELFEM_FN_CREATE_BEAM_POLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "../../linalg/lapack/fn_gesv.hpp"

namespace belfem
{


    inline void
    create_beam_poly( const real aX1, const real aF1, const real adF1dX,
                  const real aX2, const real aF2, const real adF2dX,
                  Vector< real > & aCoefficients,
                  Matrix< real > & aWork,
                  Vector< int_t >  & aPivot )
    {
        BELFEM_ASSERT( aWork.n_rows() == 4 && aWork.n_cols() == 4, "Work Matrix must be 4x4" );
        BELFEM_ASSERT( aPivot.length() >= 4, "Pivot Vector must be at least of length 4" );

        // populate the matrix
        aWork( 0, 0 ) = aX1*aX1*aX1;
        aWork( 1, 0 ) = 3.0 * aX1 * aX1;
        aWork( 2, 0 ) = aX2*aX2*aX2;
        aWork( 3, 0 ) = 3.0 * aX2 * aX2;
        aWork( 0, 1 ) = aX1*aX1;
        aWork( 1, 1 ) = 2.0 * aX1;
        aWork( 2, 1 ) = aX2*aX2;
        aWork( 3, 1 ) = 2.0 * aX2 ;
        aWork( 0, 2 ) = aX1;
        aWork( 1, 2 ) = 1.0;
        aWork( 2, 2 ) = aX2;
        aWork( 3, 2 ) = 1.0;
        aWork( 0, 3 ) = 1.0;
        aWork( 1, 3 ) = 0.0;
        aWork( 2, 3 ) = 1.0;
        aWork( 3, 3 ) = 0.0;

        // create the right hand side
        aCoefficients.set_size( 4 );

        aCoefficients( 0 ) = aF1;
        aCoefficients( 1 ) = adF1dX;
        aCoefficients( 2 ) = aF2;
        aCoefficients( 3 ) = adF2dX;

        // solve the system and return the coefficients
        gesv( aWork, aCoefficients, aPivot );
    }

    inline void
    create_beam_poly( const real aX1, const real aF1, const real adF1dX,
                  const real aX2, const real aF2, const real adF2dX,
                  Vector< real > & aCoefficients )
    {
        // create the vandermonde matrix
        Matrix< real > tVandermonde( 4, 4 );

        // allocate the pivot vector
        Vector< int_t > tPivot( 4 );

        create_beam_poly( aX1, aF1, adF1dX, aX2, aF2, adF2dX, aCoefficients, tVandermonde, tPivot );
    }
}

#endif //BELFEM_FN_CREATE_BEAM_POLY_HPP
