// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI12_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI12_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * 6th order interpolation
         *
         * G.R. Cowper :
         * Gaussian quadrature formulas for triangles
         * Numerical Methods in Engineering, vol. 7, no. 3, pp. 405–408, 1973
         * https://doi.org/10.1002/nme.1620070316
         */
        inline void
        gauss_tri12(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 12 );

            aPoints( 0, 0 ) = 0.06308901449150225;
            aPoints( 1, 0 ) = 0.06308901449150225;

            aPoints( 0, 1 ) = 0.2492867451709104;
            aPoints( 1, 1 ) = 0.2492867451709104;

            aPoints( 0, 2 ) = 0.06308901449150225;
            aPoints( 1, 2 ) = 0.8738219710169954;

            aPoints( 0, 3 ) = 0.2492867451709104;
            aPoints( 1, 3 ) = 0.5014265096581791;

            aPoints( 0, 4 ) = 0.8738219710169954;
            aPoints( 1, 4 ) = 0.06308901449150225;

            aPoints( 0, 5 ) = 0.5014265096581791;
            aPoints( 1, 5 ) = 0.2492867451709104;

            aPoints( 0, 6 ) = 0.6365024991213987;
            aPoints( 1, 6 ) = 0.3103524510337844;

            aPoints( 0, 7 ) = 0.05314504984481694;
            aPoints( 1, 7 ) = 0.6365024991213987;

            aPoints( 0, 8 ) = 0.3103524510337844;
            aPoints( 1, 8 ) = 0.05314504984481694;

            aPoints( 0, 9 ) = 0.3103524510337844;
            aPoints( 1, 9 ) = 0.6365024991213987;

            aPoints( 0, 10 ) = 0.05314504984481694;
            aPoints( 1, 10 ) = 0.3103524510337844;

            aPoints( 0, 11 ) = 0.6365024991213987;
            aPoints( 1, 11 ) = 0.05314504984481694;

            for( uint k=0; k<12; ++k )
            {
                aPoints( 2, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k );
            }

            aWeights.set_size( 12 );

            aWeights( 0 ) = 0.025422453185103444;
            aWeights( 1 ) = 0.058393137863189656;
            aWeights( 2 ) = 0.025422453185103444;
            aWeights( 3 ) = 0.058393137863189656;
            aWeights( 4 ) = 0.025422453185103444;
            aWeights( 5 ) = 0.058393137863189656;
            aWeights( 6 ) = 0.04142553780918681;
            aWeights( 7 ) = 0.04142553780918681;
            aWeights( 8 ) = 0.04142553780918681;
            aWeights( 9 ) = 0.04142553780918681;
            aWeights( 10 ) = 0.04142553780918681;
            aWeights( 11 ) = 0.04142553780918681;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI12_HPP
