// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI16_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI16_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_sum.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * 8th order interpolation
         *
         * J.N. Lyness, D. Jespersen :
         * Moderate Degree Symmetric Quadrature Rules for the Triangle
         * IMA Journal of Applied Mathematics, vol. 15, no. 1, pp. 19-32, 1975
         * https://doi.org/10.1093/imamat/15.1.19
         */
        inline void
        gauss_tri16(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 16 );

            aPoints( 0, 0 ) = 1./3.;
            aPoints( 1, 0 ) = 1./3.;

            aPoints( 0, 1 ) = 0.4592925882927229;
            aPoints( 1, 1 ) = 0.4592925882927229;

            aPoints( 0, 2 ) = 0.05054722831703103;
            aPoints( 1, 2 ) = 0.05054722831703103;

            aPoints( 0, 3 ) = 0.1705693077517601;
            aPoints( 1, 3 ) = 0.1705693077517601;

            aPoints( 0, 4 ) = 0.4592925882927229;
            aPoints( 1, 4 ) = 0.0814148234145542;

            aPoints( 0, 5 ) = 0.05054722831703103;
            aPoints( 1, 5 ) = 0.898905543365938;

            aPoints( 0, 6 ) = 0.1705693077517601;
            aPoints( 1, 6 ) = 0.6588613844964798;

            aPoints( 0, 7 ) = 0.0814148234145542;
            aPoints( 1, 7 ) = 0.4592925882927229;

            aPoints( 0, 8 ) = 0.898905543365938;
            aPoints( 1, 8 ) = 0.05054722831703103;

            aPoints( 0, 9 ) = 0.6588613844964798;
            aPoints( 1, 9 ) = 0.1705693077517601;

            aPoints( 0, 10 ) = 0.008394777409957211;
            aPoints( 1, 10 ) = 0.7284923929554041;

            aPoints( 0, 11 ) = 0.26311282963463867;
            aPoints( 1, 11 ) = 0.008394777409957211;

            aPoints( 0, 12 ) = 0.7284923929554041;
            aPoints( 1, 12 ) = 0.26311282963463867;

            aPoints( 0, 13 ) = 0.7284923929554041;
            aPoints( 1, 13 ) = 0.008394777409957211;

            aPoints( 0, 14 ) = 0.26311282963463867;
            aPoints( 1, 14 ) = 0.7284923929554041;

            aPoints( 0, 15 ) = 0.008394777409957211;
            aPoints( 1, 15 ) = 0.26311282963463867;

            for( uint k=0; k<16; ++k )
            {
                aPoints( 2, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k );
            }

            aWeights.set_size( 16 );

            aWeights( 0 ) = 0.0721578038388931;
            aWeights( 1 ) = 0.04754581713364248;
            aWeights( 2 ) = 0.016229248811599067;
            aWeights( 3 ) = 0.0516086852673592;
            aWeights( 4 ) = 0.04754581713364248;
            aWeights( 5 ) = 0.016229248811599067;
            aWeights( 6 ) = 0.0516086852673592;
            aWeights( 7 ) = 0.04754581713364248;
            aWeights( 8 ) = 0.016229248811599067;
            aWeights( 9 ) = 0.0516086852673592;
            aWeights( 10 ) = 0.013615157087217432;
            aWeights( 11 ) = 0.013615157087217432;
            aWeights( 12 ) = 0.013615157087217432;
            aWeights( 13 ) = 0.013615157087217432;
            aWeights( 14 ) = 0.013615157087217432;
            aWeights( 15 ) = 0.013615157087217432;

            aWeights( 0 )  = 0.0;
            aWeights( 0 )  = 0.5 - sum( aWeights );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI16_HPP
