// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI28_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI28_HPP

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
         * 11th order interpolation
         *
         * Hong Xiao, Zydrunas Gimbutas :
         * A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions
         * Computers & Mathematics with Applications, vol. 59, no. 2, pp. 663–676, 2010
         * https://doi.org/10.1016/j.camwa.2009.10.027
         */
        inline void
        gauss_tri28(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 28 );

            aPoints( 0, 0 ) = 1./3.;
            aPoints( 1, 0 ) = 1./3.;

            aPoints( 0, 1 ) = 0.030846895635588123;
            aPoints( 1, 1 ) = 0.030846895635588123;

            aPoints( 0, 2 ) = 0.49878016517846074;
            aPoints( 1, 2 ) = 0.49878016517846074;

            aPoints( 0, 3 ) = 0.11320782728669404;
            aPoints( 1, 3 ) = 0.11320782728669404;

            aPoints( 0, 4 ) = 0.4366550163931761;
            aPoints( 1, 4 ) = 0.4366550163931761;

            aPoints( 0, 5 ) = 0.21448345861926937;
            aPoints( 1, 5 ) = 0.21448345861926937;

            aPoints( 0, 6 ) = 0.030846895635588123;
            aPoints( 1, 6 ) = 0.9383062087288238;

            aPoints( 0, 7 ) = 0.49878016517846074;
            aPoints( 1, 7 ) = 0.0024396696430785125;

            aPoints( 0, 8 ) = 0.11320782728669404;
            aPoints( 1, 8 ) = 0.7735843454266119;

            aPoints( 0, 9 ) = 0.4366550163931761;
            aPoints( 1, 9 ) = 0.12668996721364778;

            aPoints( 0, 10 ) = 0.21448345861926937;
            aPoints( 1, 10 ) = 0.5710330827614613;

            aPoints( 0, 11 ) = 0.9383062087288238;
            aPoints( 1, 11 ) = 0.030846895635588123;

            aPoints( 0, 12 ) = 0.0024396696430785125;
            aPoints( 1, 12 ) = 0.49878016517846074;

            aPoints( 0, 13 ) = 0.7735843454266119;
            aPoints( 1, 13 ) = 0.11320782728669404;

            aPoints( 0, 14 ) = 0.12668996721364778;
            aPoints( 1, 14 ) = 0.4366550163931761;

            aPoints( 0, 15 ) = 0.5710330827614613;
            aPoints( 1, 15 ) = 0.21448345861926937;

            aPoints( 0, 16 ) = 0.014366662569555624;
            aPoints( 1, 16 ) = 0.1593036198376935;

            aPoints( 0, 17 ) = 0.04766406697215078;
            aPoints( 1, 17 ) = 0.31063121631346313;

            aPoints( 0, 18 ) = 0.8263297175927509;
            aPoints( 1, 18 ) = 0.014366662569555624;

            aPoints( 0, 19 ) = 0.6417047167143861;
            aPoints( 1, 19 ) = 0.04766406697215078;

            aPoints( 0, 20 ) = 0.1593036198376935;
            aPoints( 1, 20 ) = 0.8263297175927509;

            aPoints( 0, 21 ) = 0.31063121631346313;
            aPoints( 1, 21 ) = 0.6417047167143861;

            aPoints( 0, 22 ) = 0.1593036198376935;
            aPoints( 1, 22 ) = 0.014366662569555624;

            aPoints( 0, 23 ) = 0.31063121631346313;
            aPoints( 1, 23 ) = 0.04766406697215078;

            aPoints( 0, 24 ) = 0.8263297175927509;
            aPoints( 1, 24 ) = 0.1593036198376935;

            aPoints( 0, 25 ) = 0.6417047167143861;
            aPoints( 1, 25 ) = 0.31063121631346313;

            aPoints( 0, 26 ) = 0.014366662569555624;
            aPoints( 1, 26 ) = 0.8263297175927509;

            aPoints( 0, 27 ) = 0.04766406697215078;
            aPoints( 1, 27 ) = 0.6417047167143861;

            for( uint k=0; k<28; ++k )
            {
                aPoints( 2, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k );
            }

            aWeights.set_size( 28 );

            aWeights( 0 ) = 0.040722567354675644;
            aWeights( 1 ) = 0.006124648475353982;
            aWeights( 2 ) = 0.0062327459369406904;
            aWeights( 3 ) = 0.02006462119065416;
            aWeights( 4 ) = 0.031547436079949344;
            aWeights( 5 ) = 0.033922553871847574;
            aWeights( 6 ) = 0.006124648475353982;
            aWeights( 7 ) = 0.0062327459369406904;
            aWeights( 8 ) = 0.02006462119065416;
            aWeights( 9 ) = 0.031547436079949344;
            aWeights( 10 ) = 0.033922553871847574;
            aWeights( 11 ) = 0.006124648475353982;
            aWeights( 12 ) = 0.0062327459369406904;
            aWeights( 13 ) = 0.02006462119065416;
            aWeights( 14 ) = 0.031547436079949344;
            aWeights( 15 ) = 0.033922553871847574;
            aWeights( 16 ) = 0.007278811668904623;
            aWeights( 17 ) = 0.020321424327943236;
            aWeights( 18 ) = 0.007278811668904623;
            aWeights( 19 ) = 0.020321424327943236;
            aWeights( 20 ) = 0.007278811668904623;
            aWeights( 21 ) = 0.020321424327943236;
            aWeights( 22 ) = 0.007278811668904623;
            aWeights( 23 ) = 0.020321424327943236;
            aWeights( 24 ) = 0.007278811668904623;
            aWeights( 25 ) = 0.020321424327943236;
            aWeights( 26 ) = 0.007278811668904623;
            aWeights( 27 ) = 0.020321424327943236;

            aWeights( 0 )  = 0.0;
            aWeights( 0 )  = 0.5 - sum( aWeights );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI28_HPP
