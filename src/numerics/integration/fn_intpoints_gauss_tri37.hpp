// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI37_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI37_HPP

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
         * 13th order interpolation
         *
         * Hong Xiao, Zydrunas Gimbutas :
         * A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions
         * Computers & Mathematics with Applications, vol. 59, no. 2, pp. 663–676, 2010
         * https://doi.org/10.1016/j.camwa.2009.10.027
         */
        inline void
        gauss_tri37(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 37 );

            aPoints( 0, 0 ) = 1./3.;
            aPoints( 1, 0 ) = 1./3.;

            aPoints( 0, 1 ) = 0.4961358947410461;
            aPoints( 1, 1 ) = 0.4961358947410461;

            aPoints( 0, 2 ) = 0.4696086896534919;
            aPoints( 1, 2 ) = 0.4696086896534919;

            aPoints( 0, 3 ) = 0.23111028494908226;
            aPoints( 1, 3 ) = 0.23111028494908226;

            aPoints( 0, 4 ) = 0.4144775702790546;
            aPoints( 1, 4 ) = 0.4144775702790546;

            aPoints( 0, 5 ) = 0.11355991257213327;
            aPoints( 1, 5 ) = 0.11355991257213327;

            aPoints( 0, 6 ) = 0.024895931491216494;
            aPoints( 1, 6 ) = 0.024895931491216494;

            aPoints( 0, 7 ) = 0.4961358947410461;
            aPoints( 1, 7 ) = 0.007728210517907841;

            aPoints( 0, 8 ) = 0.4696086896534919;
            aPoints( 1, 8 ) = 0.06078262069301621;

            aPoints( 0, 9 ) = 0.23111028494908226;
            aPoints( 1, 9 ) = 0.5377794301018355;

            aPoints( 0, 10 ) = 0.4144775702790546;
            aPoints( 1, 10 ) = 0.17104485944189085;

            aPoints( 0, 11 ) = 0.11355991257213327;
            aPoints( 1, 11 ) = 0.7728801748557335;

            aPoints( 0, 12 ) = 0.024895931491216494;
            aPoints( 1, 12 ) = 0.950208137017567;

            aPoints( 0, 13 ) = 0.007728210517907841;
            aPoints( 1, 13 ) = 0.4961358947410461;

            aPoints( 0, 14 ) = 0.06078262069301621;
            aPoints( 1, 14 ) = 0.4696086896534919;

            aPoints( 0, 15 ) = 0.5377794301018355;
            aPoints( 1, 15 ) = 0.23111028494908226;

            aPoints( 0, 16 ) = 0.17104485944189085;
            aPoints( 1, 16 ) = 0.4144775702790546;

            aPoints( 0, 17 ) = 0.7728801748557335;
            aPoints( 1, 17 ) = 0.11355991257213327;

            aPoints( 0, 18 ) = 0.950208137017567;
            aPoints( 1, 18 ) = 0.024895931491216494;

            aPoints( 0, 19 ) = 0.01898800438375904;
            aPoints( 1, 19 ) = 0.2920786885766364;

            aPoints( 0, 20 ) = 0.09773603106601653;
            aPoints( 1, 20 ) = 0.26674525331035115;

            aPoints( 0, 21 ) = 0.021966344206529244;
            aPoints( 1, 21 ) = 0.1267997757838373;

            aPoints( 0, 22 ) = 0.6889333070396046;
            aPoints( 1, 22 ) = 0.01898800438375904;

            aPoints( 0, 23 ) = 0.6355187156236324;
            aPoints( 1, 23 ) = 0.09773603106601653;

            aPoints( 0, 24 ) = 0.8512338800096335;
            aPoints( 1, 24 ) = 0.021966344206529244;

            aPoints( 0, 25 ) = 0.2920786885766364;
            aPoints( 1, 25 ) = 0.6889333070396046;

            aPoints( 0, 26 ) = 0.26674525331035115;
            aPoints( 1, 26 ) = 0.6355187156236324;

            aPoints( 0, 27 ) = 0.1267997757838373;
            aPoints( 1, 27 ) = 0.8512338800096335;

            aPoints( 0, 28 ) = 0.2920786885766364;
            aPoints( 1, 28 ) = 0.01898800438375904;

            aPoints( 0, 29 ) = 0.26674525331035115;
            aPoints( 1, 29 ) = 0.09773603106601653;

            aPoints( 0, 30 ) = 0.1267997757838373;
            aPoints( 1, 30 ) = 0.021966344206529244;

            aPoints( 0, 31 ) = 0.6889333070396046;
            aPoints( 1, 31 ) = 0.2920786885766364;

            aPoints( 0, 32 ) = 0.6355187156236324;
            aPoints( 1, 32 ) = 0.26674525331035115;

            aPoints( 0, 33 ) = 0.8512338800096335;
            aPoints( 1, 33 ) = 0.1267997757838373;

            aPoints( 0, 34 ) = 0.01898800438375904;
            aPoints( 1, 34 ) = 0.6889333070396046;

            aPoints( 0, 35 ) = 0.09773603106601653;
            aPoints( 1, 35 ) = 0.6355187156236324;

            aPoints( 0, 36 ) = 0.021966344206529244;
            aPoints( 1, 36 ) = 0.8512338800096335;

            for( uint k=0; k<37; ++k )
            {
                aPoints( 2, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k );
            }

            aWeights.set_size( 37 );

            aWeights( 0 ) = 0.02581132333214541;
            aWeights( 1 ) = 0.004970738180536294;
            aWeights( 2 ) = 0.01639062080186149;
            aWeights( 3 ) = 0.023031204796389124;
            aWeights( 4 ) = 0.0234735477710776;
            aWeights( 5 ) = 0.015451548987879897;
            aWeights( 6 ) = 0.0040146998976292115;
            aWeights( 7 ) = 0.004970738180536294;
            aWeights( 8 ) = 0.01639062080186149;
            aWeights( 9 ) = 0.023031204796389124;
            aWeights( 10 ) = 0.0234735477710776;
            aWeights( 11 ) = 0.015451548987879897;
            aWeights( 12 ) = 0.0040146998976292115;
            aWeights( 13 ) = 0.004970738180536294;
            aWeights( 14 ) = 0.01639062080186149;
            aWeights( 15 ) = 0.023031204796389124;
            aWeights( 16 ) = 0.0234735477710776;
            aWeights( 17 ) = 0.015451548987879897;
            aWeights( 18 ) = 0.0040146998976292115;
            aWeights( 19 ) = 0.00906274932310044;
            aWeights( 20 ) = 0.018605980228630768;
            aWeights( 21 ) = 0.007696536341891089;
            aWeights( 22 ) = 0.00906274932310044;
            aWeights( 23 ) = 0.018605980228630768;
            aWeights( 24 ) = 0.007696536341891089;
            aWeights( 25 ) = 0.00906274932310044;
            aWeights( 26 ) = 0.018605980228630768;
            aWeights( 27 ) = 0.007696536341891089;
            aWeights( 28 ) = 0.00906274932310044;
            aWeights( 29 ) = 0.018605980228630768;
            aWeights( 30 ) = 0.007696536341891089;
            aWeights( 31 ) = 0.00906274932310044;
            aWeights( 32 ) = 0.018605980228630768;
            aWeights( 33 ) = 0.007696536341891089;
            aWeights( 34 ) = 0.00906274932310044;
            aWeights( 35 ) = 0.018605980228630768;
            aWeights( 36 ) = 0.007696536341891089;

            aWeights( 0 )  = 0.0;
            aWeights( 0 )  = 0.5 - sum( aWeights );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI37_HPP
