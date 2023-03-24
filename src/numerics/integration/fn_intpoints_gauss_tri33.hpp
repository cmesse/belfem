// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI33_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI33_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * 12th order interpolation
         *
         * Hong Xiao, Zydrunas Gimbutas :
         * A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions
         * Computers & Mathematics with Applications, vol. 59, no. 2, pp. 663–676, 2010
         * https://doi.org/10.1016/j.camwa.2009.10.027
         */
        inline void
        gauss_tri33(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 33 );

            aPoints( 0, 0 ) = 0.27146250701492614;
            aPoints( 1, 0 ) = 0.27146250701492614;

            aPoints( 0, 1 ) = 0.10925782765935432;
            aPoints( 1, 1 ) = 0.10925782765935432;

            aPoints( 0, 2 ) = 0.4401116486585931;
            aPoints( 1, 2 ) = 0.4401116486585931;

            aPoints( 0, 3 ) = 0.4882037509455415;
            aPoints( 1, 3 ) = 0.4882037509455415;

            aPoints( 0, 4 ) = 0.02464636343633564;
            aPoints( 1, 4 ) = 0.02464636343633564;

            aPoints( 0, 5 ) = 0.27146250701492614;
            aPoints( 1, 5 ) = 0.45707498597014773;

            aPoints( 0, 6 ) = 0.10925782765935432;
            aPoints( 1, 6 ) = 0.7814843446812914;

            aPoints( 0, 7 ) = 0.4401116486585931;
            aPoints( 1, 7 ) = 0.11977670268281382;

            aPoints( 0, 8 ) = 0.4882037509455415;
            aPoints( 1, 8 ) = 0.02359249810891695;

            aPoints( 0, 9 ) = 0.02464636343633564;
            aPoints( 1, 9 ) = 0.9507072731273287;

            aPoints( 0, 10 ) = 0.45707498597014773;
            aPoints( 1, 10 ) = 0.27146250701492614;

            aPoints( 0, 11 ) = 0.7814843446812914;
            aPoints( 1, 11 ) = 0.10925782765935432;

            aPoints( 0, 12 ) = 0.11977670268281382;
            aPoints( 1, 12 ) = 0.4401116486585931;

            aPoints( 0, 13 ) = 0.02359249810891695;
            aPoints( 1, 13 ) = 0.4882037509455415;

            aPoints( 0, 14 ) = 0.9507072731273287;
            aPoints( 1, 14 ) = 0.02464636343633564;

            aPoints( 0, 15 ) = 0.1162960196779266;
            aPoints( 1, 15 ) = 0.25545422863851736;

            aPoints( 0, 16 ) = 0.021382490256170623;
            aPoints( 1, 16 ) = 0.12727971723358936;

            aPoints( 0, 17 ) = 0.023034156355267166;
            aPoints( 1, 17 ) = 0.29165567973834094;

            aPoints( 0, 18 ) = 0.6282497516835561;
            aPoints( 1, 18 ) = 0.1162960196779266;

            aPoints( 0, 19 ) = 0.85133779251024;
            aPoints( 1, 19 ) = 0.021382490256170623;

            aPoints( 0, 20 ) = 0.6853101639063919;
            aPoints( 1, 20 ) = 0.023034156355267166;

            aPoints( 0, 21 ) = 0.25545422863851736;
            aPoints( 1, 21 ) = 0.6282497516835561;

            aPoints( 0, 22 ) = 0.12727971723358936;
            aPoints( 1, 22 ) = 0.85133779251024;

            aPoints( 0, 23 ) = 0.29165567973834094;
            aPoints( 1, 23 ) = 0.6853101639063919;

            aPoints( 0, 24 ) = 0.25545422863851736;
            aPoints( 1, 24 ) = 0.1162960196779266;

            aPoints( 0, 25 ) = 0.12727971723358936;
            aPoints( 1, 25 ) = 0.021382490256170623;

            aPoints( 0, 26 ) = 0.29165567973834094;
            aPoints( 1, 26 ) = 0.023034156355267166;

            aPoints( 0, 27 ) = 0.6282497516835561;
            aPoints( 1, 27 ) = 0.25545422863851736;

            aPoints( 0, 28 ) = 0.85133779251024;
            aPoints( 1, 28 ) = 0.12727971723358936;

            aPoints( 0, 29 ) = 0.6853101639063919;
            aPoints( 1, 29 ) = 0.29165567973834094;

            aPoints( 0, 30 ) = 0.1162960196779266;
            aPoints( 1, 30 ) = 0.6282497516835561;

            aPoints( 0, 31 ) = 0.021382490256170623;
            aPoints( 1, 31 ) = 0.85133779251024;

            aPoints( 0, 32 ) = 0.023034156355267166;
            aPoints( 1, 32 ) = 0.6853101639063919;

            for( uint k=0; k<33; ++k )
            {
                aPoints( 2, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k );
            }

            aWeights.set_size( 33 );

            aWeights( 0 ) = 0.03127060659795138;
            aWeights( 1 ) = 0.014243026034438775;
            aWeights( 2 ) = 0.024959167464030475;
            aWeights( 3 ) = 0.012133419040726017;
            aWeights( 4 ) = 0.0039658212549868194;
            aWeights( 5 ) = 0.03127060659795138;
            aWeights( 6 ) = 0.014243026034438775;
            aWeights( 7 ) = 0.024959167464030475;
            aWeights( 8 ) = 0.012133419040726017;
            aWeights( 9 ) = 0.0039658212549868194;
            aWeights( 10 ) = 0.03127060659795138;
            aWeights( 11 ) = 0.014243026034438775;
            aWeights( 12 ) = 0.024959167464030475;
            aWeights( 13 ) = 0.012133419040726017;
            aWeights( 14 ) = 0.0039658212549868194;
            aWeights( 15 ) = 0.021613681829707104;
            aWeights( 16 ) = 0.007541838788255721;
            aWeights( 17 ) = 0.01089179251930378;
            aWeights( 18 ) = 0.021613681829707104;
            aWeights( 19 ) = 0.007541838788255721;
            aWeights( 20 ) = 0.01089179251930378;
            aWeights( 21 ) = 0.021613681829707104;
            aWeights( 22 ) = 0.007541838788255721;
            aWeights( 23 ) = 0.01089179251930378;
            aWeights( 24 ) = 0.021613681829707104;
            aWeights( 25 ) = 0.007541838788255721;
            aWeights( 26 ) = 0.01089179251930378;
            aWeights( 27 ) = 0.021613681829707104;
            aWeights( 28 ) = 0.007541838788255721;
            aWeights( 29 ) = 0.01089179251930378;
            aWeights( 30 ) = 0.021613681829707104;
            aWeights( 31 ) = 0.007541838788255721;
            aWeights( 32 ) = 0.01089179251930378;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI33_HPP
