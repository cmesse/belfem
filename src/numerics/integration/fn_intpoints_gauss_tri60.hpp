// auto generated with Python using the quadpy library

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TRI60_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TRI60_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        /**
         * 17th order interpolation
         *
         * Hong Xiao, Zydrunas Gimbutas :
         * A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions
         * Computers & Mathematics with Applications, vol. 59, no. 2, pp. 663–676, 2010
         * https://doi.org/10.1016/j.camwa.2009.10.027
         */
        inline void
        gauss_tri60(
                Vector< real > & aWeights,
                Matrix< real > & aPoints )
        {
            aPoints.set_size( 3, 60 );

            aPoints( 0, 0 ) = 0.4171034443615992;
            aPoints( 1, 0 ) = 0.4171034443615992;

            aPoints( 0, 1 ) = 0.18035811626637066;
            aPoints( 1, 1 ) = 0.18035811626637066;

            aPoints( 0, 2 ) = 0.2857065024365867;
            aPoints( 1, 2 ) = 0.2857065024365867;

            aPoints( 0, 3 ) = 0.06665406347959701;
            aPoints( 1, 3 ) = 0.06665406347959701;

            aPoints( 0, 4 ) = 0.014755491660754072;
            aPoints( 1, 4 ) = 0.014755491660754072;

            aPoints( 0, 5 ) = 0.46559787161889027;
            aPoints( 1, 5 ) = 0.46559787161889027;

            aPoints( 0, 6 ) = 0.4171034443615992;
            aPoints( 1, 6 ) = 0.16579311127680163;

            aPoints( 0, 7 ) = 0.18035811626637066;
            aPoints( 1, 7 ) = 0.6392837674672587;

            aPoints( 0, 8 ) = 0.2857065024365867;
            aPoints( 1, 8 ) = 0.42858699512682663;

            aPoints( 0, 9 ) = 0.06665406347959701;
            aPoints( 1, 9 ) = 0.866691873040806;

            aPoints( 0, 10 ) = 0.014755491660754072;
            aPoints( 1, 10 ) = 0.9704890166784919;

            aPoints( 0, 11 ) = 0.46559787161889027;
            aPoints( 1, 11 ) = 0.06880425676221946;

            aPoints( 0, 12 ) = 0.16579311127680163;
            aPoints( 1, 12 ) = 0.4171034443615992;

            aPoints( 0, 13 ) = 0.6392837674672587;
            aPoints( 1, 13 ) = 0.18035811626637066;

            aPoints( 0, 14 ) = 0.42858699512682663;
            aPoints( 1, 14 ) = 0.2857065024365867;

            aPoints( 0, 15 ) = 0.866691873040806;
            aPoints( 1, 15 ) = 0.06665406347959701;

            aPoints( 0, 16 ) = 0.9704890166784919;
            aPoints( 1, 16 ) = 0.014755491660754072;

            aPoints( 0, 17 ) = 0.06880425676221946;
            aPoints( 1, 17 ) = 0.46559787161889027;

            aPoints( 0, 18 ) = 0.011575175903180683;
            aPoints( 1, 18 ) = 0.07250547079900238;

            aPoints( 0, 19 ) = 0.013229672760086951;
            aPoints( 1, 19 ) = 0.41547545929522905;

            aPoints( 0, 20 ) = 0.013135870834002753;
            aPoints( 1, 20 ) = 0.27179187005535477;

            aPoints( 0, 21 ) = 0.15750547792686992;
            aPoints( 1, 21 ) = 0.29921894247697034;

            aPoints( 0, 22 ) = 0.06734937786736123;
            aPoints( 1, 22 ) = 0.3062815917461865;

            aPoints( 0, 23 ) = 0.07804234056828245;
            aPoints( 1, 23 ) = 0.16872251349525944;

            aPoints( 0, 24 ) = 0.016017642362119337;
            aPoints( 1, 24 ) = 0.15919228747279268;

            aPoints( 0, 25 ) = 0.9159193532978169;
            aPoints( 1, 25 ) = 0.011575175903180683;

            aPoints( 0, 26 ) = 0.5712948679446841;
            aPoints( 1, 26 ) = 0.013229672760086951;

            aPoints( 0, 27 ) = 0.7150722591106424;
            aPoints( 1, 27 ) = 0.013135870834002753;

            aPoints( 0, 28 ) = 0.5432755795961598;
            aPoints( 1, 28 ) = 0.15750547792686992;

            aPoints( 0, 29 ) = 0.6263690303864522;
            aPoints( 1, 29 ) = 0.06734937786736123;

            aPoints( 0, 30 ) = 0.7532351459364581;
            aPoints( 1, 30 ) = 0.07804234056828245;

            aPoints( 0, 31 ) = 0.824790070165088;
            aPoints( 1, 31 ) = 0.016017642362119337;

            aPoints( 0, 32 ) = 0.07250547079900238;
            aPoints( 1, 32 ) = 0.9159193532978169;

            aPoints( 0, 33 ) = 0.41547545929522905;
            aPoints( 1, 33 ) = 0.5712948679446841;

            aPoints( 0, 34 ) = 0.27179187005535477;
            aPoints( 1, 34 ) = 0.7150722591106424;

            aPoints( 0, 35 ) = 0.29921894247697034;
            aPoints( 1, 35 ) = 0.5432755795961598;

            aPoints( 0, 36 ) = 0.3062815917461865;
            aPoints( 1, 36 ) = 0.6263690303864522;

            aPoints( 0, 37 ) = 0.16872251349525944;
            aPoints( 1, 37 ) = 0.7532351459364581;

            aPoints( 0, 38 ) = 0.15919228747279268;
            aPoints( 1, 38 ) = 0.824790070165088;

            aPoints( 0, 39 ) = 0.07250547079900238;
            aPoints( 1, 39 ) = 0.011575175903180683;

            aPoints( 0, 40 ) = 0.41547545929522905;
            aPoints( 1, 40 ) = 0.013229672760086951;

            aPoints( 0, 41 ) = 0.27179187005535477;
            aPoints( 1, 41 ) = 0.013135870834002753;

            aPoints( 0, 42 ) = 0.29921894247697034;
            aPoints( 1, 42 ) = 0.15750547792686992;

            aPoints( 0, 43 ) = 0.3062815917461865;
            aPoints( 1, 43 ) = 0.06734937786736123;

            aPoints( 0, 44 ) = 0.16872251349525944;
            aPoints( 1, 44 ) = 0.07804234056828245;

            aPoints( 0, 45 ) = 0.15919228747279268;
            aPoints( 1, 45 ) = 0.016017642362119337;

            aPoints( 0, 46 ) = 0.9159193532978169;
            aPoints( 1, 46 ) = 0.07250547079900238;

            aPoints( 0, 47 ) = 0.5712948679446841;
            aPoints( 1, 47 ) = 0.41547545929522905;

            aPoints( 0, 48 ) = 0.7150722591106424;
            aPoints( 1, 48 ) = 0.27179187005535477;

            aPoints( 0, 49 ) = 0.5432755795961598;
            aPoints( 1, 49 ) = 0.29921894247697034;

            aPoints( 0, 50 ) = 0.6263690303864522;
            aPoints( 1, 50 ) = 0.3062815917461865;

            aPoints( 0, 51 ) = 0.7532351459364581;
            aPoints( 1, 51 ) = 0.16872251349525944;

            aPoints( 0, 52 ) = 0.824790070165088;
            aPoints( 1, 52 ) = 0.15919228747279268;

            aPoints( 0, 53 ) = 0.011575175903180683;
            aPoints( 1, 53 ) = 0.9159193532978169;

            aPoints( 0, 54 ) = 0.013229672760086951;
            aPoints( 1, 54 ) = 0.5712948679446841;

            aPoints( 0, 55 ) = 0.013135870834002753;
            aPoints( 1, 55 ) = 0.7150722591106424;

            aPoints( 0, 56 ) = 0.15750547792686992;
            aPoints( 1, 56 ) = 0.5432755795961598;

            aPoints( 0, 57 ) = 0.06734937786736123;
            aPoints( 1, 57 ) = 0.6263690303864522;

            aPoints( 0, 58 ) = 0.07804234056828245;
            aPoints( 1, 58 ) = 0.7532351459364581;

            aPoints( 0, 59 ) = 0.016017642362119337;
            aPoints( 1, 59 ) = 0.824790070165088;

            for( uint k=0; k<60; ++k )
            {
                aPoints( 2, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k );
            }

            aWeights.set_size( 60 );

            aWeights( 0 ) = 0.013655463264051053;
            aWeights( 1 ) = 0.013156315294008993;
            aWeights( 2 ) = 0.01885811857639764;
            aWeights( 3 ) = 0.006229500401152722;
            aWeights( 4 ) = 0.001386943788818821;
            aWeights( 5 ) = 0.01250972547524868;
            aWeights( 6 ) = 0.013655463264051053;
            aWeights( 7 ) = 0.013156315294008993;
            aWeights( 8 ) = 0.01885811857639764;
            aWeights( 9 ) = 0.006229500401152722;
            aWeights( 10 ) = 0.001386943788818821;
            aWeights( 11 ) = 0.01250972547524868;
            aWeights( 12 ) = 0.013655463264051053;
            aWeights( 13 ) = 0.013156315294008993;
            aWeights( 14 ) = 0.01885811857639764;
            aWeights( 15 ) = 0.006229500401152722;
            aWeights( 16 ) = 0.001386943788818821;
            aWeights( 17 ) = 0.01250972547524868;
            aWeights( 18 ) = 0.002292174200867934;
            aWeights( 19 ) = 0.005199219977919768;
            aWeights( 20 ) = 0.004346107250500596;
            aWeights( 21 ) = 0.013085812967668494;
            aWeights( 22 ) = 0.011243886273345534;
            aWeights( 23 ) = 0.01027894916022726;
            aWeights( 24 ) = 0.003989150102964797;
            aWeights( 25 ) = 0.002292174200867934;
            aWeights( 26 ) = 0.005199219977919768;
            aWeights( 27 ) = 0.004346107250500596;
            aWeights( 28 ) = 0.013085812967668494;
            aWeights( 29 ) = 0.011243886273345534;
            aWeights( 30 ) = 0.01027894916022726;
            aWeights( 31 ) = 0.003989150102964797;
            aWeights( 32 ) = 0.002292174200867934;
            aWeights( 33 ) = 0.005199219977919768;
            aWeights( 34 ) = 0.004346107250500596;
            aWeights( 35 ) = 0.013085812967668494;
            aWeights( 36 ) = 0.011243886273345534;
            aWeights( 37 ) = 0.01027894916022726;
            aWeights( 38 ) = 0.003989150102964797;
            aWeights( 39 ) = 0.002292174200867934;
            aWeights( 40 ) = 0.005199219977919768;
            aWeights( 41 ) = 0.004346107250500596;
            aWeights( 42 ) = 0.013085812967668494;
            aWeights( 43 ) = 0.011243886273345534;
            aWeights( 44 ) = 0.01027894916022726;
            aWeights( 45 ) = 0.003989150102964797;
            aWeights( 46 ) = 0.002292174200867934;
            aWeights( 47 ) = 0.005199219977919768;
            aWeights( 48 ) = 0.004346107250500596;
            aWeights( 49 ) = 0.013085812967668494;
            aWeights( 50 ) = 0.011243886273345534;
            aWeights( 51 ) = 0.01027894916022726;
            aWeights( 52 ) = 0.003989150102964797;
            aWeights( 53 ) = 0.002292174200867934;
            aWeights( 54 ) = 0.005199219977919768;
            aWeights( 55 ) = 0.004346107250500596;
            aWeights( 56 ) = 0.013085812967668494;
            aWeights( 57 ) = 0.011243886273345534;
            aWeights( 58 ) = 0.01027894916022726;
            aWeights( 59 ) = 0.003989150102964797;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* end namespace belfem */

#endif  // BELFEM_FN_INTPOINTS_GAUSS_TRI60_HPP
