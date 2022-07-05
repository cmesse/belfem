//
// Created by Christian Messe on 15.07.20.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_HEX14_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_HEX14_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        // source 10.1090/S0025-5718-1958-0102176-6
        inline void
        gauss_hex14(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aWeights.set_size( 14 );
            aPoints.set_size( 3, 14 ) ;

            aPoints( 0, 0 ) =   0.7958224257542215 ;
            aPoints( 1, 0 ) =   0.0000000000000000 ;
            aPoints( 2, 0 ) =   0.0000000000000000 ;

            aPoints( 0, 1 ) =  -0.7958224257542215 ;
            aPoints( 1, 1 ) =   0.0000000000000000 ;
            aPoints( 2, 1 ) =   0.0000000000000000 ;

            aPoints( 0, 2 ) =   0.0000000000000000 ;
            aPoints( 1, 2 ) =   0.7958224257542215 ;
            aPoints( 2, 2 ) =   0.0000000000000000 ;

            aPoints( 0, 3 ) =   0.0000000000000000 ;
            aPoints( 1, 3 ) =  -0.7958224257542215 ;
            aPoints( 2, 3 ) =   0.0000000000000000 ;

            aPoints( 0, 4 ) =   0.0000000000000000 ;
            aPoints( 1, 4 ) =   0.0000000000000000 ;
            aPoints( 2, 4 ) =   0.7958224257542215 ;

            aPoints( 0, 5 ) =   0.0000000000000000 ;
            aPoints( 1, 5 ) =   0.0000000000000000 ;
            aPoints( 2, 5 ) =  -0.7958224257542215 ;

            aPoints( 0, 6 ) =   0.7587869106393281 ;
            aPoints( 1, 6 ) =   0.7587869106393281 ;
            aPoints( 2, 6 ) =   0.7587869106393281 ;

            aPoints( 0, 7 ) =   0.7587869106393281 ;
            aPoints( 1, 7 ) =   0.7587869106393281 ;
            aPoints( 2, 7 ) =  -0.7587869106393281 ;

            aPoints( 0, 8 ) =   0.7587869106393281 ;
            aPoints( 1, 8 ) =  -0.7587869106393281 ;
            aPoints( 2, 8 ) =   0.7587869106393281 ;

            aPoints( 0, 9 ) =   0.7587869106393281 ;
            aPoints( 1, 9 ) =  -0.7587869106393281 ;
            aPoints( 2, 9 ) =  -0.7587869106393281 ;

            aPoints( 0, 10 ) =  -0.7587869106393281 ;
            aPoints( 1, 10 ) =   0.7587869106393281 ;
            aPoints( 2, 10 ) =   0.7587869106393281 ;

            aPoints( 0, 11 ) =  -0.7587869106393281 ;
            aPoints( 1, 11 ) =   0.7587869106393281 ;
            aPoints( 2, 11 ) =  -0.7587869106393281 ;

            aPoints( 0, 12 ) =  -0.7587869106393281 ;
            aPoints( 1, 12 ) =  -0.7587869106393281 ;
            aPoints( 2, 12 ) =   0.7587869106393281 ;

            aPoints( 0, 13 ) =  -0.7587869106393281 ;
            aPoints( 1, 13 ) =  -0.7587869106393281 ;
            aPoints( 2, 13 ) =  -0.7587869106393281 ;

            aWeights(  0 ) = 0.8864265927977839 ;
            aWeights(  1 ) = 0.8864265927977839 ;
            aWeights(  2 ) = 0.8864265927977839 ;
            aWeights(  3 ) = 0.8864265927977839 ;
            aWeights(  4 ) = 0.8864265927977839 ;
            aWeights(  5 ) = 0.8864265927977839 ;
            aWeights(  6 ) = 0.3351800554016621 ;
            aWeights(  7 ) = 0.3351800554016621 ;
            aWeights(  8 ) = 0.3351800554016621 ;
            aWeights(  9 ) = 0.3351800554016621 ;
            aWeights( 10 ) = 0.3351800554016621 ;
            aWeights( 11 ) = 0.3351800554016621 ;
            aWeights( 12 ) = 0.3351800554016621 ;
            aWeights( 13 ) = 0.3351800554016621 ;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_HEX14_HPP
