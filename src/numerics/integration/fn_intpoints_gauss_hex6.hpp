//
// Created by Christian Messe on 15.07.20.
//

#ifndef BELFEM_FN_INITPOINTS_GAUSS_HEX6_HPP
#define BELFEM_FN_INITPOINTS_GAUSS_HEX6_HPP
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
        gauss_hex6(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aWeights.set_size( 6 );
            aPoints.set_size( 3, 6 ) ;

            aPoints( 0, 0 ) =   1.0000000000000000 ;
            aPoints( 1, 0 ) =   0.0000000000000000 ;
            aPoints( 2, 0 ) =   0.0000000000000000 ;

            aPoints( 0, 1 ) =  -1.0000000000000000 ;
            aPoints( 1, 1 ) =   0.0000000000000000 ;
            aPoints( 2, 1 ) =   0.0000000000000000 ;

            aPoints( 0, 2 ) =   0.0000000000000000 ;
            aPoints( 1, 2 ) =  1.0000000000000000 ;
            aPoints( 2, 2 ) =   0.0000000000000000 ;

            aPoints( 0, 3 ) =   0.0000000000000000 ;
            aPoints( 1, 3 ) =  -1.0000000000000000 ;
            aPoints( 2, 3 ) =   0.0000000000000000 ;

            aPoints( 0, 4 ) =   0.0000000000000000 ;
            aPoints( 1, 4 ) =   0.0000000000000000 ;
            aPoints( 2, 4 ) =   1.0000000000000000 ;

            aPoints( 0, 5 ) =   0.0000000000000000 ;
            aPoints( 1, 5 ) =   0.0000000000000000 ;
            aPoints( 2, 5 ) =  -1.0000000000000000 ;

            aWeights( 0 ) = 1.3333333333333333 ;
            aWeights( 1 ) = 1.3333333333333333 ;
            aWeights( 2 ) = 1.3333333333333333 ;
            aWeights( 3 ) = 1.3333333333333333 ;
            aWeights( 4 ) = 1.3333333333333333 ;
            aWeights( 5 ) = 1.3333333333333333 ;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */
#endif //BELFEM_FN_INITPOINTS_GAUSS_HEX6_HPP
