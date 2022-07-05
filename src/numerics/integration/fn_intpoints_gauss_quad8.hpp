//
// Created by Christian Messe on 18.07.20.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_QUAD8_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_QUAD8_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        // doi.org/10.1016/j.camwa.2015.03.017
        inline void
        gauss_quad8(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aPoints.set_size( 2, 8 ) ;
            aWeights.set_size( 8 ) ;

            aPoints( 0, 0 ) =  0.6831300510639732 ;
            aPoints( 1, 0 ) =  0.0000000000000000 ;

            aPoints( 0, 1 ) =  -0.6831300510639732 ;
            aPoints( 1, 1 ) =  0.0000000000000000 ;

            aPoints( 0, 2 ) =  0.0000000000000000 ;
            aPoints( 1, 2 ) =  0.6831300510639732 ;

            aPoints( 0, 3 ) =  0.0000000000000000 ;
            aPoints( 1, 3 ) =  -0.6831300510639732 ;

            aPoints( 0, 4 ) =  0.8819171036881969 ;
            aPoints( 1, 4 ) =  0.8819171036881969 ;

            aPoints( 0, 5 ) =  0.8819171036881969 ;
            aPoints( 1, 5 ) =  -0.8819171036881969 ;

            aPoints( 0, 6 ) =  -0.8819171036881969 ;
            aPoints( 1, 6 ) =  0.8819171036881969 ;

            aPoints( 0, 7 ) =  -0.8819171036881969 ;
            aPoints( 1, 7 ) =  -0.8819171036881969 ;

            aWeights( 0 ) = 0.8163265306122449 ;
            aWeights( 1 ) = 0.8163265306122449 ;
            aWeights( 2 ) = 0.8163265306122449 ;
            aWeights( 3 ) = 0.8163265306122449 ;
            aWeights( 4 ) = 0.1836734693877551 ;
            aWeights( 5 ) = 0.1836734693877551 ;
            aWeights( 6 ) = 0.1836734693877551 ;
            aWeights( 7 ) = 0.1836734693877551 ;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_QUAD8_HPP
