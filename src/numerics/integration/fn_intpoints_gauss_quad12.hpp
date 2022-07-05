//
// Created by Christian Messe on 18.07.20.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_QUAD12_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_QUAD12_HPP


#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
//----------------------------------------------------------------------------

        // doi.org/10.1090/S0025-5718-1958-0102176-6
        inline void
        gauss_quad12(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aPoints.set_size( 2, 12 ) ;
            aWeights.set_size( 12 ) ;

            aPoints( 0,  0 ) =  0.9258200997725514 ;
            aPoints( 1,  0 ) =  0.0000000000000000 ;

            aPoints( 0,  1 ) = -0.9258200997725514 ;
            aPoints( 1,  1 ) =  0.0000000000000000 ;

            aPoints( 0,  2 ) =  0.0000000000000000 ;
            aPoints( 1,  2 ) =  0.9258200997725514 ;

            aPoints( 0,  3 ) =  0.0000000000000000 ;
            aPoints( 1,  3 ) = -0.9258200997725514 ;

            aPoints( 0,  4 ) =  0.3805544332083157 ;
            aPoints( 1,  4 ) =  0.3805544332083157 ;

            aPoints( 0,  5 ) =  0.3805544332083157 ;
            aPoints( 1,  5 ) = -0.3805544332083157 ;

            aPoints( 0,  6 ) = -0.3805544332083157 ;
            aPoints( 1,  6 ) =  0.3805544332083157 ;

            aPoints( 0,  7 ) = -0.3805544332083157 ;
            aPoints( 1,  7 ) = -0.3805544332083157 ;

            aPoints( 0,  8 ) =  0.8059797829185987 ;
            aPoints( 1,  8 ) =  0.8059797829185987 ;

            aPoints( 0,  9 ) =  0.8059797829185987 ;
            aPoints( 1,  9 ) = -0.8059797829185987 ;

            aPoints( 0, 10 ) = -0.8059797829185987 ;
            aPoints( 1, 10 ) =  0.8059797829185987 ;

            aPoints( 0, 11 ) = -0.8059797829185987 ;
            aPoints( 1, 11 ) = -0.8059797829185987 ;

            aWeights(  0 ) = 0.2419753086419753 ;
            aWeights(  1 ) = 0.2419753086419753 ;
            aWeights(  2 ) = 0.2419753086419753 ;
            aWeights(  3 ) = 0.2419753086419753 ;
            aWeights(  4 ) = 0.5205929166673945 ;
            aWeights(  5 ) = 0.5205929166673945 ;
            aWeights(  6 ) = 0.5205929166673945 ;
            aWeights(  7 ) = 0.5205929166673945 ;
            aWeights(  8 ) = 0.2374317746906302 ;
            aWeights(  9 ) = 0.2374317746906302 ;
            aWeights( 10 ) = 0.2374317746906302 ;
            aWeights( 11 ) = 0.2374317746906302 ;
        }

//----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_QUAD12_HPP
