//
// Created by Christian Messe on 18.07.20.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_QUAD20_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_QUAD20_HPP
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
        gauss_quad20(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aPoints.set_size( 2, 20 ) ;
            aWeights.set_size( 20 ) ;

            aPoints( 0,  0 ) =  0.4889268569743691 ;
            aPoints( 1,  0 ) =  0.0000000000000000 ;

            aPoints( 0,  1 ) = -0.4889268569743691 ;
            aPoints( 1,  1 ) =  0.0000000000000000 ;

            aPoints( 0,  2 ) =  0.0000000000000000 ;
            aPoints( 1,  2 ) =  0.4889268569743691 ;

            aPoints( 0,  3 ) =  0.0000000000000000 ;
            aPoints( 1,  3 ) = -0.4889268569743691 ;

            aPoints( 0,  4 ) =  0.9396552580968377 ;
            aPoints( 1,  4 ) =  0.9396552580968377 ;

            aPoints( 0,  5 ) =  0.6908805504863439 ;
            aPoints( 1,  5 ) =  0.6908805504863439 ;

            aPoints( 0,  6 ) =  0.9396552580968377 ;
            aPoints( 1,  6 ) = -0.9396552580968377 ;

            aPoints( 0,  7 ) =  0.6908805504863439 ;
            aPoints( 1,  7 ) = -0.6908805504863439 ;

            aPoints( 0,  8 ) = -0.9396552580968377 ;
            aPoints( 1,  8 ) =  0.9396552580968377 ;

            aPoints( 0,  9 ) = -0.6908805504863439 ;
            aPoints( 1,  9 ) =  0.6908805504863439 ;

            aPoints( 0, 10 ) = -0.9396552580968377 ;
            aPoints( 1, 10 ) = -0.9396552580968377 ;

            aPoints( 0, 11 ) = -0.6908805504863439 ;
            aPoints( 1, 11 ) = -0.6908805504863439 ;

            aPoints( 0, 12 ) =  0.9186204410567222 ;
            aPoints( 1, 12 ) =  0.3448720253644036 ;

            aPoints( 0, 13 ) = -0.9186204410567222 ;
            aPoints( 1, 13 ) =  0.3448720253644036 ;

            aPoints( 0, 14 ) =  0.9186204410567222 ;
            aPoints( 1, 14 ) = -0.3448720253644036 ;

            aPoints( 0, 15 ) = -0.9186204410567222 ;
            aPoints( 1, 15 ) = -0.3448720253644036 ;

            aPoints( 0, 16 ) =  0.3448720253644036 ;
            aPoints( 1, 16 ) =  0.9186204410567222 ;

            aPoints( 0, 17 ) = -0.3448720253644036 ;
            aPoints( 1, 17 ) =  0.9186204410567222 ;

            aPoints( 0, 18 ) =  0.3448720253644036 ;
            aPoints( 1, 18 ) = -0.9186204410567222 ;

            aPoints( 0, 19 ) = -0.3448720253644036 ;
            aPoints( 1, 19 ) = -0.9186204410567222 ;

            aWeights(  0 ) = 0.4541639606867490 ;
            aWeights(  1 ) = 0.4541639606867490 ;
            aWeights(  2 ) = 0.4541639606867490 ;
            aWeights(  3 ) = 0.4541639606867490 ;
            aWeights(  4 ) = 0.0427312318657758 ;
            aWeights(  5 ) = 0.2142003609268616 ;
            aWeights(  6 ) = 0.0427312318657758 ;
            aWeights(  7 ) = 0.2142003609268616 ;
            aWeights(  8 ) = 0.0427312318657758 ;
            aWeights(  9 ) = 0.2142003609268616 ;
            aWeights( 10 ) = 0.0427312318657758 ;
            aWeights( 11 ) = 0.2142003609268616 ;
            aWeights( 12 ) = 0.1444522232603068 ;
            aWeights( 13 ) = 0.1444522232603068 ;
            aWeights( 14 ) = 0.1444522232603068 ;
            aWeights( 15 ) = 0.1444522232603068 ;
            aWeights( 16 ) = 0.1444522232603068 ;
            aWeights( 17 ) = 0.1444522232603068 ;
            aWeights( 18 ) = 0.1444522232603068 ;
            aWeights( 19 ) = 0.1444522232603068 ;
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */
#endif //BELFEM_FN_INTPOINTS_GAUSS_QUAD20_HPP
