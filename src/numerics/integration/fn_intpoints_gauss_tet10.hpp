//
// Created by Christian Messe on 15.07.20.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET10_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET10_HPP


#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
// ----------------------------------------------------------------------------
namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        // 4th order
        // source Shunn and Ham, 10.1016/j.cam.2012.03.032
        inline void
        gauss_tet10(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {

            aPoints.set_size( 4, 10 );

            aPoints( 0 , 0 ) =  0.7784952948213300 ;
            aPoints( 1 , 0 ) =  0.0738349017262234 ;
            aPoints( 2 , 0 ) =  0.0738349017262234 ;

            aPoints( 0 , 1 ) =  0.0738349017262234 ;
            aPoints( 1 , 1 ) =  0.7784952948213300 ;
            aPoints( 2 , 1 ) =  0.0738349017262234 ;

            aPoints( 0 , 2 ) =  0.0738349017262234 ;
            aPoints( 1 , 2 ) =  0.0738349017262234 ;
            aPoints( 2 , 2 ) =  0.7784952948213300 ;

            aPoints( 0 , 3 ) =  0.0738349017262234 ;
            aPoints( 1 , 3 ) =  0.0738349017262234 ;
            aPoints( 2 , 3 ) =  0.0738349017262234 ;

            aPoints( 0 , 4 ) =  0.4062443438840510 ;
            aPoints( 1 , 4 ) =  0.4062443438840510 ;
            aPoints( 2 , 4 ) =  0.0937556561159491 ;

            aPoints( 0 , 5 ) =  0.4062443438840510 ;
            aPoints( 1 , 5 ) =  0.0937556561159491 ;
            aPoints( 2 , 5 ) =  0.4062443438840510 ;

            aPoints( 0 , 6 ) =  0.4062443438840510 ;
            aPoints( 1 , 6 ) =  0.0937556561159491 ;
            aPoints( 2 , 6 ) =  0.0937556561159491 ;

            aPoints( 0 , 7 ) =  0.0937556561159491 ;
            aPoints( 1 , 7 ) =  0.4062443438840510 ;
            aPoints( 2 , 7 ) =  0.4062443438840510 ;

            aPoints( 0 , 8 ) =  0.0937556561159491 ;
            aPoints( 1 , 8 ) =  0.4062443438840510 ;
            aPoints( 2 , 8 ) =  0.0937556561159491 ;

            aPoints( 0 , 9 ) =  0.0937556561159491 ;
            aPoints( 1 , 9 ) =  0.0937556561159491 ;
            aPoints( 2 , 9 ) =  0.4062443438840510 ;

            for( uint k=0; k<10; ++k )
            {
                aPoints( 3, k ) =  1.0
                                   - aPoints( 0, k )
                                   - aPoints( 1, k )
                                   - aPoints( 2, k );
            }

            aWeights.set_size( 10 );

            aWeights( 0 ) =  0.778495294821330 ;
            aWeights( 1 ) =  0.073834901726223 ;
            aWeights( 2 ) =  0.073834901726223 ;
            aWeights( 3 ) =  0.073834901726223 ;
            aWeights( 4 ) =  0.406244343884051 ;
            aWeights( 5 ) =  0.406244343884051 ;
            aWeights( 6 ) =  0.406244343884051 ;
            aWeights( 7 ) =  0.0937556561159491 ;
            aWeights( 8 ) =  0.0937556561159491 ;
            aWeights( 9 ) =  0.0937556561159491 ;
        }

// ----------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_INTPOINTS_GAUSS_TET10_HPP
