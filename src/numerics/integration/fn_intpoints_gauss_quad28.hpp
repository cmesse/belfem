//
// Created by Christian Messe on 18.07.20.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_QUAD28_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_QUAD28_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
//----------------------------------------------------------------------------

        // doi.org/10.1016/j.camwa.2015.03.017
        inline void
        gauss_quad28(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aPoints.set_size( 2, 28 ) ;
            aWeights.set_size( 28 ) ;

            aPoints( 0,  0 ) =  0.7146178296646060 ;
            aPoints( 1,  0 ) =  0.0000000000000000 ;

            aPoints( 0,  1 ) = -0.7146178296646060 ;
            aPoints( 1,  1 ) =  0.0000000000000000 ;

            aPoints( 0,  2 ) =  0.0000000000000000 ;
            aPoints( 1,  2 ) =  0.7146178296646060 ;

            aPoints( 0,  3 ) =  0.0000000000000000 ;
            aPoints( 1,  3 ) = -0.7146178296646060 ;

            aPoints( 0,  4 ) =  0.2736572101714596 ;
            aPoints( 1,  4 ) =  0.2736572101714596 ;

            aPoints( 0,  5 ) =  0.6366039322123010 ;
            aPoints( 1,  5 ) =  0.6366039322123010 ;

            aPoints( 0,  6 ) =  0.2736572101714596 ;
            aPoints( 1,  6 ) = -0.2736572101714596 ;

            aPoints( 0,  7 ) =  0.6366039322123010 ;
            aPoints( 1,  7 ) = -0.6366039322123010 ;

            aPoints( 0,  8 ) = -0.2736572101714596 ;
            aPoints( 1,  8 ) =  0.2736572101714596 ;

            aPoints( 0,  9 ) = -0.6366039322123010 ;
            aPoints( 1,  9 ) =  0.6366039322123010 ;

            aPoints( 0, 10 ) = -0.2736572101714596 ;
            aPoints( 1, 10 ) = -0.2736572101714596 ;

            aPoints( 0, 11 ) = -0.6366039322123010 ;
            aPoints( 1, 11 ) = -0.6366039322123010 ;

            aPoints( 0, 12 ) =  0.9516303887840335 ;
            aPoints( 1, 12 ) =  0.8155654336896384 ;

            aPoints( 0, 13 ) =  0.3462072000476454 ;
            aPoints( 1, 13 ) =  0.9355678714875911 ;

            aPoints( 0, 14 ) = -0.9516303887840335 ;
            aPoints( 1, 14 ) =  0.8155654336896384 ;

            aPoints( 0, 15 ) = -0.3462072000476454 ;
            aPoints( 1, 15 ) =  0.9355678714875911 ;

            aPoints( 0, 16 ) =  0.9516303887840335 ;
            aPoints( 1, 16 ) = -0.8155654336896384 ;

            aPoints( 0, 17 ) =  0.3462072000476454 ;
            aPoints( 1, 17 ) = -0.9355678714875911 ;

            aPoints( 0, 18 ) = -0.9516303887840335 ;
            aPoints( 1, 18 ) = -0.8155654336896384 ;

            aPoints( 0, 19 ) = -0.3462072000476454 ;
            aPoints( 1, 19 ) = -0.9355678714875911 ;

            aPoints( 0, 20 ) =  0.8155654336896384 ;
            aPoints( 1, 20 ) =  0.9516303887840335 ;

            aPoints( 0, 21 ) =  0.9355678714875911 ;
            aPoints( 1, 21 ) =  0.3462072000476454 ;

            aPoints( 0, 22 ) = -0.8155654336896384 ;
            aPoints( 1, 22 ) =  0.9516303887840335 ;

            aPoints( 0, 23 ) = -0.9355678714875911 ;
            aPoints( 1, 23 ) =  0.3462072000476454 ;

            aPoints( 0, 24 ) =  0.8155654336896384 ;
            aPoints( 1, 24 ) = -0.9516303887840335 ;

            aPoints( 0, 25 ) =  0.9355678714875911 ;
            aPoints( 1, 25 ) = -0.3462072000476454 ;

            aPoints( 0, 26 ) = -0.8155654336896384 ;
            aPoints( 1, 26 ) = -0.9516303887840335 ;

            aPoints( 0, 27 ) = -0.9355678714875911 ;
            aPoints( 1, 27 ) = -0.3462072000476454 ;

            aWeights(  0 ) = 0.2174004398687120 ;
            aWeights(  1 ) = 0.2174004398687120 ;
            aWeights(  2 ) = 0.2174004398687120 ;
            aWeights(  3 ) = 0.2174004398687120 ;
            aWeights(  4 ) = 0.2772741029838511 ;
            aWeights(  5 ) = 0.2139336378782481 ;
            aWeights(  6 ) = 0.2772741029838511 ;
            aWeights(  7 ) = 0.2139336378782481 ;
            aWeights(  8 ) = 0.2772741029838511 ;
            aWeights(  9 ) = 0.2139336378782481 ;
            aWeights( 10 ) = 0.2772741029838511 ;
            aWeights( 11 ) = 0.2139336378782481 ;
            aWeights( 12 ) = 0.0440745691149831 ;
            aWeights( 13 ) = 0.1016213405196113 ;
            aWeights( 14 ) = 0.0440745691149831 ;
            aWeights( 15 ) = 0.1016213405196113 ;
            aWeights( 16 ) = 0.0440745691149831 ;
            aWeights( 17 ) = 0.1016213405196113 ;
            aWeights( 18 ) = 0.0440745691149831 ;
            aWeights( 19 ) = 0.1016213405196113 ;
            aWeights( 20 ) = 0.0440745691149831 ;
            aWeights( 21 ) = 0.1016213405196113 ;
            aWeights( 22 ) = 0.0440745691149831 ;
            aWeights( 23 ) = 0.1016213405196113 ;
            aWeights( 24 ) = 0.0440745691149831 ;
            aWeights( 25 ) = 0.1016213405196113 ;
            aWeights( 26 ) = 0.0440745691149831 ;
            aWeights( 27 ) = 0.1016213405196113 ;
        }

//----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */

#endif //BELFEM_FN_INTPOINTS_GAUSS_QUAD28_HPP
