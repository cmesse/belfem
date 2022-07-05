//
// Created by Christian Messe on 18.12.20.
//

#include <gtest/gtest.h>
#include "typedefs.hpp"

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Tensor.hpp"

#include "fn_inv.hpp"
#include "fn_invert_symmetric.hpp"
#include "fn_TR_mat_to_ten.hpp"
#include "fn_identity.hpp"
#include "fn_norm.hpp"
using namespace belfem ;

TEST( TENSOR, invert_symmetric )
{
    // example data for C/C-SiC - like

    real E1 = 60. ;
    real E2 = 65. ;
    real E3 = 20. ;
    real G23 = 8. ;
    real G13 = 8.5 ;
    real G12 = 9. ;
    real nu23 = 0.2 ;
    real nu13 = 0.25 ;
    real nu12 = 0.1 ;


    // create compliance matrix
    Matrix< real > tS( 6, 6, 0.0 );

    tS( 0, 0 ) = 1. / E1 ;
    tS( 1, 0 ) = -nu12/E1 ;
    tS( 2, 0 ) = -nu13 / E1 ;

    tS( 0, 1 ) =  -nu12/E1 ;
    tS( 1, 1 ) =  1./ E2 ;
    tS( 2, 1 ) = -nu23 / E2 ;

    tS( 0, 2 ) = -nu13 / E1 ;
    tS( 1, 2 ) = -nu23 / E2 ;
    tS( 2, 2 ) = 1. / E3 ;

    tS( 3, 3 ) = 1.0 / G23 ;
    tS( 4, 4 ) = 1.0 / G13 ;
    tS( 5, 5 ) = 1.0 / G12 ;

    // create elasticity matrix
    Matrix< real > tC( 6, 6, 0.0 );

    tC = inv( tS );

    // create a tensor from the matrix
    Tensor< real > tT( tC );

    // create an identity matrix
    Tensor< real > tI( 3, 3, 3, 3 );
    tensor::identity( tI );

    // copy the tensor
    Tensor< real > tM( 3, 3, 3, 3 );
    tM = tT ;

    // allocate work matrices
    Vector< real > tWork( 72 );

    // allokate work vector
    Vector< int > tPivot( 36 );

    // invert the tensor
    tensor::invert_symmetric( tM, tWork, tPivot );

    // create a test tensor
    Tensor< real > tJ( 3, 3, 3, 3 );

    // contract both tensors
    tJ = tT % tM ;

    // create a vector to check the error
    Vector < real > tError( 81 );

    real * tX = tI.data();
    real * tY = tJ.data();

    for( uint k=0; k<81; ++k )
    {
        tError( k ) = tX[ k ] - tY[ k ] ;
    }
    // hope this works
    EXPECT_NEAR( norm( tError ),  0.0,  1e-12 );
}