//
// Created by Christian Messe on 23.12.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Tensor.hpp"

#include "fn_fiber_polarization.hpp"
#include "fn_inv.hpp"
#include "fn_identity.hpp"
#include "fn_invert_symmetric.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

/**
 * This is a prototype to demonstrate the Mori-Tanaka
 * and the iterative Self-Consistent Method.
 *
 * To be tidied up and wrapped in functions.
 */
int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    // fiber data ( unit does not matter )
    real E1   = 235 ;
    real E2   = 20 ;
    real G12  = 18 ;
    real nu12 = 0.013 ;
    real nu23 = 0.25 ;

    real G23 = E2 / ( 2.0 * ( 1.0 + nu23 ) );

    // matrix data
    real Em = 3.5 ;
    real num = 0.35 ;

    // fiber volume fraction
    real phi = 0.6 ;

    // mixed elasticity tensor
    Tensor< real > tL0( 3, 3, 3, 3 );
    Tensor< real > tL( 3, 3, 3, 3 );

    // elasticity tensor for fibers
    Tensor< real > tLf( 3, 3, 3, 3 );
    Matrix< real > tCf( 6, 6 );

    // elasticity tensor for matrix
    Tensor< real > tLm( 3, 3, 3, 3 );
    tLf.fill_orthotropic_elasticity( E1, E2, E2, nu23, nu12, nu12, G23, G12, G12 );
    tLf.to_matrix( tCf );

    tLm.fill_isotropic_elasticity( Em, num );

    // Fiber polarization tensor
    Tensor< real > tP( 3, 3, 3, 3 );
    fiber_polarization( tCf, tP );

    // Unity tensor
    Tensor< real > tI( 3, 3, 3, 3 );
    tensor::identity( tI );

    // Help tensors
    Tensor< real > tA( 3, 3, 3, 3 );
    Tensor< real > tB( 3, 3, 3, 3 );
    Tensor< real > tTf( 3, 3, 3, 3 );
    Tensor< real > tTm( 3, 3, 3, 3 );
    Vector< real > tWork( 72 );
    Vector< int > tPivot( 36 );

    // Mori-Tanaka Method:

    // difference between fiber and matrix
    tB = tLf - tLm ;

    // compute infuluence tensor for fibers
    tP.ddot( tB , tTf );
    tTf += tI ;
    tensor::invert_symmetric( tTf, tWork, tPivot );

    tB.ddot( tTf, tL );
    tL *= phi ;
    tL += tLm ;

    Matrix< real > tC( 6, 6 );
    tL.to_matrix( tC );
    Matrix< real > tS = inv( tC );

    E1 = 1.0 / tS( 0, 0 );
    E2 = 1.0 / tS( 1, 1 );
    G12 = 1.0 / tS( 5, 5 );
    G23 = 1.0 / tS( 3, 3 );
    nu12 = -tS(1,0 ) * E1 ;


    std::cout << "mori-tanaka " << std::endl ;
    std::cout << "E1   : " << E1 << std::endl ;
    std::cout << "E2   : " << E2 << std::endl ;
    std::cout << "nu12 : " << nu12 << std::endl ;
    std::cout << "G12 : " << G12 << std::endl ;
    std::cout << "G23 : " << G23 << std::endl ;

    // self-consistent mixing
    for( uint k=0; k<10; ++k )
    {
        // shift elasticity
        tL0 = tL ;

        // influence tensor for fibers
        tB = tLf - tL0;
        tP.ddot( tB , tTf );
        tTf += tI ;
        tensor::invert_symmetric( tTf, tWork, tPivot );

        tB.ddot( tTf, tL ) ;
        tL *= phi ;

        tTm = tTf ;
        tTm *= - phi ;
        tTm += tI ;
        //tTm /= 1.0 - phi ;
        tB = tLm - tL0 ;
        tB.ddot( tTm, tA );
        //tA *= 1.0 - phi ;
        tL += tA ;
        tL += tL0 ;
        tL0 *= 0.1 ;
        tL *= 0.9 ;
        tL += tL0 ;
        //std::cout << "iteration " << k << " " << tL( 0, 0, 0, 0 ) << std::endl;
    }

    tL.to_matrix( tC );
    tS = inv( tC );

    E1 = 1.0 / tS( 0, 0 );
    E2 = 1.0 / tS( 1, 1 );
    G12 = 1.0 / tS( 5, 5 );
    G23 = 1.0 / tS( 3, 3 );
    nu12 = -tS(1,0 ) * E1 ;

    std::cout << "self consistent: " << std::endl ;
    std::cout << "E1   : " << E1 << std::endl ;
    std::cout << "E2   : " << E2 << std::endl ;
    std::cout << "nu12 : " << nu12 << std::endl ;
    std::cout << "G12 : " << G12 << std::endl ;
    std::cout << "G23 : " << G23 << std::endl ;

    // close communicator
    return gComm.finalize();
}