//
// Created by Christian Messe on 30.04.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "typedefs.hpp"
#include "cl_BS_TMatrix.hpp"
#include "cl_BS_Mapper.hpp"

#define private public
#include "cl_Gas.hpp"
#undef private

#include "cl_Progressbar.hpp"
#include "banner.hpp"
#include "cl_Timer.hpp"

#include "cl_TableGas.hpp"
#include "cl_GM_EoS_TableGas.hpp"
#include "cl_Atmosphere.hpp"

using namespace belfem;
using namespace bspline;

Communicator gComm;
Logger       gLog( 5 );


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    // create a lookup table
    //LookupTable tTable("hotair.hdf5" );
    Atmosphere tAtmos( AtmosphereType::ISA1976 );

    TableGas tAir("hotair.hdf5");

    std::cout << "done" << std::endl;
    Gas tColdAir;

    real tAltitude = 40e3;
    real tAoA = 25.0 * constant::deg;
    real tMa = 20.0;

    // test phsical properties
    real tT ;
    real tP ;
    tAtmos.compute_T_and_P( tAltitude, tT, tP );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - -
    // General Gas Properties
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - -
    std::cout << "M      " << tAir.M( tT, tP ) * 1000 << " " << tColdAir.M( tT, tP ) * 1000 << std::endl ;
    std::cout << "R      " << tAir.R( tT, tP ) << " " <<  tColdAir.R( tT, tP ) << std::endl ;
    std::cout << std::endl;
    std::cout << "h      " << tAir.h( tT, tP ) * 0.001 << " " << tColdAir.h( tT, tP ) * 0.001 << std::endl ;
    std::cout << "s      " << tAir.s( tT, tP ) << " " << tColdAir.s( tT, tP ) << std::endl ;
    std::cout << std::endl;
    std::cout << "cp     " << tAir.cp( tT, tP ) << " " << tColdAir.cp( tT, tP ) << std::endl ;
    std::cout << "cv     " << tAir.cv( tT, tP ) << " " << tColdAir.cv( tT, tP ) << std::endl ;
    std::cout << "gamma  " << tAir.gamma( tT, tP ) << " " << tColdAir.gamma( tT, tP ) << std::endl ;
    std::cout << "c      " << tAir.c( tT, tP ) << " " << tColdAir.c( tT, tP ) << std::endl ;
    std::cout << std::endl;
    std::cout << "dsdT   " << tAir.dsdT( tT, tP ) << " " << tColdAir.dsdT( tT, tP ) << std::endl ;
    std::cout << "dsdp   " << tAir.dsdp( tT, tP ) << " " << tColdAir.dsdp( tT, tP ) << std::endl ;
    std::cout << "dcpdT  " << tAir.dcpdT( tT, tP ) << " " << tColdAir.dcpdT( tT, tP ) << std::endl ;
    std::cout << std::endl;
    std::cout << "mu     " << tAir.mu( tT, tP ) * 1e6 << " " << tColdAir.mu( tT, tP ) * 1e6 << std::endl ;
    std::cout << "lambda " << tAir.lambda( tT, tP ) * 1e3 << " " << tColdAir.lambda( tT, tP ) * 1e3 << std::endl ;
    std::cout << "Pr     " << tAir.Pr( tT, tP ) << " " << tColdAir.Pr( tT, tP ) << std::endl ;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - -
    // Shocks and Expansions
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - -
    std::cout << std::endl;

    // test with cold air
    real tU0c = tColdAir.c( tT, tP ) * tMa ;
    real tTsc ;
    real tPsc ;
    real tUsc ;
    tColdAir.shock( tT, tP, tU0c, tTsc, tPsc, tUsc );

    // hot air
    real tU0 = tAir.c( tT, tP ) * tMa ;
    real tTs;
    real tPs;
    real tUs ;

    tAir.shock( tT, tP, tU0, tTs, tPs, tUs );

    std::cout << "perp. shock " << tTs << " " << tPs*1.0e-5 << " " << tUs/tAir.c( tTs, tPs ) << " | " << tTsc << " " << tPsc *1.0e-5 << " " << tUsc/tAir.c( tTsc, tPsc ) << std::endl;

    real tT1c ;
    real tP1c ;
    real tU1c ;
    real tBetac;

    tColdAir.shock( tT, tP, tU0c, tAoA, tT1c, tP1c, tU1c, tBetac );

    real tT1 ;
    real tP1 ;
    real tU1 ;
    real tBeta;
    tAir.shock( tT, tP, tU0c, tAoA, tT1, tP1, tU1, tBeta );
    std::cout << "shock " << tT1 << " " << tP1*1.0e-5 << " " << tU1/tAir.c( tT1, tP1 ) << " | " << tT1c << " " << tP1c *1.0e-5 << " " << tU1c/tColdAir.c( tT1c, tP1c ) << std::endl;

    real tT2c ;
    real tP2c ;
    real tU2c ;

    tColdAir.prandtl_meyer( tT, tP, tU0, -tAoA, tT2c, tP2c, tU2c );

    real tT2 ;
    real tP2 ;
    real tU2 ;
    tColdAir.prandtl_meyer( tT, tP, tU0, -tAoA, tT2, tP2, tU2 );
    std::cout << "expand " << tT2 << " " << tP2*1.0e-5 << " " << tU2 / tAir.c( tT2, tP2 ) << " | " << tT2c << " " << tP2c *1.0e-5 << " " << tU2c / tAir.c( tT2c, tP2c ) << std::endl;


    // std::cout <<
    return gComm.finalize();
}