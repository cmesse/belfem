//
// Created by Christian Messe on 11.05.20.
//
#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "cl_TableGas.hpp"
#include "cl_GM_EoS_TableGas.hpp"
#include "cl_Atmosphere.hpp"
#include "cl_BL_State.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_BL_StagnationPointOld.hpp"

#include "fn_sum.hpp"

using namespace belfem;
using namespace boundarylayer;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    TableGas tAir ;

    Atmosphere tAtmosphere( AtmosphereType::ISA1976 ) ;

    real tAltitude = 40e3 ;
    real tMa = 12.0 ;

    // Freestream condition
    real tT0 ;
    real tP0 ;
    real tU0 ;
    tAtmosphere.compute_T_and_P( tAltitude, tT0, tP0 );

    tU0 = tMa * tAir.c( tT0, tP0 );

    State tFreestream( tAir );
    State tStagnation( tAir );
    tFreestream.compute( tT0, tP0, tU0 );

    std::cout << "T0 " << tT0 << std::endl;
    std::cout << "p0" << tP0 << std::endl;
    std::cout << "u0 " << tU0 << std::endl;
    std::cout << "---" << std::endl;


    // StagnationPoint tStag( tAir, tFreestream, tStagnation, 0.5 );

    // tStag.compute( 600.0 );

    /*real tT1 ;
    real tP1 ;
    real tU1 ;

    tAir.shock( tT0, tP0, tU0, tT1, tP1, tU1 );

    std::cout << "T1 " << tT1 << std::endl;
    std::cout << "p1 " << tP1 << std::endl;
    std::cout << "u1 " << tU1 << std::endl;
    std::cout << "---" << std::endl;

    State tStagnation( tAir );
    real tTs;
    real tPs;


    tAir.total( tT1, tP1, tU1, tTs, tPs );

    std::cout << "Ts " << tTs << std::endl;
    std::cout << "ps " << tPs << std::endl;

    tStagnation.compute( tTs, tPs, 0.0 );

    real tHs = tStagnation.h() ; */



    return gComm.finalize();
}