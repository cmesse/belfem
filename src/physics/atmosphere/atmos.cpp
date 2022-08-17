//
// Created by Christian Messe on 26.12.19.
//

#include <iostream>

#include "typedefs.hpp"

#include "cl_Logger.hpp"
#include "cl_Communicator.hpp"

#include "cl_Vector.hpp"
#include "en_AtmosphereType.hpp"
#include "cl_Atmosphere.hpp"
#include "cl_Gas.hpp"
#include "fn_linspace.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    // create atmosphere
    Atmosphere Atmos(AtmosphereType::ISA1976 );

    index_t tNumSamples = 200;

    // heigts
    Vector< real > tHeights = linspace( 0.0, 170.0e3, tNumSamples ) ;

    real tT0;
    real tP0;
    Atmos.compute_T_and_P( 0.0, tT0, tP0 );

    real tT;
    real tP;

    for( real tH : tHeights )
    {
        Atmos.compute_T_and_P( tH, tT, tP );
        std::fprintf( stdout, "%8.3f %12.8f %12.8f\n", tH * 0.001, tT, tP );
    }
    return gComm.finalize();
}