//
// Created by Christian Messe on 02.05.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "typedefs.hpp"
#include "cl_BS_TMatrix.hpp"
#include "cl_BS_Mapper.hpp"

#include "cl_Gas.hpp"

#include "cl_Progressbar.hpp"
#include "banner.hpp"
#include "cl_Timer.hpp"

#include "cl_TableGas.hpp"
#include "cl_GM_EoS_TableGas.hpp"
#include "cl_Atmosphere.hpp"
#include "cl_GT_RefGas.hpp"
using namespace belfem;
using namespace bspline;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    Gas tCombgas( {
                    "H2",
                    "H",
                    "O",
                    "O2",
                    "OH",
                    "H2O",
                    "HO2",
                    "H2O2"},
                  {
                          0.75,
                          0.0,
                          0.0,
                          0.25,
                          0.0,
                          0.0,
                          0.0,
                          0.0
                  } );

    const Vector< real > & tY = tCombgas.mass_fractions();

    Vector< real > tY0( tCombgas.mass_fractions() );
    Vector< real > tDeltaY;

    tCombgas.h( 200.0, 20e5 );

    const Vector< real > & tX = tCombgas.molar_fractions();

    tY0 = tY ;
    real tT = 200.0;
    real tP = 10e5;
    Cell< gastables::RefGas * > & tComponents = tCombgas.components();

    real tDeltaT = 1.0 ;
    uint tCount = 0;
    while( std::abs( tDeltaT ) > 0.01 )
    {
        tCombgas.remix_to_equilibrium( tT, tP );

        tDeltaY = tY - tY0;

        tDeltaT = 0.0;
        for ( uint k = 0; k < tCombgas.number_of_components(); ++k )
        {
            tDeltaT += tCombgas.h( k, tT, tP ) * tDeltaY( k );
        }
        tDeltaT /= -tCombgas.cp( tT, tP );

        tT += tDeltaT;

        std::cout << tCount << " " << tT << " " << tCombgas.gamma( tT, tP ) << std::endl;
    }

    for ( uint k = 0; k < tCombgas.number_of_components(); ++k )
    {
        std::cout << tComponents( k )->label() << " " << tX( k ) << " " << tY( k ) << std::endl;
    }

    return gComm.finalize();
}

