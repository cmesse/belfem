//
// Created by Christian Messe on 30.12.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_OneDMapper.hpp"
#include "fn_linspace.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    Vector< real > tTargetNodes ;
    linspace( 0., 2 * constant::pi, 21 , tTargetNodes );
    Vector< real > tTargetValues ;

    Vector< real > tSourceNodes ;
    linspace( 0., 2 * constant::pi, 41 , tSourceNodes );

    Vector< real > tSourceValues( tSourceNodes.length() );

    for( index_t k=0; k<tSourceNodes.length(); ++k )
    {
        tSourceValues( k ) = std::sin( tSourceNodes( k ) );
    }

    OneDMapper tMapper( tTargetNodes, 2 );

    tMapper.project( tSourceNodes, tSourceValues, tTargetValues );


    tSourceNodes.print("X0");
    tSourceValues.print("Y0");


    tTargetNodes.print("X1");
    tTargetValues.print("Y1");

    return  gComm.finalize();
}