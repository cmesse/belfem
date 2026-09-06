/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */


#include <iostream>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "fn_linspace.hpp"

#include "cl_Gas.hpp"
#include "cl_GT_Arguments.hpp"
#include "cl_HDF5.hpp"
#include "banner.hpp"

using namespace belfem;
using namespace gastables;
using namespace gasmodels;

Communicator gComm;
Logger       gLog( 3 );

//------------------------------------------------------------------------------

void
print_help()
{
    std::cout << "gas -- evaluate the Gas class over a temperature sweep at 1 bar" << std::endl;
    std::cout << std::endl;
    std::cout << "  --gas      -g       pure gas to evaluate ( default: the air mixture )" << std::endl;
    std::cout << "  --Tmin     -a       minimum temperature in table      ( default: 200  )" << std::endl;
    std::cout << "  --Tmax     -b       maximum temperature in table      ( default: 2000 )" << std::endl;
    std::cout << "  --deltaT   -d       temperature step                  ( default: 50 )" << std::endl;
    std::cout << "  --molar    -m       print heats in J/mol instead of J/kg" << std::endl;
    std::cout << "  --version  -V       print the version banner" << std::endl;
    std::cout << std::endl;
    std::cout << "Without arguments, the default air mixture is swept from 100 to 2500 K." << std::endl;
    std::cout << "The table is also written to gasdata.h5." << std::endl;
}

//------------------------------------------------------------------------------

void
print_header( const bool aMolarFlag )
{
    if ( aMolarFlag )
    {
        printf( "%12s %12s %12s %12s %12s %12s\n",
                "T", "Cp", "H", "S", "mu", "lambda" );
        printf( "%12s %12s %12s %12s %12s %12s\n",
                "K", "J/(mol*K)", "kJ/mol", "J/(mol*K)", "muPa*s", "mW/(m*K)" );
    }
    else
    {
        printf( "%12s %12s %12s %12s %12s %12s\n",
                "T", "cp", "h", "s", "mu", "lambda" );
        printf( "%12s %12s %12s %12s %12s %12s\n",
                "K", "J/(kg*K)", "kJ/kg", "J/(kg*K)", "muPa*s", "mW/(m*K)" );
    }
}

//------------------------------------------------------------------------------

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

//------------------------------------------------------------------------------

    // create argument list, shared with the gastable executable
    gastables::Arguments tArgs( argc, argv );

    if ( tArgs.state() == State::PrintHelp )
    {
        print_help();
        return gComm.finalize();
    }
    if ( tArgs.state() == State::PrintBanner )
    {
        print_banner();
        return gComm.finalize();
    }
    // without arguments the historical demo behavior is kept:
    // the default air mixture over its own sweep
    const bool tIsAir = tArgs.state() != State::PrintTable ;

    Vector< real > tT;

    if ( tIsAir )
    {
        linspace( 100.0, 2500.0, 25, tT );
    }
    else
    {
        uint tNumSteps = ( tArgs.T_max() - tArgs.T_min() ) / tArgs.delta_T() + 1;
        linspace( tArgs.T_min(), tArgs.T_max(), tNumSteps, tT );
    }

    const uint tN = tT.length();

    // create the gas
    Gas * tGas = tIsAir ? new Gas() : new Gas( tArgs.gasname() );

    std::cout << "gas " << ( tIsAir ? "air" : tArgs.gasname() ) << std::endl;

    real tP = 1.0e5;

    // molar mass, for the molar output option
    const real tM = tGas->M( 298.15, tP );
    const real tScale = tArgs.molar() ? tM : 1.0 ;

    // Matrix with values
    Matrix< real > tData( tN, 6 );
    Vector< real > tRefpoint( 6 );

    print_header( tArgs.molar() );

    for( uint k=0; k<tN; ++k )
    {
        printf( "%12.3f %12.6f %12.6f %12.6f %12.6f %12.6f ;\n",
                tT( k ),
                tGas->cp( tT( k ) , tP  ) * tScale,
                tGas->h( tT( k ) , tP ) * tScale * 1e-3,
                tGas->s( tT( k ) , tP ) * tScale,
                tGas->mu( tT( k ) , tP ) * 1e6,
                tGas->lambda( tT( k ) , tP ) * 1e3 );

        tData( k, 0 ) = tT( k );
        tData( k, 1 ) = tGas->cp( tT( k ) , tP  );
        tData( k, 2 ) = tGas->h( tT( k ) , tP ) ;
        tData( k, 3 ) = tGas->s( tT( k ) , tP ) ;
        tData( k, 4 ) = tGas->mu( tT( k ) , tP ) ;
        tData( k, 5 ) = tGas->lambda( tT( k ) , tP ) ;

    }

    HDF5 tFile( "gasdata.h5", FileMode::NEW );
    tFile.save_data( "Table", tData );

    real tTref = 298.15 ;

    tRefpoint( 0 ) = tTref ;
    tRefpoint( 1 ) = tGas->cp(tTref , tP  );
    tRefpoint( 2 ) = tGas->h( tTref , tP ) ;
    tRefpoint( 3 ) = tGas->s( tTref , tP ) ;
    tRefpoint( 4 ) = tGas->mu( tTref , tP ) ;
    tRefpoint( 5 ) = tGas->lambda( tTref , tP ) ;

    tFile.save_data( "Refpoint", tRefpoint );

    std::cout << tTref
              << " " << tGas->cp( tTref , tP  )
              << " " << tGas->h( tTref , tP ) * 1e-3
              << " " << tGas->s( tTref, tP )
              << " " << tGas->mu( tTref , tP ) * 1e6
              << " " << tGas->lambda( tTref , tP ) * 1e3 << std::endl;

    printf( "R: %12.5f\n", tGas->R( tTref, tP ) );

    delete tGas;

    /* close before the communicator goes down: PetscFinalize() closes the
     * HDF5 library when PETSc is built with HDF5 support, and the file
     * object's destructor would then close an already dead handle */
    tFile.close();

//------------------------------------------------------------------------------

    return gComm.finalize();

}
