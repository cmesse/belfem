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
#include "cl_Timer.hpp"

#include "GT_globals.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "cl_Vector.hpp"
#include "fn_linspace.hpp"
#include "cl_GT_RefGasFactory.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_GT_Arguments.hpp"
#include "cl_Matrix.hpp"
#include "banner.hpp"

using namespace belfem;
using namespace gastables;

Communicator gComm;
Logger       gLog( 3 );

//------------------------------------------------------------------------------

void
print_usage()
{
    std::cout << "Usage:" << std::endl;
    std::cout << "gastable --gas <Gasname>" << std::endl;
    std::cout << "type gastable -h for more options." << std::endl;
}

//------------------------------------------------------------------------------

void
print_help()
{
    std::cout << "gastable" << std::endl;
    std::cout << "  --gas      -g       specify the gas to be printed" << std::endl;
    std::cout << "  --Tmin     -a       minimum temperature in table      ( default: 200  )" << std::endl;
    std::cout << "  --Tmax     -b       maximum temperature in table      ( default: 2000 )" << std::endl;
    std::cout << "  --deltaT   -d       temperature step                  ( default: 50 )" << std::endl;
    std::cout << "  --molar    -m       print heats in J/Mol instead J/kg" << std::endl;
    std::cout << "  --version  -V       print the version banner" << std::endl;
}

//------------------------------------------------------------------------------

void
calculate_values( const gastables::Arguments & aArgs, Matrix< real > & aValues )
{
    // the polynomials divide by T, so the table must start above zero
    real tTmin = std::max( aArgs.T_min(), gDeltaT );
    real tTmax = std::min( aArgs.T_max(), gTmax );

    BELFEM_ERROR( aArgs.delta_T() > 0.0, "--deltaT must be positive" );
    BELFEM_ERROR( tTmax > tTmin, "--Tmax must be larger than --Tmin" );

    uint tNumberOfSamples= ( tTmax - tTmin ) / aArgs.delta_T() + 1;

    aValues.set_size( tNumberOfSamples, 6 );

    RefGasFactory tFactory;
    RefGas * tGas = tFactory.create_refgas( aArgs.gasname() );
    // tGas->set_mode( RefGasMode::POLY );

    std::cout << "gas " << tGas->label() << std::endl;


    if ( aArgs.molar() )
    {
        real tT = tTmin;
        std::cout << "Href " << tGas->H_ref() << std::endl;
        for( uint k=0; k<tNumberOfSamples; ++k )
        {
            aValues( k, 0 ) = tT;
            aValues( k, 1 ) = tGas->Cp( tT );
            aValues( k, 2 ) = tGas->H( tT ) * 0.001;
            aValues( k, 3 ) = tGas->S( tT );
            aValues( k, 4 ) = tGas->mu( tT ) * 1e6;
            aValues( k, 5 ) = tGas->lambda( tT ) * 1e3;
            tT += aArgs.delta_T();
        }
    }
    else
    {
        real tT = tTmin;
        std::cout << "href " << tGas->h_ref() << std::endl;
        for( uint k=0; k<tNumberOfSamples; ++k )
        {
            aValues( k, 0 ) = tT;
            aValues( k, 1 ) = tGas->cp( tT );
            aValues( k, 2 ) = tGas->h( tT ) * 0.001;
            aValues( k, 3 ) = tGas->s( tT );
            aValues( k, 4 ) = tGas->mu( tT ) * 1e6;
            aValues( k, 5 ) = tGas->lambda( tT ) * 1e3;
            tT += aArgs.delta_T();
        }
    }
    delete tGas;
}

//------------------------------------------------------------------------------

void
print_values( const gastables::Arguments & aArgs, Matrix< real > & aValues )
{
    uint tN = aValues.n_rows();
    if ( aArgs.molar() )
    {
        printf( "      T         Cp          H          S         mu     lambda\n");
        printf( "      K   J/(mol*K)     kJ/mol  J/(mol*K)     µPa*s   mW/(m*K)\n");
    }
    else
    {
        printf( "      T         cp          h          s         mu     lambda\n");
        printf( "      K    J/(kg*K)      kJ/kg   J/(kg*K)     µPa*s   mW/(m*K)\n");
    }

    for( uint k=0; k<tN; ++k )
    {
        printf( " %6.0f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
        //printf( " %6.0f %13.6f %13.6f %13.6f %13.6f %13.6f\n",
                aValues( k, 0 ),
                aValues( k, 1 ),
                aValues( k, 2 ),
                aValues( k, 3 ),
                aValues( k, 4 ),
                aValues( k, 5 ) );
    }
}

//------------------------------------------------------------------------------
int
main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // create argument list
    gastables::Arguments tArgs(  argc, argv );

    // get state
    auto tState = tArgs.state();

    switch ( tState )
    {
        case  State::PrintHelp :
        {
            print_help();
            break;
        }
        case State::PrintBanner :
        {
            print_banner();
            break ;
        }
        case State::PrintTable :
        {
            Matrix< real > tValues;
            calculate_values( tArgs, tValues );
            print_values( tArgs, tValues );
            break;
        }
        default:
        {
            print_usage();
            break;
        }
    }
    return 0;
}