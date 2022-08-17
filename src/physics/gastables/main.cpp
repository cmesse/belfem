//
// Created by Christian Messe on 2019-08-18.
//

// thermostructural analysis of reusable spacecraft

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
    std::cout << "  --Tmin     -a       minumum temperature in table      ( default: 200  )" << std::endl;
    std::cout << "  --Tmax     -b       maximum temperature in table      ( default: 2000 )" << std::endl;
    std::cout << "  --deltaT   -d       temperature step                  ( default: 50 )" << std::endl;
    std::cout << "  --molar    -m       print heats in J/Mol instead J/kg" << std::endl;
}

//------------------------------------------------------------------------------

void
calculate_values( const gastables::Arguments & aArgs, Matrix< real > & aValues )
{
    real tTmin = std::max( aArgs.T_min(), 0.0 );
    real tTmax = std::min( aArgs.T_max(), gTmax );

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
            aValues( k, 2 ) = ( tGas->H( tT ) -tGas->H_ref() ) * 0.001;
            aValues( k, 3 ) = tGas->S( tT );
            aValues( k, 4 ) = tGas->mu( tT ) * 1e6;
            aValues( k, 5 ) = tGas->lambda( tT ) * 1e3;
            tT += aArgs.delta_T();
        }
    }
    else
    {
        real tT = tTmin;
        std::cout << "href " << tGas->H_ref() << std::endl;
        for( uint k=0; k<tNumberOfSamples; ++k )
        {
            aValues( k, 0 ) = tT;
            aValues( k, 1 ) = tGas->cp( tT );
            aValues( k, 2 ) = ( tGas->h( tT ) - tGas->h_ref() ) * 0.001;
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
        printf( "      K   J/(mol*K) kJ/(mol*K) J/(mol*K)        µPa   mW/(m*K)\n");
    }
    else
    {
        printf( "      T         cp          h          s         mu     lambda\n");
        printf( "      K    J/(kg*K)  kJ/(kg*K)   J/(kg*K)       µPa   mW/(m*K)\n");
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
    gComm = Communicator( argc, argv );

    // create argument list
    gastables::Arguments tArgs(  argc, argv );

    // get state
    auto tState = tArgs.state();

    switch ( tState )
    {
        case ( State::PrintHelp ):
        {
            print_help();
            break;
        }
        case( State::PrintTable ):
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