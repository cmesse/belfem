//
// Created by christian on 7/27/22.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "cl_Spline.hpp"
#include "banner.hpp"
#include "cl_Timer.hpp"
#include "stringtools.hpp"
#include "cl_MaterialFactory.hpp"
#include "fn_linspace.hpp"
#include "cl_Progressbar.hpp"

#include "cl_HDF5.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    print_banner();


    MaterialFactory tFactory ;


    string tMatLabel = "Silver" ;

    const string tDatabase = "materials.hdf5" ;

    bool tComputeCp = true ;
    bool tComputeK = true ;
    bool tComputeRho = true ;

    string tCpLabel = sprint( "%sCp", tMatLabel.c_str() ) ;
    string tKLabel = sprint( "%sK", tMatLabel.c_str() ) ;
    string tRhoLabel = sprint( "%sRho", tMatLabel.c_str() ) ;
    string tDensityLabel = sprint( "%sDensity", tMatLabel.c_str() ) ;

    // create the material
    Material * tMaterial = tFactory.create_material( tMatLabel );

    // the temperature list
    uint tNumSamples = 251 ;
    Vector< real > tT = linspace( 0., 500., tNumSamples );
    Vector< real > tValues( tNumSamples );

    // help matrix
    SpMatrix tH ;
    spline::create_helpmatrix( tT.length(), tT( 1 ), tH );

    // populate data
    if ( tComputeCp )
    {
        for( uint k=0; k<tNumSamples; ++k )
        {
            tValues( k ) = tMaterial->c( tT( k ) );
        }

        Spline tSpline( tT, tValues, tH );
        tSpline.save_to_database( tDatabase, tCpLabel );
    }

    if ( tComputeK )
    {
        for( uint k=0; k<tNumSamples; ++k )
        {
            tValues( k )  = tMaterial->lambda( tT( k ) );
            //tK( k )
            //tRho( k ) = tMaterial->rho_el( 0, tT( k ), 0 );
        }

        Spline tSpline( tT, tValues, tH );
        tSpline.save_to_database( tDatabase, tKLabel );
    }


    if ( tComputeRho )
    {
        for( uint k=0; k<tNumSamples; ++k )
        {
            tValues( k )  = tMaterial->rho_el( 0, tT( k ), 0 );
        }

        Spline tSpline( tT, tValues, tH );
        tSpline.save_to_database( tDatabase, tRhoLabel );
    }

    // also save density
    HDF5 tFile( tDatabase, FileMode::OPEN_RDWR );
    tFile.save_data( tDensityLabel, tMaterial->rho() );
    tFile.close() ;
    ./ma
    delete tMaterial ;

    return gComm.finalize();
}