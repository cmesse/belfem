//
// Created by Christian Messe on 28.04.20.
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
#include "cl_GT_RefGas.hpp"
#include "stringtools.hpp"
using namespace belfem;
using namespace bspline;

Communicator gComm;
Logger       gLog( 5 );


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    print_banner();

    Gas tAir( {
                      "N2",
                      "O2",
                      "Ar",
                      "CO2",
                      "Ne",
                      "NO",
                      "CO",
                      "O3",
                      "N2O",
                      "NO2",
                      "e-",
                      "Ar+",
                      "C",
                      "C-",
                      "C+",
                      "CO+",
                      "CO2+",
                      "N",
                      "N-",
                      "N+",
                      "N2-",
                      "N2+",
                      "N2O+",
                      "Ne+",
                      "NO+",
                      "NO2-",
                      "O",
                      "O-",
                      "O+",
                      "O2-",
                      "O2+"
              },
              {
                      0.78084,
                      0.20942,
                      0.00934,
                      0.00038182,
                      0.00001818,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0
              }
    );


       Mapper tMapper( 2, 3,
                       { 646, 251 },// 646 x 251 - 64 x 35
                       { 100.0, -4000.0 },
                       { 13000.0, 6500.0 } );

       const Matrix< real > & tGrid = tMapper.integration_grid();

       // get size of grid
       uint tGridSize = tGrid.n_cols();

       Progressbar tProgress( tGridSize );

       // Molar Masses
       Vector< real > & tM = tMapper.create_field( "M" );

       // enthalpy
       Vector< real > & tH = tMapper.create_field( "h" );

       // entropy
       Vector< real > & tS = tMapper.create_field( "s" );

       // frozen specific heat
       // Vector< real > & tCp = tMapper.create_field( "cp" );

       // viscosity
       Vector< real > & tMu = tMapper.create_field( "mu" );

       // thermal conductivity ( frozen )
       Vector< real > & tLambda = tMapper.create_field( "lambda" );

       // Diffusion enthalpy
       Vector< real > & tHd = tMapper.create_field( "hd" );

       uint tFieldCount = tMapper.number_of_fields();

       // Volume fractions for air
       Vector< real > tX0( tAir.molar_fractions());
       Vector< real > tX( tAir.molar_fractions());

       std::cout << " Computing Equilibrium properties of Air" << std::endl;

       Timer tTimer;

       // print molar fractions for visualization in exodus
       bool tMolarFlag = false;



       Cell< gastables::RefGas * > tComp = tAir.components();

       if( tMolarFlag )
       {
           for ( uint k = 0; k < tComp.size(); ++k )
           {
               //string tLabel = search_and_replace( search_and_replace( tComp( k )->label(), "+", "plus"), "-", "minus");
               std::cout << tComp( k )->label() << " " << tComp( k )->M() * 1000 << std::endl;
               tMapper.create_field(  tComp( k )->label() );
           }
      }

       // enthalpies at zero K
       Vector< real > tHf( tAir.number_of_components() );

       // compute values
       for( index_t k=0; k<tGridSize; ++k )
       {

           tProgress.step( k );

           // compute temperature
           real tT = tGrid( 0, k ) ;

           // compute pressure
           real tP = std::pow( 10, tGrid( 1, k ) * 0.001 ) ;

           tX = tX0 ;
           if( tT >= 350.0 )
               tAir.compute_equilibrium( tT, tP, tX );

           tAir.remix_R( tX );
           tAir.remix_heat();

           // tAir.remix_critical_point() ;
           const Vector< real > & tY = tAir.mass_fractions() ;

           // save molar mass
           tM( k ) = tAir.M( tT, tP ) ;

           BELFEM_ERROR( std::abs( tM(k) ) > 0, "Error for T= %f, p= %f, M = %f",
                   (float ) tT, (float) tP, (float ) tM( k ) );

           // enthalpy
           tH( k ) = tAir.h( tT, tP );

           // entropy
           tS( k ) = tAir.s( tT, tP );

           // tCp( k ) = tAir.cp( tT, tP );

           tMu( k ) = tAir.cea_mu( tT );

           tLambda( k ) = tAir.cea_lambda( tT ) ;

           tHd( k ) = 0.0;

           // compute dissociation enthalpy
           tAir.Hf( tT, tHf );

           // compute dissociation enthalpy
           for( uint i=0; i<tAir.number_of_components(); ++i )
           {
               tHd( k ) += tY( i ) * tHf( i ) / tComp( i )->M();
           }

           if( tMolarFlag )
           {
               for( uint j=0; j<tAir.number_of_components(); ++j )
               {
                   Vector< real > & tV = tMapper.field( j + tFieldCount );
                   tV( k ) = tX( j );
               }
           }
       }

       tProgress.finish();

       std::cout << "    Time for computing field: " << tTimer.stop() << std::endl;

       tMapper.compute_node_values();




       tMapper.mesh()->save( "hotair.hdf5" );

        // only if interpolation order < 3
       // tMapper.mesh()->save( "hotair.exo" );

       return gComm.finalize();
   }