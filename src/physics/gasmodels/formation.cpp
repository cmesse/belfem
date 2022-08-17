//
// Created by Christian Messe on 20.03.20.
// Note: This is a RAND Algorithm, but there seems to be a mistake somewhere
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"

#define private public
#include "cl_Gas.hpp"
#undef private

#include "cl_GT_RefGas.hpp"
#include "fn_trans.hpp"
#include "fn_sum.hpp"
#include "fn_dot.hpp"
#include "fn_gesv.hpp"
#include "fn_max.hpp"
#include "fn_min.hpp"
#include "fn_norm.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{

    // create communicator
    gComm = Communicator( argc, argv );


   Gas tAir( {
                      "N2",
                      "O2",
                      "Ar",
                      "e-",
                      "N+",
                      "NO+",
                      "O+",
                      "N",
                      "O2",
                      "O",
                      "NO",
                      "Ar+"
              },
              {
                      7.8084e-1,
                      2.0948e-1,
                      9.34e-3,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.00
              } );


    /*Gas tAir( {
                      "Ar", "CH4", "CO", "CO2", "H2", "He", "Kr",
                      "N2", "N2O", "Ne", "NO2", "O2", "O3", "Xe",
                      "O", "N", "NO" },
              { 9.34e-3,   // Ar
                2.0e-6,    // CH4
                2.5e-7,    // CO
                3.14e-4,   // CO2
                5.0e-7,    // H2
                5.24e-6,   // He
                1.14e-6,   // Kr
                7.8084e-1, // N2
                3.1e-7,    // N2O
                1.818e-5,  // Ne
                2.0e-8,    // NO2
                2.0948e-1, // O2
                8.0e-6,    // O3
                8.7e-8,    // Xe
                0.0,       // O
                0.0,       // N
                0.0        // NO
              } ); */
    tAir.print();

    // real tT = 5000.0 ;

    real tT = 300.0;
    real tP = 1.00E5 ;


    while( tT <= 6000.0 )
    {
        tAir.reset_mixture();
        tAir.remix_to_equilibrium( tT, tP );

        //const Vector< real > & tX = tAir.molar_fractions();
        //tX.print("Molar");
        //std::cout << tT << "; " << tX( 0 )  << "; " << tX( 1 ) << "; "
        //    << tX( 2 ) << "; " << tX( 3 ) << "; " << tX( 4 ) << "; "   << tX( 5 ) << "; "   << tX( 6 ) << "; "  << std::endl ;
        //std::cout << tT << " ; " << tAir.M() *1000.0 <<  "; " << tAir.R() <<std::endl;
        std::cout << tT << " ; " << tAir.M( tT, tP ) *1000.0 << " ; " << tAir.R( tT, tP ) << " ; " << tAir.cp( tT, tP ) << "; " << tAir.h( tT, tP ) << "; " << tAir.s( tT, tP )  << "; " <<  tAir.mu( tT, tP )*1e6 <<  " ; " << tAir.lambda( tT, tP )*1E3 << "; "  << tAir.gamma( tT, tP ) << "; " << tAir.Pr( tT, tP ) << std::endl;
        tT += 100 ;
    }

    return gComm.finalize();
}