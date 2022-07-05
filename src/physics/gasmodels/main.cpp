#include <iostream>


#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "fn_linspace.hpp"

#include "cl_Gas.hpp"
#include "cl_HDF5.hpp"


using namespace belfem;
using namespace gastables;
using namespace gasmodels;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

//------------------------------------------------------------------------------

    uint tN = 25 ;

    Vector< real > tT = linspace( 100.0, 2500.0, tN );

    // crate Air model
    Gas tAir;
    real tP = 1.0e5;

    // Matrix with values
    Matrix< real > tData( tN, 6 );
    Vector< real > tRefpoint( 6 );

    for( uint k=0; k<tN; ++k )
    {
        printf( "%12.3f %12.6f %12.6f %12.6f %12.6f %12.6f ;\n",
                tT( k ),
                tAir.cp( tT( k ) , tP  ),
                tAir.h( tT( k ) , tP ) * 1e-3,
                tAir.s( tT( k ) , tP ),
                tAir.mu( tT( k ) , tP ) * 1e6,
                tAir.lambda( tT( k ) , tP ) * 1e3 );

        tData( k, 0 ) = tT( k );
        tData( k, 1 ) = tAir.cp( tT( k ) , tP  );
        tData( k, 2 ) = tAir.h( tT( k ) , tP ) ;
        tData( k, 3 ) = tAir.s( tT( k ) , tP ) ;
        tData( k, 4 ) = tAir.mu( tT( k ) , tP ) ;
        tData( k, 5 ) = tAir.lambda( tT( k ) , tP ) ;

    }

    HDF5 tFile( "gasdata.h5", FileMode::NEW );
    tFile.save_data( "Table", tData );

    real tTref = 298.15 ;

    tRefpoint( 0 ) = tTref ;
    tRefpoint( 1 ) = tAir.cp(tTref , tP  );
    tRefpoint( 2 ) = tAir.h( tTref , tP ) ;
    tRefpoint( 3 ) = tAir.s( tTref , tP ) ;
    tRefpoint( 4 ) = tAir.mu( tTref , tP ) ;
    tRefpoint( 5 ) = tAir.lambda( tTref , tP ) ;

    tFile.save_data( "Refpoint", tRefpoint );

    std::cout << tTref
              << " " << tAir.cp( tTref , tP  )
              << " " << tAir.h( tTref , tP ) * 1e-3
              << " " << tAir.s( tTref, tP )
              << " " << tAir.mu( tTref , tP ) * 1e6
              << " " << tAir.lambda( tTref , tP ) * 1e3 << std::endl;

    printf( "R: %12.5f\n", tAir.R( tTref, tP ) );

//------------------------------------------------------------------------------

    return gComm.finalize();

}