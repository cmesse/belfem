//
// Created by christian on 7/27/22.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "cl_BS_Mapper.hpp"
#include "banner.hpp"
#include "cl_Timer.hpp"
#include "stringtools.hpp"
#include "cl_MaterialFactory.hpp"
#include "fn_linspace.hpp"
#include "cl_Progressbar.hpp"

using namespace belfem;
using namespace bspline;

Communicator gComm;
Logger       gLog( 3 );


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    print_banner();


    MaterialFactory tFactory ;


    Vector< uint > tRRR = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
                            125, 150, 175, 200, 225, 250, 275, 300 };


    Material * tMaterial = tFactory.create_material( "silver" );
    tMaterial->use_splines( true );

    string tFile = "materials.hdf5" ;

    for( uint r=0; r<tRRR.length(); ++r )
    {
        std::cout << "RRR " << tRRR( r ) << std::endl ;

        tMaterial->set_rrr( tRRR( r ));

        string tRhoLabel = sprint( "SilverRho_RRR%u", tRRR( r ));
        string tKLabel = sprint( "SilverK_RRR%u", tRRR( r ));

        Mapper tMapper( 2, 2,
                        { 251, 100 },
                        { 0.0, 0.0 },
                        { 500.0, 50.0 } );


        const Matrix< real > & tGrid = tMapper.integration_grid();
        Vector< real > & tK = tMapper.create_field( tKLabel );
        Vector< real > & tRho = tMapper.create_field( tRhoLabel );

        index_t tGridSize = tGrid.n_cols();

        Progressbar tProgress( tGridSize );

        // compute values
        for ( index_t k = 0; k < tGridSize; ++k )
        {
            tProgress.step( k );

            // temperature
            real tT = tGrid( 0, k );

            // magnetic field
            real tB = tGrid( 1, k );

            tK( k ) = tMaterial->lambda( tT, tB );
            tRho( k ) = tMaterial->rho_el( 0, tT, tB );
        }
        tProgress.finish();
        tMapper.compute_node_values();

        // tMapper.mesh()->save( "data.exo");
        tMapper.write_field_to_database( tRhoLabel, tFile );
        tMapper.write_field_to_database( tKLabel, tFile );

    }

    delete tMaterial ;

    return gComm.finalize();
}