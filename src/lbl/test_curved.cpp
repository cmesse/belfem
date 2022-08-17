//
// Created by Christian Messe on 05.01.22.
//

#include <iostream>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Mesh.hpp"
#include "cl_Logger.hpp"

using namespace belfem;


Communicator gComm;
Logger       gLog( 4 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    // open mesh
    Mesh tMesh( "test.msh");

    // scale mesh to meters
    tMesh.scale_mesh( 0.001 );

    tMesh.flag_curved_elements() ;

    // create a new field
    Vector< real > & tField = tMesh.create_field( "IsCurved", EntityType::ELEMENT ) ;

    for( mesh::Element * tElement : tMesh.elements() )
    {
        if( tElement->is_curved() )
        {
            tField( tElement->index() ) = 1.0 ;
            std::cout << "curved " << tElement->id() << " " << tElement->index() << std::endl ;
        }
        else
        {
            tField( tElement->index() ) = 0.0 ;
        }
    }

    tMesh.save( "test.exo");

    // close communicator
    return gComm.finalize();
}