//
// Created by Christian Messe on 2019-07-25.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Mesh_ExodusWriter.hpp"
#include "cl_Mesh.hpp"
#include "cl_Mesh_OrientationChecker.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );


    // print_banner();

    // load the mesh
    Mesh * tMesh = new Mesh( "tetsindie.msh" );

    // convert node coordinates from mm to m
    //tMesh->scale_mesh( 0.001 );

    // create an orientation checker
    for( mesh::Block * tBlock : tMesh->blocks() )
    {

        mesh::OrientationChecker * tChecker = new mesh::OrientationChecker();
        tChecker->set_element_type( tBlock->element_type() );

        for ( mesh::Element * tElement: tMesh->elements())
        {
            tChecker->process_element( tElement );
        }

        delete tChecker ;
    }

    tMesh->create_edges( false );
    tMesh->create_faces( false );

    //mesh::ExodusWriter tWriter( tMesh );
    //tMesh->save( "mesh.vtk" );
    tMesh->save( "mesh.exo" );


    delete tMesh ;

    return  gComm.finalize();
}