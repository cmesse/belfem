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
    tMesh->create_edges( false );
    tMesh->create_faces( false );

    // check orientations (optional)
    index_t tCount = tMesh->check() ;
    std::cout << "reoriented " << tCount << " elements." << std::endl ;

    mesh::Edge * tEdge = tMesh->edges()(0);
    std::cout << tEdge->number_of_faces() << std::endl ;

    //mesh::ExodusWriter tWriter( tMesh );
    //tMesh->save( "mesh.vtk" );
    tMesh->save( "mesh.exo" );


    delete tMesh ;

    return  gComm.finalize();
}