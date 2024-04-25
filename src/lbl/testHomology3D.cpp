//
// Created by Gregory Giard on 2014-03-12.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Mesh_ExodusWriter.hpp"
#include "cl_Mesh.hpp"
#include "cl_SimplicialComplex.hpp"

using namespace belfem ;
using namespace mesh ;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );


    // print_banner();

    // load the mesh
    Mesh * tMesh = new Mesh( "mesh.msh" );
    tMesh->create_edges();
    tMesh->create_faces();

    std::cout << tMesh->number_of_faces() << std::endl;
    std::cout << tMesh->number_of_dimensions() << std::endl;


    SimplicialComplex * tSComplex = new SimplicialComplex(tMesh) ;

    // convert node coordinates from mm to m
    tMesh->scale_mesh( 0.001 );

    tMesh->save( "mesh.exo" );

    delete tMesh;

    delete tSComplex ;

    return gComm.finalize();
}
