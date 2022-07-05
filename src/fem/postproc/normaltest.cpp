//
// Created by Christian Messe on 14.06.20.
//


#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Mesh_ExodusWriter.hpp"
#include "cl_Mesh.hpp"
#include "fn_Mesh_compute_surface_normals.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );


    // print_banner();

    // create a reader object
    Mesh * tMesh = new Mesh( "/tmp/mesh.msh" );
    tMesh->scale_mesh( 0.001 );

    mesh::compute_surface_normals( tMesh, { 30, 31 }, GroupType::SIDESET );

    mesh::ExodusWriter tWriter( tMesh );
    tWriter.save( "/tmp/mesh.exo" );

    delete tMesh ;

    return  gComm.finalize();
}