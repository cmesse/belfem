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
    Mesh * tMesh = new Mesh( "mesh.msh" );

    // convert node coordinates from mm to m
    tMesh->scale_mesh( 0.001 );

    /*// create a field
    Vector< real > & tElemField = tMesh->create_field( "MyElementField", EntityType::ELEMENT );
    Vector< real > & tNodeField = tMesh->create_field( "MyNodeField", EntityType::NODE );

    // create a global variable
    tMesh->create_global_variable( "Pi", 3.141 );

    for( uint k=0; k<tMesh->number_of_nodes(); ++k )
    {
        tNodeField( k ) = ( real ) k;
    }

    for( uint k=0; k<tMesh->number_of_elements(); ++k )
    {
        tElemField( k ) = ( real ) k;
    } */

    //mesh::ExodusWriter tWriter( tMesh );
    //tMesh->save( "mesh.vtk" );
    tMesh->save( "mesh.exo" );

    delete tMesh ;

    return  gComm.finalize();
}