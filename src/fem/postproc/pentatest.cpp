/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Mesh_ExodusWriter.hpp"
#include "cl_Mesh.hpp"
#include "fn_Mesh_compute_volume.hpp"
#include "fn_Mesh_compute_surface_normals.hpp"
#include "fn_Mesh_compute_surface.hpp"
#include "constants.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );



    Mesh * tMesh2 = new Mesh( "cube.msh" );
    tMesh2->scale_mesh( 0.001 );


    // compute the volume
    real tVolume = mesh::compute_volume( tMesh2, 1 );

    // compute the surface
    real tSurface = mesh::compute_surface( tMesh2, { 2 } );

    std::cout << "Volume: " << tVolume << std::endl ;
    std::cout << "Surface: " << tSurface << std::endl ;


    delete tMesh2 ;


    // create a reader object
    Mesh * tMesh = new Mesh( "sphere.msh" );
    tMesh->scale_mesh( 0.001 );

    mesh::compute_surface_normals( tMesh, { 6, 9} );

    // compute the volume
    tVolume = mesh::compute_volume( tMesh, 1 );

    // compute the surface
    tSurface = mesh::compute_surface( tMesh, { 6, 9 } );

    /*
     // create a reader object
    Mesh * tMesh = new Mesh( "pentatest.msh" );

    mesh::compute_surface_normals( tMesh, { 4, 5, 6 } );

    // compute the volume
    real tVolume = mesh::compute_volume( tMesh, { 1, 2, 3 } );

    // compute the surface
    real tSurface = mesh::compute_surface( tMesh, { 4, 5, 6 } );

     */

    std::cout << "Volume: " << tVolume * 3 / constant::pi << std::endl ;
    std::cout << "Surface: " << tSurface / constant::pi << std::endl ;

    tMesh->save( "mesh.exo");

    delete tMesh ;




    return  gComm.finalize();
}