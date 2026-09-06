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
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"

#include "cl_Mesh.hpp"

#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "../iwg/cl_IWG_StaticHeatConduction.hpp"
#include "fn_norm.hpp"
#include "fn_Mesh_compute_surface_normals.hpp"

using namespace belfem;
using namespace fem;

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    print_banner();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Load the Mesh
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // load the mesh
    Mesh * tMesh = new Mesh( "dipole.msh" );

    // assume the mesh was set in mm
    //tMesh->scale_mesh( 0.001 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Setup the problem
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // create the parameters object
    KernelParameters tParams( tMesh );

    // with the parameters object set, we create the kernel
    Kernel tKernel( &tParams );

    // create the equation object
    //IWG_StaticHeatConduction tIWG( tMesh->number_of_dimensions() );

    IWG * tIWG = tKernel.create_equation( IwgType::Poisson );

    tIWG->select_sidesets( { 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } );
    tIWG->select_block( 3 );

    // now we grab the first field of the kernel
    //Field * tField = tKernel.field( 0 );
    DofManager * tField = tKernel.create_field( tIWG );


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// User Settings and Boundary conditions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    // select the solver for this problem
#ifdef BELFEM_PETSC
    tField->set_solver( SolverType::PETSc );
#else
    tField->set_solver( gDefaultSolver );
#endif




    // select the material
    tField->block( 26 )->set_material( MaterialType::Copper ); // 26

    // set hotgas temperature
    tField->sideset( 1 )->impose_dirichlet( 1 ); // 1
    tField->sideset( 2 )->impose_dirichlet( 1 ); // 1
    tField->sideset( 3 )->impose_dirichlet( 1 ); // 1
    tField->sideset( 4 )->impose_dirichlet( 1 ); // 1
    tField->sideset( 5 )->impose_dirichlet( 1 ); // 1
    tField->sideset( 6 )->impose_dirichlet( 1 ); // 1
    tField->sideset( 7 )->impose_dirichlet( 1 ); // 1
    tField->sideset( 8 )->impose_dirichlet( 1 ); // 1

    //tField->sideset( 1 )->impose_neumann( 60e6 );
    //tField->sideset( 1 )->impose_alpha( 1e5, 800.0 );

    // set coldgas temperature
    tField->sideset( 9 )->impose_dirichlet( 0 ); // 1
    tField->sideset( 10 )->impose_dirichlet( 0 ); // 1
    tField->sideset( 11 )->impose_dirichlet( 0 ); // 1
    tField->sideset( 12 )->impose_dirichlet( 0 ); // 1

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Initialize Kernel
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    comm_barrier();

    // initialize the dofs and detect wetted sidesets
    //tField->init_dofs() ;

    // compute the surface normals of the mesh and redistribute over all procs
    //mesh::compute_surface_normals( tField->mesh(), tIWG.wetted_sidesets() );

    //tField->initialize_jacobian();
    tField->initialize() ;


    // initialize starting conditions
    tField->field_data( "phi" ).fill( 0 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Start the computation
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    tField->compute_jacobian();
    tField->solve();

    tMesh->save( "poisson.exo");

    // tidy up mesh
    delete tMesh ;

    return gComm.finalize();
}