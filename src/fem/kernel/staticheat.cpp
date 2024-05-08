//
// Created by Christian Messe on 27.07.20.
//

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
#include "cl_IWG_StaticHeatConduction.hpp"
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
    gComm = Communicator( argc, argv );

    print_banner();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Load the Mesh
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // load the mesh
    Mesh * tMesh = new Mesh( "combustor.msh" );

    // assume the mesh was set in mm
    tMesh->scale_mesh( 0.001 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Setup the problem
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // create the parameters object
    KernelParameters tParams( tMesh );

    // with the parameters object set, we create the kernel
    Kernel tKernel( &tParams );

    // create the equation object
    //IWG_StaticHeatConduction tIWG( tMesh->number_of_dimensions() );

    IWG * tIWG = tKernel.create_equation( IwgType::StaticHeatConduction );

    tIWG->select_sidesets( { 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 } );
    tIWG->select_block( 26 );

    // now we grab the first field of the kernel
    //Field * tField = tKernel.field( 0 );
    DofManager * tField = tKernel.create_field( tIWG );


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// User Settings and Boundary conditions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // initial guess
    real tTw0 = 300.0 ;


    // select the solver for this problem
#ifdef BELFEM_PETSC
    tField->set_solver( SolverType::PETSC );
#else
    tField->set_solver( gDefaultSolver );
#endif




    // select the material
    tField->block( 26 )->set_material( MaterialType::Copper ); // 26

    // set hotgas temperature
    tField->sideset( 1 )->impose_dirichlet( 600.0 ); // 1
    //tField->sideset( 1 )->impose_neumann( 60e6 );
    //tField->sideset( 1 )->impose_alpha( 1e5, 800.0 );

    // set coldgas temperature
    Vector< id_t > tChannels = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };

    for ( id_t tID: tChannels )
    {
        tField->sideset( tID )->impose_dirichlet( 100.0 );
    }

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
    tField->field_data( "T" ).fill( tTw0 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Start the computation
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real tResidual = BELFEM_REAL_MAX ;

    uint tIter = 0 ;

    tMesh->save( "combustor.exo");

    while( tResidual > 1e-7 )
    {
        tField->reset_convection() ;
        // tField->compute_surface_loads();
        tField->compute_jacobian_and_rhs();


        tField->solve();
        tResidual = tField->residual( tIter );

        if( comm_rank() == 0 )
        {

            std::cout << "    iteration:  " << tIter++ << " epsilon: " << tResidual << std::endl;
        }
    }


    tMesh->save( "combustor.exo");

    // tidy up mesh
    delete tMesh ;

    return gComm.finalize();
}