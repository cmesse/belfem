//
// Created by Christian Messe on 03.08.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Timer.hpp"

#include "cl_Mesh.hpp"

#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_IWG_TransientHeatConduction.hpp"
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
    Kernel tKernel( &tParams, true );

    // create the equation object
    IWG_TransientHeatConduction tIWG( tMesh->number_of_dimensions() );

    // now we grab the first field of the kernel
    Field * tField = tKernel.field( 0 );

    // and link it with the equation object
    tField->set_integrated_weak_governing_equation( &tIWG );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// User Settings and Boundary conditions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // initial guess
    real tTw0 = 100.0;

    // initialize starting conditions
    tField->field_data( "T" ).fill( tTw0 );
    tField->field_data( "T0" ).fill( tTw0 );

    // select the solver for this problem
#ifdef BELFEM_PETSC
    tField->set_solver( SolverType::PETSC );
#else
    tField->set_solver( gDefaultSolver );
#endif
    // select the material
    tField->block( 26 )->set_material( MaterialType::Copper ); // 26

    // set hotgas temperature
    //tField->sideset( 1 )->impose_dirichlet( 600.0 ); // 1
    //tField->sideset( 1 )->impose_neumann( 60e6 );
     tField->sideset( 1 )->impose_alpha( 1e5, 800.0 );

    // set coldgas temperature
    Vector< id_t > tChannels = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };

    for ( id_t tID: tChannels )
    {
        tField->sideset( tID )->impose_dirichlet( 100.0 );
    }


    // set the end time
    real tStartTime = 0.0 ;
    real tEndTime   = 2.0 ;

    // computation timestep
    tIWG.set_timestep( 0.05 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Initialize Kernel
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    comm_barrier();

    // initialize the dofs and detect wetted sidesets
    tField->init_dofs() ;

    tField->initialize_jacobian();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Start the computation
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // grab the rank from MPI
    const proc_t tMyRank = comm_rank();

    // initialize time circuits and flux capacitor
    uint  tTimeCount = 1 ;    // <-- for Exodus. Must begin with 1, not zero!
    real  tTime = tStartTime ;
    Timer tTimer ;


    if( tMyRank == 0 )
    {
        std::cout << "  --------------------------------------------------------------------------" << std::endl;
        std::cout << "   computing problem using " << tField->solver()->wrapper()->label() << std::endl;
        std::cout << "  --------------------------------------------------------------------------" << std::endl;
    }

    // begin time loop
    while( tTime < tEndTime )
    {
        tTimer.reset() ;

        // write step into mesh
        tMesh->set_time_step( tTimeCount );
        tMesh->time_stamp() = tTime ;

        // residual
        real tResidual = BELFEM_REAL_MAX;

        uint tIter = 0 ;

        if( tMyRank == 0 )
        {
            std::cout << std::endl  << "  --------------------------------------------------------------------------"
                      << std::endl;
            std::cout << "   time : " << tTime << " s " << std::endl ;
            std::cout << "  --------------------------------------------------------------------------" << std::endl  ;
        }


        // catch error
        while ( tResidual > 1e-6 )
        {
            // compute the nodal distribution
            tField->compute_surface_loads();

            tField->compute_jacobian_and_rhs();

            tField->solve();

            tResidual = tField->residual( tIter ) ;

            if( tMyRank == 0 )
            {

                std::cout << "    iteration:  " << tIter++ << " epsilon: " << tResidual << std::endl;
            }

        }

        if( tMyRank == 0 )
        {
            std::cout << std::endl;
        }

        // shift temperature
        tField->field_data( "T0" ) = tField->field_data( "T" );

        if( tMyRank == 0 )
        {
            Vector< real > tT = tMesh->field_data("T");

            std::cout << "    time for computing timestep:  " << tTimer.stop() * 0.001 << " s" << std::endl << std::endl;
        }

        tMesh->save( "combustor.e-s" );

        tTime += tIWG.timestep() ;
        ++tTimeCount;
    }

    // tidy up mesh
    delete tMesh ;

    return gComm.finalize();
}