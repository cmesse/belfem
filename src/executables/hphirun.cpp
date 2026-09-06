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

#include "typedefs.hpp"
#include "cl_Arguments.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_CutFactory.hpp"
#include "fn_linspace.hpp"
#include "cl_MaxwellFactory.hpp"
#include "cl_ElectricalCircuitFactory.hpp"
#include "cl_Timer.hpp"
#include "cl_FEM_Controller.hpp"
#include "globals.hpp"
#include "banner.hpp"

using namespace belfem;
using namespace fem;

Communicator gComm;

#if !defined( NDEBUG ) || defined( DEBUG )
Logger       gLog( InfoLevel::Detailed );
#else
Logger       gLog( InfoLevel::Default );
#endif

// retired reference driver, not built ( see CMakeLists.txt ); the equivalent run is:
// mpirun -np $nprocs --bind-to socket belfem
int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // parse the command line ( -v / --verbose sets the logger info level;
    // PETSc and STRUMPACK read their own flags from gComm, see the
    // sparse module documentation )
    Arguments tArguments( argc, argv );

    print_banner( "h-ɸ electromagnetic simulation");


    // load the mesh
    MaxwellFactory tMFactory( "input.conf" );

    if ( std::isnan( gTbulk ) )
    {
        if ( gComm.rank() == 0 )
            message( InfoLevel::Default, "    Note: the temperature value for a non-thermal problem was not set, fixing it at 77 K " );
        gTbulk = 77.0 ;
    }

    //Create the electrical circuit, and adding the required boundary conditions to the Maxwell Factory
    electronics::ElectricalCircuitFactory tElFactory( "input.conf", tMFactory.boundary_conditions() ) ;

    auto tKernel = tMFactory.create_magnetic_kernel() ;

    Vector< real > & tPhi = tKernel->mesh()->field_data("phi");
    tPhi.fill( 0.0 );

    //Get the current boundary conditions
    Cell < PhysicalBoundaryCondition *> tCurrentBCs = tMFactory.current_BCs() ;

    //Initialize the current vector (BELFEM_EPS to avoid singular matrices)
    Vector< real > tI = Vector< real >(tCurrentBCs.size(),BELFEM_EPS);
    reinterpret_cast< IWG_Maxwell * >( tKernel->dofmgr()->iwg() )->set_currents( tI );


    auto tControl = tMFactory.create_controller();


    //send the circuit to the controller
    tControl->set_circuit( tElFactory.circuit() ) ;

    // check if backup exists
    tControl->load_memdump( "memdump.hdf5" );

    // this is needed if we want to compute the matrix conditioning in mumps
    //tKernel->dofmgr()->solver()->set_mumps_error_analysis( MumpsErrorAnalysis::Full );
    while( tControl->time() < tControl->simulation_time() )
    {

        // also updates boundary condition values
        tControl->initialize_timestep();



        // compute currents, if tCurrentBCs.length is smaller than the incidence matrix size, it means that the remaining currents are 0.
        uint tCount = 0;
        for ( PhysicalBoundaryCondition * tBC : tCurrentBCs )
        {
            tI(tCount++) = tBC->value() ;
        }

        // impose currents
        reinterpret_cast< IWG_Maxwell * >( tKernel->dofmgr()->iwg() )->set_currents( tI );
        // the certified-exit loop lives in the controller now
        tControl->solve_coupled();
        if (tControl->reset())
        {
            continue ;
        }

        // check if user wants to save timesteps but always save last timestep
        if ( tControl->save() || tControl->time() >= tControl->simulation_time() )
        {
            tControl->finalize( true ) ;
            // save_IV must run BEFORE the exodus file is written: it refreshes
            // the I_/U_ mesh globals that save() then snapshots into the frame
            tControl->save_IV( "iv_results.csv" );
            tControl->save( "hphi_results.e-s");
            tControl->save_memdump( "memdump.hdf5" );
        }
        else
        {
            tControl->finalize( false ) ;
        }
    }

    return gComm.finalize();
}
