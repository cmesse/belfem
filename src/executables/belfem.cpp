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

#include <memory>

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
#include "cl_InputFile.hpp"
#include "cl_MaxwellFactory.hpp"
#include "cl_ThermalFactory.hpp"
#include "cl_ElectricalCircuitFactory.hpp"
#include "cl_Timer.hpp"
#include "cl_FEM_Controller.hpp"
#include "globals.hpp"

using namespace belfem;
using namespace fem;

Communicator gComm;

#if !defined( NDEBUG ) || defined( DEBUG )
Logger       gLog( InfoLevel::Detailed );
#else
Logger       gLog( InfoLevel::Default );
#endif

//------------------------------------------------------------------------------

// The deck decides the physics: an unlabeled `linear thermal` or
// `nonlinear thermal` section inside the solver block selects the coupled
// h-ɸ/T problem; without one, belfem solves the magnetic problem only.
// Every rank parses the same file and takes the same branch — a rank-split
// decision would deadlock the collective calls inside the factories.
bool
deck_requests_thermal( const InputFile & aInputFile )
{
    if ( ! aInputFile.section_exists( "solver" ) )
    {
        return false ;
    }

    const input::Section * tSolver = aInputFile.section( "solver" );

    return tSolver->section_exists( "linear thermal" )
        || tSolver->section_exists( "nonlinear thermal" );
}

//------------------------------------------------------------------------------

// -V / --version is the version pair shared by the BELFEM executables
// ( -v / --verbose is taken by the logger level, see cl_Arguments.hpp ).
// Every rank scans its own copy of the command line, so all ranks take the
// same branch and reach the same finalize()
bool
version_requested( const Arguments & aArguments )
{
    const Cell< string > & tArguments = aArguments.data();

    const index_t tNumArgs = tArguments.size();

    for ( index_t k = 1; k < tNumArgs; ++k )
    {
        if ( tArguments( k ) == "-V" || tArguments( k ) == "--version" )
        {
            return true ;
        }
    }

    return false ;
}

//------------------------------------------------------------------------------

// comment:
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

    // a version request prints the banner and exits before anything reads
    // the deck: `belfem --version` must work in a directory that has no
    // input.conf. print_banner() is rank-0 only, finalize() is collective
    if ( version_requested( tArguments ) )
    {
        print_banner();

        return gComm.finalize();
    }

    const string tFile = "input.conf" ;

    // decide the physics from the deck before any factory runs
    bool tHaveThermal ;
    bool tHaveThermalBCs ;
    {
        InputFile tInput( tFile );

        tHaveThermal    = deck_requests_thermal( tInput );

        tHaveThermalBCs = tInput.section_exists( "boundary conditions" )
            && tInput.section( "boundary conditions" )->section_exists( "thermal" );
    }

    print_banner( tHaveThermal ?
          "mixed h-ɸ/T electromagmetic-thermal simulation"
        : "h-ɸ electromagnetic simulation" );

    // a deck that declares thermal boundary conditions but no thermal solver
    // is inconsistent: refuse it rather than silently dropping the BCs
    BELFEM_ERROR( tHaveThermal || ! tHaveThermalBCs,
        "the input file defines thermal boundary conditions but no thermal solver.\n"
        "Add a 'linear thermal' or 'nonlinear thermal' section to the solver block\n"
        "to request the coupled problem, or remove the thermal boundary conditions." );

    if ( gComm.rank() == 0 )
    {
        message( InfoLevel::Default, tHaveThermal ?
              "    Thermal solver section found in the input file: solving the coupled h-ɸ/T problem."
            : "    No thermal solver section in the input file: solving the magnetic problem only." );
    }

    // load the mesh
    MaxwellFactory tMFactory( tFile );

    // reference temperature: the material laws need one even without a
    // thermal solve. Set on all ranks, after the factory has parsed an
    // explicit value and before create_magnetic_kernel() publishes the
    // temperature mesh global
    if ( std::isnan( gTbulk ) )
    {
        if ( gComm.rank() == 0 )
        {
            message( InfoLevel::Default, tHaveThermal ?
                  "    No initial temperature given in input.conf. Assuming 77 K."
                : "    No bulk temperature given in input.conf. Assuming 77 K." );
        }
        gTbulk = 77.0 ;
    }

    //Create the electrical circuit, and adding the required boundary conditions to the Maxwell Factory
    electronics::ElectricalCircuitFactory tElFactory( tFile, tMFactory.boundary_conditions() ) ;

    auto tKernel = tMFactory.create_magnetic_kernel() ;

    Vector< real > & tPhi = tKernel->mesh()->field_data("phi");
    tPhi.fill( 0.0 );

    //Get the current boundary conditions
    Cell < PhysicalBoundaryCondition *> tCurrentBCs = tMFactory.current_BCs() ;

    //Initialize the current vector (BELFEM_EPS to avoid singular matrices)
    Vector< real > tI = Vector< real >(tCurrentBCs.size(),BELFEM_EPS);
    reinterpret_cast< IWG_Maxwell * >( tKernel->dofmgr()->iwg() )->set_currents( tI );

    auto tControl = tMFactory.create_controller();

    // the thermal kernel exists only when the deck asked for it. The
    // controller is never handed a null kernel ( set_thermal_kernel
    // dereferences its argument ), and the magnetic-only run never enters
    // the thermal iteration methods. Both objects live until the end of
    // main, mirroring the lifetime hphiTrun gives them
    std::unique_ptr< ThermalFactory > tTFactory ;
    std::shared_ptr< Kernel >         tKernel2 ;

    if ( tHaveThermal )
    {
        //Create the thermal problem
        tTFactory.reset( new ThermalFactory( tFile, tKernel.get() ) );

        tKernel2 = tTFactory->create_thermal_kernel() ;

        tKernel2->mesh()->field_data("T").fill( gTbulk );

        reinterpret_cast< IWG_Timestep * >( tKernel2->dofmgr()->iwg() )->set_timestepping_method( tControl->euler_method() );

        //Add thermal kernel to controller
        tControl->set_thermal_kernel( tKernel2.get() ) ;
    }

    //send the circuit to the controller
    tControl->set_circuit( tElFactory.circuit() ) ;

    // check if backup exists
    tControl->load_memdump( "memdump.hdf5" );

    // outfile name
    string tExodusFile = sprint( "%s.e-s" , tMFactory.label().c_str() );

    // magnetic-only runs use the same loop as the fully coupled problem:
    // with no thermal kernel set, solve_coupled() skips the thermal phase.
    // The segregated loop is only reachable with a thermal kernel — its
    // thermal methods are unguarded by design
    if ( ( ! tHaveThermal ) || tControl->is_fullycoupled() )
    {
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

            tControl->solve_coupled();

            if (tControl->reset())
            {
                continue ;
            }

            // check if user wants to save timesteps but always save last timestep
            if ( tControl->save() || tControl->time() >= tControl->simulation_time() )
            {
                tControl->finalize( true ) ;

                // save IV must be called before the exodus file is written
                tControl->save_IV( "iv_results.csv" );
                tControl->save( tExodusFile );
                tControl->save_memdump( "memdump.hdf5" );
            }
            else
            {
                tControl->finalize( false ) ;
            }
        }
    }
    else //segregated case
    {
        while( tControl->time() < tControl->simulation_time() )
        {
            tControl->initialize_magnetic() ;

            // compute currents, if tCurrentBCs.length is smaller than the incidence matrix size, it means that the remaining currents are 0.
            uint tCount = 0;
            for ( PhysicalBoundaryCondition * tBC : tCurrentBCs )
            {
                tI(tCount++) = tBC->value() ;
            }

            // impose currents
            reinterpret_cast< IWG_Maxwell * >( tKernel->dofmgr()->iwg() )->set_currents( tI );

            tControl->solve_magnetic();

            if (tControl->reset())
            {
                continue ;
            }

            while ( abs(tControl->time() -  tControl->time_thermal()) > (1e-12))
            {

                tControl->initialize_thermal() ;
                tControl->solve_thermal();

                if (tControl->reset())
                {
                    break ;
                }

            }

            if (tControl->reset())
            {
                continue ;
            }

            // check if user wants to save timesteps but always save last timestep
            if ( tControl->save() || tControl->time() >= tControl->simulation_time() )
            {
                tControl->finalize( true ) ;

                tControl->save_IV( "iv_results.csv" );
                tControl->save( "hphi_results.e-s");
                tControl->save_memdump( "memdump.hdf5" );
            }
            else
            {
                tControl->finalize( false ) ;
            }

        }
    }

    return gComm.finalize();
}
