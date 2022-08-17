//
// Created by christian on 9/17/21.
//

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Timer.hpp"

#include "cl_InputFile.hpp"
#include "cl_MaxwellFactory.hpp"
#include "fn_FEM_compute_element_current.hpp"
#include "cl_Profiler.hpp"
#include "fn_FEM_compute_biot_savart.hpp"
#include "fn_FEM_compute_normb.hpp"

using namespace belfem ;
using namespace fem ;

Communicator gComm;
Logger       gLog( 4 );

//------------------------------------------------------------------------------


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    InputFile tInputFile( "input.txt" );

    MaxwellFactory * tFactory = new MaxwellFactory( tInputFile );

    Mesh * tMesh = tFactory->mesh() ;

    // get kernel
    Kernel * tKernel = tFactory->kernel() ;


    // get the magnetic field
    DofManager * tMagfield = tFactory->magnetic_field() ;

    // get the formulation
    IWG_Maxwell_Old * tFormulation = reinterpret_cast< IWG_Maxwell_Old * >( tMagfield->iwg() );

    // get the name for the outmesh
    const string tOutFile = tFactory->outfile();

    // get the name of the backup file
    const string tBackupFile = tFactory->backupfile() ;

    // get the restart flag
    const bool tRestart =  tFactory->restart();

    const real tMaxTime = tFactory->maxtime() ;

//------------------------------------------------------------------------------

    // get the nonlinear settings
    NonlinearSettings tNonlinear = tFactory->nonlinear_settings() ;

    // extra postproccessing flags
    bool tComputeBiotSavart = tFactory->compute_biot_savart() ;
    bool tComputeNormB      = tFactory->compute_normb() ;

    delete tFactory ;

//------------------------------------------------------------------------------

    // get the timestep
    uint & tTimeCount = tMesh->time_step() ;


    int tStatus = 1 ;
    if( file_exists( tBackupFile ) && ! tRestart )
    {
        tStatus = tFormulation->load( tBackupFile, tMesh  );
    }

//------------------------------------------------------------------------------


    const real tDeltaTime = tFormulation->timestep();
    real & tTime = tMesh->time_stamp() ;

//------------------------------------------------------------------------------

    // #hack
    while( tTime < tMaxTime )
    {
        // increment timestep
        tTime += tDeltaTime;

        // increment time counter
        tTimeCount++;

        // reset number of iterations
        tFormulation->shift_fields();
        tFormulation->compute_boundary_conditions( tTime );


        uint tIter = 0;

        // residual
        real tEpsilon = BELFEM_REAL_MAX;

        if ( tKernel->is_master())
        {
            std::cout << std::endl << "  --------------------------------------------------------------------------"
                      << std::endl;
            std::cout << "   time : " << tTime << " s " << std::endl;
            std::cout << "  --------------------------------------------------------------------------" << std::endl;
        }

        Timer tTimer;


        while ( tEpsilon > tNonlinear.newtonEpsilon || tIter < tNonlinear.minIter )
        {
            if ( tEpsilon > tNonlinear.picardEpsilon )
            {
                tFormulation->set_algorithm( SolverAlgorithm::Picard );
                tFormulation->set_omega( tNonlinear.picardOmega );
            }
            else
            {
                tFormulation->set_algorithm( SolverAlgorithm::NewtonRaphson );
                tFormulation->set_omega( tNonlinear.newtonOmega );
            }

            tMagfield->compute_jacobian_and_rhs();
            tMagfield->save_system("matrix.hdf5");
            tMagfield->solve();

            tEpsilon = tMagfield->residual( tIter++ );

            if ( tKernel->is_master() )
            {
                string tAlgLabel = tFormulation->algorithm() == SolverAlgorithm::Picard ? " P " : " NR";

                std::cout << "    iteration:  " << tIter << tAlgLabel << " omega " << tFormulation->omega()
                          << " log10(epsilon): " << std::round( std::log10( tEpsilon ) * 100 ) * 0.01 << std::endl;
            }

            BELFEM_ERROR( tIter < tNonlinear.maxIter, "too many iterations" );
        }
        if ( tKernel->is_master() )
        {
            gLog.message( 1, "    timestep completed in %4.2f seconds", ( float ) tTimer.stop() * 0.001 );
        }

        // postprocess
        compute_element_current( tMagfield );
        tMagfield->postprocess() ;

        // todo: tidy this up

        if ( tComputeNormB )
        {
            compute_normb( tMagfield, false );
        }


        if( tComputeBiotSavart )
        {
            compute_biot_savart( tMagfield );

            if( comm_rank() == 0 &&  tComputeNormB )
            {
                compute_normb( tMagfield, true );
            }
        }

        tMesh->save( tOutFile );
        tFormulation->save( tBackupFile );

        comm_barrier() ;
    }
    delete tKernel ;
    delete tMesh ;

//------------------------------------------------------------------------------

    // close communicator
    return gComm.finalize();
}