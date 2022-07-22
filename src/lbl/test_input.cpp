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
#include "cl_Profiler.hpp"
#include "fn_FEM_compute_normb.hpp"

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

using namespace belfem ;
using namespace fem ;

Communicator gComm;
Logger       gLog( 3 );

//------------------------------------------------------------------------------


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // initialize job
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // read the input file
    InputFile tInputFile( "input.txt" );

    // create a new factory
    MaxwellFactory * tFactory = new MaxwellFactory( tInputFile );

    // get the mesh
    Mesh * tMesh = tFactory->mesh() ;

    /*tMesh->unflag_all_nodes();
    tMesh->block( 1 )->flag_nodes();
    for( mesh::Node * tNode : tMesh->nodes() )
    {
        if( tNode->is_flagged() )
        {
            real tX = tNode->x() ;
            real tY = tNode->y() ;
            tNode->set_coords( tX, tY, 0.001 );
        }
    }

    tMesh->unflag_all_nodes();
    tMesh->block( 2 )->flag_nodes();
    for( mesh::Node * tNode : tMesh->nodes() )
    {
        if( tNode->is_flagged() )
        {
            real tX = tNode->x() ;
            real tY = tNode->y() ;
            tNode->set_coords( tX, tY, 0.00066 );
        }
    }

    tMesh->unflag_all_nodes();
    tMesh->block( 3 )->flag_nodes();
    for( mesh::Node * tNode : tMesh->nodes() )
    {
        if( tNode->is_flagged() )
        {
            real tX = tNode->x() ;
            real tY = tNode->y() ;
            tNode->set_coords( tX, tY, -0.00066 );
        }
    }


    tMesh->unflag_all_nodes();
    tMesh->block( 4 )->flag_nodes();
    for( mesh::Node * tNode : tMesh->nodes() )
    {
        if( tNode->is_flagged() )
        {
            real tX = tNode->x() ;
            real tY = tNode->y() ;
            tNode->set_coords( tX, tY, -0.001 );
        }
    }
    tMesh->save( "mesh.vtk"); */

    // get the kernel
    Kernel * tKernel = tFactory->kernel() ;

    // get the magnetic field
    DofManager * tMagfield = tFactory->magnetic_field() ;

    // get the formulation
    IWG_Maxwell * tFormulation = reinterpret_cast< IWG_Maxwell * >( tMagfield->iwg() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // some parameters
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // get the name for the outmesh
    const string tOutFile = tFactory->outfile();

    // get the name of the backup file
    const string tBackupFile = tFactory->backupfile() ;

    // get the restart flag
    const bool tRestart =  tFactory->restart();

    // time when simulation ends
    const real tMaxTime = tFactory->maxtime() ;

    // get the nonlinear settings
    NonlinearSettings tNonlinear = tFactory->nonlinear_settings() ;

    // delete the factory
    delete tFactory ;

    // get the timestep
    uint & tTimeCount = tMesh->time_step() ;

    tMesh->save( "test.exo");

    // check if we have to load a field from HDF5


    tMagfield->initialize() ; // must be called before loading mesh data
    //tFormulation->set_nan_values() ;

    if ( file_exists( tBackupFile ) && ! tRestart )
    {
        int tStatus = ( file_exists( tBackupFile ) && ! tRestart ) ?
                      tFormulation->load( tBackupFile, tMesh  ) : 1 ;

        if ( tStatus == 1 )
        {
            tTimeCount = 0 ;
        }
    }
    else
    {
        tTimeCount = 0 ;
    }

    //tMesh->save( "test.exo");

   const real tDeltaTime = tFormulation->timestep();
   real & tTime = tMesh->time_stamp() ;

   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // begin timeloop
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   real tOmegaP = tNonlinear.picardOmega ;
   real tOmegaN = tNonlinear.newtonOmega ;


    /*
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
       real tEpsilon0 ;

       if ( tKernel->is_master() )
       {
           std::cout << std::endl << "  --------------------------------------------------------------------------"
                     << std::endl;
           std::cout << "   time : " << tTime << " s " << std::endl;
           std::cout << "  --------------------------------------------------------------------------" << std::endl;
       }

       Timer tTimer;

       index_t tCountP = 0 ;
       index_t tCountN = 0 ;

       while ( tEpsilon > tNonlinear.newtonEpsilon || tIter < tNonlinear.minIter )
       {
           if ( tEpsilon > tNonlinear.picardEpsilon )
           {
               tFormulation->set_algorithm( SolverAlgorithm::Picard );
               tFormulation->set_omega( tOmegaP );
               if( tCountP++ > 0 )
               {
                   if ( tEpsilon < tEpsilon0 )
                   {

                       tOmegaP *= tCountP + std::min( std::pow( tEpsilon0 / tEpsilon, 0.25 ), 1.05 );
                       tOmegaP /= tCountP + 1 ;
                       if ( tOmegaP > 0.99 )
                       {
                           tOmegaP = 0.99;
                       }
                   }
                   else
                   {
                       tOmegaP *= tCountP;
                       tOmegaP /= tCountP + 1 ;
                   }
               }
               tCountN = 0 ;
           }
           else
           {
               tFormulation->set_algorithm( SolverAlgorithm::NewtonRaphson );
               tFormulation->set_omega( tOmegaN );
               if( tCountN++ > 1 )
               {
                   if ( tEpsilon < tEpsilon0 )
                   {

                       tOmegaN *= tCountN + std::min( std::pow( tEpsilon0 / tEpsilon, 0.25 ), 1.05 );
                       tOmegaN /= tCountN + 1 ;
                       if ( tOmegaN > 0.99 )
                       {
                           tOmegaN = 0.99;
                       }
                   }
                   else
                   {
                       tOmegaN *= tCountN ;
                       tOmegaN /= tCountN + 1 ;
                   }
               }
           }

           tMagfield->compute_jacobian_and_rhs();
           tMagfield->solve();

           tEpsilon0 = tEpsilon ;
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
       tMagfield->postprocess() ;

       // save mesh
       tMesh->save( tOutFile );
       tMesh->save( "test.exo");

       // save backup
       tFormulation->save( tBackupFile );

       comm_barrier() ;
   } */

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // tidy up
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    delete tKernel ;
    delete tMesh ;

    // close communicator
    return gComm.finalize();
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif