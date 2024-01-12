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
#include "fn_sum.hpp"
#include "fn_max.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Homology.hpp"
#include "cl_Cohomology.hpp"


//#include "fn_FEM_compute_element_current_thinshell.hpp"
#include "cl_MaxwellJob.hpp" // delete me

#include "cl_Pipette.hpp"
#include "fn_rcond.hpp"
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
    gComm = Communicator( argc, argv );

    if( comm_rank() == 0 )
    {
        print_banner( "- input file test program -" );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // initialize job
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // read the input file
    InputFile tInputFile( "input.txt" );


    // create a new factory
    MaxwellFactory * tFactory = new MaxwellFactory( tInputFile );

    // get the mesh
    Mesh * tMesh = tFactory->magnetic_mesh() ;

    // get the simplicial complex
    mesh::SimplicialComplex * tSComplex = tFactory->simplicial_complex();

    // Reduce the simplicial complex
    // start timer
    Timer tTimer;
    std::cout << "Reducing the complex ..." << std::endl;
    tSComplex->reduce_complexCCR();
    std::cout << "Complex reduced in: "<< tTimer.stop()*1e-3 << " s" << std::endl;

    std::cout << "Number of 0-chains: " << tSComplex->number_of_0simplices() <<std::endl;
    std::cout << "Number of 1-chains: " << tSComplex->number_of_1simplices() <<std::endl;
    std::cout << "Number of 2-chains: " << tSComplex->number_of_2simplices() <<std::endl;

    // Coreduce the simplicial complex
    Timer tTimer2;
    std::cout << "(Co)reducing the complex ..." << std::endl;
    tSComplex->coreduce_complexCCR();
    std::cout << "Complex (co)reduced in: "<< tTimer2.stop()*1e-3 << " s" << std::endl;

    std::cout << "Number of 0-cochains: " << tSComplex->number_of_kcosimplices(0) <<std::endl;
    std::cout << "Number of 1-cochains: " << tSComplex->number_of_kcosimplices(1) <<std::endl;
    std::cout << "Number of 2-cochains: " << tSComplex->number_of_kcosimplices(2) <<std::endl;

    //Compute the homology and cohomology generators
    Timer tTimer3;
    std::cout << "Computing Homology ..." << std::endl;
    mesh::Homology * tHomology = new mesh::Homology(tSComplex, tMesh);
    std::cout << "Homology generators computed in: "<< tTimer3.stop() << " ms" << std::endl;
    tHomology->create_kGeneratorsField(1);

    Timer tTimer4;
    std::cout << "Computing Cohomology ..." << std::endl;
    mesh::Cohomology * tCohomology = new mesh::Cohomology(tSComplex, tMesh);
    std::cout << "Cohomology generators computed in: "<< tTimer4.stop() << " ms" << std::endl;
    tCohomology->create_kGeneratorsField(1);

    // flag telling if we compute the norm of b
    bool tComputeNormB      = tFactory->compute_normb() ;

    // get the kernel
    Kernel * tKernel = tFactory->magnetic_kernel() ;


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
    NonlinearSettings tNonlinMagnetic = tFactory->nonlinear_settings() ;

    const uint tMeshDump = tFactory->meshdump() ;
    const uint tCsvDump  = tFactory->csvdump() ;


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //  thermal stuff
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // check if thermal field exists
    const bool tHaveThermal = tFactory->have_thermal() ;

    Kernel * tThermalKernel = tFactory->thermal_kernel() ;
    Mesh   * tThermalMesh   = tFactory->thermal_mesh() ;
    MaxwellMeshSynch * tSynch = tFactory->thermal_synch() ;
    DofManager * tThermalField = tFactory->thermal_field() ;
    IWG_Timestep * tFourier = tFactory->fourier() ;

    NonlinearSettings tNonlinThermal = tHaveThermal ?
                                       tFactory->nonlinear_settings( MaxwellFieldType::THERMAL )
                                                    : tNonlinMagnetic ;

    if( tHaveThermal )
    {


        tFourier->set_algorithm( SolverAlgorithm::NewtonRaphson );

    }
    // delete the factory
    delete tFactory ;

    // get the timestep
    uint & tTimeCount = tMesh->time_step() ;

    // check if we have to load a field from HDF5
    tMagfield->initialize() ; // must be called before loading mesh data

    if( tHaveThermal )
    {
        tThermalField->initialize() ;
    }

    //tFormulation->set_nan_values() ;

    if ( file_exists( tBackupFile ) && ! tRestart )
    {
        int tStatus = ( file_exists( tBackupFile ) && ! tRestart ) ?
                      tFormulation->load( tBackupFile, tMesh  ) : 1 ;

        if ( tStatus == 1 )
        {
            tTimeCount = 0 ;
        }
        else if ( tHaveThermal )
        {
            tSynch->magnetic_to_thermal_T() ;
        }
    }
    else
    {
        tTimeCount = 0 ;
    }

    // end delete me
    const real tDeltaTime = tFormulation->timestep();
    real & tTime = tMesh->time_stamp() ;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // hide fields we don't need
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /*if( comm_rank() == 0 )
    {
        Cell< string > tHide = { "SurfaceNormalsx",
                                 "SurfaceNormalsy",
                                 "SurfaceNormalsz",
                                 "az",
                                 "jz" };

        for ( const string & tLabel: tHide )
        {
            tMesh->field( tLabel )->set_write_to_file_flag( false );
        }
    }*/

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // begin timeloop
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    tMagfield->solver()->set_mumps_error_analysis( MumpsErrorAnalysis::Full );

    uint & tTimeLoop = tFormulation->time_loop();

    //if( tTimeLoop == 0 )
    //{
    //    tMesh->save( tOutFile );
    //}
    //tMesh->save( tOutFile );
    uint tTimeLoopCSV = 1 ;



    while( tTime < tMaxTime )
    {
        if( tHaveThermal )
        {
            tFourier->shift_fields() ;
        }

        if( tTimeLoopCSV++ >= tCsvDump )
        {
            tTimeLoopCSV = 1 ;
        }

        if( tTimeLoop++ >= tMeshDump )
        {
            // save backup
            /*string tString = sprint("%s.%04u", tBackupFile.c_str(), ( unsigned int ) tTimeCount );

            // save backup
            tFormulation->save( tString ); */

            // save backup
            tFormulation->save( tBackupFile );

            // todo: move this somewhere else
            if( tMesh->number_of_dimensions() == 3 )
            {
                tMagfield->synchronize_fields( { "elementJx", "elementJy", "elementJz", "elementEJ", "elementJJc", "elementRho" } );
            }
            else
            {
                tMagfield->synchronize_fields( { "elementJz", "elementEJ", "elementJJc", "elementRho" } );
            }


            tMagfield->postprocess() ;

            if( tComputeNormB )
            {
                compute_normb( tMagfield, false );
            }

            // save mesh
            real tOldTime = tTime;
            tTime *= 1000.0;
            tMesh->save( tOutFile );
            tTime = tOldTime;
            tTimeLoop = 1;
        }

        real tOmegaP = tTime == 0 ? 0.5 : tNonlinMagnetic.picardOmega ;
        real tOmegaN = tTime == 0 ? 0.1 : tNonlinMagnetic.newtonOmega ;
        //real tOmegaP = tNonlinMagnetic.picardOmega ;
        //real tOmegaN = tNonlinMagnetic.newtonOmega ;

        // increment timestep
        tTime += tDeltaTime;

        // increment time counter
        tTimeCount++;

        // reset number of iterations
        tFormulation->shift_fields();
        tFormulation->compute_boundary_conditions( tTime );

        if( tHaveThermal )
        {
            tFourier->shift_fields() ;
        }

        uint tIter = 0;

        // residual
        real tEpsilon = BELFEM_REAL_MAX;
        real tEpsilon0 ;

        real tEpsilonT = 0 ;


        if ( tKernel->is_master() )
        {
            std::cout << std::endl << "  --------------------------------------------------------------------------"
                      << std::endl;
            std::cout << "   time : " << tTime*1000.0 << " ms " << std::endl;
            std::cout << "  --------------------------------------------------------------------------" << std::endl;
        }

        Timer tTimer;

        //index_t tCountP = 0 ;
        //index_t tCountN = 0 ;


        while ( tEpsilon > tNonlinMagnetic.newtonEpsilon || tIter < tNonlinMagnetic.minIter || tEpsilonT > tNonlinThermal.newtonEpsilon )
        {
            if( tHaveThermal )
            {
                if( tEpsilonT > tNonlinThermal.picardEpsilon )
                {
                    tFourier->set_algorithm( SolverAlgorithm::Picard );
                    tFourier->set_omega( tNonlinThermal.picardOmega );
                }
                else
                {
                    tFourier->set_algorithm( SolverAlgorithm::NewtonRaphson );
                    tFourier->set_omega( tNonlinThermal.newtonOmega );
                }
            }


            if ( tEpsilon > tNonlinMagnetic.picardEpsilon )
            {
                tFormulation->set_algorithm( SolverAlgorithm::Picard );
                tFormulation->set_omega( tOmegaP );

                /*if( tCountP++ > 0 )
                {
                        tOmegaP *= std::max(std::min( std::pow( tEpsilon0 / tEpsilon, 0.25 ), 1.05 ), 0.5 );
                        if ( tOmegaP > 1.0 )
                        {
                            tOmegaP = 1.0 ;
                        }
                }
                tCountN = 0 ;
                tOmegaN = std::min( tOmegaN, tOmegaP ); */
            }
            else
            {
                //tFormulation->set_algorithm( SolverAlgorithm::Picard );
                tFormulation->set_algorithm( SolverAlgorithm::NewtonRaphson );
                tFormulation->set_omega( tOmegaN );
                /*if( tCountN++ > 1 )
                {
                    tOmegaN *= std::max(std::min( std::pow( tEpsilon0 / tEpsilon, 0.25 ), 1.05 ), 0.5 );

                    if ( tOmegaN > 1.0 )
                    {
                        tOmegaN = 1.0;
                    }
                }*/
            }

            //tFormulation->compute_thin_shell_error_2d();

            tMagfield->compute_jacobian_and_rhs();

            // synchronize ej because we need this for the quench
            tMagfield->collect_field("elementEJ");
            tMagfield->collect_field( "elementJz");

            tMagfield->solve();

            if( tHaveThermal )
            {
                tSynch->magnetic_to_thermal_b_and_ej() ;
                tThermalField->compute_jacobian_and_rhs();
                tThermalField->solve() ;
                tEpsilonT = tThermalField->residual( tIter );
                tSynch->thermal_to_magnetic_T() ;
            }

            tEpsilon0 = tEpsilon ;
            tEpsilon = tMagfield->residual( tIter++ );

            if ( tKernel->is_master() )
            {
                string tAlgLabel = tFormulation->algorithm() == SolverAlgorithm::Picard ? " P " : " NR";
                if( tHaveThermal )
                {
                    real tTmax = max( tThermalMesh->field_data( "T") );

                    std::cout << "    it:  " << tIter << tAlgLabel << " omega " << tFormulation->omega()
                              << " log10(eps): " << std::round( std::log10( tEpsilon ) * 100 ) * 0.01
                              << " log10(epsT): " << std::round( std::log10( tEpsilonT ) * 100 ) * 0.01
                              << " Tmax: " << tTmax
                              << std::endl;
                }
                else
                {
                    std::cout << "    it:  " << tIter << tAlgLabel << " omega " << tFormulation->omega()
                              << " log10(eps): " << std::round( std::log10( tEpsilon ) * 100 ) * 0.01 << std::endl;
                }
            }

            if( tIter > tNonlinMagnetic.maxIter )
            {
                if ( tKernel->is_master() )
                {
                    const Vector< real > & tRhs = tMagfield->rhs_vector() ;
                    std::cout << "    WARNING: too many iterations. Exiting this time step" << std::endl ;
                    tMagfield->print_worst_dof();
                }
                break ;
            }
            else if ( (  std::abs( std::log10( tEpsilon ) - std::log10( tEpsilon0 ) ) < 1e-4 ) &&
                      ( tEpsilon > tNonlinMagnetic.newtonEpsilon ) )
            {
                if ( tKernel->is_master() )
                {
                    //tMagfield->solver()->
                    std::cout << "    WARNING: desired convergence could not be reached" << std::endl ;
                    tMagfield->print_worst_dof();
                }
                break ;
            }
        }
        if ( tKernel->is_master() )
        {
            gLog.message( 1, "    timestep completed in %4.2f seconds", ( float ) tTimer.stop() * 0.001 );
        }


        // todo: move into postprocess routine
        if( tTimeLoopCSV == 1 )
        {
            if ( comm_rank() == 0 )
            {
                real tEI = 0;
                real tI = 0;
                mesh::Pipette tPipette;

                Vector< real > & tEJel = tMesh->field_data( "elementEJ" );
                Vector< real > & tJel = tMesh->field_data( "elementJz" );
                Vector< real > & tJJcel = tMesh->field_data( "elementJJc" );
                Vector< real > & tRhoel = tMesh->field_data( "elementRho" );

                string tTimeString = sprint( "%04u", tTimeCount );

                /* std::ofstream tTimeFile ;
                 tTimeFile.open( "time.csv", std::ios_base::app);
                 tTimeFile << tTimeCount << "," << tTime << std::endl ;
                 tTimeFile.close() ; */


                // compute the rcond function
                real tK = 0;
                real tR = 0;
                if ( tMagfield->solver()->type() == SolverType::MUMPS )
                {
                    real tC0 = tMagfield->solver()->wrapper()->get_cond0();
                    real tC1 = tMagfield->solver()->wrapper()->get_cond1();
                    tK = tC0 > tC1 ? tC0 : tC1;

                    // Sofia's truncation function
                    tR = tK / ( 1 - tK * BELFEM_EPS ) * 2.0 * BELFEM_EPS;
                }


                for ( id_t tID: tMesh->ghost_block_ids())
                {
                    string tFileName = "tape_" + std::to_string( tID ) + "_" + tTimeString + ".csv";

                    std::ofstream tTapeFile;
                    tTapeFile.open( tFileName, std::ios_base::app );
                    ////tTapeFile << "t=" << tTime << ", " << std::endl ;

                    tPipette.set_element_type( tMesh->block( tID )->element_type());
                    tTapeFile << "x,v,jjc,e*j,rho" << std::endl;

                    for ( mesh::Element * tElement: tMesh->block( tID )->elements())
                    {
                        real tV = tPipette.measure( tElement );
                        tEI += tEJel( tElement->index()) * tV;
                        tI += tJel( tElement->index()) * tV;
                        tTapeFile << 0.5 * ( tElement->node( 0 )->x() + tElement->node( 1 )->x()) << "," << tV << ","
                                  << tJJcel( tElement->index()) << "," << tEJel( tElement->index()) << ","
                                  << tRhoel( tElement->index()) << std::endl;
                    }
                    tTapeFile.close();

                }

                for ( fem::Block * tBlock: tMagfield->blocks())
                {
                    if ( tBlock->domain_type() == DomainType::Conductor && tBlock->element_type() == ElementType::TRI6 )
                    {

                        tFormulation->link_to_group( tBlock );

                        // get elements
                        Cell< fem::Element * > & tElements = tBlock->elements();

                        for ( fem::Element * tElement: tElements )
                        {
                            tI += tFormulation->compute_element_current( tElement );
                        }
                    }
                }
                tMagfield->synchronize_fields( { "elementJz" } );
                Vector< real > tAllI( comm_size(), 0.0 );
                receive( tKernel->comm_table(), tAllI );
                tAllI( 0 ) = 0;

                tI += sum( tAllI );


                //std::cout << "-----" << std::endl ;
                // std::cout  <<  tTime << ", " << tI << ", " << tEI << std::endl ;

                std::ofstream tCSV;
                tCSV.open( "acloss.csv", std::ios_base::app );
                tCSV << tTime << ", " << tI << ", " << tEI << ", " << std::log10( tK ) << ", " << std::log10( tR )
                     << std::endl;
                tCSV.close();

                if ( tMagfield->solver()->type() == SolverType::MUMPS )
                {
                    std::cout << "    cond: " << tK << " " << tR << std::endl;
                }
            }
            else
            {
                real tI = 0.0;

                for ( fem::Block * tBlock: tMagfield->blocks())
                {
                    if ( tBlock->domain_type() == DomainType::Conductor && tBlock->element_type() == ElementType::TRI6 )
                    {

                        tFormulation->link_to_group( tBlock );

                        // get elements
                        Cell< fem::Element * > & tElements = tBlock->elements();

                        for ( fem::Element * tElement: tElements )
                        {
                            tI += tFormulation->compute_element_current( tElement );
                        }
                    }
                }
                tMagfield->synchronize_fields( { "elementJz" } );

                send( 0, tI );


            }
        }

        //fem::compute_element_current_thinshell_tri3( tMagfield );
        //tMagfield->postprocess() ;

        //if( tComputeNormB )
        //{
        //    compute_normb( tMagfield, false );
        //}

        //tMesh->save( "test.exo");


        comm_barrier() ;
    }

    //fem::compute_element_current_thinshell_tri3( tMagfield );
    tMagfield->postprocess() ;

    if( tComputeNormB )
    {
        compute_normb( tMagfield, false );
    }

    tMesh->save( tOutFile );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // tidy up
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( tHaveThermal )
    {
        delete tSynch ;
        delete tThermalKernel ;
        delete tThermalMesh ;
    }

    delete tKernel ;
    delete tMesh ;

    delete tHomology;

    delete tCohomology;

    // close communicator
    return gComm.finalize();
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
