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

real
compute_criteria( Mesh * aMesh , const id_t aBlockID )
{
    // make sure that we don't run this with MPO in
    BELFEM_ERROR( comm_rank() == 0 , "This fuction is not made for parallel yet" );

    // grab field
    Vector< real > & tField = aMesh->field_exists( "Checkerboard" ) ?
                              aMesh->field_data( "Checkerboard") :
                              aMesh->create_field("Checkerboard", EntityType::ELEMENT );

    // grab jjc data
    const Vector< real > & tJJc = aMesh->field_data("elementJJc" );

    // reset data
    tField.fill( 0.0 );

    // grab elements from selected block
    Cell< mesh::Element * > & tElements = aMesh->block( aBlockID )->elements() ;

    index_t tCount = 0 ;
    real tSum = 0.0 ;

    // loop over all elements on selected block
    for( mesh::Element * tElement : tElements )
    {
        // get number of neigbors for this element
        uint tNumNeighbors = tElement->number_of_elements() ;

        // approximate value by averaging neighbor values
        real tApprox = 0 ;

        // loop over all element neighbors
        for( uint e = 0; e<tNumNeighbors; ++e )
        {
            // grab field data and add to value
            tApprox += tJJc( tElement->element( e )->index() );
        }

        // get value for this element
        real tValue = tJJc( tElement-> index() );

        // compute local average of field
        //real tAverage = ( tValue + tApprox ) / ( tNumNeighbors + 1 );

        // approximated value
        tApprox /= tNumNeighbors ;

        if( tApprox < BELFEM_EPSILON )
        {
            real tCriterion = std::abs( ( tValue - tApprox ) / tApprox );

            // save the value into the field
            tField( tElement->index() )  = tCriterion ;

            tSum += tCriterion * tCriterion ;
            ++tCount ;
        }

    }

    return std::sqrt( tSum / tCount );

}

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

    //
    const real tMaxdt = tFactory->maxdt();

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
    real tDeltaTime = tFormulation->timestep();
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

    //get the initial relaxation factors
    real tOmegaP = tNonlinMagnetic.picardOmega ;
    real tOmegaN = tNonlinMagnetic.newtonOmega ;
    real tOmegaPT = tNonlinThermal.picardOmega ;
    real tOmegaNT = tNonlinThermal.newtonOmega ;
    real tRate; //rate of the iterative method
    real tRateT; //rate of the thermal iterative method
    uint tIter0 = 0; // number of iterations of the last time step
    bool tRestartTime = false; //do we restart the time step because of poor convergence

    const uint tIterFast = tNonlinMagnetic.maxIter/2 ; //number of iteration defining fast convergence (ideally chosen by user)
    const real tTimeFastMult = 1.5 ; //Multiplication factor of time step when fast convergence (ideally chosen by user)
    const real tTimeSlowMult = 0.5 ; //Multiplication factor of time step when too many iter (ideally chosen by user)
    real tf = 1.0;

    //real tOmegaP = tTime == 0 ? 0.5 : tNonlinMagnetic.picardOmega ;
    //real tOmegaN = tTime == 0 ? 0.1 : tNonlinMagnetic.newtonOmega ;

    while( tTime < tMaxTime )
    {

        if( tHaveThermal )
        {
            tOmegaP = tNonlinMagnetic.picardOmega; //reset EM relax if thermal
            tOmegaN = tNonlinMagnetic.newtonOmega;
            //tFourier->shift_fields() ;
        }

        if( tTimeLoopCSV++ >= tCsvDump )
        {
            tTimeLoopCSV = 1 ;
        }

        if( (tTimeLoop++ >= tMeshDump) & (!tRestartTime) )
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


        // increment timestep (or restart with a smaller time step if necessary)
        if ( tRestartTime ) //restart with a smaller time step if necessary
        {
            tTime -= tDeltaTime; //go back one time step
            tTimeCount--;
            tDeltaTime *= tTimeSlowMult; //update the time step
            tFormulation->timestep()=tDeltaTime;
            //tFormulation->load( tBackupFile, tMesh  );
        }

        tTime += tDeltaTime; //increment

        // increment time counter
        tTimeCount++;

        if ( tRestartTime )
        {
            tFormulation->reset_fields( );
        }
        else
        {
            tFormulation->shift_fields();
        }

        tFormulation->compute_boundary_conditions( tTime );

        if( tHaveThermal )
        {
            if ( tRestartTime )
            {
                tFourier->reset_fields();
            }
            else
            {
                tFourier->shift_fields();
            }
        }

        tRestartTime = false ;
        uint tIter = 0;

        // residual
        real tEpsilon = BELFEM_REAL_MAX;
        real tEpsilon0 ;
        real tEpsilonT ;

        if (tHaveThermal)
        {
            tEpsilonT = BELFEM_REAL_MAX ;

        }
        else
        {
            tEpsilonT = 0 ;
        }
        real tEpsilonT0 = 0 ;


        if ( tKernel->is_master() )
        {
            std::cout << std::endl << "  --------------------------------------------------------------------------"
                      << std::endl;
            std::cout << "   time : " << tTime*1000.0 << " ms " << std::endl;
            std::cout << "  --------------------------------------------------------------------------" << std::endl;
        }

        Timer tTimer;

        bool tNewton = false; //did the algorithm go in Newton yet?
        while ( tEpsilon > tNonlinMagnetic.newtonEpsilon || tIter < tNonlinMagnetic.minIter || tEpsilonT > tNonlinThermal.newtonEpsilon )
        {
            if( tHaveThermal )
            {
                if( tEpsilonT > tNonlinThermal.picardEpsilon )
                {
                    tFourier->set_algorithm( SolverAlgorithm::Picard );
                    tFourier->set_omega( tOmegaPT );
                }
                else
                {
                    tFourier->set_algorithm( SolverAlgorithm::NewtonRaphson );
                    tFourier->set_omega( tOmegaNT );
                }
            }


            if ( tEpsilon > tNonlinMagnetic.picardEpsilon )
            {
                tFormulation->set_algorithm( SolverAlgorithm::Picard );
                //tFormulation->set_algorithm( SolverAlgorithm::NewtonRaphson );
                if (tNewton) //did the algorithm already been in Newton
                {
                    tOmegaP=tOmegaN;
                    tNewton = false;
                }
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
                if (!tNewton) //first time in Newton
                {
                    tOmegaN=std::min(tOmegaP,tOmegaN) ; //we don't want the algorithm to have a too high relax if slowly converging
                    tNewton = true;
                }
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
            tMagfield->collect_field( "elementJJc");

            tMagfield->solve();
            tIter++;

            if( tHaveThermal )
            {
                tSynch->magnetic_to_thermal_b_and_ej() ;
                tThermalField->compute_jacobian_and_rhs();
                tThermalField->solve() ;

                //is the iteration convergent or divergent for the thermal problem
                tEpsilonT0 = tEpsilonT ; //last residual
                tEpsilonT = tThermalField->residual( tIter ); //current residual
                tRateT = (std::abs(tEpsilonT)-std::abs(tEpsilonT0))/std::abs(tEpsilonT0);

                if (tRateT >= 0.0) //diverging
                {
                    //tf = 0.5;
                    //std::cout << "Multiplication factor: " << tf << std::endl;
                    if ( tEpsilonT > tNonlinThermal.picardEpsilon )
                    {
                        //tOmegaPT=std::max(tf*tOmegaPT, 0.05);
                        tOmegaPT=0.5;
                    }
                    else
                    {
                        //tOmegaNT=std::max(tf*tOmegaNT, 0.05);
                        tOmegaNT=0.5;
                    }
                }
                else //converging
                {
                    //tf = 1.1+2*0.4*std::atan(-1.0*tRateT)/constant::pi;
                    //tf = 1.5;
                    //std::cout << "Multiplication factor: " << tf << std::endl;
                    if ( tEpsilonT > tNonlinThermal.picardEpsilon )
                    {
                        //tOmegaP=std::min(tRelaxConvMult*tOmegaP, 0.5);
                        //tOmegaPT=std::min(tf*tOmegaPT, 2.0); //allow sur-relax?
                        tOmegaPT=1.0;
                    }
                    else
                    {
                        //tOmegaN=std::min(tRelaxConvMult*tOmegaN, 0.5);
                        //tOmegaNT=std::min(tf*tOmegaNT, 2.0); //allow sur-relax?
                        tOmegaNT=1.0;
                    }
                }

                tSynch->thermal_to_magnetic_T() ;
            }

            tEpsilon0 = tEpsilon ;
            tEpsilon = tMagfield->residual( tIter );
            tRate = (std::abs(tEpsilon)-std::abs(tEpsilon0))/std::abs(tEpsilon0);
            //std::cout << "Rate of convergence: " << tRate << std::endl;

            // Increase or decrease the relax factors depending on the previous and the current residual
            if (tRate >= 0.0) //diverging
            {
                tf = 0.5;
                if ( tEpsilon > tNonlinMagnetic.picardEpsilon )
                {
                    tOmegaP=std::max(tf*tOmegaP, 0.0005);
                }
                else
                {
                    tOmegaN=std::max(tf*tOmegaN, 0.0005);
                }
            }
            else //converging
            {
                tf = 1.1+2*0.4*std::atan(-1.0*tRate)/constant::pi;
                if ( tEpsilon > tNonlinMagnetic.picardEpsilon )
                {
                    tOmegaP=std::min(tf*tOmegaP, 1.0);
                }
                else
                {
                    tOmegaN=std::min(tf*tOmegaN, 1.0);
                }
            }


            if ( tKernel->is_master() )
            {


                string tAlgLabel = tFormulation->algorithm() == SolverAlgorithm::Picard ? " P " : " NR";
                if( tHaveThermal )
                {
                    string tAlgLabelT = tFourier->algorithm() == SolverAlgorithm::Picard ? " P " : " NR";
                    real tTmax = max( tThermalMesh->field_data( "T") );

                    /*std::cout << "    it:  " << tIter << tAlgLabel << " omega " << tFormulation->omega()
                              << " log10(eps): " << std::round( std::log10( tEpsilon ) * 100 ) * 0.01
                            << " log10(epsT): " << std::round( std::log10( tEpsilonT ) * 100 ) * 0.01
                            << " Tmax: " << tTmax
                            << std::endl;*/
                    std::cout << "    it:  " << tIter << " EM " << tAlgLabel << " omega " << tFormulation->omega()
                              << " log10(eps) " << std::round( std::log10( tEpsilon ) * 100 ) * 0.01 << std::endl;
                    std::cout << "           Th " << tAlgLabelT <<  " omega " << tFourier->omega()
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
                tRestartTime = true; //set flag to restart with smaller time step
                tOmegaP = tNonlinMagnetic.picardOmega ; //reset the relaxation factors
                tOmegaN = tNonlinMagnetic.newtonOmega ;
                if ( tKernel->is_master() )
                {
                    const Vector< real > & tRhs = tMagfield->rhs_vector() ;
                    //std::cout << "    WARNING: too many iterations. Exiting this time step" << std::endl ;
                    std::cout << "    WARNING: too many iterations. Exiting this time step and decreasing time step to "<<
                              tTimeSlowMult*tDeltaTime*1000 << " ms" << std::endl ;
                    tMagfield->print_worst_dof();
                }
                break ;
            }
            /*else if ( (  std::abs( std::log10( tEpsilon ) - std::log10( tEpsilon0 ) ) < 1e-4 ) &&
            ( tEpsilon > tNonlinMagnetic.newtonEpsilon ) )
            {
                tRestartTime = true; //set flag to restart with smaller time step
                tOmegaP = tNonlinMagnetic.picardOmega ; //reset the relaxation factors
                tOmegaN = tNonlinMagnetic.newtonOmega ;
                if ( tKernel->is_master() )
                {
                    //tMagfield->solver()->
                    //std::cout << "    WARNING: desired convergence could not be reached" << std::endl ;
                    std::cout << "    WARNING: desired convergence could not be reached, decreasing time step to "<<
                              tTimeSlowMult*tDeltaTime*1000 << " ms" << std::endl;
                    tMagfield->print_worst_dof();
                }
                break ;
            }*/
        }

        if ( tKernel->is_master() )
        {
            gLog.message( 1, "    timestep completed in %4.2f seconds", ( float ) tTimer.stop() * 0.001 );
        }


        // todo: move into postprocess routine
        if( (tTimeLoopCSV == 1) & (!tRestartTime))
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

                // test new checkerboarding criteria
                real tCrit = compute_criteria( tMesh, 1 );
                std::cout << "    checkerboard criterion : " << tCrit << std::endl ;

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


        //adaptive time step
        /*if(( tIter > 40) & (!tRestartTime)) //slow convergence (can be modified)
        {
            tDeltaTime = 0.75*tDeltaTime; //factor can be modified
            tFormulation->timestep()=tDeltaTime;
            std::cout << "Slow convergence, decreasing the time step to " << tDeltaTime << std::endl;
        }*/
        if (( tIter < tIterFast ) & ( !tRestartTime )) //fast convergence
        {
            tDeltaTime = std::min( tTimeFastMult * tDeltaTime, tMaxdt ); //increasing time step
            tFormulation->timestep() = tDeltaTime;
            if ( tKernel->is_master() )
            {
                std::cout << "Fast convergence, increasing the time step to " << tDeltaTime * 1000 << " ms" << std::endl;
            }
        }


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

    // close communicator
    return gComm.finalize();
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
