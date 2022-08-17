
//
// Created by christian on 8/2/21.
//

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_IWG_Maxwell_Old.hpp"
#include "cl_FEM_Kernel.hpp"
#include "fn_Mesh_compute_volume.hpp"
#include "cl_Profiler.hpp"
#include "fn_sum.hpp"

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

    print_banner();

    Mesh tMesh = Mesh("coil4.msh") ;
    tMesh.scale_mesh( 0.001 );
    tMesh.create_edges();

    tMesh.create_field("ElementJz", EntityType::ELEMENT );

//------------------------------------------------------------------------------

    KernelParameters tParams( &tMesh );

    Kernel tKernel( &tParams );

    IWG_Maxwell * tMaxwell = reinterpret_cast< IWG_Maxwell *  >
            ( tKernel.create_equation( IwgType::MAXWELL_HA2D ) );

    IWG_Maxwell * tCurrent = reinterpret_cast< IWG_Maxwell *  >
            ( tKernel.create_equation( IwgType::MAXWELL_E2J ) );

    IWG_Maxwell * tMagneticA = reinterpret_cast< IWG_Maxwell *  >
            ( tKernel.create_equation( IwgType::MAXWELL_E2B ) );

    IWG_Maxwell * tMagneticB = reinterpret_cast< IWG_Maxwell *  >
            ( tKernel.create_equation( IwgType::MAXWELL_E2B ) );

//------------------------------------------------------------------------------
// User Settings
//------------------------------------------------------------------------------

    /*Vector< id_t > tPlus( 1, 2 );
    Vector< id_t > tMinus( 1, 3 );

    // block identification
    tMaxwell->set_block( 1, DomainType::SuperConductor );
    tMaxwell->set_block( { 2, 3 }, DomainType::Coil );
    tMaxwell->set_block( 4, DomainType::Air );
    tMaxwell->set_sideset( { 1, 2, 3, 4, 5, 6, 7, 8 },  DomainType::Interface ); */

    Vector< id_t > tPlus( 1, 3 );
    Vector< id_t > tMinus( 1, 4 );

    tMaxwell->set_block( { 1, 2 }, DomainType::SuperConductor );
    tMaxwell->set_block( { 3, 4 }, DomainType::Coil );
    tMaxwell->set_block( 5, DomainType::Air );
    if( tMaxwell->type() == IwgType::MAXWELL_HA2D )
    {
        tMaxwell->set_sideset( { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 },  DomainType::Interface );
    }
    else if ( tMaxwell->type() == IwgType::MAXWELL_HPHI2D )
    {
        tMaxwell->set_sideset( { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 },  DomainType::Interface );
    }

    // settings for power law
    tMaxwell->powerlaw()->set_coefficients( 1e-5, 150, 25 );

    tMaxwell->set_euler_method( EulerMethod::BackwardImplicit );

    // current ramp
    real tImax     = 150 ;
    real tRamptime = 10 ;
    real tMaxTime = 0.001 ;
    real tDeltaTime = 1e-7 ;

    real tOmegaPicard = 0.1 ;
    real tOmegaNewtonRaphson = 0.01 ;
    real tEpsilonPicard = 1e-3 ;
    real tEpsilonNewtonRaphson = 1e-7 ;

    DofManager * tField = tKernel.create_field( tMaxwell );


    tCurrent->set_block( 1, DomainType::SuperConductor );
    DofManager * tField2 = tKernel.create_field( tCurrent );


    tMagneticA->set_block( { 1, 2 }, DomainType::SuperConductor );
    tMagneticB->set_block( { 3, 4 }, DomainType::Coil );
    tMagneticB->set_block( 5, DomainType::Air );

    DofManager * tField3 = tKernel.create_field( tMagneticA );
    DofManager * tField4 = tKernel.create_field( tMagneticB );

    // select solver library
    tField->set_solver( SolverType::PARDISO );
    //tField->solver()->set_petsc(
    //        Preconditioner::BJACOBI,
    //        KrylovMethod::GMRES,
    //        1e-9 );

    tField2->set_solver( SolverType::MUMPS );
    tField3->set_solver( SolverType::MUMPS );
    tField4->set_solver( SolverType::MUMPS );

    // return gComm.finalize() ;

    tField->initialize() ;


    // compute the size of the matrix
    uint tNumDofs = tField->jacobian()->n_cols() ;
    broadcast( 0, tNumDofs );

    if( tKernel.is_master() )
    {
        tMesh.save("mesh.exo");
    }

    // special setting for sc block
    //if ( tField->block_exists( 1 ) )
    //{
    //    tField->block( 1 )->set_integration_order( 20 );
    //}

    tMaxwell->set_euler_method(EulerMethod::BackwardImplicit );

    tMaxwell->set_timestep( tDeltaTime );

//------------------------------------------------------------------------------


    // cross section of positive coil
    real tAplus = mesh::compute_volume( tKernel.mesh(), tPlus );

    // cross section of negative coil
    real tAminus = mesh::compute_volume( tKernel.mesh(), tMinus );


//------------------------------------------------------------------------------



    uint tTimeCount = 1 ;
    real & tTime = tMesh.time_stamp() ;


    tTime = tMaxwell->timestep() ;

    //uint tPow = 1 ;
    if( tMaxwell->type() == IwgType::MAXWELL_HA2D )
    {
        tField2->compute_jacobian();
        tField3->compute_jacobian();
        tField4->compute_jacobian();
    }
    else if ( tMaxwell->type() == IwgType::MAXWELL_HPHI2D )
    {
        tField2->compute_jacobian();
    }

    real & tI = tMesh.create_global_variable( "Current");
    real & tTimeStamp = tMesh.create_global_variable( "Time");


    while( tTime < tMaxTime )
    {

        tMaxwell->shift_fields();

        uint tIter = 0 ;

        tMesh.set_time_step( tTimeCount );


        // compute imposed current density

        tTimeStamp = tTime ;
        tI = tTime < tRamptime ? tTime / tRamptime * tImax : tImax ;
        //real tI = tTime < tRamptime ? tDeltaTime / tRamptime * tImax : tImax ;
        if( tKernel.is_master() )
        {
            std::cout << std::endl  << "  --------------------------------------------------------------------------"
            << std::endl;
            std::cout << "   time : " << tTime << " s " << ",   I : " << tI << " A,  dofs : " << tNumDofs <<std::endl ;
            std::cout << "  --------------------------------------------------------------------------" << std::endl  ;
        }


        comm_barrier();

        tMaxwell->set_current_density( tPlus, tI / tAplus );
        tMaxwell->set_current_density( tMinus, -tI / tAminus );

        // compute rhs

        tField->compute_volume_loads( { 3, 4 } );


        real tResidual = BELFEM_REAL_MAX ;

        //Profiler tProf ;
        //tProf.start() ;

        real tOmega = 0 ;

        tMaxwell->set_omega( tOmegaPicard );
        tMaxwell->set_algorithm( SolverAlgorithm::Picard );

        while( tResidual > tEpsilonNewtonRaphson || tIter < 3 )
        {
            if( tMaxwell->algorithm() == SolverAlgorithm::NewtonRaphson )
            {
                if( tIter > 500 )
                {
                    tOmega = tOmegaNewtonRaphson * 0.05;
                }
                else if( tIter > 250 )
                {
                    tOmega = tOmegaNewtonRaphson * 0.25 ;
                }
                else if( tIter > 100 )
                {
                    tOmega = tOmegaNewtonRaphson*0.5;  ;
                }
                else
                {
                    tOmega = tOmegaNewtonRaphson ;
                }
            }
            else
            {
                tOmega = tOmegaPicard ;
            }
            tMaxwell->set_omega( tOmega );

            tField->compute_jacobian_and_rhs();
            //tField->save_system("matrix.hdf5");

            tField->solve();

            tResidual = tField->residual( tIter++ );



            //if( tKernel.is_master() )
            //{
            //    tField->save_system("matrix.hdf5");
            //}

            /*if( tAvgCount < tNumAvgSteps )
            {
                tAvgVec( 0 ) += tAvgCount * std::log10( tResidual );
                tAvgVec( 1 ) += std::log10( tResidual );
                ++tAvgCount;
            }
            else
            {
                real tOldGrad = tAvgGrad( 1 ) ;
                tAvgGrad = tAvgMat * tAvgVec ;

                if(  tAvgGrad( 1 ) < 4 )
                {
                    tOmega *= 0.9 * std::pow( tAvgGrad( 1 ) / tOldGrad, 0.25 );
                }

                tAvgVec.fill( 0.0 );
                tAvgCount = 0 ;
            }*/


            if( tKernel.is_master() )
            {
                string tAlgLabel = tMaxwell->algorithm() == SolverAlgorithm::Picard ?  " P " : " NR";

                std::cout << "    iteration:  " << tIter << tAlgLabel << " omega " << tOmega << " log10(epsilon): " << std::round( std::log10( tResidual ) * 100 ) * 0.01  << std::endl;
            }

            // switch algorithm
            if( tResidual < tEpsilonPicard &&  tMaxwell->algorithm() == SolverAlgorithm::Picard )
            {
                tMaxwell->set_algorithm( SolverAlgorithm::NewtonRaphson );
            }
            else if ( tResidual > tEpsilonPicard && tMaxwell->algorithm() == SolverAlgorithm::NewtonRaphson )
            {
                tMaxwell->set_algorithm( SolverAlgorithm::Picard );
            }

        }

        if( tMaxwell->type() == IwgType::MAXWELL_HA2D )
        {

            //tField2->compute_jacobian();   // todo: this is not necessary
            tField2->compute_rhs();
            tField2->solve();

            //tField3->compute_jacobian();   // todo: this is not necessary
            tField3->compute_rhs();
            tField3->solve() ;

            tField4->compute_rhs();
            tField4->solve() ;

            tMaxwell->compute_elementwise_current();
        }

        tMesh.set_time_step( tTimeCount++ );
        tMesh.save("mesh.e-s");
        tTime += tMaxwell->timestep() ;


    }


//------------------------------------------------------------------------------

    // close communicator
    return gComm.finalize();
}
