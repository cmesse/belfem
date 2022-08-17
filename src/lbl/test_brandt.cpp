//
// Created by Christian Messe on 10.01.22.
//


#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "en_IWGs.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG_Maxwell_L2_old.hpp"
#include "cl_IWG_Maxwell_HA_Tri3.hpp"
#include "cl_IWG_Maxwell_HPhi_Tri3.hpp"
#include "cl_IWG_Maxwell_HPhi_Tri6.hpp"
#include "cl_MaxwellMaterial.hpp"
#include "cl_FEM_MaxwellBoundaryConditionMagfield.hpp"
#include "fn_FEM_compute_brandt.hpp"
#include "cl_Timer.hpp"

using namespace belfem;
using namespace fem ;

Communicator gComm;
Logger       gLog( 0 );

int main( int    argc,
          char * argv[] )
          {
    // create communicator
    gComm = Communicator( argc, argv );

    // get the mesh
    Mesh tMesh( "brick.msh" );
    tMesh.scale_mesh( 0.001 );

    tMesh.create_edges( false, { 1 } );
    tMesh.create_faces( false, { 1 } );
    tMesh.flag_curved_elements();

    // user settings
    real tJc = 1000 ;
    real tEc = 0.001 ;
    real tN  = 40 ;
    real tHaHp = 0.25 ;
    real tDeltaTime = 0.0005 ;
    real tRampTime = tDeltaTime * 20 ;
    real tMaxTime  = tDeltaTime * 30 ;
    real tOmega = 0.05 ;

    // prepare the kernel
    KernelParameters tParams( tMesh );
    Kernel tKernel( &tParams );

    // get the element type
    ElementType tType = tMesh.block( 1 )->element_type() ;

    // create the IWGs
    IWG_Maxwell * tIWG = nullptr ;

    if( tType == ElementType::TRI3 )
    {
        tIWG = new IWG_Maxwell_HPhi_Tri3 ;
    }
    else if( tType == ElementType::TRI6 )
    {
        tIWG = new IWG_Maxwell_HPhi_Tri6 ;
    }
    IWG_Maxwell_L2 * tL2J = new IWG_Maxwell_L2( tType, Maxwell_L2_Mode::H2J );
    //IWG_Maxwell_L2 * tL2B = new IWG_Maxwell_L2( tType, Maxwell_L2_Mode::H2B );

    tIWG->set_blocks( { 1, 2 }, { DomainType::SuperConductor, DomainType::Air } );

    if( true )
    {
        tIWG->set_sidesets( { 1, 2, 3, 4, 5, 6, 7, 8 }, {
                DomainType::Boundary,
                DomainType::Boundary,
                DomainType::Boundary,
                DomainType::Boundary,
                DomainType::InterfaceScAir,
                DomainType::InterfaceScAir,
                DomainType::InterfaceScAir,
                DomainType::InterfaceScAir } );
    }
    else
    {
        tIWG->set_sidesets( { 5, 6, 7, 8 }, {
                DomainType::InterfaceScAir,
                DomainType::InterfaceScAir,
                DomainType::InterfaceScAir,
                DomainType::InterfaceScAir } );
    }
    // create the boundary condition
    MaxwellBoundaryConditionMagfield * tBC = new MaxwellBoundaryConditionMagfield( "Farfield" );
    //tBC->set_scaling( BoundaryConditionScaling::Constant, 0.01 );
    tBC->set_sidesets( { 1, 2, 3, 4 } );
    tIWG->add_boundary_condition( tBC );

    tL2J->set_block( 1, DomainType::SuperConductor );
    //tL2H->set_block( 2, DomainType::Air );

    // crate the fields
    DofManager * tMagfield = tKernel.create_field( tIWG );
    tMagfield->block( 1 )->set_domain_type( DomainType::SuperConductor );
    tMagfield->block( 2 )->set_domain_type( DomainType::Air );

    if( true )
    {
        tMagfield->sideset( 1 )->set_domain_type( DomainType::Boundary );
        tMagfield->sideset( 2 )->set_domain_type( DomainType::Boundary );
        tMagfield->sideset( 3 )->set_domain_type( DomainType::Boundary );
        tMagfield->sideset( 4 )->set_domain_type( DomainType::Boundary );
    }
    tMagfield->sideset( 5 )->set_domain_type( DomainType::InterfaceScAir );
    tMagfield->sideset( 6 )->set_domain_type( DomainType::InterfaceScAir );
    tMagfield->sideset( 7 )->set_domain_type( DomainType::InterfaceScAir );
    tMagfield->sideset( 8 )->set_domain_type( DomainType::InterfaceScAir );

    tIWG->set_timestep( tDeltaTime );

    tMagfield->set_solver( SolverType::UMFPACK );

    tBC->set_dof_manager( tMagfield );
    tBC->compute( 0.0 );

    // hack
    /*if( ! tMesh.field_exists( "bx") )
    {
        tMesh.create_field("bx");
        tMesh.create_field("by");
    }*/
    tMagfield->initialize();

    DofManager * tL2       = tKernel.create_field( tL2J );
    tL2->set_solver( SolverType::UMFPACK );
    tL2->initialize() ;

    // create the material
    MaxwellMaterial * tMat = new MaxwellMaterial( "superconductor" );
    tMat->set_rho_el_ej( tEc, tJc, tN );

    // set the material
    tMagfield->block( 1 )->set_material( tMat );
    tL2->block( 1 )->set_material( tMat );
    tL2->block( 1 )->set_domain_type( DomainType::SuperConductor );

    Vector< real > tB( 2, 0.0 );
    tB( 1 ) = compute_brandt( tMagfield, tHaHp );

    tB.print("Brandt");

    tBC->set( BoundaryConditionScaling::Sigmoid, tB, tRampTime, tDeltaTime );

    // begin timestep
    real & tTime = tMesh.time_stamp() ;
    uint & tTimeCount = tMesh.time_step() ;
    tTimeCount = 0 ;
    tTime = 0.0 ;

    while( tTime < tMaxTime )
    {
        // increment timestep
        tTime += tDeltaTime;
        tTimeCount++ ;

        // reset number of iterations
        tIWG->shift_fields();
        tIWG->compute_boundary_conditions( tTime );

        uint tIter = 0;

        // residual
        real tEpsilon = BELFEM_REAL_MAX;

        if ( tKernel.is_master() )
        {
          std::cout << std::endl
                    << "  --------------------------------------------------------------------------"
                    << std::endl;
          std::cout << "   time : " << tTime << " s " << std::endl;
          std::cout << "  --------------------------------------------------------------------------"
                    << std::endl;
        }

        Timer tTimer;

        tIWG->set_algorithm( SolverAlgorithm::Picard );
        tIWG->set_omega( tOmega );



        while ( tEpsilon > 1e-6 || tIter < 3 )
        {
          tBC->compute( tTime );
          tMagfield->compute_jacobian_and_rhs();
          tMagfield->solve();
          tEpsilon = tMagfield->residual( tIter++ );

          if ( tKernel.is_master() )
          {
              string tAlgLabel = tIWG->algorithm() == SolverAlgorithm::Picard ? " P " : " NR";

              std::cout << "    iteration:  " << tIter << tAlgLabel << " omega " << tIWG->omega()
                        << " log10(epsilon): " << std::round( std::log10( tEpsilon ) * 100 ) * 0.01
                        << std::endl;
          }
          BELFEM_ERROR( tIter < 5000, "too many iterations" );
        }
        tL2->compute_jacobian();
        tL2->compute_rhs() ;
        tL2->solve();
        //tL2->save_system("test.hdf5");
        if ( tKernel.is_master() )
        {
          gLog.message( 1, "    timestep completed in %4.2f seconds", ( float ) tTimer.stop() * 0.001 );
        }

        tMesh.save( "hphi.e-s");
        //exit( 0 );
    }


    // close communicator
    return gComm.finalize();
}