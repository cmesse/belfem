//
// Created by Christian Messe on 03.03.22.
//

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
    Mesh tMesh( "phitest.msh" );
    tMesh.scale_mesh( 0.001 );

    real tDeltaTime = 0.0001 ;
    real tMaxTime = 10 * tDeltaTime ;

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
    //IWG_Maxwell_L2 * tL2J = new IWG_Maxwell_L2( tType, Maxwell_L2_Mode::H2J );
    //IWG_Maxwell_L2 * tL2B = new IWG_Maxwell_L2( tType, Maxwell_L2_Mode::H2B );

    tIWG->set_blocks( { 1 }, { DomainType::Air } );

    // create the boundary condition
    MaxwellBoundaryConditionMagfield * tBC = new MaxwellBoundaryConditionMagfield( "Farfield" );
    tBC->set_scaling( BoundaryConditionScaling::Constant, 0.01 );
    tBC->set_sidesets( { 1, 2, 3, 4 } );
    tIWG->add_boundary_condition( tBC );

    //tL2H->set_block( 2, DomainType::Air );

    // crate the fields
    DofManager * tMagfield = tKernel.create_field( tIWG );
    tMagfield->block( 1 )->set_domain_type( DomainType::Air );

    tIWG->set_timestep( tDeltaTime );

    tMagfield->set_solver( SolverType::PARDISO );

    tBC->set_dof_manager( tMagfield );

    // hack
    /*if( ! tMesh.field_exists( "bx") )
    {
        tMesh.create_field("bx");
        tMesh.create_field("by");
    }*/
    tMagfield->initialize();

    Vector< real > tB( 2, 0 );
    tB(1) = 1.0 ;

    tBC->set( BoundaryConditionScaling::Constant, tB );

    // begin timestep
    real & tTime = tMesh.time_stamp() ;
    uint & tTimeCount = tMesh.time_step() ;
    tTime = 0.0 ;

    while( tTime < tMaxTime )
    {
        // increment timestep
        tTime += tDeltaTime;

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
        tIWG->set_omega( 0.1 );
        tBC->compute( tTime );
        while ( tEpsilon > 1e-6 )
        {

            tMagfield->compute_jacobian_and_rhs();
            //tMagfield->save_system( "matrix.hdf5" );
            tMagfield->solve();

            tEpsilon = tMagfield->residual( tIter++ );

            if ( tKernel.is_master() )
            {
                string tAlgLabel = tIWG->algorithm() == SolverAlgorithm::Picard ? " P " : " NR";

                std::cout << "    iteration:  " << tIter << tAlgLabel << " omega " << tIWG->omega()
                          << " log10(epsilon): " << std::round( std::log10( tEpsilon ) * 100 ) * 0.01
                          << std::endl;
            }
            BELFEM_ERROR( tIter < 1000, "too many iterations" );
        }

        if ( tKernel.is_master() )
        {
            gLog.message( 1, "    timestep completed in %4.2f seconds", ( float ) tTimer.stop() * 0.001 );
        }

        tMesh.save( "mesh.exo");
        tTimeCount++ ;
    }


    // close communicator
    return gComm.finalize();
}