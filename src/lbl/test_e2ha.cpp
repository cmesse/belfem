//
// Created by christian on 8/5/21.
//

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG_Maxwell_Old.hpp"
#include "fn_Mesh_compute_volume.hpp"
#include "cl_Element_Factory.hpp"
using namespace belfem ;
using namespace fem ;

Communicator gComm;
Logger       gLog( 3 );

//-----------------------------------------------------------------------------

Mesh *
create_1st_order_mesh()
{
    Mesh * aMesh = new Mesh();

    // - - - - - - - - - - - - - - - - -
    // STEP 1: create nodes on mesh
    // - - - - - - - - - - - - - - - - -

    // grab node container of mesh
    Cell< mesh::Node * > & tNodes = aMesh->nodes() ;
    tNodes.set_size( 6, nullptr );

    real tX ;
    real tY = 0.0 ;

    index_t tCount = 0 ;

    for( uint k=0; k<3; ++k )
    {
        tX = 0.0 ;
        for( uint i=0; i<2; ++i )
        {
            tNodes( tCount ) = new mesh::Node( tCount+1, tX, tY );
            ++tCount ;
            tX += 2.0 ;
        }
        tY += 2.0 ;
    }

    aMesh->set_number_of_dimensions( 2 );

    // - - - - - - - - - - - - - - - - -
    // STEP 2: superconducting block
    // - - - - - - - - - - - - - - - - -
    // crate an element factory
    mesh::ElementFactory tFactory ;


    mesh::Block * tBlock1 = new mesh::Block( 1, 1 );
    mesh::Element * tEl1 = tFactory.create_lagrange_element( ElementType::TRI3, 1 );
    tEl1->insert_node( tNodes( 3 ), 0 );
    tEl1->insert_node( tNodes( 5 ), 1 );
    tEl1->insert_node( tNodes( 4 ), 2 );
    tBlock1->insert_element( tEl1 );

    // - - - - - - - - - - - - - - - - -
    // STEP 3: coil block
    // - - - - - - - - - - - - - - - - -
    mesh::Block * tBlock2 = new mesh::Block( 2, 1 );
    mesh::Element * tEl2 = tFactory.create_lagrange_element( ElementType::TRI3, 2 );
    tEl2->insert_node( tNodes( 0 ), 0 );
    tEl2->insert_node( tNodes( 1 ), 1 );
    tEl2->insert_node( tNodes( 2 ), 2 );
    tBlock2->insert_element( tEl2 );

    // - - - - - - - - - - - - - - - - -
    // STEP 4: air block
    // - - - - - - - - - - - - - - - - -
    mesh::Block * tBlock3 = new mesh::Block( 3, 2 );
    mesh::Element * tEl3 = tFactory.create_lagrange_element( ElementType::TRI3, 3 );
    tEl3->insert_node( tNodes( 1 ), 0 );
    tEl3->insert_node( tNodes( 3 ), 1 );
    tEl3->insert_node( tNodes( 2 ), 2 );
    tBlock3->insert_element( tEl3 );

    mesh::Element * tEl4 = tFactory.create_lagrange_element( ElementType::TRI3, 4 );
    tEl4->insert_node( tNodes( 2 ), 0 );
    tEl4->insert_node( tNodes( 3 ), 1 );
    tEl4->insert_node( tNodes( 4 ), 2 );
    tBlock3->insert_element( tEl4 );

    // - - - - - - - - - - - - - - - - -
    // STEP 5: Interface
    // - - - - - - - - - - - - - - - - -

    mesh::Element * tEl5 = tFactory.create_lagrange_element( ElementType::LINE2, 5 );
    tEl5->insert_node( tNodes( 3 ), 0 );
    tEl5->insert_node( tNodes( 4 ), 1 );

    mesh::SideSet * tInterface = new mesh::SideSet( 1, 1 );
    mesh::Facet * tFacet1 = new mesh::Facet( tEl5 );

    tInterface->insert_facet( tFacet1 );

    // - - - - - - - - - - - - - - - - - - -
    // STEP 6: finalize mesh
    // - - - - - - - - - - - - - - - - - -

    Cell< mesh::Block * > & tBlocks = aMesh->blocks() ;
    tBlocks.set_size( 3, nullptr );
    tBlocks( 0 ) = tBlock1 ;
    tBlocks( 1 ) = tBlock2 ;
    tBlocks( 2 ) = tBlock3 ;

    aMesh->sidesets().set_size( 1, tInterface );
    aMesh->finalize() ;

    return aMesh ;
}

Mesh *
create_2nd_order_mesh()
{
    Mesh * aMesh = new Mesh();

    // - - - - - - - - - - - - - - - - -
    // STEP 1: create nodes on mesh
    // - - - - - - - - - - - - - - - - -

    // grab node container of mesh
    Cell< mesh::Node * > & tNodes = aMesh->nodes() ;
    tNodes.set_size( 15, nullptr );

    real tX ;
    real tY = 0.0 ;

    index_t tCount = 0 ;

    for( uint k=0; k<5; ++k )
    {
        tX = 0.0 ;
        for( uint i=0; i<3; ++i )
        {
            tNodes( tCount ) = new mesh::Node( tCount+1, tX, tY );
            ++tCount ;
            tX += 1.0 ;
        }
        tY += 1.0 ;
    }

    aMesh->set_number_of_dimensions( 2 );

    // - - - - - - - - - - - - - - - - -
    // STEP 2: superconducting block
    // - - - - - - - - - - - - - - - - -
    // crate an element factory
    mesh::ElementFactory tFactory ;


    mesh::Block * tBlock1 = new mesh::Block( 1, 1 );
    mesh::Element * tEl1 = tFactory.create_lagrange_element( ElementType::TRI6, 1 );
    tEl1->insert_node( tNodes( 12 ), 0 );
    tEl1->insert_node( tNodes( 8 ), 1 );
    tEl1->insert_node( tNodes( 14 ), 2 );
    tEl1->insert_node( tNodes( 10 ), 3 );
    tEl1->insert_node( tNodes( 11 ), 4 );
    tEl1->insert_node( tNodes( 13 ), 5 );
    tBlock1->insert_element( tEl1 );

    // - - - - - - - - - - - - - - - - -
    // STEP 3: coil block
    // - - - - - - - - - - - - - - - - -
    mesh::Block * tBlock2 = new mesh::Block( 2, 1 );
    mesh::Element * tEl2 = tFactory.create_lagrange_element( ElementType::TRI6, 2 );
    tEl2->insert_node( tNodes( 0 ), 0 );
    tEl2->insert_node( tNodes( 2 ), 1 );
    tEl2->insert_node( tNodes( 6 ), 2 );
    tEl2->insert_node( tNodes( 1 ), 3 );
    tEl2->insert_node( tNodes( 4 ), 4 );
    tEl2->insert_node( tNodes( 3 ), 5 );
    tBlock2->insert_element( tEl2 );

    // - - - - - - - - - - - - - - - - -
    // STEP 4: air block
    // - - - - - - - - - - - - - - - - -
    mesh::Block * tBlock3 = new mesh::Block( 3, 2 );
    mesh::Element * tEl3 = tFactory.create_lagrange_element( ElementType::TRI6, 3 );
    tEl3->insert_node( tNodes( 6 ), 0 );
    tEl3->insert_node( tNodes( 8 ), 1 );
    tEl3->insert_node( tNodes( 12 ), 2 );
    tEl3->insert_node( tNodes( 7 ), 3 );
    tEl3->insert_node( tNodes( 10 ), 4 );
    tEl3->insert_node( tNodes( 9 ), 5 );
    tBlock3->insert_element( tEl3 );

    mesh::Element * tEl4 = tFactory.create_lagrange_element( ElementType::TRI6, 4 );
    tEl4->insert_node( tNodes( 6 ), 0 );
    tEl4->insert_node( tNodes( 2 ), 1 );
    tEl4->insert_node( tNodes( 8 ), 2 );
    tEl4->insert_node( tNodes( 4 ), 3 );
    tEl4->insert_node( tNodes( 5 ), 4 );
    tEl4->insert_node( tNodes( 8 ), 5 );
    tBlock3->insert_element( tEl4 );

    // - - - - - - - - - - - - - - - - -
    // STEP 5: Interface
    // - - - - - - - - - - - - - - - - -

    mesh::Element * tEl5 = tFactory.create_lagrange_element( ElementType::LINE3, 5 );
    tEl5->insert_node( tNodes( 12 ), 0 );
    tEl5->insert_node( tNodes( 8 ), 1 );
    tEl5->insert_node( tNodes( 10 ), 2 );

    mesh::SideSet * tInterface = new mesh::SideSet( 1, 1 );
    mesh::Facet * tFacet1 = new mesh::Facet( tEl5 );

    tInterface->insert_facet( tFacet1 );

    // - - - - - - - - - - - - - - - - - - -
    // STEP 6: finalize mesh
    // - - - - - - - - - - - - - - - - - -

    Cell< mesh::Block * > & tBlocks = aMesh->blocks() ;
    tBlocks.set_size( 3, nullptr );
    tBlocks( 0 ) = tBlock1 ;
    tBlocks( 1 ) = tBlock2 ;
    tBlocks( 2 ) = tBlock3 ;

    aMesh->sidesets().set_size( 1, tInterface );
    aMesh->finalize() ;

    return aMesh ;
}

//------------------------------------------------------------------------------

void
write_test_data_to_mesh( Mesh * aMesh )
{
    // write dummy data into mesh
    Vector< real > & tHx = aMesh->create_field("hx");
    Vector< real > & tHy = aMesh->create_field("hy");
    tHx(  8 ) = 50 ;
    tHx( 10 ) = 25 ;
    tHx( 11 ) = 62.5 ;
    tHx( 13 ) = 25 ;
    tHx( 14 ) = 75 ;
    tHy = 2.0 * tHx ;

    Vector< real > & tH0x = aMesh->create_field("h0x");
    Vector< real > & tH0y = aMesh->create_field("h0y");
    tH0x = 0.75 * tHx ;
    tH0y = 0.75 * tHy ;
}

//------------------------------------------------------------------------------
int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    Mesh * tMesh = create_1st_order_mesh() ;

    //write_test_data_to_mesh( tMesh );
    tMesh->create_edges( true );

    //------------------------------------------------------------------------------

    KernelParameters tParams( tMesh );

    tParams.enforce_linear( Vector< uint > ( 1, 1 ) );

    Kernel tKernel( &tParams );


    IWG_Maxwell * tMaxwell = reinterpret_cast< IWG_Maxwell  *  >
            ( tKernel.create_equation( IwgType::MAXWELL_HPHI2D ) );

    //------------------------------------------------------------------------------
    // User Settings
    //------------------------------------------------------------------------------

    Vector< id_t > tPlus( 1, 3 );
    //Vector< id_t > tMinus( 1, 3 );

    // block identification
    tMaxwell->set_block( 1, DomainType::Conductor );
    tMaxwell->set_block( 3, DomainType::Coil );
    tMaxwell->set_block( 2, DomainType::Air );
    tMaxwell->set_sideset( 1,  DomainType::Interface );

    // settings for power law
    tMaxwell->powerlaw()->set_coefficients( 1e-5, 150, 5 );
    tMaxwell->set_euler_method( EulerMethod::BackwardImplicit );

    // current ramp
    //real tImax     = 150 ;
    //real tRamptime = 10 ;

    //real tDeltaTime = 0.05 ;
    //real tMaxTime = tDeltaTime ;
    //real tRelax = 0.5 ;

    DofManager * tField = tKernel.create_field( tMaxwell );

    // select solver library
    tField->set_solver( SolverType::UMFPACK );
    //tField->solver()->set_petsc(Preconditioner::ASM,
    //                            KrylovMethod::GMRES,
    //                            1e-8 );


    // get sideset
    SideSet * tSideSet = tField->sideset( 1 );
    Element * tElement = tSideSet->element( 5 );

    tField->initialize() ;
    tMesh->save("mesh.exo");

    Matrix< real > tK( 7, 7 );
    Vector< real > tB( 7 );
    tMaxwell->link_to_group( tSideSet );
    tMaxwell->print_dofs( tElement );
    tMaxwell->compute_jacobian_and_rhs( tElement, tK, tB );

    tK.print("K");

    exit( 0 );

    tMaxwell->set_current_density( Vector<id_t>(1, 2 ) , 100 );
    tField->compute_volume_loads( Vector<id_t>(1, 2 ) );
    tField->compute_jacobian_and_rhs();

    std::cout << "===DOFS=== " << std::endl ;

    for( Dof * tDof : tField->dofs() )
    {
        std::cout << tDof->index() << " " << tDof->id() << " " << tDof->is_node() << " " << tDof->mesh_basis()->id() << " " << tDof->dof_index_on_field() << std::endl ;
    }

    tField->solve();
    //exit( 0 );

    std::cout << "R2 " << tField->residual( 0 ) << std::endl ;

    /*tMesh->save("mesh.exo");

    // special setting for sc block
    //tField->block( 1 )->set_integration_order( 20 );
    tMaxwell->set_euler_method(EulerMethod::BackwardImplicit );

    tMaxwell->set_timestep( tDeltaTime );
    tMaxwell->set_omega( tRelax );

    //------------------------------------------------------------------------------


    // cross section of positive coil
    real tAplus = mesh::compute_volume( tKernel.mesh(), tPlus );

    // cross section of negative coil
    //real tAminus = mesh::compute_volume( &tMesh, tMinus );


    //------------------------------------------------------------------------------



    uint tTimeCount = 1 ;
    real & tTime = tMesh->time_stamp() ;


    tTime = tMaxwell->timestep() ;

    //uint tPow = 1 ;

    while( tTime < 2*tMaxTime )
    {
        // shift fields
        if( tKernel.is_master() )
        {
            Vector< real > & tAz = tMesh->field_data("az");
            Vector< real > & tHx = tMesh->field_data("hx");
            Vector< real > & tHy = tMesh->field_data("hy");

            Vector< real > & tA0z = tMesh->field_data("a0z");
            Vector< real > & tH0x = tMesh->field_data("h0x");
            Vector< real > & tH0y = tMesh->field_data("h0y");

            tA0z = tAz ;
            tH0x = tHx ;
            tH0y = tHy ;
        }
        comm_barrier() ;
        tField->distribute_fields({"a0z", "h0x","h0y"});

        uint tIter = 0 ;

        tMesh->set_time_step( tTimeCount );



        // compute imposed current density
        real tI = tTime < tRamptime ? tTime / tRamptime * tImax : tImax ;

        if( tKernel.is_master() )
        {
            std::cout << std::endl  << "  --------------------------------------------------------------------------"
            << std::endl;
            std::cout << "   time : " << tTime << " s " << ",   I : " << tI << " A,   n :"<< tMaxwell->power_law_exponent() <<std::endl ;
            std::cout << "  --------------------------------------------------------------------------" << std::endl  ;
        }



        tMaxwell->set_current_density( tPlus, tI / tAplus );
        //tMaxwell->set_current_density( tMinus, -tI / tAminus );

        // compute rhs
        tField->compute_volume_loads( { 2, 3 } );

        real tResidual = BELFEM_REAL_MAX ;
        while( tResidual > 1e-5 )
        {
            if( tIter > 200 )
            {
                tMaxwell->set_omega( 0.05 );
            }
            else if( tIter > 100 )
            {
                tMaxwell->set_omega( 0.1 );
            }
            else if( tIter > 50 )
            {
                tMaxwell->set_omega( 0.2 );
            }
            else
            {
                tMaxwell->set_omega( tRelax );
            }

            if( tKernel.is_master() && comm_size() == 1 )
            {
                tMaxwell->print_dofs( tField->block( 1 )->element( 1 ) );
            }
            if( ! tKernel.is_master() )
            {
                tMaxwell->print_dofs( tField->block( 1 )->element( 1 ) );
            }


            tField->compute_jacobian_and_rhs();
            tField->solve();

            if( tKernel.is_master() )
            {
                tField->jacobian()->print();
            }
            // return gComm.finalize();
            tResidual = tField->residual( tIter );
            if( tKernel.is_master() )
            {
                std::cout << "    iteration:  " << tIter++ << " log10(epsilon): " << std::round( std::log10( tResidual ) * 10 ) * 0.1 << std::endl;
            }
        }

        tMesh->set_time_step( tTimeCount++ );
        tMesh->save("mesh.e-s");
        tTime += tMaxwell->timestep() ;

        //if( tPow < 20 )
        //{
        //    tMaxwell->set_power_law_exponent( tPow++ );
        // }
    } */


    delete tMesh ;

    return gComm.finalize();
}
