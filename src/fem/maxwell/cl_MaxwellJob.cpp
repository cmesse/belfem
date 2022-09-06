//
// Created by Christian Messe on 13.01.22.
//
#include "commtools.hpp"
#include "cl_MaxwellJob.hpp"
#include "assert.hpp"
#include "cl_HDF5.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_IWG_Maxwell_HPhi_Tri3.hpp"
#include "cl_IWG_Maxwell_HPhi_Tri6.hpp"
//#include "cl_IWG_Maxwell_HPhi_Tet4.hpp"
//#include "cl_IWG_Maxwell_HPhi_Tet10.hpp"
#include "cl_MaxwellMaterial.hpp"
#include "fn_r2.hpp"

namespace belfem
{
    namespace fem
    {
        MaxwellJob::MaxwellJob() :
        mMyRank( comm_rank() ),
        mCommSize( comm_size() )
        {

        }

//------------------------------------------------------------------------------

        MaxwellJob::MaxwellJob( Mesh * aMesh, Kernel * aKernel ) :
        mMyRank( comm_rank() ),
        mCommSize( comm_size() ),
        mMesh( aMesh ),
        mKernel( aKernel )
        {

        }

//------------------------------------------------------------------------------

        MaxwellJob::~MaxwellJob()
        {

            for( Material * tMaterial : mMaterials )
            {
                // check if material is owned by the kernel
                if( ! tMaterial->is_flagged() )
                {
                    delete tMaterial ;
                }
            }

            if( mKernel != nullptr )
            {
                // we must release the parameters since we destroy them in
                // this object
                mKernel->claim_parameter_ownership( false );
                delete mKernel ;
            }

            if( mParams != nullptr ) delete mParams ;
            if( mMesh != nullptr )   delete mMesh ;

        }

//------------------------------------------------------------------------------

        void
        MaxwellJob::initialize_test( const IwgType aIwgType, const bool aCurved )
        {
            // create the mesh and the equation
            switch( aIwgType )
            {
                case( IwgType::MAXWELL_HPHI_TRI3 ) :
                {
                    mLabel = "hphi_tri3";

                    this->create_test_mesh_tri3() ;

                    mEquation = new IWG_Maxwell_HPhi_Tri3 ;

                    break ;
                }
                case( IwgType::MAXWELL_HPHI_TRI6 ) :
                {

                    mLabel = aCurved ? "hphi_tri6_curved" : "hphi_tri6_straight" ;

                    this->create_test_mesh_tri6( aCurved ) ;

                    mEquation = new IWG_Maxwell_HPhi_Tri6 ;

                    break ;
                }
                /*case( IwgType::MAXWELL_HPHI_TET4 ) :
                {
                    mLabel = "hphi_tet4";

                    this->create_test_mesh_tet4() ;

                    mEquation = new IWG_Maxwell_HPhi_Tet4 ;

                    break ;
                }
                case( IwgType::MAXWELL_HPHI_TET10 ) :
                {
                    mLabel = aCurved ? "hphi_tet10_curved" : "hphi_tet10_straight" ;

                    this->create_test_mesh_tet10( aCurved ) ;

                    mEquation = new IWG_Maxwell_HPhi_Tet10 ;

                    break ;
                }*/
                default :
                {
                    BELFEM_ERROR( false, "No test implemented for this type" );
                }
            }

            // define the topology
            Vector< id_t > tBlockIDs = { 1, 2, 3, 4 };
            Cell< DomainType > tBlockTypes = {
                    DomainType::SuperConductor,
                    DomainType::FerroMagnetic,
                    DomainType::Coil,
                    DomainType::Air };

            mEquation->set_blocks( tBlockIDs, tBlockTypes );


            Vector< id_t > tSideSetIDs = { 1, 2, 3 };
            Cell< DomainType > tSideSetTypes = {
                    DomainType::InterfaceScAir,
                    DomainType::InterfaceScFm,
                    DomainType::InterfaceFmAir };
            mEquation->set_sidesets( tSideSetIDs, tSideSetTypes );



            // create the kernel parameters
            mParams = new fem::KernelParameters( mMesh );

            // create the kernel
            mKernel = new fem::Kernel( mParams );

            // make sure that Kernel owns and destroys the parameters
            mKernel->claim_parameter_ownership();

            // add equation to kernel
            mKernel->add_equation( mEquation );

            // create dof manager for magnetic field
            mMagneticField = mKernel->create_field( mEquation );

            // create data arrays on mesh
            mMagneticField->create_fields( mEquation );

            // more settings
            mEquation->set_field( mMagneticField );
            mEquation->set_timestep( 0.0001 );
            mEquation->set_euler_method( 1.0 ); // <-- must be 1.0

            mMaterials.set_size( 2, nullptr );

            // create materials
            mMaterials( 0 ) = new MaxwellMaterial( "superconductor" );

            // set an arbitrary prime number so that it is not zero
            reinterpret_cast< MaxwellMaterial * >( mMaterials( 0 ))->set_rho_el_const( 0.1 );
            reinterpret_cast< MaxwellMaterial * >( mMaterials( 0 ))->set_mu_r( 1000 );

            // set the material type for the superconductor ( it will be deleted by the Dof Manager)
            mMagneticField->block( 1 )->set_material( mMaterials( 0 ) );

            // create the material type for the ferro
            mMaterials( 1 ) = new MaxwellMaterial( "ferro" );

            // set an arbitrary prime number so that it is not zero
            reinterpret_cast< MaxwellMaterial * >( mMaterials( 1 ))->set_nu_s( 3.0 * constant::nu0 );

            // set the material type for the ferro
            mMagneticField->block( 2 )->set_material( mMaterials( 1 ) );

            // we must also set the block and sideset types.
            mMagneticField->block( 1 )->set_domain_type( DomainType::SuperConductor );
            mMagneticField->block( 2 )->set_domain_type( DomainType::FerroMagnetic );
            mMagneticField->block( 3 )->set_domain_type( DomainType::Coil );
            mMagneticField->block( 4 )->set_domain_type( DomainType::Air );

            mMagneticField->sideset( 1 )->set_domain_type( DomainType::InterfaceScAir );
            mMagneticField->sideset( 2 )->set_domain_type( DomainType::InterfaceScFm );
            mMagneticField->sideset( 3 )->set_domain_type( DomainType::InterfaceFmAir );

        }

//------------------------------------------------------------------------------

        void
        MaxwellJob::run_test(  const string & aReferenceFilePath, Vector< real > & aResult, const bool aPrint )
        {

            // container with test results
            aResult.set_size( 6, BELFEM_QUIET_NAN );

            // load the reference file
            HDF5 tFile( aReferenceFilePath, FileMode::OPEN_RDONLY );

            // select the group
            tFile.select_group( mLabel );

            // matrices
            Matrix < real > tJacobian ;
            Vector < real > tRHS ;

            // reference matrices
            Matrix< real > tK ;
            Matrix< real > tM ;
            Matrix< real > tJ ;

            fem::Element * tElement = nullptr ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // superconducting block
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // check superconductor
            mEquation->link_to_group( mMagneticField->block( 1 ) );

            // get first element on block
            tElement = mMagneticField->block( 1 )->elements()( 0 );

            // print element connectivity
            if ( aPrint )  mEquation->print_dofs( tElement );

            tJacobian.set_size( tElement->number_of_dofs(), tElement->number_of_dofs() );
            tRHS.set_size( tElement->number_of_dofs() );

            mEquation->compute_jacobian_and_rhs( tElement, tJacobian, tRHS );

            // load the reference solution
            tFile.load_data( "K_sc", tK );
            tFile.load_data( "M_sc", tM );
            tJ = tM + mEquation->timestep() * tK ;


            // compare results
            aResult( 0 ) = r2( tJacobian, tJ );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // ferromagnetic block
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // link to ferro block
            mEquation->link_to_group( mMagneticField->block( 2 ) );

            // get first element on block
            tElement = mMagneticField->block( 2 )->elements()( 0 );

            // print element connectivity
            if ( aPrint ) mEquation->print_dofs( tElement );

            tJacobian.set_size( tElement->number_of_dofs(), tElement->number_of_dofs() );
            tRHS.set_size( tElement->number_of_dofs() );

            mEquation->compute_jacobian_and_rhs( tElement, tJacobian, tRHS );

            // load the reference solution
            tFile.load_data( "K_fm", tK );
            tFile.load_data( "M_fm", tM );
            tJ = tM + mEquation->timestep() * tK ;

            // compare results
            aResult( 1 ) = r2( tJacobian, tJ );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // air block
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // check air
            mEquation->link_to_group( mMagneticField->block( 4 ) );

            // get first element on block
            tElement = mMagneticField->block( 4 )->elements()( 0 );

            if ( aPrint ) mEquation->print_dofs( tElement );

            tJacobian.set_size( tElement->number_of_dofs(), tElement->number_of_dofs() );
            tRHS.set_size( tElement->number_of_dofs() );

            mEquation->compute_jacobian_and_rhs( tElement, tJacobian, tRHS );

            // load the reference solution
            tFile.load_data( "K_air", tK );
            tFile.load_data( "M_air", tM );
            tJ = tM + mEquation->timestep() * tK ;

            aResult( 2 ) = r2( tJacobian, tJ );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // interface sc-air
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // check air
            mEquation->link_to_group( mMagneticField->sideset( 1 ) );
            tElement = mMagneticField->sideset( 1 )->elements()( 0 );

            // print element connectivity
            if ( aPrint ) mEquation->print_dofs( tElement );

            tJacobian.set_size( tElement->number_of_dofs(), tElement->number_of_dofs() );
            tRHS.set_size( tElement->number_of_dofs() );

            mEquation->compute_jacobian_and_rhs( tElement, tJacobian, tRHS );

            // load reference solution
            tFile.load_data( "K_scair", tK );
            tFile.load_data( "M_scair", tM );
            tJ = tM + mEquation->timestep() * tK ;

            aResult( 3 ) = r2( tJacobian, tJ );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // interface sc-fm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


            // check air
            mEquation->link_to_group( mMagneticField->sideset( 2 ) );
            tElement = mMagneticField->sideset( 2 )->elements()( 0 );

            // print element connectivity
            if ( aPrint ) mEquation->print_dofs( tElement );

            tJacobian.set_size( tElement->number_of_dofs(), tElement->number_of_dofs() );
            tRHS.set_size( tElement->number_of_dofs() );

            mEquation->compute_jacobian_and_rhs( tElement, tJacobian, tRHS );


            // load reference solution
            tFile.load_data( "K_scfm", tK );
            tFile.load_data( "M_scfm", tM );

            tJ = tM + mEquation->timestep() * tK ;

            aResult( 4 ) = r2( tJacobian, tJ );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // interface fm-air
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mEquation->link_to_group( mMagneticField->sideset( 3 ) );
            tElement = mMagneticField->sideset( 3 )->elements()( 0 );

            // print element connectivity
            if ( aPrint ) mEquation->print_dofs( tElement );

            tJacobian.set_size( tElement->number_of_dofs(), tElement->number_of_dofs() );
            tRHS.set_size( tElement->number_of_dofs() );

            mEquation->compute_jacobian_and_rhs( tElement, tJacobian, tRHS );

            // load reference solution
            tFile.load_data( "K_fmair", tK );
            tFile.load_data( "M_fmair", tM );
            tJ = tM + mEquation->timestep() * tK ;

            aResult( 5 ) = r2( tJacobian, tJ );

            if ( aPrint )  aResult.print("result");
        }

//------------------------------------------------------------------------------
// private
//------------------------------------------------------------------------------

        void
        MaxwellJob::create_test_mesh_tri3( )
        {

            // make a new mesh
            mMesh = new Mesh( 0 );
            mMesh->set_number_of_dimensions( 2 );

            // get the nodes
            Cell< mesh::Node * > & tNodes = mMesh->nodes() ;

            // define coordinates
            Vector< real > tX = { -1, 1, 1, -1, 1, -1, 0  };
            Vector< real > tY = { -1, -1, 1, 1, 0, 0, 0 };

            tX *= 0.5 ;
            tY *= 0.5 ;

            index_t tNumNodes = 7 ;

            tNodes.set_size( tNumNodes, nullptr );

            // create nodes
            for( uint k=0; k<tNumNodes; ++k )
            {
                tNodes( k ) = new mesh::Node( k+1, tX( k ), tY( k ) );
            }

            // create the element factory

            mesh::ElementFactory tFactory ;

            // container for new element
            mesh::Element * tElement ;

            id_t tID = 1 ;

            Cell< mesh::Block * > & tBlocks = mMesh->blocks() ;

            // - - - - - - - - - - superconductor block
            mesh::Block * tScBlock = new mesh::Block( 1,  0 );
            tBlocks.push( tScBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++ );
            tElement->insert_node( tNodes(  3 ), 0 );
            tElement->insert_node( tNodes(  5 ), 1 );
            tElement->insert_node( tNodes(  6 ), 2 );
            tScBlock->elements().push( tElement );

            // - - - - - - - - - - ferro block
            mesh::Block * tFerroBlock = new mesh::Block( 2,  0 );
            tBlocks.push( tFerroBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++ );
            tElement->insert_node( tNodes(  6 ), 0 );
            tElement->insert_node( tNodes(  5 ), 1 );
            tElement->insert_node( tNodes(  0 ), 2 );
            tFerroBlock->elements().push( tElement );

            // - - - - - - - - - - coil block
            mesh::Block * tCoilBlock = new mesh::Block( 3,  0 );
            tBlocks.push( tCoilBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++ );
            tElement->insert_node( tNodes(  6 ), 0 );
            tElement->insert_node( tNodes(  1 ), 1 );
            tElement->insert_node( tNodes(  4 ), 2 );
            tCoilBlock->elements().push( tElement );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++ );
            tElement->insert_node( tNodes(  6 ), 0 );
            tElement->insert_node( tNodes(  4 ), 1 );
            tElement->insert_node( tNodes(  2 ), 2 );
            tCoilBlock->elements().push( tElement );

            // - - - - - - - - - - air block
            mesh::Block * tAirBlock = new mesh::Block( 4,  0 );
            tBlocks.push( tAirBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++ );
            tElement->insert_node( tNodes(  3 ), 0 );
            tElement->insert_node( tNodes(  6 ), 1 );
            tElement->insert_node( tNodes(  2 ), 2 );
            tAirBlock->elements().push( tElement );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++ );
            tElement->insert_node( tNodes(  1 ), 0 );
            tElement->insert_node( tNodes(  6 ), 1 );
            tElement->insert_node( tNodes(  0 ), 2 );
            tAirBlock->elements().push( tElement );


            // --- sidesets
            Cell< mesh::SideSet * > & tSideSets = mMesh->sidesets() ;

            // interface sc/air
            mesh::SideSet * tScAir = new mesh::SideSet( 1, 0 );
            tSideSets.push( tScAir );

            tElement = tFactory.create_lagrange_element( ElementType::LINE2 , tID++ );
            tElement->insert_node( tNodes( 6 ), 0 );
            tElement->insert_node( tNodes( 3 ), 1 );
            tScAir->facets().push( new mesh::Facet( tElement ) );

            // interface sc/ferro
            mesh::SideSet * tScFerro = new mesh::SideSet( 2, 0 );
            tSideSets.push( tScFerro );

            tElement = tFactory.create_lagrange_element( ElementType::LINE2 , tID++  );
            tElement->insert_node( tNodes( 5 ), 0 );
            tElement->insert_node( tNodes( 6 ), 1 );
            tScFerro->facets().push( new mesh::Facet( tElement ) );


            // interface ferro/air
            mesh::SideSet * tFerroAir = new mesh::SideSet( 3, 0 );
            tSideSets.push( tFerroAir );

            tElement = tFactory.create_lagrange_element( ElementType::LINE2 , tID++ );
            tElement->insert_node( tNodes( 0 ), 0 );
            tElement->insert_node( tNodes( 6 ), 1 );
            tFerroAir->facets().push( new mesh::Facet( tElement ) );

            mMesh->finalize();
            mMesh->create_edges( false, { 1 }, { 1, 2 } );
        }

//------------------------------------------------------------------------------

        void
        MaxwellJob::create_test_mesh_tri6( const bool aCurved )
        {

            // make a new mesh
            mMesh = new Mesh( 0 );
            mMesh->set_number_of_dimensions( 2 );

            // get the nodes
            Cell< mesh::Node * > & tNodes = mMesh->nodes() ;

            // define coordinates
            Vector< real > tX = { -1, 1, 1, -1, 1, -1, 0, 0, 1, 1, 0, -1, -1, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5 };
            Vector< real > tY = { -1, -1, 1, 1, 0, 0, 0, -1, -0.5, 0.5, 1, 0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0, 0 };

            // shift points if mesh is curved
            if( aCurved )
            {
                tX( 13 ) = -0.382 ;
                tY( 13 ) = -0.618 ;
                tX( 16 ) = -0.382;
                tY( 16 ) =  0.618 ;
            }

            tX *= 0.5 ;
            tY *= 0.5 ;

            index_t tNumNodes = 19 ;

            tNodes.set_size( tNumNodes, nullptr );

            // create nodes
            for( uint k=0; k<tNumNodes; ++k )
            {
                tNodes( k ) = new mesh::Node( k+1, tX( k ), tY( k ) );
            }

            // create the element factory
            mesh::ElementFactory tFactory ;

            // container for new element
            mesh::Element * tElement ;

            id_t tID = 1 ;

            Cell< mesh::Block * > & tBlocks = mMesh->blocks() ;

            // - - - - - - - - - - superconductor block
            mesh::Block * tScBlock = new mesh::Block( 1,  0 );
            tBlocks.push( tScBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++ );
            tElement->insert_node( tNodes(   3 ), 0 );
            tElement->insert_node( tNodes(   5 ), 1 );
            tElement->insert_node( tNodes(   6 ), 2 );
            tElement->insert_node( tNodes(  11 ), 3 );
            tElement->insert_node( tNodes(  17 ), 4 );
            tElement->insert_node( tNodes(  16 ), 5 );
            tScBlock->elements().push( tElement );

            // - - - - - - - - - - ferro block
            mesh::Block * tFerroBlock = new mesh::Block( 2,  0 );
            tBlocks.push( tFerroBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++ );
            tElement->insert_node( tNodes(   6 ), 0 );
            tElement->insert_node( tNodes(   5 ), 1 );
            tElement->insert_node( tNodes(   0 ), 2 );
            tElement->insert_node( tNodes(  17 ), 3 );
            tElement->insert_node( tNodes(  12 ), 4 );
            tElement->insert_node( tNodes(  13 ), 5 );
            tFerroBlock->elements().push( tElement );

            // - - - - - - - - - - coil block
            mesh::Block * tCoilBlock = new mesh::Block( 3,  0 );
            tBlocks.push( tCoilBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++ );
            tElement->insert_node( tNodes(   6 ), 0 );
            tElement->insert_node( tNodes(   1 ), 1 );
            tElement->insert_node( tNodes(   4 ), 2 );
            tElement->insert_node( tNodes(  14 ), 3 );
            tElement->insert_node( tNodes(  18 ), 4 );
            tElement->insert_node( tNodes(  18 ), 5 );
            tCoilBlock->elements().push( tElement );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++ );
            tElement->insert_node( tNodes(   6 ), 0 );
            tElement->insert_node( tNodes(   4 ), 1 );
            tElement->insert_node( tNodes(   2 ), 2 );
            tElement->insert_node( tNodes(  18 ), 3 );
            tElement->insert_node( tNodes(  9 ), 4 );
            tElement->insert_node( tNodes(  15 ), 5 );
            tCoilBlock->elements().push( tElement );

            // - - - - - - - - - - air block
            mesh::Block * tAirBlock = new mesh::Block( 4,  0 );
            tBlocks.push( tAirBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++ );
            tElement->insert_node( tNodes(   3 ), 0 );
            tElement->insert_node( tNodes(   6 ), 1 );
            tElement->insert_node( tNodes(   2 ), 2 );
            tElement->insert_node( tNodes(  16 ), 3 );
            tElement->insert_node( tNodes(  15 ), 4 );
            tElement->insert_node( tNodes(  10 ), 5 );
            tAirBlock->elements().push( tElement );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++ );
            tElement->insert_node( tNodes(   1 ), 0 );
            tElement->insert_node( tNodes(   6 ), 1 );
            tElement->insert_node( tNodes(   0 ), 2 );
            tElement->insert_node( tNodes(  14 ), 3 );
            tElement->insert_node( tNodes(  13 ), 4 );
            tElement->insert_node( tNodes(   7 ), 5 );
            tAirBlock->elements().push( tElement );


            // --- sidesets
            Cell< mesh::SideSet * > & tSideSets = mMesh->sidesets() ;

            // interface sc/air
            mesh::SideSet * tScAir = new mesh::SideSet( 1, 0 );
            tSideSets.push( tScAir );

            tElement = tFactory.create_lagrange_element( ElementType::LINE3 , tID++  );
            tElement->insert_node( tNodes(  6 ), 0 );
            tElement->insert_node( tNodes(  3 ), 1 );
            tElement->insert_node( tNodes( 16 ), 2 );
            tScAir->facets().push( new mesh::Facet( tElement ) );

            // interface sc/ferro
            mesh::SideSet * tScFerro = new mesh::SideSet( 2, 0 );
            tSideSets.push( tScFerro );

            tElement = tFactory.create_lagrange_element( ElementType::LINE3 , tID++  );
            tElement->insert_node( tNodes(  5 ), 0 );
            tElement->insert_node( tNodes(  6 ), 1 );
            tElement->insert_node( tNodes( 17 ), 2 );
            tScFerro->facets().push( new mesh::Facet( tElement ) );


            // interface ferro/air
            mesh::SideSet * tFerroAir = new mesh::SideSet( 3, 0 );
            tSideSets.push( tFerroAir );

            tElement = tFactory.create_lagrange_element( ElementType::LINE3 , tID++  );
            tElement->insert_node( tNodes(  0 ), 0 );
            tElement->insert_node( tNodes(  6 ), 1 );
            tElement->insert_node( tNodes( 13 ), 2 );
            tFerroAir->facets().push( new mesh::Facet( tElement ) );

            mMesh->finalize();
            mMesh->create_edges( false, { 1 }, { 1, 2 } );
            mMesh->create_faces( false, { 1 }, { 1, 2 } );
            mMesh->flag_curved_elements();
        }

//------------------------------------------------------------------------------

        void
        MaxwellJob::create_test_mesh_tet4()
        {

            // make a new mesh
            mMesh = new Mesh( 0 );
            mMesh->set_number_of_dimensions( 3 );

            // get the nodes
            Cell< mesh::Node * > & tNodes = mMesh->nodes() ;

            Vector< real > tX = { -1., 1., 1., -1., 1., -1., 0., 0. };
            Vector< real > tY = { -1., -1., 1., 1., 0., 0., 0., 0 };
            Vector< real > tZ = {  0., 0., 0., 0., 0., 0., -0.5, 0.5 };

            tX *= 0.5 ;
            tY *= 0.5 ;
            tZ *= 0.5 ;

            index_t tNumNodes = 8 ;

            tNodes.set_size( tNumNodes, nullptr );

            // create nodes
            for( uint k=0; k<tNumNodes; ++k )
            {
                tNodes( k ) = new mesh::Node( k+1, tX( k ), tY( k ), tZ( k ) );
            }


            // create the element factory
            mesh::ElementFactory tFactory ;

            // container for new element
            mesh::Element * tElement ;

            id_t tID = 1 ;

            Cell< mesh::Block * > & tBlocks = mMesh->blocks() ;

            // - - - - - - - - - - superconductor block
            mesh::Block * tScBlock = new mesh::Block( 1,  0 );
            tBlocks.push( tScBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TET4 , tID++ );
            tElement->insert_node( tNodes(   3 ), 0 );
            tElement->insert_node( tNodes(   5 ), 1 );
            tElement->insert_node( tNodes(   6 ), 2 );
            tElement->insert_node( tNodes(   7 ), 3 );
            tScBlock->elements().push( tElement );

            // - - - - - - - - - - ferromagnetic block
            mesh::Block * tFeBlock = new mesh::Block( 2,  0 );
            tBlocks.push( tFeBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TET4 , tID++ );
            tElement->insert_node( tNodes(   0 ), 0 );
            tElement->insert_node( tNodes(   6 ), 1 );
            tElement->insert_node( tNodes(   5 ), 2 );
            tElement->insert_node( tNodes(   7 ), 3 );
            tFeBlock->elements().push( tElement );

            // - - - - - - - - - - coil block
            mesh::Block * tCoilBlock = new mesh::Block( 3,  0 );
            tBlocks.push( tCoilBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TET4 , tID++ );
            tElement->insert_node( tNodes(   6 ), 0 );
            tElement->insert_node( tNodes(   1 ), 1 );
            tElement->insert_node( tNodes(   4 ), 2 );
            tElement->insert_node( tNodes(   7 ), 3 );
            tCoilBlock->elements().push( tElement );

            tElement = tFactory.create_lagrange_element( ElementType::TET4 , tID++ );
            tElement->insert_node( tNodes(   6 ), 0 );
            tElement->insert_node( tNodes(   4 ), 1 );
            tElement->insert_node( tNodes(   2 ), 2 );
            tElement->insert_node( tNodes(   7 ), 3 );
            tCoilBlock->elements().push( tElement );

            // - - - - - - - - - - air block
            mesh::Block * tAirBlock = new mesh::Block( 4,  0 );
            tBlocks.push( tAirBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TET4 , tID++ );
            tElement->insert_node( tNodes(   3 ), 0 );
            tElement->insert_node( tNodes(   6 ), 1 );
            tElement->insert_node( tNodes(   2 ), 2 );
            tElement->insert_node( tNodes(   7 ), 3 );
            tAirBlock->elements().push( tElement );

            tElement = tFactory.create_lagrange_element( ElementType::TET4 , tID++ );
            tElement->insert_node( tNodes(   0 ), 0 );
            tElement->insert_node( tNodes(   1 ), 1 );
            tElement->insert_node( tNodes(   6 ), 2 );
            tElement->insert_node( tNodes(   7 ), 3 );
            tAirBlock->elements().push( tElement );

            Cell< mesh::SideSet * > & tSideSets = mMesh->sidesets() ;

            // interface sc/air
            mesh::SideSet * tScAir = new mesh::SideSet( 1, 0 );
            tSideSets.push( tScAir );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++ );
            tElement->insert_node( tNodes( 6 ), 0 );
            tElement->insert_node( tNodes( 3 ), 1 );
            tElement->insert_node( tNodes( 7 ), 2 );
            tScAir->facets().push( new mesh::Facet( tElement ) );

            // interface sc/ferro
            mesh::SideSet * tScFerro = new mesh::SideSet( 2, 0 );
            tSideSets.push( tScFerro );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++  );
            tElement->insert_node( tNodes( 5 ), 0 );
            tElement->insert_node( tNodes( 6 ), 1 );
            tElement->insert_node( tNodes( 7 ), 2 );
            tScFerro->facets().push( new mesh::Facet( tElement ) );


            // interface ferro/air
            mesh::SideSet * tFerroAir = new mesh::SideSet( 3, 0 );
            tSideSets.push( tFerroAir );

            tElement = tFactory.create_lagrange_element( ElementType::TRI3 , tID++ );
            tElement->insert_node( tNodes( 0 ), 0 );
            tElement->insert_node( tNodes( 6 ), 1 );
            tElement->insert_node( tNodes( 7 ), 2 );
            tFerroAir->facets().push( new mesh::Facet( tElement ) );


            // - - - - - - - - - - finalize mesh
            mMesh->finalize();
            mMesh->create_edges( false, { 1 }, { 1, 2 } );
            mMesh->flag_curved_elements();

        }

//------------------------------------------------------------------------------

        void
        MaxwellJob::create_test_mesh_tet10( const bool aCurved  )
        {
            // make a new mesh
            mMesh = new Mesh( 0 );
            mMesh->set_number_of_dimensions( 3 );

            // get the nodes
            Cell< mesh::Node * > & tNodes = mMesh->nodes() ;

            Vector< real > tX = { -1., 1., 1., -1., 1., -1., 0., 0, 0., 1., 1., 0., -1., -1., -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, 0. };
            Vector< real > tY = { -1., -1., 1., 1., 0., 0, 0., 0, -1., -0.5, 0.5, 1., 0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0., 0, 0.5, 0.5, 0.5, 0.5, 0., 0, 0. };
            Vector< real > tZ = {  0., 0, 0., 0, 0., 0, -0.5, 0.5, 0., 0, 0., 0, 0., 0, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25, 0. };

            // shift points if mesh is curved
            if( aCurved )
            {
                tX( 14 ) = -0.382 ;
                tY( 14 ) = -0.618 ;
                tX( 15 ) = -0.382 ;
                tY( 15 ) = -0.618 ;
                tX( 22 ) = -0.382;
                tY( 22 ) =  0.618 ;
                tX( 23 ) = -0.382;
                tY( 23 ) =  0.618 ;
                tY( 24 ) = -0.1 ;
                tY( 25 ) =  0.1 ;
            }

            tX *= 0.5 ;
            tY *= 0.5 ;
            tZ *= 0.5 ;

            index_t tNumNodes = 27 ;

            tNodes.set_size( tNumNodes, nullptr );

            // create nodes
            for( uint k=0; k<tNumNodes; ++k )
            {
                tNodes( k ) = new mesh::Node( k+1, tX( k ), tY( k ), tZ( k ) );
            }

            // create the element factory
            mesh::ElementFactory tFactory ;

            // container for new element
            mesh::Element * tElement ;

            id_t tID = 1 ;

            Cell< mesh::Block * > & tBlocks = mMesh->blocks() ;

            // - - - - - - - - - - superconductor block
            mesh::Block * tScBlock = new mesh::Block( 1,  0 );
            tBlocks.push( tScBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TET10 , tID++ );
            tElement->insert_node( tNodes(    3 ), 0 );
            tElement->insert_node( tNodes(    5 ), 1 );
            tElement->insert_node( tNodes(    6 ), 2 );
            tElement->insert_node( tNodes(    7 ), 3 );
            tElement->insert_node( tNodes(   12 ), 4 );
            tElement->insert_node( tNodes(   24 ), 5 );
            tElement->insert_node( tNodes(   22 ), 6 );
            tElement->insert_node( tNodes(   23 ), 7 );
            tElement->insert_node( tNodes(   25 ), 8 );
            tElement->insert_node( tNodes(   26 ), 9 );
            tScBlock->elements().push( tElement );

            // - - - - - - - - - - ferro block
            mesh::Block * tFeBlock = new mesh::Block( 2,  0 );
            tBlocks.push( tFeBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TET10 , tID++ );
            tElement->insert_node( tNodes(    0 ), 0 );
            tElement->insert_node( tNodes(    6 ), 1 );
            tElement->insert_node( tNodes(    5 ), 2 );
            tElement->insert_node( tNodes(    7 ), 3 );
            tElement->insert_node( tNodes(   14 ), 4 );
            tElement->insert_node( tNodes(   24 ), 5 );
            tElement->insert_node( tNodes(   13 ), 6 );
            tElement->insert_node( tNodes(   15 ), 7 );
            tElement->insert_node( tNodes(   26 ), 8 );
            tElement->insert_node( tNodes(   25 ), 9 );

            tFeBlock->elements().push( tElement );

            // - - - - - - - - - - ferro block
            mesh::Block * tCoilBlock = new mesh::Block( 3,  0 );
            tBlocks.push( tCoilBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TET10 , tID++ );
            tElement->insert_node( tNodes(    6 ), 0 );
            tElement->insert_node( tNodes(    1 ), 1 );
            tElement->insert_node( tNodes(    4 ), 2 );
            tElement->insert_node( tNodes(    7 ), 3 );
            tElement->insert_node( tNodes(   16 ), 4 );
            tElement->insert_node( tNodes(    9 ), 5 );
            tElement->insert_node( tNodes(   18 ), 6 );
            tElement->insert_node( tNodes(   26 ), 7 );
            tElement->insert_node( tNodes(   17 ), 8 );
            tElement->insert_node( tNodes(   19 ), 9 );

            tCoilBlock->elements().push( tElement );

            tElement = tFactory.create_lagrange_element( ElementType::TET10 , tID++ );
            tElement->insert_node( tNodes(    6 ), 0 );
            tElement->insert_node( tNodes(    4 ), 1 );
            tElement->insert_node( tNodes(    2 ), 2 );
            tElement->insert_node( tNodes(    7 ), 3 );
            tElement->insert_node( tNodes(   18 ), 4 );
            tElement->insert_node( tNodes(   10 ), 5 );
            tElement->insert_node( tNodes(   20 ), 6 );
            tElement->insert_node( tNodes(   26 ), 7 );
            tElement->insert_node( tNodes(   19 ), 8 );
            tElement->insert_node( tNodes(   21 ), 9 );

            tCoilBlock->elements().push( tElement );

            // - - - - - - - - - - air block
            mesh::Block * tAirBlock = new mesh::Block( 4,  0 );
            tBlocks.push( tAirBlock );

            tElement = tFactory.create_lagrange_element( ElementType::TET10 , tID++ );
            tElement->insert_node( tNodes(    3 ), 0 );
            tElement->insert_node( tNodes(    6 ), 1 );
            tElement->insert_node( tNodes(    2 ), 2 );
            tElement->insert_node( tNodes(    7 ), 3 );
            tElement->insert_node( tNodes(   22 ), 4 );
            tElement->insert_node( tNodes(   20 ), 5 );
            tElement->insert_node( tNodes(   11 ), 6 );
            tElement->insert_node( tNodes(   23 ), 7 );
            tElement->insert_node( tNodes(   26 ), 8 );
            tElement->insert_node( tNodes(   21 ), 9 );
            tAirBlock->elements().push( tElement );

            tElement = tFactory.create_lagrange_element( ElementType::TET10 , tID++ );
            tElement->insert_node( tNodes(    0 ), 0 );
            tElement->insert_node( tNodes(    1 ), 1 );
            tElement->insert_node( tNodes(    6 ), 2 );
            tElement->insert_node( tNodes(    7 ), 3 );
            tElement->insert_node( tNodes(    8 ), 4 );
            tElement->insert_node( tNodes(   16 ), 5 );
            tElement->insert_node( tNodes(   14 ), 6 );
            tElement->insert_node( tNodes(   15 ), 7 );
            tElement->insert_node( tNodes(   17 ), 8 );
            tElement->insert_node( tNodes(   26 ), 9 );
            tAirBlock->elements().push( tElement );

            // - - - - - - - - - sidesets
            Cell< mesh::SideSet * > & tSideSets = mMesh->sidesets() ;

            // interface sc/air
            mesh::SideSet * tScAir = new mesh::SideSet( 1, 0 );
            tSideSets.push( tScAir );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++ );
            tElement->insert_node( tNodes( 6 ), 0 );
            tElement->insert_node( tNodes( 3 ), 1 );
            tElement->insert_node( tNodes( 7 ), 2 );
            tElement->insert_node( tNodes( 22 ), 3 );
            tElement->insert_node( tNodes( 23 ), 4 );
            tElement->insert_node( tNodes( 26 ), 5 );

            tScAir->facets().push( new mesh::Facet( tElement ) );

            // interface sc/ferro
            mesh::SideSet * tScFerro = new mesh::SideSet( 2, 0 );
            tSideSets.push( tScFerro );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++  );
            tElement->insert_node( tNodes( 5 ), 0 );
            tElement->insert_node( tNodes( 6 ), 1 );
            tElement->insert_node( tNodes( 7 ), 2 );
            tElement->insert_node( tNodes( 24 ), 3 );
            tElement->insert_node( tNodes( 26 ), 4 );
            tElement->insert_node( tNodes( 25 ), 5 );
            tScFerro->facets().push( new mesh::Facet( tElement ) );


            // interface ferro/air
            mesh::SideSet * tFerroAir = new mesh::SideSet( 3, 0 );
            tSideSets.push( tFerroAir );

            tElement = tFactory.create_lagrange_element( ElementType::TRI6 , tID++ );
            tElement->insert_node( tNodes( 0 ), 0 );
            tElement->insert_node( tNodes( 6 ), 1 );
            tElement->insert_node( tNodes( 7 ), 2 );
            tElement->insert_node( tNodes( 14 ), 3 );
            tElement->insert_node( tNodes( 26 ), 4 );
            tElement->insert_node( tNodes( 15 ), 5 );
            tFerroAir->facets().push( new mesh::Facet( tElement ) );


            // - - - - - - - - - - finalize mesh
            mMesh->finalize();
            mMesh->create_edges( false, { 1 }, { 1, 2 } );
            mMesh->create_faces( false,{ 1 }, { 1, 2 } );
            mMesh->flag_curved_elements();
        }

//------------------------------------------------------------------------------
    }
}