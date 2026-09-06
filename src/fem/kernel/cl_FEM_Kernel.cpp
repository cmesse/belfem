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

#include "commtools.hpp"
#include "stringtools.hpp"



#include "cl_FEM_Kernel.hpp"



#include "cl_Element_Factory.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_Graph_Vertex.hpp"
#include "fn_Graph_METIS.hpp"
#include "cl_IwgFactory.hpp"
#include "cl_Logger.hpp"
#include "cl_Map.hpp"
#include "cl_MeshChecker.hpp"
#include "cl_Mesh_Distributor.hpp"
#include "cl_Pipette.hpp"
#include "cl_Queue.hpp"
#include "cl_Timer.hpp"
#include "cl_Mesh_ConnectivityCalculator.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "fn_max.hpp"
#include "fn_sum.hpp"



namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Kernel::Kernel( KernelParameters * aKernelParameters ) :
            mCommRank( comm_rank() ),
            mCommSize( comm_size() ),
            mParams( aKernelParameters ),
            mMesh( aKernelParameters->mesh() ),
            mFieldOffset( aKernelParameters->mesh()->number_of_fields() )
        {
            // set the communication table
            if( mCommRank == 0 )
            {
                // check the mesh for flipped nodes and reorient them if needed
                // this will throw an error if edges or faces already exist
                // and the checker has not been called yet. This is by design.
                // run the checker BEFORE creating edges or faces!
                if ( ! mMesh->mesh_checker_flag() )
                {
                    MeshChecker tCheck( this->mesh() );
                }

                if ( mMesh->edges_exist() )
                {
                    for ( mesh::Block * tBlock : mMesh->blocks() )
                    {
                        if ( tBlock->number_of_elements() == 0 ) continue;

                        if ( tBlock->elements()( 0 )->has_edges() )
                        {
                            tBlock->set_edges_flag();
                        }
                    }
                }

                if ( mMesh->faces_exist() )
                {
                    for ( mesh::Block * tBlock : mMesh->blocks() )
                    {
                        if ( tBlock->number_of_elements() == 0 ) continue;
                        if ( tBlock->elements()( 0 )->has_faces() )
                        {
                            tBlock->set_faces_flag();
                        }
                    }
                }

                // set the flag for curved elements
                // for meshes on other procs, this is done in
                // receive_submesh
                mMesh->flag_curved_elements() ;

                // flag telling that this mesh is linked to a kernel
                mMesh->set_kernel_flag() ;
            }

            comm_barrier();

            // check if a Kernel has already been created
            if ( aKernelParameters->kernel() == nullptr )
            {

                // partition the mesh if in parallel mode
                if ( mCommSize > 1 && mCommRank == 0 )
                {

                    this->partition_mesh() ;
                }

                if ( mCommRank == 0 )
                {
                    mMesh->collect_hanging_basis() ;
                    mMesh->expand_hanging_basis_sources() ;
                }

                this->distribute_mesh() ;

                // having the element volumes is useful for example in the postproc
                this->compute_element_volumes();
            }
            else
            {
                this->distribute_mesh() ;
            }

            comm_barrier() ;
        }

//------------------------------------------------------------------------------

        Kernel::~Kernel()
        {
            for ( DofManager * tField: mDofManagers )
            {
                delete tField;
            }

            for ( IWG * tIWG: mIWGs )
            {
                delete tIWG;
            }

            if ( mOwnCommTables )
            {
                for ( auto tTable : mCommTables )
                {
                    delete tTable;
                }
            }

            if ( mOwnMesh )
            {
                delete mMesh;
            }
            if ( mOwnSubmesh )
            {
                delete mSubMesh;
            }

            for( Material * tMaterial : mMaterials )
            {
                delete tMaterial;
            }

            for( PhysicalBoundaryCondition * tBC : mBoundaryConditions )
            {
                delete tBC;
            }

            if( mOwnParameters )
            {
                delete mParams ;
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::partition_mesh()
        {
            Timer tTimer;

            mMesh->unflag_all_elements() ;
            mMesh->unflag_all_facets() ;
            mMesh->unflag_all_nodes();

            message( InfoLevel::Verbose, "Partitionig mesh ... ");

            // reset the element flags and indices
            for ( mesh::Element * tElement : mMesh->elements() )
            {
                tElement->set_index( gNoIndex );
                tElement->set_owner( mCommSize );
            }

            if ( mMesh->thin_shells().size() > 0 )
            {
                for ( mesh::Facet * tFacet : mMesh->facets() )
                {
                    tFacet->set_index( gNoIndex );
                    tFacet->set_owner( mCommSize );
                }
            }

            index_t tCount = 0 ;

            Cell< index_t > tBlockIndices ;

            // the following lines determine which blocks are considered as solid
            Map< id_t, index_t > tBlockMap ;
            for ( mesh::Block * tBlock : mMesh->blocks() )
            {
                tBlockMap[ tBlock->id() ] = tCount++ ;
            }
            DynamicBitset * tBitset = new DynamicBitset( tCount );

            BELFEM_ASSERT( mParams->selected_blocks().length() > 0, "No blocks selected for this Kernel. Don't know what to partition" );

            for ( id_t tID : mParams->selected_blocks() )
            {
                tBitset->set( tBlockMap[ tID ] );
            }

            for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
            {
                for ( mesh::Block * tBlock : tShell->blocks() )
                {
                    tBitset->reset( tBlockMap( tBlock->id() ) );
                }
            }

            tBitset->where( tBlockIndices );

            delete tBitset ;
            tBlockMap.clear() ;

            // now we can count the elements
            tCount = 0 ;
            for ( auto b : tBlockIndices )
            {
                tCount += mMesh->blocks()(b)->number_of_elements() ;
            }


            // populate the element graph
            Graph tElementGraph( tCount, nullptr);
            tCount = 0 ;
            Cell< index_t > tIndices ;
            for ( index_t b : tBlockIndices )
            {
                Cell< mesh::Element * > & tElements = mMesh->blocks()(b)->elements() ;
                for ( mesh::Element * tElement : tElements )
                {
                    tElement->reset_vertex_container() ;
                    tElement->set_index( tCount );
                    tElementGraph( tCount++ ) = tElement ;
                }
            }
            tBitset = new DynamicBitset( tCount ) ;

            for ( graph::Vertex * tVertex : tElementGraph )
            {
                mesh::Element * tElement = reinterpret_cast< mesh::Element * >( tVertex );
                tBitset->reset();

                 // search for neighbors
                for ( uint k=0; k<tElement->number_of_corner_nodes(); ++k )
                {
                    mesh::Node * tOrg = tElement->node( k )->original() ;
                    for ( uint e=0; e<tOrg->number_of_elements(); ++e )
                    {
                        mesh::Element * tNeighbor = tOrg->element( e ) ;
                        if ( tNeighbor->index() != gNoIndex )
                        {
                            tBitset->set( tNeighbor->index() );
                        }
                    }
                    for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                    {
                        mesh::Node * tDup = tOrg->duplicate( d ) ;
                        for ( uint e=0; e<tDup->number_of_elements(); ++e )
                        {
                            mesh::Element * tNeighbor = tDup->element( e ) ;
                            if ( tNeighbor->index() != gNoIndex )
                            {
                                tBitset->set( tNeighbor->index() );
                            }
                        }
                    }
                }

                tBitset->reset( tElement->index() );
                tBitset->where( tIndices );
                tElement->init_vertex_container( tIndices.size() );
                for ( index_t e : tIndices )
                {
                    tElement->insert_vertex( tElementGraph( e ) );
                }
            }
            delete tBitset ;
            graph::metis_partition( tElementGraph, mCommSize, false );
            Graph tFacetGraph ;
            if ( mMesh->thin_shells().size() > 0 )
            {
                // populate the facet graph
                tCount = 0 ;
                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    tCount += tShell->facets().size() ;
                }
                tFacetGraph.set_size( tCount, nullptr );

                tBitset = new DynamicBitset( tCount ) ;
                tCount = 0 ;
                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    for ( mesh::Facet * tFacet : tShell->facets() )
                    {
                        tFacet->set_index( tCount );
                        tFacetGraph( tCount++ ) = tFacet ;
                    }
                }
                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    for ( mesh::Facet * tFacet : tShell->facets() )
                    {
                        tBitset->reset();

                        for ( uint k=0; k<tFacet->number_of_corner_nodes(); ++k )
                        {
                            mesh::Node * tOrg = tFacet->node( k )->original() ;
                            for ( uint f=0; f<tOrg->number_of_facets(); ++f )
                            {
                                mesh::Facet * tNeighbor = tOrg->facet( f ) ;
                                if ( tNeighbor->index() != gNoIndex )
                                {
                                    tBitset->set( tNeighbor->index() );
                                }
                            }
                            for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                            {
                                mesh::Node * tDup = tOrg->duplicate( d ) ;
                                for ( uint f=0; f<tDup->number_of_facets(); ++f )
                                {
                                    mesh::Facet * tNeighbor = tDup->facet( f ) ;
                                    if ( tNeighbor->index() != gNoIndex )
                                    {
                                        tBitset->set( tNeighbor->index() );
                                    }
                                }
                            }
                        }
                        tBitset->reset( tFacet->index() );
                        tBitset->where( tIndices );
                        tFacet->init_vertex_container( tIndices.size() );
                        for ( index_t f : tIndices )
                        {
                            tFacet->insert_vertex( tFacetGraph( f ) );
                        }
                    }
                }

                delete tBitset ;
                graph::metis_partition( tFacetGraph, mCommSize, false );

                // smooth facet partition to reduce checkerboarding via majority voting
                for ( uint s=0; s<3; ++s )
                {
                    tCount = 0 ;
                    for ( graph::Vertex * tVertex : tFacetGraph )
                    {
                        // count owner IDs among neighbors
                        Vector< index_t > tOwnerCount( mCommSize, 0 );
                        for ( uint n=0; n<tVertex->number_of_vertices(); ++n )
                        {
                            proc_t tNeighborOwner = tVertex->vertex( n )->owner() ;
                            if ( tNeighborOwner < mCommSize )
                            {
                                ++tOwnerCount( tNeighborOwner ) ;
                            }
                        }

                        // find majority owner among neighbors
                        proc_t tMajorityOwner = tVertex->owner() ;
                        uint tMaxCount = 0 ;
                        for ( proc_t p=0; p<mCommSize; ++p )
                        {
                            if ( tOwnerCount( p ) > tMaxCount )
                            {
                                tMaxCount = tOwnerCount( p ) ;
                                tMajorityOwner = p ;
                            }
                        }

                        // only change if clear majority (> 50% of neighbors)
                        if ( tMaxCount > tVertex->number_of_vertices() / 2 && tMajorityOwner != tVertex->owner() )
                        {
                            tVertex->set_owner( tMajorityOwner ) ;
                            ++tCount ;
                        }
                    }
                    if ( tCount == 0 ) break ;
                }

                for ( graph::Vertex * tVertex : tFacetGraph )
                {
                    mesh::Facet * tFacet = reinterpret_cast< mesh::Facet * >( tVertex );
                    tFacet->master()->set_owner( tFacet->owner() );
                    tFacet->slave()->set_owner( tFacet->owner() );
                    tFacet->master()->flag();
                    tFacet->slave()->flag();
                    tFacet->flag();
                }
            }

            index_t tSweep = 1 ;
            while ( tSweep > 0 )
            {
                tSweep = 0 ;

                // enforce consistent master-slave-facet owners
                for ( graph::Vertex * tVertex : tFacetGraph )
                {
                    mesh::Facet * tFacet = reinterpret_cast< mesh::Facet * >( tVertex );
                    proc_t tOwner = std::min( tFacet->owner(), tFacet->master()->owner() ) ;
                    tOwner = std::min(tOwner, tFacet->slave()->owner() ) ;

                    if ( tOwner != tFacet->owner() || tOwner != tFacet->master()->owner() || tOwner != tFacet->slave()->owner() )
                    {
                        tFacet->set_owner( tOwner ) ;
                        tFacet->master()->set_owner( tOwner ) ;
                        tFacet->slave()->set_owner( tOwner ) ;
                        ++tSweep ;
                    }
                }
            }

            Vector< index_t > tNumElements( mCommSize, 0 );
            for ( graph::Vertex * tVertex : tElementGraph )
            {
                ++tNumElements( tVertex->owner() );
            }

            if ( mMesh->thin_shells().size() > 0 )
            {
                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    Cell< mesh::Facet * > & tShellFacets = tShell->facets() ;

                    for ( mesh::Block* tBlock : tShell->blocks() )
                    {
                        Cell< mesh::Element * > & tLayer = tBlock->elements() ;

                        index_t n = tLayer.size() ;
                        for ( index_t e=0; e<n; ++e )
                        {
                            tLayer( e )->set_owner( tShellFacets( e )->owner() );
                        }
                    }
                }

                // rule: side connector elements inherit ownership from the
                // master of their recovery facet, which is the layer block
                // element the wall spans ( its owner was set above ); the
                // recovery facet follows the same owner. The facet always
                // carries the id of the wall element plus one
                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    for ( mesh::Block * tBlock : tShell->side_connector_blocks() )
                    {
                        for ( mesh::Element * tElement : tBlock->elements() )
                        {
                            mesh::Facet * tFacet = mMesh->facet( tElement->id() + 1 );
                            tFacet->set_owner( tFacet->master()->owner() );
                            tElement->set_owner( tFacet->owner() );
                        }
                    }
                }
#ifdef DEBUG
                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    Cell< mesh::Facet * > & tShellFacets = tShell->facets() ;
                    for ( mesh::Facet * tFacet : tShellFacets )
                    {
                        BELFEM_ASSERT( tFacet->owner() < mCommSize, "facet %d has owner %d, which is larger than the number of processors %d", tFacet->index(), tFacet->owner(), mCommSize );
                    }
                    for ( mesh::Block* tBlock : tShell->blocks() )
                    {
                        for ( mesh::Element * tElement : tBlock->elements() )
                        {
                            BELFEM_ASSERT( tElement->owner() < mCommSize, "element %d has owner %d, which is larger than the number of processors %d", tElement->index(), tElement->owner(), mCommSize );
                        }
                    }
                    for ( mesh::Block* tBlock : tShell->side_connector_blocks() )
                    {
                        for ( mesh::Element * tElement : tBlock->elements() )
                        {
                            BELFEM_ASSERT( tElement->owner() < mCommSize, "side connector element %d has owner %d, which is larger than the number of processors %d", tElement->index(), tElement->owner(), mCommSize );
                        }
                    }
                }
#endif
                Vector< index_t > tNumFacets( mCommSize, 0 );
                for ( graph::Vertex * tFacet : tFacetGraph )
                {
                    ++tNumFacets( tFacet->owner() ) ;
                }
            }


            mMesh->update_element_indices();
            mMesh->update_facet_indices();

            mMesh->update_ownerships();

            // fix ownerships of non-thin-shell facets
            mMesh->unflag_all_facets() ;
            for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
            {
                Cell< mesh::Facet * > & tShellFacets = tShell->facets() ;
                for ( mesh::Facet * tFacet : tShellFacets )
                {
                    tFacet->flag() ;
                }
            }
            for ( mesh::Facet * tFacet : mMesh->facets() )
            {
                if ( tFacet->is_flagged() ) continue;
                if ( tFacet->has_master() )
                {
                    tFacet->set_owner( tFacet->master()->owner() );
                }
                else if ( tFacet->has_slave() )
                {
                    tFacet->set_owner( tFacet->slave()->owner() );
                }
            }

            // mMesh->save( "metis.exo" );

            message( InfoLevel::Verbose, "   ... time for partitioniong: %u ms ",  ( unsigned int )  tTimer.stop() );
        }

//------------------------------------------------------------------------------

        IWG *
        Kernel::create_equation(
                const IwgType aEquationType,
                const ModelDimensionality aModelDimensionality,
                const Vector< id_t > aBlocks,
                const Vector< id_t > aSideSets )
        {

            // create the equation object
            IWG * aIWG = nullptr ;

            if( is_maxwell( aEquationType ) )
            {
                BELFEM_ERROR( false, "the kernel can't create maxwell type equations by itself. Use the MaxwellFactory!" );
            }
            else
            {
                IwgFactory tFactory( mMesh );

                // make sure that the dimensionality has been set, otherwise assume 2D or 3D from mesh dimension
                ModelDimensionality tModelDimensionality =
                    aModelDimensionality == ModelDimensionality::UNDEFINED ?
                        ( mMesh->number_of_dimensions() == 3 ? ModelDimensionality::ThreeD : ModelDimensionality::TwoD ) :
                        aModelDimensionality ;

                aIWG = tFactory.create_iwg( aEquationType, tModelDimensionality );

                if( aBlocks.length() > 0 )
                {
                    // select the specified blocks for the IWG
                    aIWG->select_blocks( aBlocks );
                }
                else
                {
                    // select all blocks
                    aIWG->select_blocks( tFactory.all_block_ids() );
                }

                if ( aSideSets.length() > 0 )
                {
                    aIWG->select_sidesets( aSideSets );
                }

                // deliberately not selecting all sidesets here
            }

            // add to containers
            mIWGs.push( aIWG );

            // return the field
            return aIWG ;
        }

//------------------------------------------------------------------------------

        void
        Kernel::add_equation( IWG * aEquation )
        {
            mIWGs.push( aEquation );
        }

//------------------------------------------------------------------------------

        DofManager *
        Kernel::create_field( IWG * aEquation )
        {
            aEquation->initialize() ;
            // create the dof manager
            DofManager * aField = new DofManager( this, mDofManagers.size() ) ;

            // check the sanity of the mesh
            if( mCommRank == 0 )
            {
                int tStatus = aEquation->check_mesh( aField->mesh(), 0) ;

                BELFEM_ERROR( tStatus == 0 , " mIWG->check_mesh() failed with error code %d",
                         tStatus );
            }

            // wait for other procs
            comm_barrier() ;

            // link IWG to field
            aField->set_equation( aEquation );

            mDofManagers.push( aField );

            return aField ;
        }

//------------------------------------------------------------------------------

        const  KernelParameters *
        Kernel::params()
        {
            return mParams;
        }

//------------------------------------------------------------------------------

        void
        Kernel::distribute_mesh()
        {
            if ( mCommSize < 2 )
            {
                if ( mParams->kernel() == nullptr )
                {
                    Cell< mesh::Element * > & tElements = mMesh->elements() ;
                    for ( mesh::Element * tElement : tElements )
                    {
                        tElement->reset_vertex_container() ;
                    }

                    mesh::ConnectivityCalculator tCalculator( mMesh ) ;
                    tCalculator.connect_elements_to_elements() ;
                    tCalculator.connect_facets_to_facets() ;
                }
                return ;
            }

            if ( mParams->kernel() == nullptr )
            {

                if ( mCommRank == 0 )
                {
                    mesh::ConnectivityCalculator tCalculator( mMesh ) ;
                    tCalculator.connect_elements_to_elements() ;
                    tCalculator.connect_facets_to_facets() ;
                }

                mesh::Distributor tDistributor( mMesh );

                tDistributor.run();

                Cell< mesh::Element * > & tElements = mMesh->elements() ;
                for ( mesh::Element * tElement : tElements )
                {
                    tElement->reset_vertex_container() ;
                }

                if ( mCommRank == 0 )
                {
                    // get communication tables before they go out of scope
                    // the comm tables are const pointers, so we can't use an std::move here
                    Cell< mesh::CommTable * > & tTables = tDistributor.tables();
                    mCommTables.vector_data().assign( tTables.begin(), tTables.end() );
                    tTables.clear(); // distributor destructor mustn't delete the tables

                    mOwnCommTables = true ;

                }
                else
                {
                    // grab submesh before it goes out of scope
                    mSubMesh = tDistributor.partial_mesh();

                    mesh::ConnectivityCalculator tCalculator( mSubMesh ) ;
                    tCalculator.connect_elements_to_elements() ;
                    tCalculator.connect_facets_to_facets() ;

                    mOwnSubmesh = true ;
                }

                //mesh::ConnectivityCalculator tCalculator( this->mesh() ) ;
                //tCalculator.connect_elements_to_elements() ;
                //tCalculator.connect_facets_to_facets() ;
            }
            else
            {
                if ( mCommRank == 0 )
                {

                    mCommTables.set_size( mCommSize, nullptr );
                    for ( proc_t p=1; p<mCommSize; ++p )
                    {
                        mCommTables( p ) = mParams->kernel()->comm_table( p );
                    }
                }
                else
                {
                    mSubMesh = mMesh ;
                    mMesh = new Mesh( mSubMesh->number_of_dimensions() );
                }
            }
        }

//------------------------------------------------------------------------------

        DofManager *
        Kernel::dofmgr( const uint aIndex )
        {
            return mDofManagers( aIndex );
        }

//------------------------------------------------------------------------------

        Material *
        Kernel::material( const string & aLabel )
        {
            BELFEM_ERROR(  mMaterialMap.key_exists( aLabel ),
                              "Could not find material %s in kernel", aLabel.c_str() );


            return mMaterialMap( aLabel );

        }

//------------------------------------------------------------------------------

        void
        Kernel::add_material( const string aLabel, Material * aMaterial )
        {
            if ( ! mMaterialMap.key_exists( aLabel ) )
            {
                mMaterialMap[aLabel] = aMaterial ;
            }

            // add material to array
            mMaterials.push( aMaterial );
        }

//------------------------------------------------------------------------------


        void
        Kernel::add_boundary_condition( belfem::fem::PhysicalBoundaryCondition * aBC )
        {
            // add material to array
            mBoundaryConditions.push( aBC );
        }

//------------------------------------------------------------------------------

        void
        Kernel::compute_boundary_conditions( real aTime )
        {
            for( PhysicalBoundaryCondition * tBC : mBoundaryConditions )
            {
                tBC->impose_bc( aTime ) ;
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::claim_parameter_ownership( const bool aFlag )
        {
            mOwnParameters = aFlag ;
        }

//------------------------------------------------------------------------------

        void
        Kernel::collect_elements(
                const Vector< id_t >    & aBlockIDs,
                const Vector< id_t >    & aSideSetIDs,
                Cell< mesh::Element * > & aElements )
        {
            // count elements
            index_t tCount = 0 ;

            // loop over all block IDs
            for( id_t b: aBlockIDs )
            {
                if( mMesh->block_exists( b ) )
                {
                    tCount += mMesh->block( b )->number_of_elements() ;
                }
            }

            // loop over all sideset IDs
            for( id_t s: aSideSetIDs )
            {
                if( mMesh->sideset_exists( s ) )
                {
                    tCount += mMesh->sideset( s )->number_of_facets() ;
                }
            }

            // allocate memory
            aElements.set_size( tCount, nullptr );

            // reset counter
            tCount = 0 ;

            // grab elements
            for( id_t b: aBlockIDs )
            {
                if( mMesh->block_exists( b ) )
                {
                    Cell< mesh::Element * > & tElements = mMesh->block( b )->elements() ;

                    for ( mesh::Element * tElement : tElements )
                    {
                        aElements( tCount++ ) = tElement ;
                    }
                }
            }

            // grab facets
            for( id_t s: aSideSetIDs )
            {
                if( mMesh->sideset_exists( s ) )
                {
                    Cell< mesh::Facet * > & tFacets = mMesh->sideset( s )->facets() ;

                    for ( mesh::Facet * tFacet : tFacets )
                    {
                        aElements( tCount++ ) = tFacet->element() ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::compute_element_volumes()
        {
            // get the mesh
            Mesh * tMesh = this->mesh() ;
            BELFEM_ASSERT( tMesh != nullptr, "No mesh specified." );

            // create a new field on the mesh
            Vector< real > & tVolumes = tMesh->field_exists( "_Volumes" ) ?
                tMesh->field_data( "_Volumes" ) : tMesh->create_field( "_Volumes", EntityType::ELEMENT ) ;

            // make sure that the volumes field is properly sized
            tVolumes.set_size( tMesh->number_of_elements(), BELFEM_QUIET_NAN );

            // we don't want to write this field
            tMesh->field( "_Volumes" )->set_write_to_file_flag( false );

            mesh::Pipette  tPip ;

            // reset the counter
            index_t tCount = 0 ;

            index_t tNumMyElements = 0 ;
            index_t tNumNotMyElements = 0 ;

            // first, each proc computes its own elements
            // loop over all blocks
            for( mesh::Block * tBlock : tMesh->blocks() )
            {
                if ( tBlock->domain_type() == DomainType::ThinShell || tBlock->domain_type() == DomainType::Buffer ) continue ;

                // set the element type
                tPip.set_element_type( tBlock->element_type() );

                // grab the elements
                Cell< mesh::Element * > & tElements = tBlock->elements() ;


                for( mesh::Element * tElement : tElements )
                {
                    // check if this proc owns the element
                    if( tElement->owner() == mCommRank )
                    {
                        // compute the element volume
                        tVolumes( tElement->index() ) = tPip.measure( tElement );

                        // increment negative counter
                        if( tVolumes( tElement->index() ) < 0 )
                        {
                            ++tCount ;
                        }

                        ++tNumMyElements ;
                    }
                    else if ( tElement->owner() < mCommSize )
                    {
                        ++tNumNotMyElements ;
                    }
                }
            }

            // Thin-shell layer volumes are facet area x layer thickness.
            // The correspondence element <-> facet is NOT positional on a
            // worker: the distributed block holds a subset of the layer
            // elements while the shell's facet container does not shrink the
            // same way ( observed 40100 elements vs 40145 facets, np = 4 ).
            // What IS invariant on every rank is the ID arithmetic: within a
            // layer block the elements were created facet by facet with
            // consecutive ids ( ThinShellFactory::create_elements_on_blocks_* ),
            // so master facet ordinal = element id - first element id of the
            // block. Rank 0 owns the full mesh, measures the facets in master
            // order, and shares areas + first-ids; workers never touch their
            // own facet subset at all.
            //
            // The exchange is ONE flat payload for all shells, done before
            // the shell loop. A worker only holds the shells it has elements
            // of ( ProtoMesh::create_thinshells drops the others ), so any
            // communication inside a per-shell loop pairs rank 0's shell k
            // with a different shell on such a worker. Shells and blocks are
            // therefore looked up by id, never by position.
            //
            // header layout, per shell, in rank-0 shell order:
            //   shell id, area offset, facet count, block count,
            //   ( block id, first element id ) x block count
            Vector< real >    tAllSurfaces ;
            Vector< index_t > tHeader ;

            if ( mCommRank == 0 )
            {
                index_t tNumSurfaces = 0 ;
                index_t tHeaderSize  = 0 ;

                for ( mesh::ThinShell * tShell : tMesh->thin_shells() )
                {
                    tNumSurfaces += tShell->facets().size() ;
                    tHeaderSize  += 4 + 2 * tShell->blocks().size() ;
                }

                tAllSurfaces.set_size( tNumSurfaces, BELFEM_QUIET_NAN );
                tHeader.set_size( tHeaderSize, 0 );

                index_t f = 0 ;
                index_t h = 0 ;

                for ( mesh::ThinShell * tShell : tMesh->thin_shells() )
                {
                    // shell ids are sideset ids, unique by construction;
                    // the id-keyed lookup below must never see a duplicate
                    for ( index_t i = 0; i < h; i += 4 + 2 * tHeader( i + 3 ) )
                    {
                        BELFEM_ERROR( tHeader( i ) != tShell->id(),
                            "Thin shell id %u occurs twice on the mesh",
                            ( unsigned int ) tShell->id() );
                    }

                    tHeader( h++ ) = tShell->id() ;
                    tHeader( h++ ) = f ;
                    tHeader( h++ ) = tShell->facets().size() ;
                    tHeader( h++ ) = tShell->blocks().size() ;

                    tPip.set_facet_type( tShell->element_type() );

                    for ( mesh::Facet * tFacet : tShell->facets() )
                    {
                        tAllSurfaces( f++ ) = tPip.measure( tFacet );
                    }

                    for ( mesh::Block * tBlock : tShell->blocks() )
                    {
                        BELFEM_ERROR( tBlock->elements().size() > 0
                                   && tBlock->elements().size() == tShell->facets().size(),
                            "Element count of thin shell block %u does not match the facet count of its shell ( %lu vs %lu )",
                            ( unsigned int ) tBlock->id(),
                            ( long unsigned int ) tBlock->elements().size(),
                            ( long unsigned int ) tShell->facets().size() );

                        tHeader( h++ ) = tBlock->id() ;
                        tHeader( h++ ) = tBlock->elements()( 0 )->id() ;
                    }
                }

                if ( mCommSize > 1 )
                {
                    // share/receive is NOT collective: rank guard required.
                    // both payloads are variable-size, so both go chunked,
                    // in this order
                    share( tAllSurfaces );
                    share( tHeader );
                }
            }
            else
            {
                // unconditional, also on a rank without any thin shell:
                // rank 0 sends both payloads to every worker
                receive( tAllSurfaces );
                receive( tHeader );
            }

            for ( mesh::ThinShell * tShell : tMesh->thin_shells() )
            {
                // locate the shell's header record by id
                index_t h = 0 ;
                while ( h + 4 <= tHeader.length() && tHeader( h ) != tShell->id() )
                {
                    h += 4 + 2 * tHeader( h + 3 );
                }

                BELFEM_ERROR( h + 4 <= tHeader.length(),
                    "Thin shell %u is not in the volume header received from rank 0",
                    ( unsigned int ) tShell->id() );

                const index_t tOffset    = tHeader( h + 1 );
                const index_t tNumFacets = tHeader( h + 2 );
                const index_t tNumBlocks = tHeader( h + 3 );

                BELFEM_ERROR( h + 4 + 2 * tNumBlocks <= tHeader.length()
                           && tOffset <= tAllSurfaces.length()
                           && tNumFacets <= tAllSurfaces.length() - tOffset,
                    "Corrupt volume header for thin shell %u",
                    ( unsigned int ) tShell->id() );

                for ( mesh::Block * tBlock : tShell->blocks() )
                {
                    BELFEM_ASSERT( ! std::isnan( tBlock->thickness() ) && std::abs( tBlock->thickness() ) > 0,
                        "Thin shell thickness for block %u not set.", ( unsigned int ) tBlock->id() );

                    // first element id of this block, looked up by block id
                    id_t tFirstId = gNoID ;
                    for ( index_t b = 0; b < tNumBlocks; ++b )
                    {
                        if ( tHeader( h + 4 + 2 * b ) == tBlock->id() )
                        {
                            tFirstId = tHeader( h + 5 + 2 * b );
                            break ;
                        }
                    }

                    BELFEM_ERROR( tFirstId != gNoID,
                        "Block %u of thin shell %u is not in the volume header received from rank 0",
                        ( unsigned int ) tBlock->id(),
                        ( unsigned int ) tShell->id() );

                    for( mesh::Element * tElement : tBlock->elements() )
                    {
                        if ( tElement->owner() == mCommRank )
                        {
                            // master facet ordinal via the id arithmetic
                            const index_t f = tElement->id() - tFirstId ;

                            BELFEM_ASSERT( f < tNumFacets,
                                "Element %lu of thin shell block %u maps to facet ordinal %lu, but the shell has only %lu facets",
                                ( long unsigned int ) tElement->id(),
                                ( unsigned int ) tBlock->id(),
                                ( long unsigned int ) f,
                                ( long unsigned int ) tNumFacets );

                             real tVol = tAllSurfaces( tOffset + f ) * tBlock->thickness() ;
                             tVolumes( tElement->index() ) = tVol ;

                            // increment negative counter
                             if ( tVol < 0 )
                             {
                                ++tCount ;
                             }

                            ++tNumMyElements ;
                        }
                        else
                        {
                            ++tNumNotMyElements ;
                        }
                    }
                }
            }
            comm_barrier() ;

            // this routine checks that all elements are positive
            if( mCommRank == 0 )
            {
                Vector< index_t > tAllCount( mCommSize, 0 );
                collect( tAllCount, tCount );

                tCount = sum( tAllCount );

                BELFEM_ERROR(  tCount == 0, "Detected %lu elements with negative volume on mesh",
                               ( long unsigned int ) tCount );

            }
            else
            {
                // send my counter to master
                send( tCount );
            }


            // next, each proc sends its computations to the master
            if ( mCommRank == 0 )
            {
                Cell< Vector< id_t > > tAllIDs( mCommSize, {} );
                Cell< Vector< real > > tAllVolumes( mCommSize, {} );
                Cell< Vector< id_t > > tAllOtherIDs( mCommSize, {} );
                Cell< Vector< real > > tAllOtherVolumes( mCommSize, {} );

                collect( tAllIDs );
                collect( tAllVolumes );
                collect( tAllOtherIDs );

                for ( proc_t p=1; p<mCommSize; ++p )
                {
                    Vector< id_t > & tProcIDs = tAllIDs( p );
                    Vector< real > & tProcVolumes = tAllVolumes( p );
                    tCount = 0 ;
                    for ( id_t tID : tProcIDs )
                    {
                        tVolumes( tMesh->element( tID )->index() ) = tProcVolumes( tCount++ );
                    }
                }

                for ( proc_t p=1; p<mCommSize; ++p )
                {
                    Vector< id_t > & tProcIDs = tAllOtherIDs( p );
                    Vector< real > & tProcVolumes = tAllOtherVolumes( p );

                    // each rank expects its own not-my-element count, which is
                    // the length of the ID list it sent — not rank 0's count
                    tProcVolumes.set_size( tProcIDs.length() );
                    tCount = 0 ;
                    for ( id_t tID : tProcIDs )
                    {
                        tProcVolumes( tCount++ ) = tVolumes( tMesh->element( tID )->index() );
                    }
                }
                comm_barrier() ;
                distribute( tAllOtherVolumes );

            }
            else
            {
                Vector< id_t > tMyIDs( tNumMyElements );
                Vector< real > tMyVolumes( tNumMyElements );
                Vector< id_t > tNotMyIDs( tNumNotMyElements );
                Vector< real > tNotMyVolumes( tNumNotMyElements );

                index_t tCountA = 0 ;
                index_t tCountB = 0 ;
                for ( mesh::Element * tElement: tMesh->elements() )
                {
                    if ( tElement->owner() == mCommRank )
                    {
                        tMyIDs( tCountA ) = tElement->id() ;
                        tMyVolumes( tCountA ) = tVolumes( tElement->index() ) ;
                        ++tCountA ;
                    }
                    else if ( tElement->owner() < mCommSize )
                    {
                        tNotMyIDs( tCountB++ ) = tElement->id() ;
                    }
                }

                send( tMyIDs ) ;
                send( tMyVolumes ) ;
                send( tNotMyIDs ) ;

                comm_barrier() ;
                receive( tNotMyVolumes );

                tCount = 0 ;
                for ( id_t tID : tNotMyIDs )
                {
                    tVolumes( tMesh->element( tID )->index() ) = tNotMyVolumes( tCount++ );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::set_controller( Controller * aController )
        {
            mController = aController ;
        }

//------------------------------------------------------------------------------

        Controller *
        Kernel::controller()
        {
            // always-active: an uninitialized read here was a real defect.
            // Callers that may legitimately run without a Controller
            // ask has_controller() first
            BELFEM_ERROR( mController != nullptr, "Controller is not set for this kernel." );
            return mController ;
        }

//------------------------------------------------------------------------------

    }
}
