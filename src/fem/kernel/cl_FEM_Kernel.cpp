//
// Created by Christian Messe on 24.10.19.
//
#include "cl_FEM_Kernel.hpp"

#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_max.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "stringtools.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Graph_Vertex.hpp"
#include "fn_Graph_symcrm.hpp"
#include "cl_MaterialFactory.hpp"
#include "cl_IwgFactory.hpp"
#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Kernel::Kernel( KernelParameters * aKernelParameters, bool aLegacyMode ) :
            mParams( aKernelParameters ),
            mMesh( aKernelParameters->mesh() ),
            mMyRank( comm_rank() ),
            mMasterRank( mMesh->master() ),
            mNumberOfProcs( aKernelParameters->number_of_procs() ),
            mFieldOffset( aKernelParameters->mesh()->number_of_fields() )
        {

            this->activate_procs();

            // set the communication table
            if( mMyRank == mMasterRank )
            {
                mCommTable = aKernelParameters->selected_procs();

                BELFEM_ERROR( mCommTable( 0 ) == mMasterRank,
                        "First proc in comm table must be master" );

                // create communication map
                uint n = mCommTable.length();
                for( uint k=0; k<n; ++k )
                {
                    mCommMap[ mCommTable( k ) ] = k;
                }

                // flag telling that this mesh is linked to a kernel
                mMesh->set_kernel_flag() ;
            }
            else
            {
                mCommTable.set_size( 1, mMesh->master() );
            }

            // rearrange nodes
            // this algorithm is slow and probably not needed anyhow
            // causes an error on c4.mesh
            // this->symrcm();

            // partition the mesh if in parallel mode
            if ( mNumberOfProcs > 1 && mMyRank == mMesh->master() )
            {
                if( mParams->auto_partition() )
                {
                    if ( mParams->selected_blocks().length() > 0 )
                    {
                        mMesh->partition( mNumberOfProcs, mParams->selected_blocks());

                        mMesh->unflag_all_elements();
                        mMesh->unflag_all_facets();
                        if ( mMesh->ghost_block_ids().length() > 0 )
                        {
                            mMesh->partition( mNumberOfProcs, mMesh->ghost_block_ids(), true, false );

                            Cell< mesh::Element * > & tLayer0 = mMesh->block( mMesh->ghost_block_ids()( 0 ) )->elements() ;

                            index_t tNumElemsPerLayer = tLayer0.size();


                            Vector< proc_t > tOwners( tNumElemsPerLayer, mNumberOfProcs );

                            // loop over all thin shell sidesets
                            for ( id_t tID: mMesh->ghost_sideset_ids())
                            {
                                // grab facets
                                Cell< mesh::Facet * > & tFacets = mMesh->sideset( tID )->facets();

                                index_t tCount = 0;
                                for ( mesh::Facet * tFacet: tFacets )
                                {
                                    // get owner
                                    proc_t tOwner = tLayer0( tCount )->owner();
                                    tOwners( tCount++ ) = tOwner ;

                                    tFacet->set_owner( tOwner );
                                    tFacet->master()->set_owner( tOwner );
                                    tFacet->slave()->set_owner( tOwner );
                                    tFacet->master()->flag();
                                    tFacet->slave()->flag();
                                    tFacet->flag();
                                }
                            }

                            // fix ownership of other layers
                            uint tNumLayers = mMesh->ghost_block_ids().length() ;
                            for( uint l=1; l<tNumLayers; ++l )
                            {
                                Cell< mesh::Element * > & tLayer
                                    = mMesh->block( mMesh->ghost_block_ids()( l ) )->elements() ;
                                index_t tCount = 0;
                                for( mesh::Element * tElement : tLayer )
                                {
                                    tElement->set_owner( tOwners( tCount++ ) );
                                }
                            }

                            // now we need to fix the ownerships of the original facets
                            // that describe the thin shell, but are not part of the ghost
                            // block

                            // grab all facet on mesh
                            Cell< mesh::Facet * > & tFacets = mMesh->facets();

                            // and now let's make sure that the ownerships are correct
                            for ( mesh::Facet * tFacet: tFacets )
                            {
                                if ( !tFacet->is_flagged())
                                {
                                    if ( tFacet->has_master())
                                    {
                                        tFacet->set_owner( tFacet->master()->owner());
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        mMesh->partition( mNumberOfProcs );
                    }
                }

                // allocate memory for node table
                mNodeTable.set_size( mNumberOfProcs, {} );
                mElementTable.set_size( mNumberOfProcs, {} );

                // allocate memory for edge table
                mEdgeTable.set_size( mNumberOfProcs, {} );

                // allocate memory for face table
                mFaceTable.set_size( mNumberOfProcs, {} );

                // for lambda dofs
                mFacetTable.set_size( mNumberOfProcs, {} );
                mConnectorTable.set_size( mNumberOfProcs, {} );
            }

            this->distribute_mesh() ;

            if( aLegacyMode )
            {
                this->create_fields();
            }

            comm_barrier() ;
        }

//------------------------------------------------------------------------------

        Kernel::~Kernel()
        {
            for ( Field * tField: mFields )
            {
                delete tField;
            }

            for ( DofManager * tField: mDofManagers )
            {
                delete tField;
            }

            for ( IWG * tIWG: mIWGs )
            {
                delete tIWG;
            }

            if( mMyRank != mMasterRank )
            {
                delete mSubMesh;
            }

            for( Material * tMaterial : mMaterials )
            {
                delete tMaterial;
            }

            if( mOwnParameters )
            {
                delete mParams ;
            }
        }

//------------------------------------------------------------------------------

        IWG *
        Kernel::create_equation(
                const IwgType aEquationType, const Vector< id_t > aBlocks )
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

                aIWG = tFactory.create_iwg( aEquationType );

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
            int tStatus = aEquation->check_mesh( aField->mesh(), this->master() ) ;

            BELFEM_ERROR( tStatus == 0 , " mIWG->check_mesh() failed with error code %d",
                         tStatus );

            // link IWG to field
            aField->set_equation( aEquation );

            mDofManagers.push( aField );

            // link mesh of field to IWG
            // aEquation->link_to_group( aField->block(aEquation->selected_blocks()(0)));

            return aField ;
        }

//------------------------------------------------------------------------------

        proc_t
        Kernel::master()
        {
            return mMasterRank;
        }

//------------------------------------------------------------------------------

        const  KernelParameters *
        Kernel::params()
        {
            return mParams;
        }

//------------------------------------------------------------------------------

        const Vector< proc_t > &
        Kernel::comm_table() const
        {
            return mCommTable;
        }

//------------------------------------------------------------------------------

        proc_t
        Kernel::comm_table( const index_t aIndex ) const
        {
            return mCommTable( aIndex );
        }

//------------------------------------------------------------------------------

        void
        Kernel::activate_procs()
        {
            proc_t tNumberOfProcs = comm_size();
            mMyCommIndex = tNumberOfProcs;

            if ( tNumberOfProcs > 1 )
            {
                Vector <uint> tFlags;
                Vector <proc_t> tCommunicationList;
                const Vector< proc_t > & tCommTable = mParams->selected_procs();
                uint tFlag = 0;
                if ( mMyRank == mMasterRank )
                {
                    // fill an array with flags
                    tFlags.set_size( tNumberOfProcs, 0 );
                    tCommunicationList.set_size( tNumberOfProcs );

                    for ( proc_t k = 0; k < tNumberOfProcs; ++k )
                    {
                        tCommunicationList( k ) = k;
                    }

                    for ( proc_t p = 0; p < mNumberOfProcs; ++p )
                    {
                        tFlags( tCommTable( p ) ) = 1;

                        if( tCommTable( p ) == mMyRank )
                        {
                            mMyCommIndex = p;
                        }
                    }

                    send( tCommunicationList, tFlags );
                }
                else
                {
                    receive( mMasterRank, tFlag );

                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::distribute_mesh()
        {
            if( mNumberOfProcs > 1 )
            {
                // wait for other procs
                comm_barrier();

                if ( mMyRank == mMasterRank )
                {
                    message( 4, " Distributing Mesh to %u procs...",
                             ( unsigned int ) mNumberOfProcs );

                    // start the timer
                    Timer tTimer;

                    // loop over all procs
                    for( proc_t p=0; p<mNumberOfProcs; ++p )
                    {
                        this->send_submesh( mCommTable( p ) );
                    }


                    // merge tables in order to master mesh synch
                    this->merge_facet_and_connector_tables();

                    message( 4, "    ... time for distributing           : %u ms\n",
                             ( unsigned int ) tTimer.stop() );

                }
                else
                {
                    this->receive_submesh();
                }
            }

            // tidy up memory
            this->delete_maps();
        }

//------------------------------------------------------------------------------

        void
        Kernel::send_submesh( const proc_t & aTarget )
        {
            if ( aTarget == mMyRank ) return;

            BELFEM_ASSERT( mMyRank == mMasterRank,
                           "send_submesh() must only be called by maste proc only" );

            // get ref to all elements on mesh
            Cell< mesh::Element * > & tElements = mMesh->elements();
            Cell< mesh::Facet * > & tFacets = mMesh->facets();

            // unflag all nodes and elements on mesh
            mMesh->unflag_all_nodes();
            mMesh->unflag_all_elements();
            mMesh->unflag_all_edges();
            mMesh->unflag_all_faces();
            mMesh->unflag_all_facets();
            mMesh->unflag_all_vertices();

            // loop over all elements
            for ( mesh::Element * tElement: tElements )
            {
                // flag element if it belongs to this proc
                if ( tElement->owner() == aTarget )
                {
                    tElement->flag();
                }
            }

            for ( mesh::Facet * tFacet: tFacets )
            {
                if ( tFacet->owner() == aTarget )
                {
                    tFacet->flag();
                }
            }

            // now also flag all element neighbors in order to create the aura
            for ( mesh::Element * tElement: tElements )
            {
                if ( tElement->is_flagged())
                {
                    // get number of neighbors
                    index_t tNumNeighbors = tElement->number_of_elements();

                    // flag neighbors
                    for ( index_t e = 0; e < tNumNeighbors; ++e )
                    {
                        tElement->element( e )->flag();
                    }
                }
            }

            // flag nodes of flagged elements
            for ( mesh::Element * tElement : tElements )
            {
                if( tElement->is_flagged() )
                {
                    // flag all nodes
                    tElement->flag_nodes();
                }
            }

            // flag edges, if they exist
            if( mMesh->edges_exist() )
            {
                mMesh->unflag_all_edges() ;

                // flag edges of flagged elements
                for ( mesh::Element * tElement : tElements )
                {
                    if( tElement->is_flagged() )
                    {
                        // flag all edges
                        tElement->flag_edges() ;
                    }
                }

                for( mesh::Facet * tFacet : tFacets )
                {
                    if( tFacet->is_flagged() )
                    {
                        tFacet->element()->flag_edges() ;
                    }
                }
            }

            // flag faces, if they exist
            if( mMesh->faces_exist() )
            {
                mMesh->unflag_all_faces() ;

                // flag faces of flagged elements
                for ( mesh::Element * tElement : tElements )
                {
                    if( tElement->is_flagged() )
                    {
                        // flag all edges
                        tElement->flag_faces() ;
                    }
                }
                for( mesh::Facet * tFacet : tFacets )
                {
                    if( tFacet->is_flagged() )
                    {
                        tFacet->element()->flag_faces() ;
                    }
                }
            }


            // flag facets for submesh
            const Vector< id_t > & tGhostSideSetIDs = mMesh->ghost_sideset_ids() ;

            for( id_t tID : tGhostSideSetIDs )
            {
                // grab facets
                Cell< mesh::Facet* > & tFacets = mMesh->sideset( tID )->facets() ;

                // loop over all facets
                for( mesh::Facet * tFacet : tFacets )
                {
                    // check owner
                    if( tFacet->owner() == aTarget )
                    {
                        tFacet->flag() ;
                        tFacet->flag_nodes() ;
                    }
                }

                if( mMesh->edges_exist() )
                {
                    for( mesh::Facet * tFacet : tFacets )
                    {
                        // check owner
                        if( tFacet->owner() == aTarget )
                        {
                            tFacet->element()->flag_edges() ;
                        }
                    }
                }

                if( mMesh->faces_exist() && mMesh->number_of_dimensions() == 3 )
                {
                    for( mesh::Facet * tFacet : tFacets )
                    {
                        // check owner
                        if( tFacet->owner() == aTarget )
                        {
                            tFacet->face()->flag() ;
                        }
                    }
                }
            }

            Vector< index_t > tNumEntities( 8 );
            tNumEntities( 0 ) = mMesh->number_of_dimensions() ;
            tNumEntities( 1 ) = mMesh->number_of_nodes() ;
            tNumEntities( 2 ) = mMesh->number_of_elements() ;
            tNumEntities( 3 ) = mMesh->number_of_edges() ;
            tNumEntities( 4 ) = mMesh->number_of_faces() ;
            tNumEntities( 5 ) = mMesh->number_of_facets() ;
            tNumEntities( 6 ) = mMesh->number_of_connectors() ;
            tNumEntities( 7 ) = mMesh->vertices().size() ;


            send( aTarget, tNumEntities );


            this->send_nodes( aTarget );
            this->send_elements( aTarget );

            if ( tNumEntities( 3 ) > 0 ) this->send_edges( aTarget );
            if ( tNumEntities( 4 ) > 0 ) this->send_faces( aTarget );
            if ( tNumEntities( 5 ) > 0 ) this->send_facets( aTarget, false );
            if ( tNumEntities( 6 ) > 0 ) this->send_facets( aTarget, true );

            if( tNumEntities( 7 ) > 0 )
            {
                // get ref to all vertices on mesh
                Cell< mesh::Element * > & tVertices = mMesh->vertices();

                if ( tVertices.size() > 0 )
                {
                    mMesh->unflag_all_nodes();

                    // only flag nodes that belong to the target
                    for ( mesh::Element * tElement: tElements )
                    {
                        if ( tElement->is_flagged() && tElement->owner() == aTarget )
                        {
                            // flag all nodes
                            tElement->flag_nodes();
                        }
                    }


                    // flag used vertices
                    for ( mesh::Element * tVertex: tVertices )
                    {
                        if ( tVertex->node( 0 )->is_flagged())
                        {
                            tVertex->flag();
                        }
                    }
                }

                this->send_vertices( aTarget );
            }

            if( tNumEntities( 3 ) > 0 )
            {
                // we can't send const data with the current version of
                // the MPI wrapper
                Vector< id_t > tNedelecBlocks( mMesh->nedelec_blocks());
                Vector< id_t > tNedelecSideSets( mMesh->nedelec_sidesets());

                // send nedelec data
                send( aTarget, tNedelecBlocks );
                send( aTarget, tNedelecSideSets );
            }

            mMesh->distribute_ghost_sidesets( aTarget, mMyRank );
        }

//------------------------------------------------------------------------------

        void
        Kernel::receive_submesh()
        {
            // receive dimension
            Vector< index_t > tNumEntities( 8 );
            receive( mMasterRank, tNumEntities );

            uint tDimension = tNumEntities( 0 );

            // create a submesh
            mSubMesh = new Mesh( tDimension, mMyRank );

            // flag telling that this mesh is linked to a kernel
            mSubMesh->set_kernel_flag();

            // get nodes
            this->receive_nodes();

            // get elements
            this->receive_elements();

            if ( tNumEntities( 3 ) > 0 ) this->receive_edges();

            if ( tNumEntities( 4 ) > 0 ) this->receive_faces();

            // get facets
            if ( tNumEntities( 5 ) > 0 ) this->receive_facets( false );

            // get connectors
            if ( tNumEntities( 6 ) > 0 ) this->receive_facets( true );

            mSubMesh->set_facets_are_linked_flag();

            if ( tNumEntities( 7 ) > 0 ) this->receive_vertices();

            // blocks that contain edge elements
            Vector< id_t > tNedelecBlocks;
            // sidesets that contain edge elements
            Vector< id_t > tNedelecSideSets;

            if ( tNumEntities( 3 ) > 0 )
            {
                receive( mMasterRank, tNedelecBlocks );
                receive( mMasterRank, tNedelecSideSets );
            }

            mSubMesh->set_number_of_partitions( mNumberOfProcs ) ;

            mSubMesh->finalize();

            if( tNumEntities( 3 ) > 0 )
            {
                // collect elements with edge data
                Cell< mesh::Element * > tElements;

                this->collect_elements( tNedelecBlocks, tNedelecSideSets, tElements );

                mSubMesh->finalize_edges( tElements );
            }

            if( tNumEntities( 4 ) > 0 ) mSubMesh->finalize_faces() ;

            mSubMesh->distribute_ghost_sidesets( mMyRank, this->master() );

        }

//------------------------------------------------------------------------------

        void
        Kernel::send_nodes( const proc_t & aTarget )
        {
            Cell< mesh::Node    * > & tNodes = mMesh->nodes();

            // check if cuts exist
            if( mMesh->number_of_connectors() > 0 )
            {
                Cell< mesh::Facet * > & tConnectors = mMesh->connectors() ;

                for( mesh::Facet * tConnector : tConnectors )
                {
                    if( tConnector->node( 0 )->is_flagged() )
                    {
                        tConnector->node( 1 )->flag() ;
                    }
                    else if ( tConnector->node( 1 )->is_flagged() )
                    {
                        tConnector->node( 0 )->flag() ;
                    }
                }
            }

            for( mesh::Facet * tFacet : mMesh->facets() )
            {
                for( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    if ( tFacet->node( k )->is_flagged() )
                    {
                        tFacet->flag_nodes() ;
                        continue;
                    }
                }
            }

            // also flag duplicates
            for( mesh::Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    for( uint k=0 ; k<tNode->number_of_duplicates(); ++k )
                    {
                        tNode->duplicate( k )->flag() ;
                    }
                }
            }

            // count nodes
            index_t tCount = 0;
            for( mesh::Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    ++tCount;
                }
            }

            // allocate data containers
            Vector< id_t >   tIDs( tCount );
            Vector< proc_t > tOwners( tCount );
            Vector< real >   tX( tCount );
            Vector< real >   tY( tCount );
            Vector< real >   tZ( tCount );

            Vector< index_t > & tIndices = mNodeTable( mCommMap( aTarget ) );
            tIndices.set_size( tCount );

            // populate data
            tCount = 0;

            for( mesh::Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    tIDs( tCount ) = tNode->id();
                    tOwners( tCount )  = tNode->owner();
                    tX( tCount ) = tNode->x();
                    tY( tCount ) = tNode->y();
                    tZ( tCount ) = tNode->z();

                    tIndices( tCount ) = tNode->index();
                    ++tCount;
                }
            }


            // send data
            send( aTarget, tIDs );
            send( aTarget, tOwners );
            send( aTarget, tX );
            send( aTarget, tY );
            send( aTarget, tZ );

            /*for( uint s=0; s<2; ++s )
            {
                // count nodes that are cut
                tCount = 0;

                // get the cut table
                Matrix< id_t > & tNodeTable = s==0 ? mMesh->node_cut_table() : mMesh->node_tape_table() ;

                // number of nodes in table
                index_t tNumNodes = tNodeTable.n_cols();

                for ( index_t k = 0; k < tNumNodes; ++k )
                {
                    // grab node
                    mesh::Node * tNode = mMesh->node( tNodeTable( 0, k ));
                    if ( tNode->is_flagged())
                    {
                        ++tCount;
                    }
                }

                Matrix< id_t > tNodeSubTable( 2, tCount );
                tCount = 0;
                for ( index_t k = 0; k < tNumNodes; ++k )
                {
                    // grab node
                    mesh::Node * tNode = mMesh->node( tNodeTable( 0, k ));
                    if ( tNode->is_flagged())
                    {
                        tNodeSubTable( 0, tCount ) = tNodeTable( 0, k );
                        tNodeSubTable( 1, tCount ) = tNodeTable( 1, k );
                        ++tCount;
                    }
                }


                send( aTarget, tNodeSubTable );
            }*/



            // check if nodes have been set
            if( tCount > 0 )
            {
                // count duplicates
                tCount = 0 ;
                for( mesh::Node * tNode : tNodes )
                {
                    if( tNode->is_flagged() )
                    {
                        tCount += tNode->number_of_duplicates() + 1 ;
                    }
                }

                // populate duplicates
                tIDs.set_size( tCount, 0 );
                tCount = 0;
                for ( mesh::Node * tNode: tNodes )
                {
                    if ( tNode->is_flagged())
                    {
                        tIDs( tCount++ ) = tNode->number_of_duplicates();

                        for ( uint d = 0; d < tNode->number_of_duplicates(); ++d )
                        {
                            tIDs( tCount++ ) = tNode->duplicate( d )->id();
                        }
                    }
                }

                send( aTarget, tIDs );
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::send_edges( const proc_t & aTarget )
        {
            Cell< mesh::Edge    * > & tEdges = mMesh->edges();

            index_t tNumberOfEdges = tEdges.size() ;

            // send number of edges to other procs
            broadcast( mMasterRank, tNumberOfEdges );

            if (  mMesh->edges_exist() )
            {
                // count edges
                index_t tEdgeCount = 0 ;
                index_t tNodeCount = 0 ;


                for ( mesh::Edge * tEdge : tEdges )
                {
                    if ( tEdge->is_flagged() )
                    {
                        ++tEdgeCount ;
                        tNodeCount += tEdge->number_of_nodes() ;
                    }
                }

                // edge IDs
                Vector< id_t > tEdgeIDs( tEdgeCount );

                // edge Owners
                Vector< proc_t > tEdgeOwners( tEdgeCount );

                // number of nodes per edge
                Vector< uint > tNumNodes( tEdgeCount );

                // node IDs
                Vector< id_t > tNodeIDs( tNodeCount );


                Vector< index_t > & tIndices = mEdgeTable( mCommMap( aTarget ) );
                tIndices.set_size( tEdgeCount );

                // reset counters
                tEdgeCount = 0;
                tNodeCount = 0 ;

                // populate vectors
                for ( mesh::Edge * tEdge : tEdges )
                {
                    if ( tEdge->is_flagged() )
                    {

                        // set value for table
                        tIndices( tEdgeCount )    = tEdge->index() ;

                        tEdgeIDs( tEdgeCount )    = tEdge->id() ;
                        tEdgeOwners( tEdgeCount ) = tEdge->owner() ;
                        tNumNodes( tEdgeCount++ ) = tEdge->number_of_nodes() ;

                        for( uint k=0; k<tEdge->number_of_nodes(); ++k )
                        {
                            tNodeIDs( tNodeCount++ ) = tEdge->node( k )->id() ;
                        }
                    }
                }

                // reset edge counter
                tEdgeCount = 0 ;

                // get element container
                Cell< mesh::Element * > & tElements = mMesh->elements() ;

                for( mesh::Element * tElement : tElements )
                {
                    if( tElement->is_flagged() && tElement->has_edges() )
                    {
                        ++tEdgeCount;
                        tEdgeCount += tElement->number_of_edges() ;
                    }
                }

                // number of edges per element
                Vector< id_t > tEdgesPerElements( tEdgeCount );

                // reset counters
                tEdgeCount = 0 ;

                for( mesh::Element * tElement : tElements )
                {
                    if( tElement->is_flagged() && tElement->has_edges() )
                    {
                        tEdgesPerElements( tEdgeCount++ ) = tElement->id();

                        uint tN = tElement->number_of_edges() ;
                        for( uint e=0; e<tN; ++e )
                        {
                            tEdgesPerElements( tEdgeCount++ ) = tElement->edge( e )->id() ;
                        }
                    }
                }

                // send data to target
                send( aTarget, tEdgeIDs );
                send( aTarget, tEdgeOwners );
                send( aTarget, tNumNodes );
                send( aTarget, tNodeIDs );
                send( aTarget, tEdgesPerElements );

            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::receive_edges()
        {
            // check if edges exist
            index_t tNumberOfEdges = 0 ;
            broadcast( mMasterRank, tNumberOfEdges );

            if( tNumberOfEdges > 0 )
            {
                // edge IDs
                Vector< id_t > tEdgeIDs ;
                receive( mMasterRank, tEdgeIDs );

                // edge owners
                Vector< proc_t > tEdgeOwners ;
                receive( mMasterRank, tEdgeOwners );

                // number of nodes per edge
                Vector< uint > tNumNodes ;
                receive( mMasterRank, tNumNodes );

                // node IDs
                Vector< id_t > tNodeIDs ;
                receive( mMasterRank, tNodeIDs );

                // edge ids per element
                Vector< id_t > tEdgesPerElements ;
                receive( mMasterRank, tEdgesPerElements );

                // initialize counters
                index_t tCount = 0 ;

                // update number of edges
                tNumberOfEdges = tEdgeIDs.length() ;

                // grab container with edges
                Cell< mesh::Edge    * > & tEdges = mSubMesh->edges();

                // allocate memory
                tEdges.set_size( tNumberOfEdges, nullptr );

                for( index_t e=0; e<tNumberOfEdges; ++e )
                {
                    // create new edge
                    mesh::Edge * tEdge = new mesh::Edge();

                    // set basic information
                    tEdge->set_id( tEdgeIDs( e ) );
                    tEdge->set_index( e );
                    tEdge->set_owner( tEdgeOwners( e ) );

                    // link edge with nodes
                    tEdge->allocate_node_container( tNumNodes( e ) );


                    for( uint k=0; k<tNumNodes( e ); ++k )
                    {
                        tEdge->insert_node( mNodeMap[ tNodeIDs( tCount++) ], k );
                    }


                    // copy edge into container
                    tEdges( e ) = tEdge ;

                    // copy edge into map
                    mEdgeMap[ tEdge->id() ] = tEdge ;
                }

                index_t tVectorLength = tEdgesPerElements.length() ;

                // reset counter
                tCount = 0;

                while ( tCount < tVectorLength )
                {
                    // grab element
                    mesh::Element * tElement = mElementMap( tEdgesPerElements( tCount++ ) );

                    // allocate memory
                    tElement->allocate_edge_container();

                    // get number of edges from element
                    uint tNumEdges = tElement->number_of_edges();

                    for ( uint k = 0; k < tNumEdges; ++k )
                    {
                        // link edge
                        tElement->insert_edge( mEdgeMap( tEdgesPerElements( tCount++ ) ), k );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::send_faces( const proc_t & aTarget )
        {
            Cell< mesh::Face * > & tFaces = mMesh->faces();

            index_t tNumberOfFaces = tFaces.size();

            // send number of edges to other procs
            broadcast( mMasterRank, tNumberOfFaces );

            if ( mMesh->faces_exist() )
            {
                // count faces
                index_t tCount = 0;

                for ( mesh::Face * tFace: tFaces )
                {
                    if ( tFace->is_flagged() )
                    {
                        ++tCount;
                    }
                }

                // face IDs
                Vector< id_t > tFaceIDs( tCount );
                Vector< id_t > tMasterIDs( tCount );

                Vector< index_t > & tIndices = mFaceTable( mCommMap( aTarget ) );
                tIndices.set_size( tCount );

                // reset counters
                tCount = 0;

                if( mMesh->number_of_dimensions() == 2 )
                {
                    for ( mesh::Face * tFace: tFaces )
                    {
                        if ( tFace->is_flagged() )
                        {
                            tFaceIDs( tCount ) = tFace->id();

                            tMasterIDs( tCount ) = tFace->master()->id();

                            // set value for table
                            tIndices( tCount++ ) = tFace->index();
                        }
                    }

                    send( aTarget, tFaceIDs );
                    send( aTarget, tMasterIDs );

                }
                else if( mMesh->number_of_dimensions() == 3 )
                {
                    Vector< id_t > tSlaveIDs( tCount );

                    Vector< uint > tIndexOnMaster( tCount );
                    Vector< uint > tIndexOnSlave( tCount );


                    for ( mesh::Face * tFace: tFaces )
                    {
                        if ( tFace->is_flagged() )
                        {
                            tFaceIDs( tCount ) = tFace->id();

                            tMasterIDs( tCount ) = tFace->master()->id();
                            tIndexOnMaster( tCount ) = tFace->index_on_master();

                            if ( tFace->slave() != nullptr )
                            {
                                tSlaveIDs( tCount ) = tFace->slave()->id();
                                tIndexOnSlave( tCount ) = tFace->index_on_slave();
                            }
                            else
                            {
                                tSlaveIDs( tCount ) = gNoID;
                                tIndexOnSlave( tCount ) = gNoIndex;
                            }

                            // set value for table
                            tIndices( tCount++ ) = tFace->index();
                        }
                    }

                    send( aTarget, tFaceIDs );
                    send( aTarget, tMasterIDs );
                    send( aTarget, tIndexOnMaster );
                    send( aTarget, tSlaveIDs );
                    send( aTarget, tIndexOnSlave );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::receive_faces()
        {
            // check if edges exist
            index_t tNumberOfFaces = 0 ;
            broadcast( mMasterRank, tNumberOfFaces );

            if( tNumberOfFaces > 0 )
            {
                Cell< mesh::Face * > & tFaces = mSubMesh->faces();

                Vector< id_t > tFaceIDs ;
                receive( mMasterRank, tFaceIDs );

                Vector< id_t > tMasterIDs ;
                receive( mMasterRank, tMasterIDs );

                // get the number of faces
                tNumberOfFaces = tFaceIDs.length() ;

                tFaces.set_size( tNumberOfFaces, nullptr );

                if( mMesh->number_of_dimensions() == 2 )
                {
                    // create faces
                    for( index_t f=0; f<tNumberOfFaces; ++f )
                    {
                        // get master element
                        mesh::Element * tParent = mElementMap( tMasterIDs( f ) );

                        mesh::Face * tFace = new mesh::Face( tParent );

                        tFace->set_id( tFaceIDs( f ) );
                        tFace->set_index( f );
                        tFace->set_owner( tParent->owner() );

                        // copy face into container
                        tFaces( f ) = tFace ;

                        // copy face into map
                        mFaceMap[ tFaceIDs( f ) ] = tFace ;

                        // link face to parent
                        tParent->allocate_face_container() ;
                        tParent->insert_face( tFace, 0 );

                    }
                }
                else if( mMesh->number_of_dimensions() == 3 )
                {
                    Vector< uint > tIndexOnMaster ;
                    receive( mMasterRank, tIndexOnMaster );

                    Vector< id_t > tSlaveIDs ;
                    receive( mMasterRank, tSlaveIDs );

                    Vector< uint > tIndexOnSlave ;
                    receive( mMasterRank, tIndexOnSlave );

                    // create faces
                    for( index_t f=0; f<tNumberOfFaces; ++f )
                    {
                        // get master element
                        mesh::Element * tMaster =  mElementMap( tMasterIDs( f ) );

                        // get slave element
                        mesh::Element * tSlave = tSlaveIDs( f ) == gNoID ? nullptr : mElementMap( tSlaveIDs( f ) );

                        mesh::Face * tFace = new mesh::Face( tMaster, tIndexOnMaster( f ), tSlave, tIndexOnSlave( f ) );

                        tFace->set_id( tFaceIDs( f ) );
                        tFace->set_index( f );
                        tFace->set_owner( tMaster->owner() );

                        // copy face into container
                        tFaces( f ) = tFace ;

                        // copy face into map
                        mFaceMap[ tFaceIDs( f ) ] = tFace ;

                        // link face with master
                        if( ! tMaster->has_faces() )
                        {
                            tMaster->allocate_face_container() ;
                        }
                        tMaster->insert_face( tFace, tIndexOnMaster( f ) );

                        // link face with slave if it exists
                        if( tSlave != nullptr )
                        {
                            if( ! tSlave->has_faces() )
                            {
                                tSlave->allocate_face_container() ;
                                tSlave->insert_face( tFace, tIndexOnSlave( f ) );
                            }
                        }

                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::receive_nodes()
        {
            // data containers
            Vector< id_t >   tIDs;
            Vector< proc_t > tOwners;
            Vector< real >   tX;
            Vector< real >   tY;
            Vector< real >   tZ;

            Cell< mesh::Node    * > & tNodes = mSubMesh->nodes();

            // get number of nodes
            receive( mMasterRank, tIDs );
            receive( mMasterRank, tOwners );
            receive( mMasterRank, tX );
            receive( mMasterRank, tY );
            receive( mMasterRank, tZ );
            //receive( mMasterRank, mSubMesh->node_cut_table() );
            //receive( mMasterRank, mSubMesh->node_tape_table() );

            index_t tNumberOfNodes = tIDs.length();

            if( tNumberOfNodes > 0 )
            {
                tNodes.set_size( tNumberOfNodes, nullptr );

                for ( index_t k = 0; k < tNumberOfNodes; ++k )
                {
                    mesh::Node * tNode = new mesh::Node(
                            tIDs( k ), tX( k ), tY( k ), tZ( k ));

                    tNode->set_owner( tOwners( k ));

                    mNodeMap[ tIDs( k ) ] = tNode;

                    tNodes( k ) = tNode;
                }

                // receive duplicates
                receive( mMasterRank, tIDs );
                index_t tCount = 0;
                for ( mesh::Node * tNode: tNodes )
                {
                    // get number of duplicates
                    uint tN = tIDs( tCount++ );
                    if( tN > 0 )
                    {
                        tNode->allocate_duplicate_container( tN );
                        for( uint d=0; d<tN; ++d )
                        {
                            tNode->add_duplicate( mNodeMap( tIDs( tCount++ ) ) );
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::send_elements( const proc_t & aTarget )
        {
            // loop over all blocks
            uint tNumBlocks = mMesh->number_of_blocks();

            // send number of blocks
            send( aTarget, tNumBlocks );

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            index_t tAllElementCount = 0 ;

            for( uint b=0; b<tNumBlocks; ++b )
            {
                Cell< mesh::Element * > & tElements = tBlocks( b )->elements();

                // count elements
                index_t tElementCount = 0;
                for( mesh::Element * tElement : tElements )
                {
                    if( tElement->is_flagged() )
                    {
                        ++tElementCount;
                    }
                }
                tAllElementCount += tElementCount ;

                // send number of elements
                send( aTarget, tElementCount );

                if( tElementCount > 0 )
                {
                    ElementType tElementType = tBlocks( b )->element_type();
                    uint tNumNodesPerElement = mesh::number_of_nodes( tElementType );

                    // containers
                    Vector<id_t> tIDs( tElementCount );
                    Vector<uint> tGeometryTags( tElementCount );
                    Vector<uint> tPhysicalTags( tElementCount );
                    Vector<proc_t> tOwners( tElementCount );
                    Vector<id_t> tNodes( tNumNodesPerElement * tElementCount );

                    // populate data
                    tElementCount = 0;
                    index_t tNodeCount = 0;

                    for ( mesh::Element * tElement : tElements )
                    {
                        if ( tElement->is_flagged() )
                        {
                            tIDs( tElementCount ) = tElement->id();
                            tGeometryTags( tElementCount ) = tElement->geometry_tag();
                            tPhysicalTags( tElementCount ) = tElement->physical_tag();
                            tOwners( tElementCount ) = tElement->owner();

                            ++tElementCount;

                            for ( uint k = 0; k < tNumNodesPerElement; ++k )
                            {
                                tNodes( tNodeCount++ ) = tElement->node( k )->id();
                            }
                        }
                    }

                    // send block id
                    uint tID = tBlocks( b )->id();
                    send( aTarget, tID );

                    // send element type
                    int tType = static_cast< int >( tElementType );
                    send( aTarget, tType );

                    // send IDs
                    send( aTarget, tIDs );

                    // send Geometry tags
                    send( aTarget, tGeometryTags );

                    // send physical tags
                    send( aTarget, tPhysicalTags );

                    // send Owners
                    send( aTarget, tOwners );

                    // send node ids
                    send( aTarget, tNodes );
                }
            }

            Vector< index_t > & tIndices = mElementTable( aTarget );
            tIndices.set_size( tAllElementCount );
            tAllElementCount = 0 ;

            Cell< mesh::Element * > & tElements = mMesh->elements();
            for( mesh::Element * tElement : tElements )
            {
                if( tElement->is_flagged() )
                {
                    tIndices( tAllElementCount++ ) = tElement->index() ;
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Kernel::send_vertices( const proc_t & aTarget )
        {

            Cell< mesh::Element * > & tVertices = mMesh->vertices();

            // count elements
            index_t tVertexCount = 0;
            for ( mesh::Element * tVertex : tVertices )
            {
                if ( tVertex->is_flagged() )
                {
                    ++tVertexCount;
                }
            }

            // send number of vertices
            send( aTarget, tVertexCount );

            if ( tVertexCount > 0 )
            {
                // containers
                Vector< id_t > tIDs( tVertexCount );
                Vector< uint > tGeometryTags( tVertexCount );
                Vector< uint > tPhysicalTags( tVertexCount );
                Vector< proc_t > tOwners( tVertexCount );
                Vector< id_t > tNodes( tVertexCount );

                // populate data
                tVertexCount = 0;

                for ( mesh::Element * tVertex : tVertices )
                {
                    if ( tVertex->is_flagged())
                    {
                        tIDs( tVertexCount ) = tVertex->id();
                        tGeometryTags( tVertexCount ) = tVertex->geometry_tag();
                        tPhysicalTags( tVertexCount ) = tVertex->physical_tag();
                        tOwners( tVertexCount ) = tVertex->owner();
                        tNodes( tVertexCount ) = tVertex->node( 0 )->id();
                        ++tVertexCount;
                    }
                }

                // send IDs
                send( aTarget, tIDs );

                // send Geometry tags
                send( aTarget, tGeometryTags );

                // send physical tags
                send( aTarget, tPhysicalTags );

                // send Owners
                send( aTarget, tOwners );

                // send node ids
                send( aTarget, tNodes );
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::receive_elements()
        {
            // get block container
            Cell< mesh::Block * > tBlocks;

            // get number of blocks
            uint tNumBlocks = 1;
            receive( mMasterRank, tNumBlocks );

            tBlocks.set_size( tNumBlocks, nullptr );

            // the factory
            mesh::ElementFactory tFactory;

            uint tBlockCount = 0 ;

            // loop over all blocks
            for( uint b=0; b<tNumBlocks; ++b )
            {
                // get number of elements
                index_t tNumberOfElements = 0;
                receive( mMasterRank, tNumberOfElements );

                if( tNumberOfElements > 0 )
                {
                    // get id
                    uint tID = 0;
                    receive( mMasterRank, tID );

                    // get block type
                    int tType = 0 ;
                    receive( mMasterRank, tType );
                    ElementType tElementType = static_cast< ElementType >( tType );

                    // number of nodes per element
                    uint tNumNodesPerElement = mesh::number_of_nodes( tElementType );

                    // Element IDs
                    Vector<id_t> tIDs;
                    receive( mMasterRank, tIDs );

                    // Geometry Tags
                    Vector<uint> tGeometryTags;
                    receive( mMasterRank, tGeometryTags );

                    // Physical Tags
                    Vector<uint> tPhysicalTags;
                    receive( mMasterRank, tPhysicalTags );

                    // owners
                    Vector<proc_t> tOwners;
                    receive( mMasterRank, tOwners );

                    // Node IDs
                    Vector<id_t> tNodes;
                    receive( mMasterRank, tNodes );

                    mesh::Block * tBlock = new mesh::Block( tID, tNumberOfElements );

                    // initialize counter
                    index_t tCount = 0;

                    // loop over all elements
                    for ( index_t e = 0; e < tNumberOfElements; ++e )
                    {
                        // create a new element
                        mesh::Element * tElement = tFactory.create_lagrange_element(
                                tElementType, tIDs( e ));

                        // set tags
                        tElement->set_geometry_tag( tGeometryTags( e ));
                        tElement->set_physical_tag( tPhysicalTags( e ));
                        tElement->set_owner( tOwners( e ));

                        for ( uint k = 0; k < tNumNodesPerElement; ++k )
                        {
                            tElement->insert_node( mNodeMap( tNodes( tCount++ ) ), k );
                        }

                        mElementMap[ tElement->id() ] = tElement;

                        // add element to block
                        tBlock->insert_element( tElement );
                    }

                    // add block to container
                    tBlocks( tBlockCount++ ) = tBlock;
                }
            }

            // get block list on mesh
            Cell< mesh::Block * > & tMeshBlocks = mSubMesh->blocks();

            tMeshBlocks.set_size( tBlockCount, nullptr );

            // add non-empty blocks to submesh
            for( uint b=0; b<tBlockCount; ++b )
            {
                tMeshBlocks( b ) = tBlocks( b );
            }

            tBlocks.clear();
        }

//------------------------------------------------------------------------------

        void
        Kernel::receive_vertices()
        {

            // get number of elements
            index_t tNumberOfVertices = 0;
            receive( mMasterRank, tNumberOfVertices );

            if( tNumberOfVertices > 0 )
            {
                // the factory
                mesh::ElementFactory tFactory;


                // Element IDs
                Vector<id_t> tIDs;
                receive( mMasterRank, tIDs );

                // Geometry Tags
                Vector<uint> tGeometryTags;
                receive( mMasterRank, tGeometryTags );

                // Physical Tags
                Vector<uint> tPhysicalTags;
                receive( mMasterRank, tPhysicalTags );

                // owners
                Vector<proc_t> tOwners;
                receive( mMasterRank, tOwners );

                // Node IDs
                Vector<id_t> tNodes;
                receive( mMasterRank, tNodes );

                Cell< mesh::Element * > & tVertices = mSubMesh->vertices();

                tVertices.set_size( tNumberOfVertices, nullptr );


                // loop over all elements
                for ( index_t e = 0; e < tNumberOfVertices; ++e )
                {
                    // create a new element
                    mesh::Element * tVertex = tFactory.create_lagrange_element(
                            ElementType::VERTEX, tIDs( e ));

                    // set tags
                    tVertex->set_geometry_tag( tGeometryTags( e ));
                    tVertex->set_physical_tag( tPhysicalTags( e ));
                    tVertex->set_owner( tOwners( e ));
                    tVertex->insert_node( mNodeMap( tNodes ( e ) ), 0 );

                    mVertexMap[ tVertex->id() ] = tVertex;
                    tVertices( e ) = tVertex;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::send_facets( const proc_t & aTarget, const bool aConnectorSwitch )
        {
            // loop over all blockst
            uint tNumSideSets = aConnectorSwitch ?
                                mMesh->cuts().size() :
                                mMesh->number_of_sidesets() ;

            // send number of blocks
            send( aTarget, tNumSideSets );


            Cell< mesh::SideSet * > & tSideSets =
                    aConnectorSwitch ? mMesh->cuts() : mMesh->sidesets();

            // table with all indices for this proc
            // we need this to do the communication correctly
            Vector< index_t > & tAllIndices =
                    aConnectorSwitch ?
                    mConnectorTable( mCommMap( aTarget ) ) :
                    mFacetTable( mCommMap( aTarget ) );

            // get mesh dimension
            uint tDimension = mMesh->number_of_dimensions() ;

            // count all indices
            // additional counters

            index_t tFacetCount = 0 ;
            index_t tElementHasEdgesCount  = 0 ;
            index_t tElementHasFacesCount  = 0 ;
            index_t tEdgeCount = 0 ;
            index_t tFaceCount = 0 ;

            Vector< id_t > tEdgeIDs;
            Vector< id_t > tFaceIDs;
            Vector< id_t > tElementsWithEdges;
            Vector< id_t > tElementsWithFaces;
            Matrix< uint > tTags ;

            // count size of facet table
            if( aConnectorSwitch )
            {
                // connector mode
                for( uint s=0; s<tNumSideSets; ++s )
                {
                    for ( mesh::Facet * tFacet: tSideSets( s )->facets())
                    {
                        if ( tFacet->owner() == aTarget )
                        {
                            ++tFacetCount;
                        }
                    }
                }
            }
            else
            {
                // facet mode
                for( uint s=0; s<tNumSideSets; ++s )
                {
                    for ( mesh::Facet * tFacet: tSideSets( s )->facets() )
                    {
                        // grab element
                        mesh::Element * tElement = tFacet->element();

                        if ( tFacet->owner() == aTarget )
                        {
                            if ( tElement->has_edges() )
                            {
                                ++tElementHasEdgesCount;
                                tEdgeCount += tElement->number_of_edges();
                            }
                            if ( tElement->has_faces() && tDimension == 3 )
                            {
                                ++tElementHasFacesCount;
                                tFaceCount += tElement->number_of_faces();
                            }
                            ++tFacetCount;
                        }
                    }
                } // end loop over all sidesets

                // set sized for containers
                if ( tElementHasEdgesCount > 0 )
                {
                    tElementsWithEdges.set_size( tElementHasEdgesCount, 0 );
                    tEdgeIDs.set_size( tEdgeCount, 0 );
                    tElementHasEdgesCount  = 0 ;
                    tEdgeCount = 0 ;
                }
                if ( tElementHasFacesCount > 0 )
                {
                    tElementsWithFaces.set_size( tElementHasFacesCount, 0 );
                    tFaceIDs.set_size( tFaceCount, 0 );
                    tElementHasFacesCount  = 0 ;
                    tFaceCount = 0 ;
                }
            }

            // set tags
            tTags.set_size( 2, tFacetCount );

            tAllIndices.set_size( tFacetCount );
            tFacetCount = 0 ;


            for( uint s=0; s<tNumSideSets; ++s )
            {
                // flag facet if it touches target
                // reset counter
                index_t tCount = 0;

                // unflag all facets
                for( mesh::Facet * tFacet : tSideSets( s )->facets()  )
                {
                    tFacet->element()->unflag();
                }

                tCount=0;

                if( aConnectorSwitch )
                {
                    // connector mode
                    for ( mesh::Facet * tFacet: tSideSets( s )->facets() )
                    {
                        if( tFacet->owner() == aTarget )
                        {
                            tFacet->element()->flag() ;
                            ++tCount;
                        }
                    }
                }
                else
                {
                    // facet mode
                    for ( mesh::Facet * tFacet: tSideSets( s )->facets() )
                    {
                        if ( tFacet->owner() == aTarget )
                        {
                            // grab element
                            mesh::Element * tElement = tFacet->element();

                            // flag element
                            tElement->flag();
                            if ( tElement->has_edges() )
                            {
                                tElementsWithEdges( tElementHasEdgesCount++ ) = tFacet->id();

                                for ( uint e = 0; e < tElement->number_of_edges(); ++e )
                                {
                                    tEdgeIDs( tEdgeCount++ ) = tElement->edge( e )->id();
                                }
                            }
                            if ( tElement->has_faces() && tDimension == 3 )
                            {
                                tElementsWithFaces( tElementHasFacesCount++ ) = tFacet->id();

                                for ( uint f = 0; f < tElement->number_of_faces(); ++f )
                                {
                                    tFaceIDs( tFaceCount++ ) = tElement->face( f )->id();
                                }
                            }

                            ++tCount;
                        }
                    }
                }

                Vector< index_t > tIndices( tCount );

                // collect facets
                Cell< mesh::Facet * > tFacets( tCount, nullptr );
                tCount = 0;
                index_t tNodeCount = 0 ;

                for( mesh::Facet * tFacet : tSideSets( s )->facets() )
                {
                    if( tFacet->element()->is_flagged() )
                    {
                        // the vector from the current sideset or cut
                        tIndices( tCount ) = tFacet->index();

                        // add this index to the communication vector
                        tAllIndices( tFacetCount++ ) = tFacet->index() ;
                        tFacets( tCount ) = tFacet;
                        tFacet->element()->unflag();
                        tNodeCount += tFacet->element()->number_of_nodes() ;

                        // copy the tags
                        tTags( 0, tCount ) = tFacet->element()->geometry_tag() ;
                        tTags( 1, tCount ) = tFacet->element()->physical_tag() ;

                        // increment the counter
                        ++tCount ;
                    }
                }

                // send counter
                send( aTarget, tCount );

                if( tCount > 0 )
                {
                    // send the ID of this sideset
                    id_t tID = tSideSets( s )->id();
                    send( aTarget, tID );

                    // initialize id vectors
                    Vector<id_t> tIDs( tCount );

                    // populate IDs
                    for ( index_t k = 0; k < tCount; ++k )
                    {
                        tIDs( k ) = tFacets( k )->element()->id();
                    }

                    // communicate IDs
                    send( aTarget, tIDs );

                    // communicate tags
                    send( aTarget, tTags );

                    if( aConnectorSwitch )
                    {
                        // collect IDs from first node
                        for ( index_t k = 0; k < tCount; ++k )
                        {
                            tIDs( k ) = tFacets( k )->node( 0 )->id() ;
                        }
                        send( aTarget, tIDs );

                        // collect IDs from second node
                        for ( index_t k = 0; k < tCount; ++k )
                        {
                            tIDs( k ) = tFacets( k )->node( 1 )->id() ;
                        }
                        send( aTarget, tIDs );

                        // collect block IDS
                        for ( index_t k = 0; k < tCount; ++k )
                        {
                            tIDs( k ) = tFacets( k )->element()->block_id() ;
                        }
                        send( aTarget, tIDs );
                    }
                    else
                    {
                        // collect master IDs
                        for ( index_t k = 0; k < tCount; ++k )
                        {
                            tIDs( k ) = tFacets( k )->master()->id();
                        }

                        send( aTarget, tIDs );

                        // collect slave IDs
                        for ( index_t k = 0; k < tCount; ++k )
                        {
                            if ( tFacets( k )->has_slave() )
                            {
                                tIDs( k ) = tFacets( k )->slave()->id();
                            }
                            else
                            {
                                tIDs( k ) = 0;
                            }
                        }

                        send( aTarget, tIDs );

                        Vector< uint > tFacetIndices( tCount );

                        // collect master indices
                        for ( index_t k = 0; k < tCount; ++k )
                        {
                            tFacetIndices( k ) = tFacets( k )->master_index();
                        }
                        send( aTarget, tFacetIndices );

                        // collect slave indices
                        for ( index_t k = 0; k < tCount; ++k )
                        {
                            tFacetIndices( k ) = tFacets( k )->slave_index();
                        }
                        send( aTarget, tFacetIndices );

                        // collect node ids
                        tIDs.set_size( tNodeCount );

                        tNodeCount = 0 ;
                        for ( index_t k = 0; k < tCount; ++k )
                        {
                            mesh::Element * tElement = tFacets( k )->element() ;
                            for( uint i=0; i<tElement->number_of_nodes(); ++i )
                            {
                                tIDs( tNodeCount++ ) = tElement->node( i )->id() ;
                            }
                        }
                        send( aTarget, tIDs );

                    }
                }
            } // end loop over sidesets

            if( ! aConnectorSwitch )
            {

                send( aTarget, tElementsWithEdges );
                send( aTarget, tEdgeIDs );
                send( aTarget, tElementsWithFaces );
                send( aTarget, tFaceIDs );
            }
        }


//------------------------------------------------------------------------------

        void
        Kernel::receive_facets( const bool aConnectorSwitch )
        {
            // get sideset container
            Cell< mesh::SideSet * > tSideSets ;

            // get number of blocks
            uint tNumSideSets = 0;
            receive( mMasterRank, tNumSideSets );
            tSideSets.set_size( tNumSideSets, nullptr );

            // sideset counter
            uint tCount = 0;

            // element factory
            mesh::ElementFactory tFactory;

            // loop over all sidesets
            for( uint s=0; s<tNumSideSets; ++s )
            {
                // get number of facets
                index_t tNumFacets = 0;
                receive( mMasterRank, tNumFacets );

                if( tNumFacets > 0 )
                {
                    // get the ID of this sideset
                    id_t tID;
                    receive( mMasterRank, tID );

                    // IDs for the elements
                    Vector< id_t > tIDs;
                    receive( mMasterRank, tIDs );

                    // id for tags
                    Matrix< uint > tTags ;
                    receive( mMasterRank, tTags );

                    // Master IDs  or first node ID
                    Vector< id_t > tMasteIDs;

                    // Slave IDs ( 0 if no slave ) or second node ID
                    Vector< id_t > tSlaveIDs;

                    // block ids, needed for connector elements
                    Vector< id_t > tBlockIDs ;

                    // Master Facet Indices
                    Vector< uint > tMasterIndices;

                    // Slave Facet Indices
                    Vector< uint > tSlaveIndices;

                    // facets that have elements
                    Vector< id_t > tElementsWithEdges ;

                    // edges on facets
                    Vector< id_t > tEdgeIDs ;

                    // facets that have faces
                    Vector< id_t > tElementsWithFaces ;

                    // faces on facets
                    Vector< id_t > tFaceIDs ;

                    // node ids
                    Vector< id_t > tNodeIDs ;

                    receive( mMasterRank, tMasteIDs );

                    receive( mMasterRank, tSlaveIDs );

                    if( aConnectorSwitch )
                    {
                        receive( mMasterRank, tBlockIDs );
                    }
                    else
                    {
                        receive( mMasterRank, tMasterIndices );

                        receive( mMasterRank, tSlaveIndices );

                        receive( mMasterRank, tNodeIDs );
                    }

                    // create a new sideset
                    mesh::SideSet * tSideSet
                        = new mesh::SideSet( tID, tNumFacets );

                    // loop over facets
                    if( aConnectorSwitch )
                    {
                        for ( index_t f = 0; f < tNumFacets; ++f )
                        {
                            // create a new element
                            mesh::Element * tElement =
                                    tFactory.create_lagrange_element( ElementType::LINE2, tIDs( f ) );

                            tElement->insert_node( mNodeMap( tMasteIDs( f ) ), 0 );
                            tElement->insert_node( mNodeMap( tSlaveIDs( f ) ), 1 );

                            // set the tags of the element
                            tElement->set_geometry_tag( tTags( 0, f ) );
                            tElement->set_physical_tag( tTags( 1, f ) );

                            // create a facet
                            mesh::Facet * tFacet = new mesh::Facet( tElement );

                            // claim ownership
                            tFacet->set_owner( mMyRank );

                            // set block id
                            tFacet->element()->set_block_id( tBlockIDs( f ) );

                            // save facet in map
                            //mConnectorMap[ tFacet->id() ] = tFacet;

                            // add facet to cut
                            tSideSet->insert_facet( tFacet );
                        }
                    }
                    else
                    {
                        index_t tNodeCount = 0 ;

                        for ( index_t f = 0; f < tNumFacets; ++f )
                        {
                            // facet master
                            mesh::Element * tMaster = mElementMap( tMasteIDs( f ));

                            // create a new element
                            mesh::Element * tElement =
                                    tFactory.create_lagrange_element( mesh::element_type_of_facet(
                                            tMaster->type(), tMasterIndices( f )), tIDs( f ));

                            // note: this is nice but fails if we introduce thin shells
                            //       therefore, we have to communicate the node coordinates
                            // Cell< mesh::Node * > tNodes;
                            // tMaster->get_nodes_of_facet( tMasterIndices( f ), tNodes );

                            // link nodes
                            //uint tNumNodes = tNodes.size();
                            for ( uint k = 0; k < tElement->number_of_nodes(); ++k )
                            {
                                tElement->insert_node( mNodeMap( tNodeIDs( tNodeCount++ ) ), k );
                            }


                            // set the tags of the element
                            tElement->set_geometry_tag( tTags( 0, f ) );
                            tElement->set_physical_tag( tTags( 1, f ) );

                            // create a facet
                            mesh::Facet * tFacet = new mesh::Facet( tElement );

                            // claim ownership
                            tFacet->set_owner( mMyRank );

                            // link facet to master, but do not relink nodes
                            tFacet->set_master( tMaster, tMasterIndices( f ), false );

                            // save facet in map
                            mFacetMap[ tFacet->id() ] = tFacet;

                            BELFEM_ASSERT( tFacet->master_index() < tMaster->number_of_facets(),
                                          "error in assigning facet id %lu to element %lu",
                                          ( long unsigned int ) tFacet->id(),
                                          ( long unsigned int ) tMaster->id());

                            // test if a slave exists
                            if ( tSlaveIDs( f ) > 0 )
                            {
                                // link facet to slave
                                tFacet->set_slave( mElementMap( tSlaveIDs( f )),
                                                   tSlaveIndices( f ));
                            }

                            // add facet to sideset
                            tSideSet->insert_facet( tFacet );
                        }
                    }

                    // add sideset to list
                    tSideSets( tCount++ ) = tSideSet;
                }
            }

            // add sidesets to mesh
            Cell< mesh::SideSet * > & tSets = aConnectorSwitch ?
                    mSubMesh->cuts() : mSubMesh->sidesets();

            tSets.set_size( tCount, nullptr );

            for( uint s=0; s<tCount; ++s )
            {
                tSets( s ) = tSideSets( s );
            }



            if( ! aConnectorSwitch )
            {
                Vector< id_t > tEdgeIDs ;
                Vector< id_t > tFaceIDs ;
                Vector< id_t > tElementsWithEdges ;
                Vector< id_t > tElementsWithFaces ;

                receive( mMasterRank, tElementsWithEdges );
                receive( mMasterRank, tEdgeIDs );
                receive( mMasterRank, tElementsWithFaces );
                receive( mMasterRank, tFaceIDs );

                // initialize counter
                index_t tCount = 0 ;

                // link facets with edges
                for( id_t tID : tElementsWithEdges )
                {
                    // grab element
                    mesh::Element * tElement = mFacetMap[ tID ]->element();

                    tElement->allocate_edge_container() ;

                    // link edges
                    for( uint e=0; e<tElement->number_of_edges(); ++e )
                    {
                        tElement->insert_edge( mEdgeMap[ tEdgeIDs( tCount++ ) ], e );
                    }
                }

                // reset counter
                tCount = 0 ;

                // link facets with faces
                for( id_t tID : tElementsWithFaces )
                {
                    // grab element
                    mesh::Element * tElement = mFacetMap[ tID ]->element();

                    tElement->allocate_face_container();

                    // link faces
                    for ( uint f = 0; f < tElement->number_of_edges(); ++f )
                    {
                        tElement->insert_face( mFaceMap[ tFaceIDs( tCount++ ) ], f );
                    }
                }
            } // end if face mode
        }

//------------------------------------------------------------------------------

        void
        Kernel::create_fields()
        {
            if ( mMyRank == mMasterRank )
            {
                Timer tTimer;

                // send field dimension to each proc
                Vector< uint > tNumDOFsPerNode = mParams->num_dofs_per_node();
                Vector< uint > tNumDOFsPerEdge = mParams->num_dofs_per_edge();
                Vector< uint > tEnforceLinear  = mParams->linear_enforcement_flags();

                if ( tNumDOFsPerNode.length() == 1 )
                {
                    message( 4, " Creating 1 field ..." );
                }
                else
                {
                    message( 4, " Creating %u fields ...",
                             ( unsigned int ) tNumDOFsPerNode.length());
                }

                Cell< Vector< uint > > tNodeFieldDimensions( mNumberOfProcs, tNumDOFsPerNode );
                send( mCommTable, tNodeFieldDimensions );

                Cell< Vector< uint > > tEdgeFieldDimensions( mNumberOfProcs, tNumDOFsPerEdge );
                send( mCommTable, tEdgeFieldDimensions );

                Cell< Vector< uint > > tEnforceLinearFlags( mNumberOfProcs, tEnforceLinear );
                send( mCommTable, tEnforceLinearFlags );

                const Vector< uint > & tBlockIndices = mParams->block_indices();

                Cell< mesh::Block * > & tBlocks = mMesh->blocks();

                // prepare block IDs
                Vector< id_t > tBlockIDs( tBlockIndices.length() );
                for ( uint b = 0; b < tBlockIndices.length(); ++b )
                {
                    tBlockIDs( b ) = tBlocks( tBlockIndices( b ) )->id();
                }

                // send block IDs to other procs
                // note: we assume that all blocks contribute to all fields here
                Cell< Vector< id_t > > tAllBlocksIDs( mNumberOfProcs, tBlockIDs );
                send( mCommTable, tAllBlocksIDs );

                // get number of fields
                uint tNumberOfFields = tNumDOFsPerNode.length();

                // populate the fields
                mFields.set_size( tNumberOfFields, nullptr );

                for ( uint f = 0; f < tNumberOfFields; ++f )
                {
                    bool tEnforceLinearFlag = false ;
                    if( f < tEnforceLinear.length() )
                    {
                        if ( tEnforceLinear( f ) != 0 )
                        {
                            tEnforceLinearFlag = true ;
                        }
                    }

                    mFields( f ) = new Field(
                            this, tBlockIDs, tNumDOFsPerNode( f ), tNumDOFsPerEdge( f ),
                            tEnforceLinearFlag );
                }

                message( 4, "    ... time for field creation                 : %u ms\n",
                         ( unsigned int ) tTimer.stop());
            }
            else
            {
                Vector< uint > tNumDOFsPerNode;
                receive( mMasterRank, tNumDOFsPerNode );

                Vector< uint > tNumDOFsPerEdge;
                receive( mMasterRank, tNumDOFsPerEdge );

                Vector< uint > tEnforceLinear;
                receive( mMasterRank, tEnforceLinear );

                Vector< id_t > tRecBlockIDs;
                receive( mMasterRank, tRecBlockIDs );

                // count existing blocks
                uint tCount = 0;

                for ( id_t tID : tRecBlockIDs )
                {
                    if ( mSubMesh->block_exists( tID ))
                    {
                        ++tCount;
                    }
                }

                // generate container with blocks that exist on this proc
                Vector< id_t > tBlockIDs( tCount );

                tCount = 0;
                for ( id_t tID : tRecBlockIDs )
                {
                    if ( mSubMesh->block_exists( tID ))
                    {
                        tBlockIDs( tCount++ ) = tID;
                    }
                }

                // get number of fields
                uint tNumberOfFields = tNumDOFsPerNode.length();

                // populate the fields
                mFields.set_size( tNumberOfFields, nullptr );

                for ( uint f = 0; f < tNumberOfFields; ++f )
                {
                    bool tEnforceLinearFlag = false ;

                    if( f < tEnforceLinear.length() )
                    {
                        if ( tEnforceLinear( f ) != 0 )
                        {
                            tEnforceLinearFlag = true ;
                        }
                    }

                    mFields( f ) = new Field(
                            this,
                            tBlockIDs,
                            tNumDOFsPerNode( f ),
                            tNumDOFsPerEdge( f ),
                            tEnforceLinearFlag );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::merge_facet_and_connector_tables()
        {
            for( proc_t p=0; p<mNumberOfProcs; ++p )
            {
                // get the facet table
                Vector< index_t > & tFacets = mFacetTable( p );

                // get the connector table
                Vector< index_t > & tConnectors = mConnectorTable( p );

                // create backup copy
                Vector< index_t > tTemp ( tFacets );

                // reallocate facets
                tFacets.set_size( tTemp.length() + tConnectors.length() );

                // merge data
                index_t tCount = 0 ;

                for( index_t tIndex : tTemp )
                {
                    tFacets( tCount++ ) = tIndex ;
                }
                for( index_t tIndex : tConnectors )
                {
                    tFacets( tCount++ ) = tIndex ;
                }
            }

            // delete connector table
            mConnectorTable.clear() ;
        }

//------------------------------------------------------------------------------
        void
        Kernel::delete_maps()
        {
            mNodeMap.clear() ;
            mEdgeMap.clear() ;
            mFacetMap.clear() ;
            mVertexMap.clear() ;
            mElementMap.clear();
        }

//------------------------------------------------------------------------------

        void
        Kernel::symrcm()
        {
            if( mMyRank == mMesh->master() )
            {
                Cell< graph::Vertex * > & tVertices
                    = reinterpret_cast< Cell< graph::Vertex * > & >
                            ( mMesh->nodes() );

                graph::symrcm( tVertices );
            }
        }

//------------------------------------------------------------------------------

        Field *
        Kernel::field( const uint aIndex )
        {
            return mFields( aIndex );
        }

//------------------------------------------------------------------------------

        DofManager *
        Kernel::dofmgr( const uint aIndex )
        {
            return mDofManagers( aIndex );
        }

//------------------------------------------------------------------------------

        Material *
        Kernel::get_material( const MaterialType aMaterial )
        {
            // test if material exists
            if( mMaterialMap.key_exists( aMaterial ) )
            {
                return mMaterialMap( aMaterial );
            }
            else
            {
                // create a factory
                MaterialFactory tFactory;

                // create the object
                Material * aMat = tFactory.create_material( aMaterial );

                // add material to array
                mMaterials.push( aMat );

                // add material to map
                mMaterialMap[ aMaterial ] = aMat;

                // return the new object
                return aMat;
            }
        }

//------------------------------------------------------------------------------

        void
        Kernel::add_material( Material * aMaterial )
        {
            // add material to array
            mMaterials.push( aMaterial );
        }

//------------------------------------------------------------------------------

        const Vector< index_t > &
        Kernel::node_table( const uint aProcIndex )
        {
            BELFEM_ASSERT( mCommTable( aProcIndex ) != mMasterRank,
                    "Node table is not specified for master." );

            BELFEM_ASSERT( aProcIndex < mCommTable.length(),
                          "Illegal proc index." );

            BELFEM_ASSERT( mMyRank == mMasterRank,
                    "Only master may call the node table." );

            return mNodeTable( aProcIndex );
        }

//------------------------------------------------------------------------------

        const Vector< index_t > &
        Kernel::element_table( const uint aProcIndex )
        {
            BELFEM_ASSERT( mCommTable( aProcIndex ) != mMasterRank,
                           "Element table is not specified for master." );

            BELFEM_ASSERT( aProcIndex < mCommTable.length(),
                           "Illegal proc index." );

            BELFEM_ASSERT( mMyRank == mMasterRank,
                           "Only master may call the element table." );

            return mElementTable( aProcIndex );
        }

//------------------------------------------------------------------------------

        const Vector< index_t > &
        Kernel::edge_table( const uint aProcIndex )
        {
            BELFEM_ASSERT( mCommTable( aProcIndex ) != mMasterRank,
                          "Edge table is not specified for master." );

            BELFEM_ASSERT( aProcIndex < mCommTable.length(),
                          "Illegal proc index." );

            BELFEM_ASSERT( mMyRank == mMasterRank,
                          "Only master may call the edge table." );

            return mEdgeTable( aProcIndex );
        }

//------------------------------------------------------------------------------

        const Vector< index_t > &
        Kernel::face_table( const uint aProcIndex )
        {
            BELFEM_ASSERT( mCommTable( aProcIndex ) != mMasterRank,
                          "Face table is not specified for master." );

            BELFEM_ASSERT( aProcIndex < mCommTable.length(),
                          "Illegal proc index." );

            BELFEM_ASSERT( mMyRank == mMasterRank,
                          "Only master may call the edge table." );

            return mFaceTable( aProcIndex );
        }

//------------------------------------------------------------------------------

        const Vector< index_t > &
        Kernel::facet_table( const uint aProcIndex )
        {
            BELFEM_ASSERT( mCommTable( aProcIndex ) != mMasterRank,
                          "Edge table is not specified for master." );

            BELFEM_ASSERT( aProcIndex < mCommTable.length(),
                          "Illegal proc index." );

            BELFEM_ASSERT( mMyRank == mMasterRank,
                          "Only master may call the edge table." );

            return mFacetTable( aProcIndex );
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
    }
}
