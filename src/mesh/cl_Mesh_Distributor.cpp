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

#include "cl_Mesh_Distributor.hpp"

#include "assert.hpp"
#include "commtools.hpp"
#include "fn_entity_type.hpp"

namespace belfem
{
    namespace mesh
    {

//------------------------------------------------------------------------------

        Distributor::Distributor( Mesh * aMesh ) :
            mCommRank( comm_rank() ),
            mCommSize( comm_size() ),
            mMesh( mCommRank == 0 ? aMesh : nullptr )
        {
            BELFEM_ERROR( ( mCommRank == 0 && aMesh != nullptr ) || mCommRank != 0,
                "The root proc must not pass a null pointer to this constructor");

            mOwnMesh = mCommRank > 0 ;
        }

//------------------------------------------------------------------------------

        Distributor::~Distributor()
        {
            this->delete_bitsets();

            this->delete_node_data();
            this->delete_edge_data();
            this->delete_face_data();
            this->delete_element_data();
            this->delete_facet_data();

            this->delete_element_extra() ;
            this->delete_facet_extra() ;

            this->delete_vertex_data();
            this->delete_t_matrices();
            this->delete_tables();

            this->delete_mesh_data();

            if ( mOwnMesh && mMesh != nullptr )
            {
                delete mMesh ;
                mMesh = nullptr ;
            }

            // if the tables are populated, they will be destroyed
            // usually, the kernel would move the Cell and claim ownership
            for ( CommTable * tTable : mTables )
            {
                delete tTable ;
            }
        }

//-------------------------------------------------------------------------------

        void
        Distributor::run()
        {
            if ( mCommSize < 2 )return ;

            this->count_entities();

            if ( mCommRank == 0 )
            {
                //mMesh->set_node_owners();

                this->create_bitsets();

                // create the containers
                mNodeData.set_size( mCommSize, nullptr );
                mEdgeData.set_size( mCommSize, nullptr );
                mFaceData.set_size( mCommSize, nullptr );
                mElementData.set_size( mCommSize, nullptr );

                mFacetData.set_size( mCommSize, nullptr );
                mVertexData.set_size( mCommSize, nullptr );
                mControlPointData.set_size( mCommSize, nullptr );

                mElementExtra.set_size( mCommSize, nullptr );
                mFacetExtra.set_size( mCommSize, nullptr );
                mTables.set_size( mCommSize, nullptr );

                mTMatrices.set_size( mCommSize, nullptr );

                // create the container data
                for ( proc_t p=1; p<mCommSize; ++p )
                {
                    mNodeData( p )     = new proto::NodeData ;
                    mEdgeData( p )     = new proto::EdgeData ;
                    mFaceData( p )     = new proto::FaceData ;
                    mElementData( p )  = new proto::ElementData ;
                    mFacetData( p )    = new proto::FacetData ;
                    mVertexData( p )   = new proto::ElementData ;
                    mControlPointData( p ) = new proto::ControlPointData ;
                    mElementExtra( p ) = new proto::ElementExtra ;
                    mFacetExtra( p )   = new proto::FacetExtra ;
                    mTMatrices( p )    = new proto::TMatrixData ;
                    mTables( p )       = new CommTable ;
                }

                for ( proc_t p=1; p<mCommSize; ++p )
                {
                    this->select_entities( p );
                    this->populate_t_matrices( p );
                    this->populate_node_data( p );
                    this->populate_edge_data( p );
                    this->populate_face_data( p );
                    this->populate_element_data( p );
                    this->populate_facet_data( p );
                    this->populate_facet_extra( p );
                    this->populate_vertex_data( p );
                    this->populate_control_point_data( p );
                    this->flag_curved_elements( p );
                    this->flag_curved_facets( p );
                }
            }

            // wait for other proces
            comm_barrier();


            // create meshes on other procs
            if ( mCommRank != 0 )
            {
                mMesh = new Mesh( mNumberOfDimensions, mCommRank );
                mOwnMesh = true ;
            }

            if ( mCommRank == 0 )
            {
                this->send_block_data();
                this->send_sideset_data();
                this->send_thinshell_data();
                this->send_node_data() ;
                comm_barrier();

                this->send_element_data();
                comm_barrier();

                this->send_edge_data();
                comm_barrier();

                this->send_face_data();
                comm_barrier();

                this->send_facet_data();
                comm_barrier();

                this->send_vertex_data();
                comm_barrier();

                this->send_control_point_data();
                comm_barrier();

                this->send_facet_extra();
                comm_barrier();

                this->send_t_matrices();
                comm_barrier();
            }
            else
            {
                mProtoMesh = new ProtoMesh( mMesh );

                this->receive_block_data();
                this->receive_sideset_data();
                this->receive_thinshell_data();

                this->receive_node_data();
                mProtoMesh->create_nodes() ;
                comm_barrier();

                this->receive_element_data();
                mProtoMesh->create_elements() ;
                comm_barrier();

                this->receive_edge_data();
                mProtoMesh->create_edges() ;
                comm_barrier();

                this->receive_face_data();
                mProtoMesh->create_faces() ;
                comm_barrier();

                this->receive_facet_data();
                mProtoMesh->create_facets() ;
                comm_barrier();

                this->receive_vertex_data();
                mProtoMesh->create_vertices();
                comm_barrier();

                this->receive_control_point_data();
                mProtoMesh->create_control_points();
                comm_barrier();

                this->receive_facet_extra();
                mProtoMesh->create_facet_extra();
                comm_barrier();

                this->receive_t_matrices();
                comm_barrier();

            }

            if ( mCommRank == 0 )
            {
                for ( proc_t p=1; p<mCommSize; ++p )
                {
                    this->populate_element_extra( p );
                    comm_barrier();
                    this->send_element_extra( p );
                }
            }
            else
            {
                for ( proc_t p=1; p<mCommSize; ++p )
                {
                    comm_barrier();
                    if ( p == mCommRank ) this->receive_element_extra();
                }
                mProtoMesh->create_element_extra() ;
            }

            if ( mCommRank > 0 )
            {

                mProtoMesh->edge();
                mProtoMesh->create_sidesets();
                mProtoMesh->create_thinshells();
                mMesh->finalize();
                mProtoMesh->create_t_matrices();
            }

            comm_barrier();

        }

//-------------------------------------------------------------------------------


        void
        Distributor::count_entities()
        {
            if ( mCommRank == 0 )
            {
                mNumberOfDimensions       = mMesh->number_of_dimensions();
                mMaxElementOrder          = mMesh->max_element_order();
                mNumberOfAllNodes         = mMesh->nodes().size() ;
                mNumberOfAbstractNodes    = mMesh->abstract_nodes().size() ;
                mNumberOfAllEdges         = mMesh->edges().size() ;
                mNumberOfAllFaces         = mMesh->faces().size() ;
                mNumberOfAllElements      = mMesh->elements().size() ;
                mNumberOfAllFacets        = mMesh->facets().size() ;
                mNumberOfAllVertices      = mMesh->vertices().size() ;
                mNumberOfAllControlPoints = mMesh->control_points().size() ;
                mNumberOfAllBlocks        = mMesh->blocks().size() ;
                mNumberOfAllSideSets      = mMesh->sidesets().size() ;
                mNumberOfAllThinShells    = mMesh->thin_shells().size() ;
                mNumberOfAllTMatrices     = mMesh->number_of_hanging_basis();
            }

            broadcast( mNumberOfAllEntities, 0, 14 );
        }

        void
        Distributor::select_entities( const proc_t aTarget )
        {
            BELFEM_ASSERT( mCommRank == 0,
                       "select_entities() must only be called by master proc only" );

            mMesh->unflag_all_elements() ;
            mMesh->unflag_all_facets() ;

            this->reset_bitsets();

            // get entity containers
            Cell< Node     * >     & tNodes         = mMesh->nodes();
            //Cell< Edge     * >     & tEdges         = mMesh->edges();
            //Cell< Face     * >     & tFaces         = mMesh->faces();
            Cell< Element  * >     & tElements      = mMesh->elements();
            Cell< Facet    * >     & tFacets        = mMesh->facets();
            Cell< Element  * >     & tVertices      = mMesh->vertices();

            BELFEM_ERROR( tNodes.size() > 0, "No nodes in the mesh" );
            BELFEM_ERROR( tElements.size() > 0, "No elements in the mesh" );

            // select owned elements and facets
            for ( Element * tElement : tElements )
            {
                if ( tElement->owner() == aTarget )
                {
                    mElementBitset->set( tElement->index() );
                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        Node * tNode = tElement->node( k ) ;
                        this->select_sources( tNode );
                        mNodeBitset->set( tNode->index() );
                    }
                }
            }

            for ( Facet * tFacet : tFacets )
            {
                if ( tFacet->owner() == aTarget )
                {
                    mFacetBitset->set( tFacet->index() );

                    // a periodic slave facet hangs on its master facet ( a FACET-typed
                    // source set in Periodicity::set_entity_dependencies ), so the facet's
                    // own sources must be ghosted too - not just its node sources - or the
                    // master facet is missing on this rank when create_t_matrices runs
                    this->select_sources( tFacet );

                    for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                    {
                        Node * tNode = tFacet->node( k ) ;
                        this->select_sources( tNode );
                        mNodeBitset->set( tNode->index() );
                    }
                }
            }

            for ( Element * tVertex : tVertices )
            {
                if ( tVertex->owner() == aTarget )
                {
                    mVertexBitset->set( tVertex->index() );
                    for ( uint k=0; k<tVertex->number_of_nodes(); ++k )
                    {
                        this->select_sources( tVertex->node( k ) );
                    }
                }
            }
            // up until now, only nodes connected to owned elements and facets are set

            // definitley set node owned by this proc
            for ( Node * tNode : tNodes )
            {
                if ( tNode->owner() == aTarget )
                {
                    mNodeBitset->set( tNode->index() );
                }
            }

            // now, make sure that duplicates are set
            Cell< index_t > tNodeIndices ;
            mNodeBitset->where( tNodeIndices );
            for ( index_t k: tNodeIndices )
            {
                Node * tOrg = tNodes( k )->original() ;

                // the original itself must be in the table: the duplicate
                // records are keyed on it, and an element-less original
                // ( orphaned seam node ) is selected by no other rule
                mNodeBitset->set( tOrg->index() );

                for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                {
                    mNodeBitset->set( tOrg->duplicate( d )->index() );
                }
            }

            // for the post processors, we also need to flag connected nodes to this set
            mNodeBitset->where( tNodeIndices );

            for ( index_t k: tNodeIndices )
            {
                Node * tNode = tNodes( k ) ;

                for ( uint i=0; i<tNode->number_of_nodes(); ++i )
                {
                    mNodeBitset->set( tNode->node( i )->index() );
                }
            }

            // we also select the abstract nodes
            for ( Node * tNode : mMesh->abstract_nodes() )
            {
                mNodeBitset->set( tNode->index() );
            }

            // ( the former pre-aura edge/face select_sources loops were dead code:
            //   they gated on tElement->is_flagged() before any element was flagged,
            //   so hanging edge/face/facet sources were never ghosted. Source-closure
            //   now happens in the joint fixpoint below. )

            // After having the sources of the owned elements and facets, we create the aura
            mNodeBitset->where( tNodeIndices );

            for ( index_t k: tNodeIndices )
            {
                Node * tOrg = tNodes( k )->original() ;
                for ( uint e=0; e<tOrg->number_of_elements(); ++e )
                {
                    mElementBitset->set( tOrg->element( e )->index() );
                }
                for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                {
                    Node * tDup = tOrg->duplicate( d ) ;
                    for ( uint e=0; e<tDup->number_of_elements(); ++e )
                    {
                        mElementBitset->set( tDup->element( e )->index() );
                    }
                }
            }


            Cell< index_t > tFacetIndices ;
            if ( mMesh->number_of_facets() > 0 )
            {
                mFacetBitset->where( tFacetIndices );
                for ( index_t f: tFacetIndices )
                {
                    Facet * tFacet = tFacets( f ) ;
                    for ( uint n=0; n<tFacet->number_of_facets(); ++n )
                    {
                        mFacetBitset->set( tFacet->facet( n )->index() );
                    }
                }

                // pick master and slave elements of all facets
                mFacetBitset->where( tFacetIndices );
                for ( index_t f: tFacetIndices )
                {
                    Facet * tFacet = tFacets( f ) ;

                    if ( tFacet->has_master() )
                    {
                        mElementBitset->set( tFacet->master()->index() );
                    }
                    if ( tFacet->has_slave() )
                    {
                        mElementBitset->set( tFacet->slave()->index() );
                    }
                    for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                    {
                        mNodeBitset->set( tFacet->node( k )->index() );
                    }
                }
            }

            // flag the nodes of these elements and facets
            Cell< index_t > tElementIndices ;
            mElementBitset->where( tElementIndices );
            for ( index_t e: tElementIndices )
            {
                Element * tElement = tElements( e ) ;
                tElement->flag();
                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    mNodeBitset->set( tElement->node( k )->index() );
                }
            }

            // with the node originals selected, we also make sure that the duplicates are set
            // and vice versa

            if ( mMesh->edges_exist() )
            {
                for ( Block * tBlock : mMesh->blocks() )
                {
                    if ( tBlock->has_edges() )
                    {
                        for ( Element * tElement : tBlock->elements() )
                        {
                            if ( tElement->is_flagged() )
                            {
                                for ( uint e=0; e<tElement->number_of_edges(); ++e )
                                {
                                    Edge * tEdge = tElement->edge( e ) ;
                                    mEdgeBitset->set( tEdge->index() );
                                    for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                                    {
                                        mNodeBitset->set( tEdge->node( k )->index() );
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if ( mMesh->faces_exist() )
            {
                for ( Block * tBlock : mMesh->blocks() )
                {
                    if ( tBlock->has_faces() )
                    {
                        for ( Element * tElement : tBlock->elements() )
                        {
                            if ( tElement->is_flagged() )
                            {
                                for ( uint f=0; f<tElement->number_of_faces(); ++f )
                                {
                                    Face * tFace = tElement->face( f ) ;
                                    mFaceBitset->set( tFace->index() );
                                    for ( uint k=0; k<tFace->number_of_nodes(); ++k )
                                    {
                                        mNodeBitset->set( tFace->node( k )->index() );
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // 5b: source-closure fixpoint.
            // Ghost the SOURCES of every selected entity ( node, edge, face, facet ),
            // transitively, until no bitset grows. A source only needs to exist as a
            // basis ( plus its own sub-nodes, which select_sources adds ), so this does
            // NOT expand the geometric aura -- the node->element->edge/face halo was
            // already grown one layer above and must stay one layer, or the closure
            // would walk the whole mesh. The old loop closed node sources only; this
            // also closes the hanging edge/face/facet sources that were missed.
            index_t tPrevTotal = 0 ;
            while ( true )
            {
                mNodeBitset->where( tNodeIndices );
                for ( index_t k : tNodeIndices )
                {
                    Node * tNode = tNodes( k )->original();
                    for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                    {
                        Node * tDup = tNode->duplicate( d );
                        mNodeBitset->set( tDup->index() );
                        this->select_sources( tDup );
                    }
                    this->select_sources( tNode );
                }
                if ( mEdgeBitset != nullptr )
                {
                    mEdgeBitset->where( tNodeIndices );
                    for ( index_t k : tNodeIndices ) this->select_sources( mMesh->edges()( k ) );
                }
                if ( mFaceBitset != nullptr )
                {
                    mFaceBitset->where( tNodeIndices );
                    for ( index_t k : tNodeIndices ) this->select_sources( mMesh->faces()( k ) );
                }
                if ( mFacetBitset != nullptr )
                {
                    mFacetBitset->where( tNodeIndices );
                    for ( index_t k : tNodeIndices ) this->select_sources( tFacets( k ) );
                }

                // converge when no bitset grew ( bits are only ever set -> terminates ).
                // control points are summed too so the loop cannot stop early if a
                // select_sources call adds a control-point source ( closing the sources
                // OF control points / elements is deferred -- no current mesh has them ).
                index_t tTotal = mNodeBitset->count();
                if ( mEdgeBitset != nullptr ) tTotal += mEdgeBitset->count();
                if ( mFaceBitset != nullptr ) tTotal += mFaceBitset->count();
                if ( mFacetBitset != nullptr ) tTotal += mFacetBitset->count();
                if ( mControlPointBitset != nullptr ) tTotal += mControlPointBitset->count();
                if ( tTotal == tPrevTotal ) break;
                tPrevTotal = tTotal;
            }

#if !defined( NDEBUG ) || defined( DEBUG )
            // Debug-only source-closure invariant ( plan 5b ): every selected entity that
            // carries sources must have each of its sources selected in the matching
            // bitset. Guards against a regression of the source-closure gap that broke
            // parallel periodic/thin-shell distribution. Compiled out in release.
            {
                auto tSrcInSet = [&]( Basis * aSrc ) -> bool
                {
                    switch ( aSrc->entity_type() )
                    {
                        case EntityType::NODE :
                            return mNodeBitset->test( reinterpret_cast< Node * >( aSrc )->original()->index() );
                        case EntityType::EDGE :
                            return mEdgeBitset != nullptr && mEdgeBitset->test( reinterpret_cast< Edge * >( aSrc )->index() );
                        case EntityType::FACE :
                            return mFaceBitset != nullptr && mFaceBitset->test( reinterpret_cast< Face * >( aSrc )->index() );
                        case EntityType::FACET :
                            return mFacetBitset != nullptr && mFacetBitset->test( reinterpret_cast< Facet * >( aSrc )->index() );
                        case EntityType::CELL :
                        case EntityType::ELEMENT :
                            return mElementBitset != nullptr && mElementBitset->test( reinterpret_cast< Element * >( aSrc )->index() );
                        case EntityType::CONTROLPOINT :
                            return mControlPointBitset != nullptr && mControlPointBitset->test( reinterpret_cast< ControlPoint * >( aSrc )->index() );
                        default :
                            return true ;
                    }
                };

                auto tCheck = [&]( Basis * aEntity, const char * aKind )
                {
                    for ( uint s=0; s<aEntity->number_of_sources(); ++s )
                    {
                        Basis * tSrc = aEntity->source( s ) ;
                        BELFEM_ERROR( tSrcInSet( tSrc ),
                            "select_entities( target %lu ): %s %lu has source %u ( %s %lu ) not selected on this rank - source-closure regression",
                            ( long unsigned int ) aTarget, aKind, ( long unsigned int ) aEntity->id(),
                            ( unsigned int ) s, to_string( tSrc->entity_type() ).c_str(),
                            ( long unsigned int ) tSrc->id() );
                    }
                };

                Cell< index_t > tWhere ;
                mNodeBitset->where( tWhere ) ;
                for ( index_t i : tWhere ) tCheck( mMesh->nodes()( i ), "node" ) ;
                if ( mEdgeBitset != nullptr )
                {
                    mEdgeBitset->where( tWhere ) ;
                    for ( index_t i : tWhere ) tCheck( mMesh->edges()( i ), "edge" ) ;
                }
                if ( mFaceBitset != nullptr )
                {
                    mFaceBitset->where( tWhere ) ;
                    for ( index_t i : tWhere ) tCheck( mMesh->faces()( i ), "face" ) ;
                }
                if ( mFacetBitset != nullptr )
                {
                    mFacetBitset->where( tWhere ) ;
                    for ( index_t i : tWhere ) tCheck( mMesh->facets()( i ), "facet" ) ;
                }
                if ( mElementBitset != nullptr )
                {
                    mElementBitset->where( tWhere ) ;
                    for ( index_t i : tWhere ) tCheck( mMesh->elements()( i ), "element" ) ;
                }
                if ( mControlPointBitset != nullptr )
                {
                    mControlPointBitset->where( tWhere ) ;
                    for ( index_t i : tWhere ) tCheck( mMesh->control_points()( i ), "controlpoint" ) ;
                }
            }
#endif
        }


        void
        Distributor::populate_t_matrices( const proc_t aTarget )
        {
            if ( mNumberOfAllTMatrices == 0 ) return ;

            Cell< index_t > tIndices ;

            index_t tCount = 0 ;

            // check for nodes
            index_t tNumNodes = 0 ;
            Cell< mesh::Node * > & tNodes = mMesh->nodes();
            mNodeBitset->where( tIndices );
            for ( index_t i: tIndices )
            {
                if ( tNodes( i )->number_of_sources() > 0 )
                {
                    tCount += tNodes( i )->number_of_sources() ;
                    ++tNumNodes ;
                }
            }

            // check for edges
            index_t tNumEdges = 0 ;
            if ( mEdgeBitset != nullptr )
            {
                Cell< mesh::Edge * > & tEdges = mMesh->edges();
                mEdgeBitset->where( tIndices );
                for ( index_t i: tIndices )
                {
                    if ( tEdges( i )->number_of_sources() > 0 )
                    {
                        tCount += tEdges( i )->number_of_sources() ;
                        ++tNumEdges ;
                    }
                }
            }

            // check for faces
            index_t tNumFaces = 0 ;
            if ( mFaceBitset != nullptr )
            {
                Cell< mesh::Face * > & tFaces = mMesh->faces();
                mFaceBitset->where( tIndices );
                for ( index_t i: tIndices )
                {
                    if ( tFaces( i )->number_of_sources() > 0 )
                    {
                        tCount += tFaces( i )->number_of_sources() ;
                        ++tNumFaces ;
                    }
                }
            }

            // elements
            index_t tNumElements = 0 ;
            Cell< mesh::Element * > & tElements = mMesh->elements();
            mElementBitset->where( tIndices );
            for ( index_t i: tIndices )
            {
                if ( tElements( i )->number_of_sources() > 0 )
                {
                    tCount += tElements( i )->number_of_sources() ;
                    ++tNumElements ;
                }
            }

            // facets
            index_t tNumFacets = 0 ;
            if ( mMesh->number_of_facets() > 0 )
            {
                Cell< mesh::Facet * > & tFacets = mMesh->facets();
                mFacetBitset->where( tIndices );
                for ( index_t i: tIndices )
                {
                    if ( tFacets( i )->number_of_sources() > 0 )
                    {
                        tCount += tFacets( i )->number_of_sources() ;
                        ++tNumFacets ;
                    }
                }
            }

            // control points
            index_t tNumControlPoints = 0 ;
            if ( mMesh->number_of_control_points() > 0 )
            {
                Cell< mesh::ControlPoint * > & tControlPoints = mMesh->control_points();
                mControlPointBitset->where( tIndices );
                for ( index_t i: tIndices )
                {
                    if ( tControlPoints( i )->number_of_sources() > 0 )
                    {
                        tCount += tControlPoints( i )->number_of_sources() ;
                        ++tNumControlPoints ;
                    }
                }
            }

            // get the data containers
            Cell< index_t > & tNumTargets  = mTMatrices( aTarget )->mNumTargets ;
            Cell< uint >    & tCounters    = mTMatrices( aTarget )->mCounters ;
            Cell< id_t >    & tSourceIDs   = mTMatrices( aTarget )->mSourceIDs ;
            Cell< id_t >    & tTargetIDs   = mTMatrices( aTarget )->mTargetIDs ;
            Cell< uchar >   & tTypes       = mTMatrices( aTarget )->mTypes ;
            Cell< real >    & tWeights     = mTMatrices( aTarget )->mWeights ;

            tNumTargets.set_size( static_cast< uint >( EntityType::UNDEFINED ), 0 );
            tNumTargets( static_cast< uint >(EntityType::NODE )) = tNumNodes ;
            tNumTargets( static_cast< uint >(EntityType::EDGE )) = tNumEdges ;
            tNumTargets( static_cast< uint >(EntityType::FACE )) = tNumFaces ;
            tNumTargets( static_cast< uint >(EntityType::ELEMENT )) = tNumElements ;
            tNumTargets( static_cast< uint >(EntityType::FACET )) = tNumFacets ;
            tNumTargets( static_cast< uint >(EntityType::CONTROLPOINT )) = tNumControlPoints ;

            tTargetIDs.set_size( mNumberOfAllTMatrices, gNoID );
            tCounters.set_size( tTargetIDs.size(), 0 );

            tSourceIDs.set_size( tCount, gNoID );
            tTypes.set_size( tCount, 0 );
            tWeights.set_size( tCount, 0.0 );

            mNodeBitset->where( tIndices );
            index_t tPivot = 0 ;
            tCount = 0 ;
            for ( index_t i: tIndices )
            {
                Node * tNode = tNodes( i ) ;
                uint s = tNode->number_of_sources() ;
                if ( s > 0 )
                {
                    tTargetIDs( tPivot ) = tNode->id() ;
                    tCounters( tPivot++ ) = s ;

                    for ( uint k=0; k<s; ++k )
                    {
                        Basis * tSource = tNode->source( k ) ;
                        tSourceIDs( tCount ) = tSource->id() ;
                        tTypes( tCount )     = static_cast< uchar >( tSource->entity_type() );
                        tWeights( tCount++ ) = tNode->weight( k ) ;
                    }
                }
            }
            if ( mEdgeBitset != nullptr )
            {
                Cell< Edge * > & tEdges = mMesh->edges();
                mEdgeBitset->where( tIndices );
                for ( index_t i: tIndices )
                {
                    Edge * tEdge = tEdges( i ) ;
                    uint s = tEdge->number_of_sources() ;
                    if ( s > 0 )
                    {
                        tTargetIDs( tPivot ) = tEdge->id() ;
                        tCounters( tPivot++ ) = s ;
                        for ( uint k=0; k<s; ++k )
                        {
                            Basis * tSource      = tEdge->source( k ) ;
                            tSourceIDs( tCount ) = tSource->id() ;
                            tTypes( tCount )     = static_cast< uchar >( tSource->entity_type() );
                            tWeights( tCount++ ) = tEdge->weight( k ) ;
                        }
                    }
                }
            }
            if ( mFaceBitset != nullptr )
            {
                Cell< Face * > & tFaces = mMesh->faces();
                mFaceBitset->where( tIndices );
                for ( index_t i: tIndices )
                {
                    Face * tFace = tFaces( i ) ;
                    uint s = tFace->number_of_sources() ;
                    if ( s > 0 )
                    {
                        tTargetIDs( tPivot )  = tFace->id() ;
                        tCounters( tPivot++ ) = s ;
                        for ( uint k=0; k<s; ++k )
                        {
                            Basis * tSource     = tFace->source( k ) ;
                            tSourceIDs( tCount ) = tSource->id() ;
                            tTypes( tCount )     = static_cast< uchar >( tSource->entity_type() );
                            tWeights( tCount++ ) = tFace->weight( k ) ;
                        }
                    }
                }
            }

            mElementBitset->where( tIndices );
            for ( index_t i: tIndices )
            {
                Element * tElement = tElements( i ) ;
                uint s = tElement->number_of_sources() ;
                if ( s > 0 )
                {
                    tTargetIDs( tPivot ) = tElement->id() ;
                    tCounters( tPivot++ ) = s ;
                    for ( uint k=0; k<s; ++k )
                    {
                        Basis * tSource      = tElement->source( k ) ;
                        tSourceIDs( tCount ) = tSource->id() ;
                        tTypes( tCount )     = static_cast< uchar >( tSource->entity_type() );
                        tWeights( tCount++ ) = tElement->weight( k ) ;
                    }
                }
            }

            if ( mFacetBitset != nullptr )
            {
                mFacetBitset->where( tIndices );
                Cell< Facet * > & tFacets = mMesh->facets();
                for ( index_t i : tIndices )
                {
                    Facet * tFacet = tFacets( i ) ;
                    uint s = tFacet->number_of_sources() ;
                    if ( s > 0 )
                    {
                        tTargetIDs( tPivot ) = tFacet->id() ;
                        tCounters( tPivot++ ) = s ;
                        for ( uint k=0; k<s; ++k )
                        {
                            Basis * tSource      = tFacet->source( k ) ;
                            tSourceIDs( tCount ) = tSource->id() ;
                            tTypes( tCount )     = static_cast< uchar >( tSource->entity_type() );
                            tWeights( tCount++ ) = tFacet->weight( k ) ;
                        }
                    }
                }
            }

            if ( mControlPointBitset != nullptr )
            {
                mControlPointBitset->where( tIndices );
                Cell< ControlPoint * > & tControlPoints = mMesh->control_points();
                for ( index_t i : tIndices )
                {
                    ControlPoint * tControlPoint = tControlPoints( i ) ;
                    uint s = tControlPoint->number_of_sources() ;
                    if ( s > 0 )
                    {
                        tTargetIDs( tPivot ) = tControlPoint->id() ;
                        tCounters( tPivot++ ) = s ;
                        for ( uint k=0; k<s; ++k )
                        {
                            Basis * tSource      = tControlPoint->source( k ) ;
                            tSourceIDs( tCount ) = tSource->id() ;
                            tTypes( tCount )     = static_cast< uchar >( tSource->entity_type() );
                            tWeights( tCount++ ) = tControlPoint->weight( k ) ;
                        }
                    }
                }
            }
        }

        void
        Distributor::select_sources( Basis * aBasis )
        {
            for ( uint s=0; s<aBasis->number_of_sources(); ++s )
            {
                switch ( aBasis->source( s )->entity_type() )
                {
                    case EntityType::NODE:
                    {
                        Node * tNode = reinterpret_cast< Node * >( aBasis->source( s ) )->original();
                        for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                        {
                            mNodeBitset->set( tNode->duplicate( d )->index() );
                        }
                        mNodeBitset->set( tNode->index() );

                        break;
                    }
                    case EntityType::EDGE:
                    {
                        Edge * tEdge = reinterpret_cast< Edge * >( aBasis->source( s ) );
                        mEdgeBitset->set( tEdge->index() );
                        for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                        {
                            Node * tNode = tEdge->node( k )->original();
                            for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                            {
                                mNodeBitset->set( tNode->duplicate( d )->index() );
                            };
                            mNodeBitset->set( tNode->index() );
                        }
                        break ;
                    }
                    case EntityType::FACE :
                    {
                        Face * tFace = reinterpret_cast< Face * >( aBasis->source( s ) );
                        mFaceBitset->set( tFace->index() );
                        for ( uint k=0; k<tFace->number_of_nodes(); ++k )
                        {
                            Node * tNode = tFace->node( k )->original();
                            for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                            {
                                mNodeBitset->set( tNode->duplicate( d )->index() );
                            }
                            mNodeBitset->set( tNode->index() );
                        }
                        break;
                    }
                    // Element::entity_type() returns CELL, so the CELL label is
                    // the one an element source actually arrives with; ELEMENT
                    // is kept for sources that carry the legacy label. Same
                    // pairing as the source-selected test above.
                    case EntityType::CELL :
                    case EntityType::ELEMENT :
                    {
                        Element * tElement = reinterpret_cast< Element * >( aBasis->source( s ) );
                        mElementBitset->set( tElement->index() );
                        for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                        {
                            Node * tNode = tElement->node( k )->original();
                            for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                            {
                                mNodeBitset->set( tNode->duplicate( d )->index() );
                            }
                            mNodeBitset->set( tNode->index() );
                        }
                        break ;
                    }
                    case EntityType::FACET :
                    {
                        Facet * tFacet = reinterpret_cast< Facet * >( aBasis->source( s ) );
                        mFacetBitset->set( tFacet->index() );
                        for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                        {
                            Node * tNode = tFacet->node( k )->original();
                            for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                            {
                                mNodeBitset->set( tNode->duplicate( d )->index() );
                            }
                            mNodeBitset->set( tNode->index() );
                        }
                        break;
                    }
                    case EntityType::CONTROLPOINT :
                    {
                        ControlPoint * tControlPoint = reinterpret_cast< ControlPoint * >( aBasis->source( s ) );
                        mControlPointBitset->set( tControlPoint->index() );
                        break;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid source type for source %u : %s",
                                    ( unsigned int ) s, to_string( aBasis->source( s )->entity_type() ).c_str() );
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::populate_node_data( const proc_t aTarget )
        {
            if ( mNumberOfAllNodes == 0 ) return ;

            // get indices
            Cell< index_t > & tTable    = mTables( aTarget )->nodes() ;

            mNodeBitset->where( tTable );

            // get container
            Cell< Node * > & tNodes = mMesh->nodes();

            // get number of nodes that are to be sent
            index_t tCount = tTable.size();

            // get node data
            Cell< id_t >    & tIDs      = mNodeData( aTarget )->mIDs;

            Cell< proc_t >  & tOwners   = mNodeData( aTarget )->mOwners;
            Matrix< real >  & tX        = mNodeData( aTarget )->mCoords ;



            Cell< id_t >    & tDup      = mNodeData( aTarget )->mDuplicateData ;

            // reserve memory
            tIDs.set_size( tCount, gNoID );
            tOwners.set_size( tCount, gNoOwner );

            tX.set_size( 3, tCount );

            tCount = 0 ;
            index_t tDupCount = 0 ;
            index_t tMemCount = 1;

            // populate main data
            for ( index_t i: tTable )
            {
                Node * tNode = tNodes( i ) ;

                tIDs( tCount ) = tNode->id();
                tOwners( tCount ) = tNode->owner();

                tX( 0, tCount ) = tNode->x();
                tX( 1, tCount ) = tNode->y();
                tX( 2, tCount ) = tNode->z();

                ++tCount;

                if ( tNode->number_of_duplicates() > 0 )
                {
                    ++tDupCount ;
                    tMemCount += tNode->number_of_duplicates() + 2 ;
                }
            }

            // populate duplicate data
            tCount = 0 ;
            tDup.set_size( tMemCount, gNoID );

            tDup( tCount ++ ) = tDupCount ;

            for ( index_t i: tTable )
            {
                Node * tNode = tNodes( i ) ;

                if ( tNode->number_of_duplicates() > 0 )
                {
                    tDup( tCount++ ) = tNode->id();
                    tDup( tCount++ ) = tNode->number_of_duplicates();
                    for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                    {
                        tDup( tCount++ ) = tNode->duplicate( d )->id();
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::populate_element_data( const proc_t aTarget )
        {
            if ( mNumberOfAllElements == 0 ) return ;

            // get indices
            Cell< index_t > & tTable    = mTables( aTarget )->elements() ;
            mElementBitset->where( tTable );

            // get container
            Cell< Element * > & tElements = mMesh->elements();

            // sort the table by block INDEX (position in mMesh->blocks()), so the
            // elements arrive grouped in the same order create_blocks() iterates the
            // block data — independent of how the block ids are numbered.
            Map< id_t, index_t > tBlockIndex ;
            for ( Block * tBlock : mMesh->blocks() )
            {
                tBlockIndex[ tBlock->id() ] = tBlock->index() ;
            }

            Cell< std::pair< index_t, index_t > > tOrder ;
            tOrder.reserve( tTable.size() );

            for ( index_t e : tTable )
            {
                tOrder.push( std::make_pair( e, tBlockIndex( tElements( e )->block_id() ) ) );
            }
            sort( tOrder.begin(), tOrder.end() ,
                  []( const std::pair< index_t, index_t > & a, const std::pair< index_t, index_t > & b )
                  {
                      return a.second < b.second ;
                  } );

            index_t tCount = 0 ;
            for ( auto tPair : tOrder )
            {
                tTable( tCount++ ) = tPair.first ;
            }
            tOrder.clear();

            // get element data
            Cell< id_t >    & tIDs      = mElementData( aTarget )->mIDs ;
            Cell< proc_t >  & tOwners   = mElementData( aTarget )->mOwners ;
            Cell< uchar >   & tTypes    = mElementData( aTarget )->mTypes ;
            Cell< uint >    & tGeo      = mElementData( aTarget )->mGeometryTags ;
            Cell< uint >    & tPhys     = mElementData( aTarget )->mPhysicalTags ;
            Cell< id_t >    & tTopo     = mElementData( aTarget )->mTopology ;

            // allocate memory
            tIDs.set_size( tCount, gNoID );
            tOwners.set_size( tCount, gNoOwner );
            tTypes.set_size( tCount, 0 );
            tGeo.set_size( tCount, 0 );
            tPhys.set_size( tCount, 0 );

            index_t tTopoCount = 0;

            for ( index_t e: tTable )
            {
                tTopoCount += tElements( e )->number_of_nodes();
            }

            tTopo.set_size( tTopoCount, gNoID );

            tCount = 0 ;
            tTopoCount = 0 ;


            for ( index_t e: tTable )
            {
                Element * tElement = tElements( e ) ;
                tIDs( tCount )     = tElement->id();
                tOwners( tCount )  = tElement->owner();
                tTypes( tCount )   = static_cast< uchar >( tElement->type() );
                tGeo( tCount )     = tElement->geometry_tag();
                tPhys( tCount )    = tElement->physical_tag();

                ++tCount;

                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    tTopo( tTopoCount++ ) = tElement->node( k )->id();
                }
            }
        }

//-----------------------------------------------------------------------------

         void
         Distributor::populate_element_extra( const proc_t aTarget )
         {
            if ( mNumberOfAllElements == 0 ) return ;

            Cell< index_t > & tTable = mTables( aTarget )->elements();
            Cell< Element * > & tElements = mMesh->elements();

            Cell< id_t >  & tNData = mElementExtra( aTarget )->mNeighborData ;
            Cell< id_t >  & tEData = mElementExtra( aTarget )->mEdgeData ;
            Cell< id_t >  & tFData = mElementExtra( aTarget )->mFaceData ;
            Cell< uchar > & tDData = mElementExtra( aTarget )->mEdgeDirections ;

            // count memory needs
            index_t n = 1 ; // neighbor counter
            index_t e = 1 ; // edge counter
            index_t f = 1 ; // face counter
            index_t d = 0 ; // direction counter

            index_t cn = 0 ;
            index_t ce = 0 ;
            index_t cf = 0 ;

            for ( index_t t : tTable )
            {
                Element * tElement = tElements( t );

                if ( tElement->owner() == aTarget )
                {
                    ++cn ;
                    n += tElement->number_of_elements() + 2 ;
                }

                if ( tElement->has_edges() )
                {
                    ++ce ;
                    d += tElement->number_of_edges() ;
                    e += tElement->number_of_edges() + 1 ;
                }

                if ( tElement->has_faces() )
                {
                    ++cf ;
                    f += tElement->number_of_faces() + 1 ;
                }
            }

            tNData.set_size( n, gNoID );
            tEData.set_size( e, gNoID );
            tDData.set_size( d, 0 );

            tFData.set_size( f, gNoID );

            tNData( 0 ) = cn ;
            tEData( 0 ) = ce ;
            tFData( 0 ) = cf ;

            n = 1 ;
            e = 1 ;
            f = 1 ;
            d = 0 ;

            for ( index_t t : tTable )
            {
                Element * tElement = tElements( t );

                if ( tElement->owner() == aTarget )
                {
                    tNData( n++ ) = tElement->id();
                    tNData( n++ ) = tElement->number_of_elements();
                    for ( uint i=0; i<tElement->number_of_elements(); ++i )
                    {
                        tNData( n++ ) = tElement->element( i )->id();
                    }
                }

                if ( tElement->has_edges() )
                {
                    tEData( e++ ) = tElement->id();
                    for ( uint i=0; i<tElement->number_of_edges(); ++i )
                    {
                        tEData( e++ ) = tElement->edge( i )->id();
                        tDData( d++ ) = tElement->edge_direction( i ) ? 1 : 0 ;
                    }
                }

                if ( tElement->has_faces() )
                {
                    tFData( f++ ) = tElement->id();
                    for ( uint i=0; i<tElement->number_of_faces(); ++i )
                    {
                        tFData( f++ ) = tElement->face( i )->id();
                    }
                }
            }
         }

//-----------------------------------------------------------------------------

        void
        Distributor::populate_facet_extra( const proc_t aTarget )
        {
            if ( mNumberOfAllFacets == 0 ) return ;

            Cell< index_t > & tTable = mTables( aTarget )->facets();
            Cell< Facet * > & tFacets = mMesh->facets();

            Cell< id_t > & tNData = mFacetExtra( aTarget )->mNeighborData ;

            // count memory needs
            index_t tNumFacets = 0 ;
            index_t tCount = 1 ;

            for ( index_t t : tTable )
            {
                Facet * tFacet = tFacets( t );

                if ( tFacet->owner() == aTarget )
                {
                    ++tNumFacets ;
                    tCount += tFacet->number_of_facets() + 2 ;
                }
            }

            tNData.set_size( tCount, gNoID );
            tCount = 0 ;

            tNData( tCount++ ) = tNumFacets ;

            for ( index_t t : tTable )
            {
                Facet * tFacet = tFacets( t );

                if ( tFacet->owner() == aTarget )
                {
                    tNData( tCount++ ) = tFacet->id();
                    tNData( tCount++ ) = tFacet->number_of_facets();
                    for ( uint i=0; i<tFacet->number_of_facets(); ++i )
                    {
                        tNData( tCount++ ) = tFacet->facet( i )->id();
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::populate_vertex_data( const proc_t aTarget )
        {
            if ( mNumberOfAllVertices == 0 ) return ;

            Cell< index_t >  & tTable    = mTables( aTarget )->vertices();
            Cell< Element * >  & tVertices = mMesh->vertices();

            mVertexBitset->where( tTable );

            Cell< id_t >   & tIdata   = mVertexData( aTarget )->mIDs ;
            Cell< proc_t > & tOdata   = mVertexData( aTarget )->mOwners ;
            Cell< uint >   & tGeo     = mVertexData( aTarget )->mGeometryTags ;
            Cell< uint >   & tPhys    = mVertexData( aTarget )->mPhysicalTags ;
            Cell< id_t >   & tTopo    = mVertexData( aTarget )->mTopology ;

            index_t tCount = tTable.size() ;
            tIdata.set_size( tCount, gNoID );
            tOdata.set_size( tCount, gNoOwner );
            tTopo.set_size( tCount, gNoID );
            tGeo.set_size( tCount, 0 );
            tPhys.set_size( tCount, 0 );

            tCount = 0 ;

            for ( index_t v : tTable )
            {
                Element * tVertex = tVertices( v );

                tIdata( tCount ) = tVertex->id();
                tOdata( tCount ) = tVertex->owner();
                tGeo( tCount )   = tVertex->geometry_tag();
                tPhys( tCount )  = tVertex->physical_tag();

                // note that vertices only have one node,
                // but the ID might differ!
                tTopo( tCount )  = tVertex->node( 0 )->id();
                ++tCount;
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::populate_control_point_data( const proc_t aTarget )
        {
            if ( mNumberOfAllControlPoints == 0 ) return ;

            Cell< index_t >  & tTable    = mTables( aTarget )->control_points();

            Cell< ControlPoint * >  & tControlPoints = mMesh->control_points();

            mControlPointBitset->where( tTable );

            // if no control points are assigned to this rank, exit early
            if ( tTable.size() == 0 ) return ;

            Cell< id_t >   & tIdata   = mControlPointData( aTarget )->mIDs ;
            Cell< proc_t > & tOdata   = mControlPointData( aTarget )->mOwners ;
            Matrix< real > & tCoords  = mControlPointData( aTarget )->mCoords ;

            index_t tCount = tTable.size() ;
            tIdata.set_size( tCount, gNoID );
            tOdata.set_size( tCount, gNoOwner );
            tCoords.set_size( 3, tCount );

            tCount = 0 ;

            for ( id_t p : tTable )
            {
                ControlPoint * tPoint = tControlPoints( p );

                tIdata( tCount ) = tPoint->id();
                tOdata( tCount ) = tPoint->owner();
                tCoords( 0, tCount ) = tPoint->x();
                tCoords( 1, tCount ) = tPoint->y();
                tCoords( 2, tCount ) = tPoint->z();
                ++tCount;
            }

            Cell< id_t > & tTopo = mControlPointData( aTarget )->mElementTopology ;

            // compute size - only for elements belonging to this target rank
            Cell< index_t > & tElementTable = mTables( aTarget )->elements();
            Cell< Element * > & tElements = mMesh->elements();

            tCount = 0 ;
            for ( index_t e : tElementTable )
            {
                tCount += 2 + tElements( e )->number_of_control_points();
            }
            tTopo.set_size( tCount, gNoID );
            tCount = 0 ;
            for ( index_t e : tElementTable )
            {
                Element * tElement = tElements( e );
                tTopo( tCount++ ) = tElement->id();
                tTopo( tCount++ ) = tElement->number_of_control_points();
                for ( uint i=0; i<tElement->number_of_control_points(); ++i )
                {
                    tTopo( tCount++ ) = tElement->control_point( i )->id();
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::populate_edge_data( const proc_t aTarget )
        {
            if ( mNumberOfAllEdges == 0 ) return ;

            Cell< index_t > & tTable    = mTables( aTarget )->edges() ;
            mEdgeBitset->where( tTable );

            // get container
            Cell< Edge * > & tEdges = mMesh->edges();

            Cell< id_t >    & tIDs      = mEdgeData( aTarget )->mIDs ;
            Cell< proc_t >  & tOwners   = mEdgeData( aTarget )->mOwners ;
            Cell< id_t >    & tTopo     = mEdgeData( aTarget )->mTopology ;



            index_t tCount = tTable.size() ;
            index_t tTopoCount = 0 ;

            for ( index_t d: tTable )
            {
                tTopoCount += tEdges( d )->number_of_nodes() + 1 ;
            }

            tIDs.set_size( tCount, gNoID );
            tOwners.set_size( tCount, gNoOwner );
            tTopo.set_size( tTopoCount, gNoID );

            tCount = 0 ;
            tTopoCount = 0;

            for ( index_t i: tTable )
            {
                Edge * tEdge = tEdges( i ) ;

                tIDs( tCount )     = tEdge->id();
                tOwners( tCount )  = tEdge->owner();

                ++tCount;

                tTopo( tTopoCount++ ) = tEdge->number_of_nodes();

                for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                {
                    tTopo( tTopoCount++ ) = tEdge->node( k )->id();
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::populate_face_data( const proc_t aTarget )
        {
            if ( mNumberOfAllFaces == 0 ) return ;

            // get indices

            Cell< index_t > & tTable    = mTables( aTarget )->faces() ;
            mFaceBitset->where( tTable );

            // get container
            Cell< Face * > & tFaces = mMesh->faces();

            Cell< id_t >    & tIDs       = mFaceData( aTarget )->mIDs ;
            Cell< proc_t >  & tOwners    = mFaceData( aTarget )->mOwners ;

            Cell< id_t >    & tMasterIDs = mFaceData( aTarget )->mMasterIDs ;
            Cell< id_t >    & tSlaveIDs  = mFaceData( aTarget )->mSlaveIDs ;

            Cell< uchar >   & tIndicesOnMaster = mFaceData( aTarget )->mIndicesOnMaster ;
            Cell< uchar >   & tIndicesOnSlave = mFaceData( aTarget )->mIndicesOnSlave ;
            Cell< uchar >   & tOrientationsOnSlave = mFaceData( aTarget )->mOrientationsOnSlave ;

            index_t tCount = tTable.size() ;

            tIDs.set_size( tCount, gNoID );
            tOwners.set_size( tCount, gNoOwner );

            tMasterIDs.set_size( tCount, gNoID );

            if ( mMesh->number_of_dimensions() == 3 )
            {
                tSlaveIDs.set_size( tCount, gNoID );
                tIndicesOnMaster.set_size( tCount, BELFEM_UCHAR_MAX );
                tIndicesOnSlave.set_size( tCount, BELFEM_UCHAR_MAX );
                tOrientationsOnSlave.set_size( tCount, BELFEM_UCHAR_MAX );
            }

            tCount = 0 ;

            if ( mMesh->number_of_dimensions() == 2 )
            {
                for ( index_t i: tTable )
                {
                    Face * tFace = tFaces( i ) ;

                    tIDs( tCount )     = tFace->id();
                    tOwners( tCount )  = tFace->owner();
                    tMasterIDs( tCount ) = tFace->master()->id();

                    ++tCount;
                }
            }
            else
            {
                for ( index_t i: tTable )
                {
                    Face * tFace = tFaces( i ) ;

                    tIDs( tCount )     = tFace->id();
                    tOwners( tCount )  = tFace->owner();

                    if ( tFace->master() != nullptr )
                    {
                        tIndicesOnMaster( tCount ) = tFace->index_on_master();
                        if ( mElementBitset->test( tFace->master()->index() ) )
                        {
                            tMasterIDs( tCount ) = tFace->master()->id();

                        }

                    }

                    if ( tFace->slave() != nullptr )
                    {
                        tIndicesOnSlave( tCount )      = tFace->index_on_slave() ;
                        tOrientationsOnSlave( tCount ) = tFace->orientation_on_slave() ;
                        if ( mElementBitset->test( tFace->slave()->index() ) )
                        {
                            tSlaveIDs( tCount )            = tFace->slave()->id() ;
                        }
                    }

                    ++tCount;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::populate_facet_data( const proc_t aTarget )
        {
            if ( mNumberOfAllFacets == 0 ) return ;

            Cell< index_t > & tTable  = mTables( aTarget )->facets() ;
            mFacetBitset->where( tTable );

            // get container
            Cell< Facet * > & tFacets = mMesh->facets();

            // sort the table by sideset INDEX (position in mMesh->sidesets()), so the
            // facets arrive grouped in the same order create_sidesets() iterates the
            // sideset data — independent of how the sideset ids are numbered.
            Map< id_t, index_t > tSideSetIndex ;
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                tSideSetIndex[ tSideSet->id() ] = tSideSet->index() ;
            }

            Cell< std::pair< index_t, index_t > > tOrder ;
            tOrder.reserve( tTable.size() );

            for ( index_t e : tTable )
            {
                tOrder.push( std::make_pair( e, tSideSetIndex( tFacets( e )->sideset_id() ) ) );
            }
            sort( tOrder.begin(), tOrder.end() ,
                  []( const std::pair< index_t, index_t > & a, const std::pair< index_t, index_t > & b )
                  {
                      return a.second < b.second ;
                  } );

            index_t tCount = 0 ;
            for ( auto tPair : tOrder )
            {
                tTable( tCount++ ) = tPair.first ;
            }
            tOrder.clear();

            Cell< id_t >    & tIDs      = mFacetData( aTarget )->mIDs ;
            Cell< proc_t >  & tOwners   = mFacetData( aTarget )->mOwners ;
            Cell< uchar >   & tTypes    = mFacetData( aTarget )->mTypes ;

            Cell< id_t >    & tMasterIDs = mFacetData( aTarget )->mMasterIDs ;
            Cell< id_t >    & tSlaveIDs  = mFacetData( aTarget )->mSlaveIDs ;

            Cell< uchar >   & tIndicesOnMaster = mFacetData( aTarget )->mIndicesOnMaster ;
            Cell< uchar >   & tIndicesOnSlave = mFacetData( aTarget )->mIndicesOnSlave ;
            Cell< uchar >   & tOrientationOnSlave = mFacetData( aTarget )->mOrientationsOnSlave ;

            Cell< uint >    & tGeo = mFacetData( aTarget )->mGeometryTags ;
            Cell< uint >    & tPhys = mFacetData( aTarget )->mPhysicalTags ;
            Cell< id_t >    & tTopo = mFacetData( aTarget )->mTopology ;

            tIDs.set_size( tCount, gNoID );
            tOwners.set_size( tCount, gNoOwner );
            tTypes.set_size( tCount, 0 );

            tMasterIDs.set_size( tCount, gNoID );
            tSlaveIDs.set_size( tCount, gNoID );

            tIndicesOnMaster.set_size( tCount, BELFEM_UCHAR_MAX );
            tIndicesOnSlave.set_size( tCount, BELFEM_UCHAR_MAX );
            tOrientationOnSlave.set_size( tCount, BELFEM_UCHAR_MAX );

            tGeo.set_size( tCount, 0 );
            tPhys.set_size( tCount, 0 );

            index_t tTopoCount = 0;

            for ( index_t f: tTable )
            {
                tTopoCount += tFacets( f )->number_of_nodes() ;
            }

            tCount = 0 ;
            tTopo.set_size( tTopoCount, gNoID );
            tTopoCount = 0;

            for ( index_t f : tTable )
            {
                Facet * tFacet = tFacets( f );

                tIDs( tCount )     = tFacet->id();
                tOwners( tCount )  = tFacet->owner();
                tTypes( tCount ) = static_cast< uchar >( tFacet->element()->type() );

                tGeo( tCount ) = tFacet->element()->geometry_tag();
                tPhys( tCount ) = tFacet->element()->physical_tag();


                if ( tFacet->master() != nullptr )
                {
                    tMasterIDs( tCount ) = tFacet->master()->id();
                    tIndicesOnMaster( tCount ) = tFacet->index_on_master();
                }

                if ( tFacet->slave() != nullptr )
                {
                    tSlaveIDs( tCount ) = tFacet->slave()->id() ;
                    tIndicesOnSlave( tCount ) = tFacet->index_on_slave() ;
                    tOrientationOnSlave( tCount ) = tFacet->orientation_on_slave() ;
                }

                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    tTopo( tTopoCount++ ) = tFacet->node( k )->id();
                }
                ++tCount;
            }
        }

        void
        Distributor::send_thinshell_data()
        {
            Cell< id_t >   tIdata ;
            Cell< real >   tRdata ;
            Cell< string > tSdata ;

            index_t i = 1 ;
            index_t j = 0 ;

            // count memory
            for ( ThinShell * tShell : mMesh->thin_shells() )
            {
                i += tShell->blocks().size() + 3 ;
                j += tShell->blocks().size();
            }

            tIdata.set_size( i, gNoID );
            tRdata.set_size( j, 0 );
            tSdata.set_size( j, "" );

            i = 0 ;
            j = 0 ;

            tIdata( i++ ) = mMesh->thin_shells().size();

            for ( ThinShell * tShell : mMesh->thin_shells() )
            {
                tIdata( i++ ) = tShell->id();
                tIdata( i++ ) = tShell->ghost_id();
                uint n = tShell->blocks().size();
                tIdata( i++ ) = n ;
                for ( uint b=0; b<n; b++ )
                {
                    tIdata( i++ ) = tShell->blocks()( b )->id();
                    tRdata( j )   = tShell->thicknesses()( b );
                    tSdata( j++ ) = tShell->materials()( b );
                }
            }

            comm_barrier();
            broadcast( tIdata );
            broadcast( tRdata );
            broadcast( tSdata );
        }

        void
        Distributor::receive_thinshell_data()
        {
            Cell< id_t >   tIdata ;
            Cell< real >   tRdata ;
            Cell< string > tSdata ;

            comm_barrier();

            broadcast( tIdata );
            broadcast( tRdata );
            broadcast( tSdata );

            index_t i = 0 ;
            index_t j = 0 ;


            mProtoMesh->thin_shell_data().set_size( tIdata( i++ ), proto::ThinShellData() );

            for ( proto::ThinShellData & tData : mProtoMesh->thin_shell_data() )
            {
                // get the sideset
                tData.mSideSetID = tIdata( i++ );

                // get the sideset
                tData.mGhostSideSetID = tIdata( i++ );

                // get the number of blocks
                uint n = tIdata( i++ );

                tData.mBlocksIDs.set_size( n, gNoID );
                tData.mThicknesses.set_size( n, BELFEM_QUIET_NAN );
                tData.mMaterials.set_size( n, "" );

                for ( uint b=0; b<n; b++ )
                {
                    tData.mBlocksIDs( b ) = tIdata( i++ );
                    tData.mThicknesses( b ) = tRdata( j );
                    tData.mMaterials( b ) = tSdata( j++ );
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::send_block_data()
        {
            if ( mNumberOfAllBlocks == 0 ) return ;

            index_t tCount = mMesh->number_of_blocks();
            Cell< id_t > tIDs ;
            tIDs.set_size( tCount, gNoID );

            Cell< uchar > tEtypes( tCount, 0 );
            Cell< uchar > tDtypes( tCount, 0 );

            // physical thickness ( thin-shell layers and side connector
            // walls; NaN otherwise ) — EF_HEX8TB reads the wall width from
            // the block, so the value must reach the other procs
            Cell< real > tThickness( tCount, BELFEM_QUIET_NAN );

            tCount = 0 ;
            for ( Block * tBlock: mMesh->blocks() )
            {
                tIDs( tCount ) = tBlock->id();
                tEtypes( tCount ) = static_cast< uchar >( tBlock->element_type() );
                tDtypes( tCount ) = static_cast< uchar >( tBlock->domain_type() );
                tThickness( tCount ) = tBlock->thickness();
                ++tCount;
            }

            share( tIDs );
            share( tEtypes );
            share( tDtypes );
            share( tThickness );
        }

        void
        Distributor::receive_block_data()
        {
            if ( mNumberOfAllBlocks == 0 ) return ;

            Cell< id_t > tIDs ;
            Cell< uchar > tEtypes ;
            Cell< uchar > tDtypes ;
            Cell< real > tThickness ;

            receive( tIDs );
            receive( tEtypes );
            receive( tDtypes );
            receive( tThickness );

            uint tNumBlocks = tIDs.size() ;

            mProtoMesh->block_data().set_size( tNumBlocks, proto::GroupData() );

            for ( uint b=0; b<tNumBlocks; ++b )
            {
                mProtoMesh->block_data( b ).mID = tIDs( b );
                mProtoMesh->block_data( b ).mElementType = static_cast< ElementType >( tEtypes( b ) );
                mProtoMesh->block_data( b ).mDomainType = static_cast< DomainType >( tDtypes( b ) );
                mProtoMesh->block_data( b ).mThickness = tThickness( b );
            }
        }

        void
        Distributor::send_sideset_data()
        {
            if ( mNumberOfAllSideSets == 0 ) return ;

            index_t tCount = mMesh->number_of_sidesets();
            Cell< id_t > tIDs ;
            tIDs.set_size( tCount, gNoID );
            Cell< uchar > tEtypes( tCount, 0 );
            Cell< uchar > tDtypes( tCount, 0 );

            tCount = 0 ;
            for ( SideSet * tSideSet: mMesh->sidesets() )
            {
                tIDs( tCount ) = tSideSet->id();
                tEtypes( tCount ) = static_cast< uchar >( tSideSet->element_type() );
                tDtypes( tCount ) = static_cast< uchar >( tSideSet->domain_type() );
                ++tCount;
            }

            broadcast( tIDs );
            broadcast( tEtypes );
            broadcast( tDtypes );
        }

        void
        Distributor::receive_sideset_data()
        {
            if ( mNumberOfAllSideSets == 0 ) return ;

             Cell< id_t > tIDs ;
             Cell< uchar > tEtypes ;
             Cell< uchar > tDtypes ;

             broadcast( tIDs );
             broadcast( tEtypes );
             broadcast( tDtypes );

             uint tNumSideSets = tIDs.size() ;

             mProtoMesh->sideset_data().set_size( tNumSideSets, proto::GroupData() );

             for ( uint s=0; s<tNumSideSets; ++s )
             {
                 mProtoMesh->sideset_data( s ).mID = tIDs( s );
                 mProtoMesh->sideset_data( s ).mElementType = static_cast< ElementType >( tEtypes( s ) );
                 mProtoMesh->sideset_data( s ).mDomainType = static_cast< DomainType >( tDtypes( s ) );
             }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::send_node_data()
        {
            if ( mNumberOfAllNodes == 0 ) return ;

            Cell< Cell< id_t > >   tIData( mCommSize, {} );
            Cell< Cell< proc_t > > tOData( mCommSize, {} );
            Cell< Matrix< real > > tRData( mCommSize, {} );
            Cell< Cell< id_t > >   tDData( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tIData( p ) = std::move( mNodeData( p )->mIDs );
                tOData( p ) = std::move( mNodeData( p )->mOwners );
                tRData( p ) = std::move( mNodeData( p )->mCoords );
                tDData( p ) = std::move( mNodeData( p )->mDuplicateData );
            }

            distribute( tIData );
            distribute( tOData );
            distribute( tRData );
            distribute( tDData );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::receive_node_data()
        {
            if ( mNumberOfAllNodes == 0 ) return ;

            receive( mProtoMesh->node_data()->mIDs );
            receive( mProtoMesh->node_data()->mOwners );
            receive( mProtoMesh->node_data()->mCoords );
            receive( mProtoMesh->node_data()->mDuplicateData );
        }
        
//-----------------------------------------------------------------------------

        void
        Distributor::send_element_data()
        {
            if ( mNumberOfAllElements == 0 ) return ;

            Cell< Cell< id_t > >   tIData( mCommSize, {} );
            Cell< Cell< proc_t > > tOData( mCommSize, {} );
            Cell< Cell< uchar > >  tTData( mCommSize, {} );
            Cell< Cell< uint > >   tGData( mCommSize, {} );
            Cell< Cell< uint > >   tPData( mCommSize, {} );
            Cell< Cell< id_t > >   tYData( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tIData( p ) = std::move( mElementData( p )->mIDs );
                tOData( p ) = std::move( mElementData( p )->mOwners );
                tTData( p ) = std::move( mElementData( p )->mTypes );
                tGData( p ) = std::move( mElementData( p )->mGeometryTags );
                tPData( p ) = std::move( mElementData( p )->mPhysicalTags );
                tYData( p ) = std::move( mElementData( p )->mTopology );
            }

            distribute( tIData );
            distribute( tOData );
            distribute( tTData );
            distribute( tGData );
            distribute( tPData );
            distribute( tYData );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::receive_element_data()
        {
            if ( mNumberOfAllElements == 0 ) return ;

            receive( mProtoMesh->element_data()->mIDs );
            receive( mProtoMesh->element_data()->mOwners );
            receive( mProtoMesh->element_data()->mTypes );
            receive( mProtoMesh->element_data()->mGeometryTags );
            receive( mProtoMesh->element_data()->mPhysicalTags );
            receive( mProtoMesh->element_data()->mTopology );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::send_edge_data()
        {
            if ( mNumberOfAllEdges == 0 ) return ;

            Cell< Cell< id_t > >   tIData( mCommSize, {} );
            Cell< Cell< proc_t > > tOData( mCommSize, {} );
            Cell< Cell< uint > >   tYData( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tIData( p ) = std::move( mEdgeData( p )->mIDs );
                tOData( p ) = std::move( mEdgeData( p )->mOwners );
                tYData( p ) = std::move( mEdgeData( p )->mTopology );
            }

            distribute( tIData );
            distribute( tOData );
            distribute( tYData );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::receive_edge_data()
        {
            if ( mNumberOfAllEdges == 0 ) return ;

            receive( mProtoMesh->edge_data()->mIDs );
            receive( mProtoMesh->edge_data()->mOwners );
            receive( mProtoMesh->edge_data()->mTopology );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::send_face_data()
        {
            if ( mNumberOfAllFaces == 0 ) return ;

            Cell< Cell< id_t > >   tIData( mCommSize, {} );
            Cell< Cell< proc_t > > tOData( mCommSize, {} );


            Cell< Cell< id_t > >   tMData( mCommSize, {} );
            Cell< Cell< id_t > >   tSData( mCommSize, {} );

            Cell< Cell< uchar > >   tJData( mCommSize, {} );
            Cell< Cell< uchar > >   tKData( mCommSize, {} );
            Cell< Cell< uchar > >   tLData( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tIData( p ) = std::move( mFaceData( p )->mIDs );
                tOData( p ) = std::move( mFaceData( p )->mOwners );
                tMData( p ) = std::move( mFaceData( p )->mMasterIDs );
                tSData( p ) = std::move( mFaceData( p )->mSlaveIDs );
                tJData( p ) = std::move( mFaceData( p )->mIndicesOnMaster );
                tKData( p ) = std::move( mFaceData( p )->mIndicesOnSlave );
                tLData( p ) = std::move( mFaceData( p )->mOrientationsOnSlave );
            }

            distribute( tIData );
            distribute( tOData );
            distribute( tMData );
            distribute( tSData );
            distribute( tJData );
            distribute( tKData );
            distribute( tLData );

        }

        void
        Distributor::receive_face_data()
        {
            if ( mNumberOfAllFaces == 0 ) return ;

            receive( mProtoMesh->face_data()->mIDs );
            receive( mProtoMesh->face_data()->mOwners );
            receive( mProtoMesh->face_data()->mMasterIDs );
            receive( mProtoMesh->face_data()->mSlaveIDs );
            receive( mProtoMesh->face_data()->mIndicesOnMaster );
            receive( mProtoMesh->face_data()->mIndicesOnSlave );
            receive( mProtoMesh->face_data()->mOrientationsOnSlave );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::send_facet_data()
        {
            if ( mNumberOfAllFacets == 0 ) return ;

            Cell< Cell< id_t > >   tIData( mCommSize, {} );
            Cell< Cell< proc_t > > tOData( mCommSize, {} );


            Cell< Cell< id_t > >   tMData( mCommSize, {} );
            Cell< Cell< id_t > >   tSData( mCommSize, {} );

            Cell< Cell< uchar > >  tJData( mCommSize, {} );
            Cell< Cell< uchar > >  tKData( mCommSize, {} );
            Cell< Cell< uchar > >  tLData( mCommSize, {} );

            Cell< Cell< uchar > >  tTData( mCommSize, {} );
            Cell< Cell< uint > >   tGData( mCommSize, {} );
            Cell< Cell< uint > >   tPData( mCommSize, {} );

            Cell< Cell< id_t > >   tYData( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tIData( p ) = std::move( mFacetData( p )->mIDs );
                tOData( p ) = std::move( mFacetData( p )->mOwners );
                tMData( p ) = std::move( mFacetData( p )->mMasterIDs );
                tSData( p ) = std::move( mFacetData( p )->mSlaveIDs );
                tJData( p ) = std::move( mFacetData( p )->mIndicesOnMaster );
                tKData( p ) = std::move( mFacetData( p )->mIndicesOnSlave );
                tLData( p ) = std::move( mFacetData( p )->mOrientationsOnSlave );
                tTData( p ) = std::move( mFacetData( p )->mTypes );
                tGData( p ) = std::move( mFacetData( p )->mGeometryTags );
                tPData( p ) = std::move( mFacetData( p )->mPhysicalTags );
                tYData( p ) = std::move( mFacetData( p )->mTopology );
            }

            distribute( tIData );
            distribute( tOData );

            distribute( tMData );
            distribute( tSData );

            distribute( tJData );
            distribute( tKData );
            distribute( tLData );

            distribute( tTData );
            distribute( tGData );
            distribute( tPData );

            distribute( tYData );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::receive_facet_data()
        {
            if ( mNumberOfAllFacets == 0 ) return ;

            receive( mProtoMesh->facet_data()->mIDs );
            receive( mProtoMesh->facet_data()->mOwners );

            receive( mProtoMesh->facet_data()->mMasterIDs );
            receive( mProtoMesh->facet_data()->mSlaveIDs );

            receive( mProtoMesh->facet_data()->mIndicesOnMaster );
            receive( mProtoMesh->facet_data()->mIndicesOnSlave );
            receive( mProtoMesh->facet_data()->mOrientationsOnSlave );

            receive( mProtoMesh->facet_data()->mTypes );
            receive( mProtoMesh->facet_data()->mGeometryTags );
            receive( mProtoMesh->facet_data()->mPhysicalTags );

            receive( mProtoMesh->facet_data()->mTopology );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::send_element_extra( const proc_t aTarget )
        {
            if ( mNumberOfAllElements == 0 ) return ;

            if ( aTarget == 0 )  // sending everything in parallel
            {
                Cell< Cell< id_t > >  tNData( mCommSize, {} );
                Cell< Cell< id_t > >  tEData( mCommSize, {} );
                Cell< Cell< id_t > >  tFData( mCommSize, {} );
                Cell< Cell< uchar > > tDData( mCommSize, {} );
                Cell< Cell< id_t > >  tCData( mCommSize, {} );

                for ( proc_t p=1; p<mCommSize; ++p )
                {
                    tNData( p ) = std::move( mElementExtra( p )->mNeighborData );
                    tEData( p ) = std::move( mElementExtra( p )->mEdgeData );
                    tFData( p ) = std::move( mElementExtra( p )->mFaceData );
                    tDData( p ) = std::move( mElementExtra( p )->mEdgeDirections );
                    tCData( p ) = std::move( mElementExtra( p )->mCurvedElementIDs );
                }

                distribute( tNData );
                distribute( tEData );
                distribute( tFData );
                distribute( tDData );

                if ( mMaxElementOrder == 1 ) return ;
                distribute( tCData );
            }
            else  // sending everything sequential
            {
                Cell< id_t > tNData = std::move( mElementExtra( aTarget )->mNeighborData );
                Cell< id_t > tEData = std::move( mElementExtra( aTarget )->mEdgeData );
                Cell< id_t > tFData = std::move( mElementExtra( aTarget )->mFaceData );
                Cell< uchar > tDData = std::move( mElementExtra( aTarget )->mEdgeDirections );
                Cell< id_t > tCData = std::move( mElementExtra( aTarget )->mCurvedElementIDs );

                send( tNData , aTarget );
                send( tEData, aTarget );
                send( tFData, aTarget );
                send( tDData, aTarget );

                if ( mMaxElementOrder == 1 ) return ;
                send( tCData, aTarget );
            }
        }

//-----------------------------------------------------------------------------

        void
        Distributor::receive_element_extra()
        {
            if ( mNumberOfAllElements == 0 ) return ;

            receive( mProtoMesh->element_extra()->mNeighborData );
            receive( mProtoMesh->element_extra()->mEdgeData );
            receive( mProtoMesh->element_extra()->mFaceData );
            receive( mProtoMesh->element_extra()->mEdgeDirections );

            if ( mMaxElementOrder == 1 ) return ;

            receive( mProtoMesh->element_extra()->mCurvedElementIDs );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::send_facet_extra()
        {
            if ( mNumberOfAllFacets == 0 ) return ;

            Cell< Cell< id_t > > tData( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tData( p ) = std::move( mFacetExtra( p )->mNeighborData );
            }
            distribute( tData );

            if ( mMaxElementOrder == 1 ) return ;

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tData( p ).clear();
                tData( p ) = std::move( mFacetExtra( p )->mCurvedFacetIDs );
            }
            distribute( tData );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::receive_facet_extra()
        {
            if ( mNumberOfAllFacets == 0 ) return ;

            receive( mProtoMesh->facet_extra()->mNeighborData );

            if ( mMaxElementOrder == 1 ) return ;

            receive( mProtoMesh->facet_extra()->mCurvedFacetIDs );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::send_vertex_data()
        {
            if ( mNumberOfAllVertices == 0 ) return ;

            Cell< Cell< id_t >  > tIdata( mCommSize, {} );
            Cell< Cell< proc_t > > tOdata( mCommSize, {} );
            Cell< Cell< uint > >  tGeo( mCommSize, {} );
            Cell< Cell< uint > >  tPhys( mCommSize, {} );
            Cell< Cell< id_t > >  tTopo( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tIdata( p ) = std::move( mVertexData( p )->mIDs );
                tOdata( p ) = std::move( mVertexData( p )->mOwners );
                tGeo( p )   = std::move( mVertexData( p )->mGeometryTags );
                tPhys( p )  = std::move( mVertexData( p )->mPhysicalTags );
                tTopo( p )  = std::move( mVertexData( p )->mTopology );
            }

            distribute( tIdata );
            distribute( tOdata );
            distribute( tGeo );
            distribute( tPhys );
            distribute( tTopo );
        }

        void
        Distributor::receive_vertex_data()
        {
            if ( mNumberOfAllVertices == 0 ) return ;

            receive( mProtoMesh->vertex_data()->mIDs );
            receive( mProtoMesh->vertex_data()->mOwners );
            receive( mProtoMesh->vertex_data()->mGeometryTags );
            receive( mProtoMesh->vertex_data()->mPhysicalTags );
            receive( mProtoMesh->vertex_data()->mTopology );
        }

        void
        Distributor::send_control_point_data()
        {
            if ( mNumberOfAllControlPoints == 0 ) return ;

            Cell< Cell< id_t >  >  tIdata( mCommSize, {} );
            Cell< Cell< proc_t > > tOdata( mCommSize, {} );
            Cell< Matrix< real > > tCoords( mCommSize, {} );
            Cell< Cell< id_t > >   tTopo( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tIdata( p ) = std::move( mControlPointData( p )->mIDs );
                tOdata( p ) = std::move( mControlPointData( p )->mOwners );
                tCoords( p ) = std::move( mControlPointData( p )->mCoords );
                tTopo( p )  = std::move( mControlPointData( p )->mElementTopology );
            }
            distribute( tIdata );
            distribute( tOdata );
            distribute( tCoords );
            distribute( tTopo );
        }

        void
        Distributor::send_t_matrices()
        {
            if ( mNumberOfAllTMatrices == 0 ) return ;

            Cell< Cell< index_t >  >  tNumTargets( mCommSize, {} );
            Cell< Cell< id_t >  >  tTargetIDs( mCommSize, {} );
            Cell< Cell< uint >  >  tCounters( mCommSize, {} );
            Cell< Cell< id_t >  >  tSourceIDs( mCommSize, {} );
            Cell< Cell< uchar >  > tTypes( mCommSize, {} );
            Cell< Cell< real > >   tWeights( mCommSize, {} );

            for ( proc_t p=1; p<mCommSize; ++p )
            {
                tNumTargets( p ) = std::move( mTMatrices( p )->mNumTargets );
                tTargetIDs( p ) = std::move( mTMatrices( p )->mTargetIDs );
                tCounters( p )  = std::move( mTMatrices( p )->mCounters );
                tSourceIDs( p ) = std::move( mTMatrices( p )->mSourceIDs );
                tTypes( p )     = std::move( mTMatrices( p )->mTypes );
                tWeights( p )   = std::move( mTMatrices( p )->mWeights );
            }
            comm_barrier() ;
            distribute( tNumTargets );
            distribute( tTargetIDs );
            distribute( tCounters );
            distribute( tSourceIDs );
            distribute( tTypes );
            distribute( tWeights );
        }

        void
        Distributor::receive_t_matrices()
        {
            if ( mNumberOfAllTMatrices == 0 ) return ;
            comm_barrier() ;
            receive( mProtoMesh->t_matrix_data()->mNumTargets );
            receive( mProtoMesh->t_matrix_data()->mTargetIDs );
            receive( mProtoMesh->t_matrix_data()->mCounters );
            receive( mProtoMesh->t_matrix_data()->mSourceIDs );
            receive( mProtoMesh->t_matrix_data()->mTypes );
            receive( mProtoMesh->t_matrix_data()->mWeights );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::receive_control_point_data()
        {
            if ( mNumberOfAllControlPoints == 0 ) return ;

            receive( mProtoMesh->control_point_data()->mIDs );
            receive( mProtoMesh->control_point_data()->mOwners );
            receive( mProtoMesh->control_point_data()->mCoords );
            receive( mProtoMesh->control_point_data()->mElementTopology );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::flag_curved_elements(  proc_t aTarget )
        {
            if ( mMaxElementOrder == 1 ) return ;

            const Cell< index_t > & tIndices = mTables( aTarget )->elements();

            mElementBitset->reset();

            Cell< Element * > & tElements = mMesh->elements();

            for ( index_t e : tIndices )
            {
                if ( tElements( e )->is_curved() )
                {
                    mElementBitset->set( e );
                }
            }

            Cell< index_t > tCurvedIndices ;
            mElementBitset->where( tCurvedIndices );
            mElementExtra( aTarget )->mCurvedElementIDs.set_size( tCurvedIndices.size() );

            index_t tCount =  0 ;
            for ( index_t e : tCurvedIndices )
            {
                mElementExtra( aTarget )->mCurvedElementIDs( tCount++ ) = tElements( e )->id();
            }

        }

        void
        Distributor::flag_curved_facets(  proc_t aTarget )
        {
            if ( mMaxElementOrder == 1 ) return ;
            if ( mFacetBitset == nullptr ) return ;

            const Cell< index_t > & tIndices = mTables( aTarget )->facets();

            mFacetBitset->reset();

            Cell< Facet * > & tFacets = mMesh->facets();

            for ( index_t f : tIndices )
            {
                if ( tFacets( f )->is_curved() )
                {
                    mFacetBitset->set( f );
                }
            }

            mFacetBitset->where( mFacetExtra( aTarget )->mCurvedFacetIDs );
        }

//-----------------------------------------------------------------------------

        void
        Distributor::create_bitsets()
        {
            this->delete_bitsets();

            if ( mMesh->number_of_nodes() > 0 )
            {
                mNodeBitset    = new DynamicBitset( mMesh->number_of_nodes() );
            }
            if ( mMesh->number_of_edges() > 0 )
            {
                mEdgeBitset    = new DynamicBitset( mMesh->number_of_edges() );
            }

            if ( mMesh->number_of_faces() > 0 )
            {
                mFaceBitset    = new DynamicBitset( mMesh->number_of_faces() );
            }

            if ( mMesh->number_of_elements() > 0 )
            {
                mElementBitset = new DynamicBitset( mMesh->number_of_elements() );
            }

            if ( mMesh->number_of_facets() > 0 )
            {
                mFacetBitset    = new DynamicBitset( mMesh->number_of_facets() );
            }

            if ( mMesh->vertices().size() > 0 )
            {
                mVertexBitset = new DynamicBitset( mMesh->vertices().size() );
            }
            if ( mMesh->number_of_control_points() > 0 )
            {
                mControlPointBitset = new DynamicBitset( mMesh->number_of_control_points() );
            }
        }

        void
        Distributor::delete_bitsets()
        {
            if ( mNodeBitset != nullptr )
            {
                delete mNodeBitset ;
                mNodeBitset = nullptr ;
            }

            if ( mEdgeBitset != nullptr )
            {
                delete mEdgeBitset ;
                mEdgeBitset = nullptr ;
            }

            if ( mFaceBitset != nullptr )
            {
                delete mFaceBitset ;
                mFaceBitset = nullptr ;
            }

            if ( mElementBitset != nullptr )
            {
                delete mElementBitset ;
                mElementBitset = nullptr ;
            }

            if ( mFacetBitset != nullptr )
            {
                delete mFacetBitset ;
                mFacetBitset = nullptr ;
            }

            if ( mVertexBitset != nullptr )
            {
                delete mVertexBitset ;
                mVertexBitset = nullptr ;
            }

            if ( mControlPointBitset != nullptr )
            {
                delete mControlPointBitset ;
                mControlPointBitset = nullptr ;
            }
        }

        void
        Distributor::reset_bitsets()
        {
            if ( mNodeBitset != nullptr )
            {
                mNodeBitset->reset() ;
            }

            if ( mEdgeBitset != nullptr )
            {
                mEdgeBitset->reset() ;
            }

            if ( mFaceBitset != nullptr )
            {
                mFaceBitset->reset() ;
            }

            if ( mElementBitset != nullptr )
            {
                mElementBitset->reset() ;
            }

            if ( mFacetBitset != nullptr )
            {
                mFacetBitset->reset() ;
            }

            if ( mVertexBitset != nullptr )
            {
                mVertexBitset->reset() ;
            }

            if ( mControlPointBitset != nullptr )
            {
                mControlPointBitset->reset() ;
            }
        }

        void
        Distributor::delete_node_data()
        {
            for ( auto * tData: mNodeData )
            {
                delete tData ;
            }
            mNodeData.clear();
        }

        void
        Distributor::delete_element_data()
        {
            for ( auto * tData: mElementData )
            {
                delete tData ;
            }
            mElementData.clear();
        }

        void
        Distributor::delete_element_extra()
        {
            for ( auto * tData: mElementExtra )
            {
                delete tData ;
            }
            mElementExtra.clear();
        }

        void
        Distributor::delete_facet_extra()
        {
            for ( auto * tData: mFacetExtra )
            {
                delete tData ;
            }
            mFacetExtra.clear();
        }

        void
        Distributor::delete_edge_data()
        {
            for ( auto * tData: mEdgeData )
            {
                delete tData ;
            }
            mEdgeData.clear();
        }

        void
        Distributor::delete_face_data()
        {
            for ( auto * tData: mFaceData )
            {
                delete tData ;
            }
            mFaceData.clear();
        }

        void
        Distributor::delete_facet_data()
        {
            for ( auto * tData: mFacetData )
            {
                delete tData ;
            }
            mFacetData.clear();
        }

        void
        Distributor::delete_vertex_data()
        {
            for ( auto * tData: mVertexData )
            {
                delete tData ;
            }
            mVertexData.clear();
        }

        void
        Distributor::delete_tables()
        {

            for ( auto * tTable: mTables )
            {
                delete tTable ;
            }
            mTables.clear();
        }

        void
        Distributor::delete_mesh_data()
        {
            if ( mProtoMesh != nullptr )
            {
                delete mProtoMesh ;
                mProtoMesh = nullptr ;
            }
        }

        void
        Distributor::delete_t_matrices()
        {
            for ( auto * tMatrix: mTMatrices )
            {
                delete tMatrix ;
            }
        }

//-------------------------------------------------------------------------------

        Mesh *
        Distributor::partial_mesh()
        {
            BELFEM_ERROR( mCommRank > 0, "Distributor::partial_mesh() must not be called by root proc");

            mOwnMesh = false ;
            return mMesh ;
        }

 //-------------------------------------------------------------------------------
    }
}