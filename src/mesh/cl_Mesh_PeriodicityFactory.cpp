/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <map>
#include <set>
#include <vector>

#include "cl_Logger.hpp"
#include "cl_Mesh_PeriodicityFactory.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_sort.hpp"
#include "op_Node_Index.hpp"

namespace belfem
{
    namespace mesh
    {
        PeriodicityFactory::PeriodicityFactory( Mesh * aMesh, ProtoMesh * aProtoMesh ) :
            mMesh( aMesh ),
            mProtoMesh( aProtoMesh )

        {

        }

        PeriodicityFactory::~PeriodicityFactory()
        {
            if ( mNodeBitset != nullptr ) delete mNodeBitset ;

            if ( mEdgeBitset != nullptr ) delete mEdgeBitset ;

            if ( mFaceBitset != nullptr ) delete mFaceBitset ;
        }

        void
        PeriodicityFactory::set_master_plane( const Vector< id_t > & aPointIDs )
        {
            BELFEM_ERROR( aPointIDs.length() == 3, "Master plane must have 3 points" );

            this->set_master_plane( aPointIDs( 0 ), aPointIDs( 1 ), aPointIDs( 2 ) );
        }

        void
        PeriodicityFactory::set_slave_plane( const Vector< id_t > & aPointIDs )
        {
            BELFEM_ERROR( aPointIDs.length() == 3, "Slave plane must have 3 points" );
            this->set_slave_plane( aPointIDs( 0 ), aPointIDs( 1 ), aPointIDs( 2 ) );
        }

        void
        PeriodicityFactory::set_master_plane( Cell< Node * > & aNodes )
        {
            BELFEM_ERROR( aNodes.size() == 3, "Master plane must have 3 nodes" );

            mSourcePlane = aNodes ;

            this->compute_transformation_matrix(
                mSourcePlane( 0 ),
                 mSourcePlane( 1 ),
                mSourcePlane( 2 ),
                mSourceDistance,
                mSourceTransform,
                mSourceHesse );

        }

        void
        PeriodicityFactory::set_slave_plane( Cell< Node * > & aNodes )
        {
            BELFEM_ERROR( aNodes.size() == 3, "Slave plane must have 3 nodes" );

            mTargetPlane = aNodes ;

            this->compute_transformation_matrix(
             mTargetPlane( 0 ),
              mTargetPlane( 1 ),
             mTargetPlane( 2 ),
             mTargetDistance,
             mTargetTransform,
             mTargetHesse );

        }

        void
        PeriodicityFactory::set_master_plane( const id_t A, const id_t B, const id_t C )
        {
            Node * tA = mMesh->vertex( A )->node( 0 );
            Node * tB = mMesh->vertex( B )->node( 0 );
            Node * tC = mMesh->vertex( C )->node( 0 );

            mSourcePlane = { tA, tB, tC };

            this->compute_transformation_matrix( tA, tB, tC,
                mSourceDistance,
                mSourceTransform,
                 mSourceHesse );

        }

        void
        PeriodicityFactory::set_slave_plane( const id_t A, const id_t B, const id_t C )
        {
            Node * tA = mMesh->vertex( A )->node( 0 );
            Node * tB = mMesh->vertex( B )->node( 0 );
            Node * tC = mMesh->vertex( C )->node( 0 );

            mTargetPlane = { tA, tB, tC };

            this->compute_transformation_matrix( tA, tB, tC,
                mTargetDistance,
                mTargetTransform,
                mTargetHesse );
        }

        Periodicity *
        PeriodicityFactory::create_periodicity()
        {
            BELFEM_ERROR( mSourcePlane.size() == 3, "Master plane must have 3 nodes" );
            BELFEM_ERROR( mTargetPlane.size() == 3, "Slave plane must have 3 nodes" );

            Periodicity * aPeriodicity = new Periodicity( mMesh ) ;

            aPeriodicity->master_plane() = mSourcePlane;
            aPeriodicity->slave_plane()  = mTargetPlane;
            aPeriodicity->master_hesse() = mSourceHesse;
            aPeriodicity->slave_hesse()  = mTargetHesse;

            this->update_periodicity( aPeriodicity );
            return aPeriodicity ;
        }

        proto::PeriodicityData *
        PeriodicityFactory::to_proto( Periodicity * aPeriodicity )
        {
            proto::PeriodicityData * aData = new proto::PeriodicityData() ;

            index_t tCount = 0 ;
            aData->mMasterPlane.set_size( aPeriodicity->master_plane().size() );
            for ( Node * tNode : aPeriodicity->master_plane() )
            {
                aData->mMasterPlane( tCount++ ) = tNode->id();
            }

            tCount = 0 ;
            aData->mSlavePlane.set_size( aPeriodicity->slave_plane().size() );
            for ( Node * tNode : aPeriodicity->slave_plane() )
            {
                aData->mSlavePlane( tCount++ ) = tNode->id();
            }

            tCount = 0 ;
            aData->mMasterFacets.set_size( aPeriodicity->master_facets().size() );
            for ( Facet * tFacet : aPeriodicity->master_facets() )
            {
                aData->mMasterFacets( tCount++ ) = tFacet->id();
            }

            tCount = 0 ;
            aData->mSlaveFacets.set_size( aPeriodicity->slave_facets().size() );
            for ( Facet * tFacet : aPeriodicity->slave_facets() )
            {
                aData->mSlaveFacets( tCount++ ) = tFacet->id();
            }

            tCount = 0 ;
            aData->mMasterNodes.set_size( aPeriodicity->master_nodes().size() );
            for ( Node * tNode : aPeriodicity->master_nodes() )
            {
                aData->mMasterNodes( tCount++ ) = tNode->id();
            }

            tCount = 0 ;
            aData->mSlaveNodes.set_size( aPeriodicity->slave_nodes().size() );
            for ( Node * tNode : aPeriodicity->slave_nodes() )
            {
                aData->mSlaveNodes( tCount++ ) = tNode->id();
            }

            if ( mMesh->edges_exist() )
            {
                tCount = 0 ;
                aData->mMasterEdges.set_size( aPeriodicity->master_edges().size() );
                for ( Edge * tEdge : aPeriodicity->master_edges() )
                {
                    aData->mMasterEdges( tCount++ ) = tEdge->id();
                }

                tCount = 0 ;
                aData->mSlaveEdges.set_size( aPeriodicity->slave_edges().size() );
                for ( Edge * tEdge : aPeriodicity->slave_edges() )
                {
                    aData->mSlaveEdges( tCount++ ) = tEdge->id();
                }
            }
            if ( mMesh->faces_exist() )
            {
                tCount = 0 ;
                aData->mMasterFaces.set_size( aPeriodicity->master_faces().size() );
                for ( Face * tFace : aPeriodicity->master_faces() )
                {
                    aData->mMasterFaces( tCount++ ) = tFace->id();
                }

                tCount = 0 ;
                aData->mSlaveFaces.set_size( aPeriodicity->slave_faces().size() );
                for ( Face * tFace : aPeriodicity->slave_faces() )
                {
                    aData->mSlaveFaces( tCount++ ) = tFace->id();
                }
            }


            return aData ;

        }

        Periodicity *
        PeriodicityFactory::from_proto( proto::PeriodicityData * aData, const bool aCrosslinkEntities )
        {
            // the planes are needed to recompute the hesse forms; node pairs
            // without plane data indicate a corrupted or incomplete proto
            BELFEM_ERROR( aData->mMasterPlane.size() == 3,
                "proto periodicity has invalid master plane data ( %lu points, need 3 )",
                ( long unsigned int ) aData->mMasterPlane.size() );
            BELFEM_ERROR( aData->mSlavePlane.size() == 3,
                "proto periodicity has invalid slave plane data ( %lu points, need 3 )",
                ( long unsigned int ) aData->mSlavePlane.size() );

            BELFEM_ASSERT( mProtoMesh != nullptr, "no proto mesh assigned to PeriodicityFactory" );

            Periodicity * aPeriodicity = new Periodicity( mMesh ) ;

            index_t tCount = 0 ;
            aPeriodicity->mMasterPlane.set_size( aData->mMasterPlane.size() );
            for ( id_t tID : aData->mMasterPlane )
            {
                aPeriodicity->mMasterPlane( tCount++ ) = mProtoMesh->node( tID );
            }

            this->create_hesse(
                aPeriodicity->mMasterPlane,
                aPeriodicity->master_hesse() );

            tCount = 0 ;
            aPeriodicity->mSlavePlane.set_size( aData->mSlavePlane.size() );
            for ( id_t tID : aData->mSlavePlane )
            {
                aPeriodicity->mSlavePlane( tCount++ ) = mProtoMesh->node( tID );
            }
            this->create_hesse(
                aPeriodicity->mSlavePlane,
                aPeriodicity->slave_hesse() );

            tCount = 0 ;
            aPeriodicity->mMasterFacets.set_size( aData->mMasterFacets.size() );
            for ( id_t tID : aData->mMasterFacets )
            {
                aPeriodicity->mMasterFacets( tCount++ ) = mProtoMesh->facet( tID );
            }

            tCount = 0 ;
            aPeriodicity->mSlaveFacets.set_size( aData->mSlaveFacets.size() );
            for ( id_t tID : aData->mSlaveFacets )
            {
                aPeriodicity->mSlaveFacets( tCount++ ) = mProtoMesh->facet( tID );
            }

            tCount = 0 ;
            aPeriodicity->mMasterNodes.set_size( aData->mMasterNodes.size() );
            for ( id_t tID : aData->mMasterNodes )
            {
                aPeriodicity->mMasterNodes( tCount++ ) = mProtoMesh->node( tID );
            }

            tCount = 0 ;
            aPeriodicity->mSlaveNodes.set_size( aData->mSlaveNodes.size() );
            for ( id_t tID : aData->mSlaveNodes )
            {
                aPeriodicity->mSlaveNodes( tCount++ ) = mProtoMesh->node( tID );
            }

            if ( mMesh->edges_exist() )
            {
                tCount = 0 ;
                aPeriodicity->mMasterEdges.set_size( aData->mMasterEdges.size() );
                for ( id_t tID : aData->mMasterEdges )
                {
                    aPeriodicity->mMasterEdges( tCount++ ) = mProtoMesh->edge( tID );
                }

                tCount = 0 ;
                aPeriodicity->mSlaveEdges.set_size( aData->mSlaveEdges.size() );
                for ( id_t tID : aData->mSlaveEdges )
                {
                    aPeriodicity->mSlaveEdges( tCount++ ) = mProtoMesh->edge( tID );
                }
            }

            if ( mMesh->faces_exist() )
            {
                tCount = 0 ;
                aPeriodicity->mMasterFaces.set_size( aData->mMasterFaces.size() );
                for ( id_t tID : aData->mMasterFaces )
                {
                    aPeriodicity->mMasterFaces( tCount++ ) = mProtoMesh->face( tID );
                }

                tCount = 0 ;
                aPeriodicity->mSlaveFaces.set_size( aData->mSlaveFaces.size() );
                for ( id_t tID : aData->mSlaveFaces )
                {
                    aPeriodicity->mSlaveFaces( tCount++ ) = mProtoMesh->face( tID );
                }
            }

            if ( aCrosslinkEntities ) this->crosslink( aPeriodicity );

            // a proto periodicity is received final: a later update() must not
            // geometrically re-match a mesh that contains coincident duplicates
            aPeriodicity->mNodePairsRestored = true ;

            return aPeriodicity ;
        }

        void
        PeriodicityFactory::update_periodicity( Periodicity * aPeriodicity )
        {
            BELFEM_ERROR( ! std::isnan( mTargetDistance ), "Slave points are not set" );
            BELFEM_ERROR( ! std::isnan( mSourceDistance ), "Master points are not set" );

            aPeriodicity->reset_nodes();
            aPeriodicity->reset_edges();
            aPeriodicity->reset_faces();
            aPeriodicity->reset_facets();

            mMesh->update_node_indices();
            mMesh->update_edge_indices() ;
            mMesh->update_face_indices();
            mMesh->update_facet_indices() ;

            BELFEM_ASSERT( mMesh->is_finalized() , "mesh is not finalized" );

            this->map_facets(
                aPeriodicity->master_facets(),
                aPeriodicity->slave_facets() );

            this->reset_bitsets( mMesh->number_of_nodes(),
                                 mMesh->number_of_edges(),
                                 mMesh->number_of_faces() );

            this->collect_nodes(
                aPeriodicity->master_facets(),
                aPeriodicity->master_nodes() );

            this->collect_nodes(
                aPeriodicity->slave_facets(),
                aPeriodicity->slave_nodes() );

            if ( aPeriodicity->has_node_pair_backup() )
            {
                aPeriodicity->restore_node_pairs();
                aPeriodicity->clear_node_pair_backup();
                aPeriodicity->flag( EntityType::NODE );
            }
            else
            {
                // once the backup has been consumed, the mesh contains coincident
                // duplicate nodes and a geometric re-match would silently mis-pair them
                BELFEM_ERROR( ! aPeriodicity->node_pairs_restored(),
                    "Periodicity::update() called again after the node pair backup was consumed" );

                this->collect_nodes(
                    aPeriodicity->master_facets(),
                    aPeriodicity->master_nodes() );

                this->collect_nodes(
                    aPeriodicity->slave_facets(),
                    aPeriodicity->slave_nodes() );

                if ( this->match_nodes(
                    aPeriodicity->master_facets(),
                    aPeriodicity->slave_facets(),
                    aPeriodicity->master_nodes(),
                    aPeriodicity->slave_nodes() ) )
                {
                    aPeriodicity->flag( EntityType::NODE );
                }
            }

            if ( mMesh->edges_exist() )
            {
                this->collect_edges(
                    aPeriodicity->master_facets(),
                    aPeriodicity->master_edges() );

                this->collect_edges(
                    aPeriodicity->slave_facets(),
                    aPeriodicity->slave_edges() );

                // a periodic plane that crosses an edge-carrying block must
                // yield edges; a silent zero-collection leaves the conductor
                // dofs unconstrained across the seam
                if ( aPeriodicity->master_edges().size() == 0 )
                {
                    for ( Facet * tFacet : aPeriodicity->master_facets() )
                    {
                        BELFEM_ERROR( ! tFacet->has_master() || ! tFacet->master()->has_edges(),
                            "periodic plane intersects edge-carrying elements, but no periodic edges were collected ( facet %lu )",
                            ( long unsigned int ) tFacet->id() );
                    }
                }

                if ( this->match_edges(
                    aPeriodicity->master_facets(),
                    aPeriodicity->slave_facets(),
                    aPeriodicity->master_nodes(),
                    aPeriodicity->master_edges(),
                    aPeriodicity->slave_nodes(),
                    aPeriodicity->slave_edges() ) )
                {
                    aPeriodicity->flag( EntityType::EDGE );
                }

                mMesh->compute_edge_directions();
            }

            // also matches faces if they exist
            if ( this->match_facets_and_faces(
                aPeriodicity->master_nodes(),
                aPeriodicity->master_facets(),
                aPeriodicity->master_faces(),
                aPeriodicity->slave_nodes(),
                aPeriodicity->slave_facets(),
                aPeriodicity->slave_faces() ) )
            {
                aPeriodicity->flag( EntityType::FACET );

                if ( aPeriodicity->master_faces().size() > 0 )
                {
                    aPeriodicity->flag( EntityType::FACE );
                }
            }

            mMesh->update_node_indices();
            mMesh->update_edge_indices() ;
            mMesh->update_face_indices();
            mMesh->update_facet_indices() ;

            if ( aPeriodicity->is_flagged( EntityType::FACE ) )
            {
                this->fix_face_slaves( aPeriodicity );
            }

        }

        void
        PeriodicityFactory::map_facets(
                Cell< Facet * >   & aSourceFacets,
                Cell< Facet * >   & aTargetFacets )
        {
            Cell< SideSet * > tSourceSideSets ;
            this->select_sidesets(
                mSourceTransform,
                mSourceDistance,
                tSourceSideSets );

            Cell< Node * > tSourceNodes ;
            this->create_temporary_nodes(
                mSourceTransform,
                mSourceDistance,
                tSourceSideSets,
                tSourceNodes );

            Cell< SideSet * > tTargetSideSets ;
            this->select_sidesets(
                mTargetTransform,
                mTargetDistance,
                tTargetSideSets );

            Cell< Node * > tTargetNodes ;
            this->create_temporary_nodes(
                mTargetTransform,
                mTargetDistance,
                tTargetSideSets,
                tTargetNodes );

            Map< id_t, Facet * > tSourceMap ;
            for ( SideSet * tSideSet : tSourceSideSets )
            {
                for ( Facet * tFacet : tSideSet->facets() )
                {
                    tSourceMap[ tFacet->id() ] = tFacet ;
                }
            }
            Map< id_t, Facet * > tTargetMap ;
            for ( SideSet * tSideSet : tTargetSideSets )
            {
                for ( Facet * tFacet : tSideSet->facets() )
                {
                    tTargetMap[ tFacet->id() ] = tFacet ;
                }
            }
            BELFEM_ERROR( tSourceNodes.size() == tTargetNodes.size(), "Number of facets does not match" );

            // creating the tree
            Node * tRoot = this->create_kdtree( tTargetNodes );

            // now we rearrange the nodes in the target container
            index_t tCount = 0 ;
            for ( Node * tNode : tTargetNodes )
            {
                tNode->set_index( tCount++ );
            }

            for ( Node * tSource : tSourceNodes )
            {
                real tBestDist = BELFEM_REAL_MAX ;
                Node * tTarget = this->find_closest_node(
                    tRoot, tSource->x(), tSource->y(), false, nullptr, tBestDist ) ;

                // tBestDist is a plain distance, so compare against the mesh tolerance directly
                BELFEM_ERROR( tBestDist < BELFEM_MESH_EPSILON ,
                    "Could not find a corresponding periodic facet for facet %lu (master: %lu), closest candidate: %lu (master: %lu). This does not seem to be a proper periodic mesh.",
                    ( long unsigned int ) tSource->id(),
                    ( long unsigned int ) mMesh->facet( tSource->id())->master()->id(),
                    ( long unsigned int ) tTarget->id(),
                    ( long unsigned int ) mMesh->facet( tTarget->id())->master()->id()
                    );

                tSource->set_index( tTarget->index() );
            }
            sort( tSourceNodes.begin(), tSourceNodes.end(), []( Node * a, Node * b ) { return a->index() < b->index(); } );

            // with the facet mapping in place, we can now set the periodicities
            aSourceFacets.set_size( tCount, nullptr );
            aTargetFacets.set_size( tCount, nullptr );

            for ( index_t k=0; k<tCount; ++k )
            {
                Facet * A = tSourceMap( tSourceNodes( k )->id() );
                Facet * B = tTargetMap( tTargetNodes( k )->id() );
                A->set_index( k );
                B->set_index( k );
                aSourceFacets( k ) = A ;
                aTargetFacets( k ) = B ;
            }

            // the temporary nodes have done their job, we can delete them now
            for ( Node * tNode : tSourceNodes )
            {
                delete tNode ;
            }
            for ( Node * tNode : tTargetNodes )
            {
                delete tNode ;
            }
        }

        void
        PeriodicityFactory::collect_nodes(
               Cell< Facet * >   & aFacets,
               Cell< Node * >    & aNodes )
        {
            BELFEM_ASSERT( mMesh->is_finalized(), "Mesh is not finalized" );

            mNodeBitset->reset();

            for ( Facet * tFacet : aFacets )
            {
                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    mNodeBitset->set( tFacet->node( k )->original()->index() );
                }
            }
            Cell< index_t > tIndices ;
            mNodeBitset->where( tIndices );

            Cell< Node * > & tNodes = mMesh->nodes();

            aNodes.set_size( tIndices.size(), nullptr );
            index_t tCount = 0 ;
            for ( index_t tIndex : tIndices )
            {
                aNodes( tCount++ ) = tNodes( tIndex );
            }
        }


        void
        PeriodicityFactory::collect_edges(
                    Cell< Facet * >   & aFacets,
                    Cell< Edge * >    & aEdges )
        {
            if ( ! mMesh->edges_exist() )
            {
                aEdges.clear() ;
                return ;
            }

            BELFEM_ASSERT( mMesh->edges_are_finalized(), "Mesh edges are not finalized" );
            mEdgeBitset->reset() ;
            Cell< Edge * > tMasterEdges ;
            for ( Facet * tFacet : aFacets )
            {
                Element * tElement = tFacet->element();

                if ( tElement->has_edges() )
                {
                    for ( uint e=0; e<tElement->number_of_edges(); ++e )
                    {
                        mEdgeBitset->set( tElement->edge( e )->index() );
                    }
                }
                else if ( tFacet->has_master() && tFacet->master()->has_edges() )
                {
                    // bulk conductor caps: the surface wrapper carries no edge
                    // container, so we reconstruct from the master element.
                    // enumeration and intrinsic direction may differ from the
                    // partner side; match_edges() pairs by endpoint originals
                    // and aligns direction by node swap, so both are safe here
                    tFacet->master()->get_edges_of_facet( tFacet->index_on_master(), tMasterEdges );
                    for ( Edge * tEdge : tMasterEdges )
                    {
                        mEdgeBitset->set( tEdge->index() );
                    }
                }
            }

            Cell< index_t > tIndices ;
            mEdgeBitset->where( tIndices );

            aEdges.set_size( tIndices.size(), nullptr );
            index_t tCount = 0 ;
            Cell< Edge * > & tEdges = mMesh->edges();
            for ( index_t tIndex : tIndices )
            {
                aEdges( tCount++ ) = tEdges( tIndex );
            }
        }

        bool
        PeriodicityFactory::match_nodes(
                   Cell< Facet * >   & aSourceFacets,
                   Cell< Facet * >   & aTargetFacets,
                   Cell< Node * >    & aSourceNodes,
                   Cell< Node * >    & aTargetNodes )
        {

            if ( aSourceNodes.size() == 0 )
            {
                aTargetNodes.clear() ;
                return false ;
            }

            BELFEM_ASSERT( aSourceFacets.size() == aTargetFacets.size(),
                "number of nodes on source and target facets do not match ( %lu vs %lu )",
                ( long unsigned int ) aSourceFacets.size(),
                ( long unsigned int ) aTargetFacets.size() );


            for ( Node * tNode : aSourceNodes )
            {
                tNode->unflag();
                tNode->set_index( gNoIndex );
            }
            for ( Node * tNode : aTargetNodes )
            {
                tNode->unflag();
                tNode->set_index( gNoIndex );
            }
            BELFEM_ASSERT( aSourceNodes.size() == aTargetNodes.size(),
                            "number of nodes on source and target nodes do not match ( %lu vs %lu ) ",
                            ( long unsigned int ) aSourceNodes.size(),
                            ( long unsigned int ) aTargetNodes.size() );

            index_t n = aTargetNodes.size();
            for ( index_t k=0; k<n; ++k )
            {
                aTargetNodes( k ) = nullptr ;
            }


            // now we check the mapping
            Vector< real > Xp( 3 );
            Vector< real > Xq( 3 );
            Vector< real > Yp( 3 );
            Vector< real > Yq( 3 );
            Vector< real > D( 3 );

            index_t m = aSourceFacets.size();

            index_t tCount = 0 ;

            for ( index_t f=0; f<m; ++f )
            {
                Facet * A = aSourceFacets( f );
                Facet * B = aTargetFacets( f );

                for ( uint i=0; i<A->number_of_nodes(); ++i )
                {
                    Node * P = A->node( i );

                    if ( P->is_flagged() ) continue ;

                    P->get_coords( Xp );
                    Yp = mSourceTransform * Xp ;

                    for ( uint j=0; j<B->number_of_nodes(); ++j )
                    {
                        Node * Q = B->node( j );
                        if ( Q->is_flagged() ) continue ;

                        Q->get_coords( Xq );
                        Yq = mTargetTransform * Xq ;

                        D = Yq - Yp ;
                        D(2) = 0.0 ;

                        if ( norm(D) < BELFEM_MESH_EPSILON )
                        {
                            P->set_periodic( Q );
                            Q->set_periodic( P );
                            P->flag();
                            Q->flag();
                            P->set_index( tCount );
                            Q->set_index( tCount );
                            aTargetNodes( tCount++ ) = Q ;
                        }
                    }
                }
            }

            BELFEM_ASSERT( tCount == n, "Failed to match all nodes" );

            sort( aSourceNodes, opNodeIndex );

            return true ;
        }

        void
        PeriodicityFactory::create_edge_map( Cell< Node * > & aNodes, Cell< Edge * > & aEdges, Map< key_t, Edge * > & aEdgeMap )
        {
            key_t N = aNodes.size();
            index_t tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNode->set_index( tCount++ );
            }

            for ( Edge * tEdge : aEdges )
            {
                key_t A = tEdge->node(0)->original()->index();
                key_t B = tEdge->node(1)->original()->index();

                // integrity: original-identity keying must not produce a self-loop
                // ( both endpoints collapse to the same original ) nor collapse two
                // distinct periodic edges onto one key ( would silently mis-pair the
                // periodic DOFs ). These indicate an unhealthy mesh; fail loudly.
                BELFEM_ERROR( A != B,
                    "degenerate periodic edge ( index %lu ): both endpoints collapse to original node-index %lu",
                    ( long unsigned int ) tEdge->index(), ( long unsigned int ) A );

                key_t tKey = B < A ? A * N + B : B * N + A ;

                BELFEM_ERROR( ! aEdgeMap.key_exists( tKey ) || aEdgeMap( tKey ) == tEdge,
                    "two distinct periodic edges ( indices %lu and %lu ) collapse to one original-node-index key - unhealthy mesh",
                    ( long unsigned int ) tEdge->index(),
                    ( long unsigned int ) ( aEdgeMap.key_exists( tKey ) ? aEdgeMap( tKey )->index() : 0 ) );

                aEdgeMap[ tKey ] = tEdge ;
            }
        }

        void
        PeriodicityFactory::create_facet_map( Cell< Node * > & aNodes, Cell< Facet * > & aFacets, Map< key128_t, Facet * > & aFacetMap )
        {
#if !defined( NDEBUG ) || defined( DEBUG )
            // testing the integrity of the node and facet data
            for ( Facet * tFacet : aFacets )
            {
                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    tFacet->node( k )->original()->unflag( 2 );
                }
            }
            for ( Node * tNode : aNodes )
            {
                tNode->flag( 2 );
            }
            for ( Facet * tFacet : aFacets )
            {
                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    BELFEM_ERROR( tFacet->node( k )->original()->is_flagged( 2 ), "Node %lu is not flagged", ( long unsigned int ) tFacet->node( k )->id() );
                }
            }
            for ( Facet * tFacet : aFacets )
            {
                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    tFacet->node( k )->original()->unflag( 2 );
                }
            }
            for ( Node * tNode : aNodes )
            {
                tNode->unflag( 2 );
            }
#endif

            index_t tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNode->set_index( tCount++ );
            }

            key128_t n = aNodes.size();
            aFacetMap.clear();

            if ( mMesh->number_of_dimensions() == 2 )
            {
                for ( Facet * tFacet : aFacets )
                {

                    key128_t a = tFacet->node( 0 )->original()->index();
                    key128_t b = tFacet->node( 1 )->original()->index();
                    key128_t tKey = a > b ? a * n + b : b * n + a ;
                    aFacetMap[ tKey ] = tFacet ;
                }
            }
            else if ( mMesh->number_of_dimensions() == 3 )
            {
                Cell< key128_t > tIndices ;
                for ( Facet * tFacet : aFacets )
                {
                    uint m = tFacet->number_of_corner_nodes();
                    tIndices.set_size( m );

                    for ( uint k=0; k<m; ++k )
                    {
                        tIndices( k ) = tFacet->node( k )->original()->index();
                    }

                    sort( tIndices );
                    key128_t tKey = ( tIndices( 2 ) * n + tIndices( 1 ) ) * n + tIndices( 0 );
                    aFacetMap[ tKey ] = tFacet ;
                }
            }
        }

        // Facet-mediated periodic edge matching ( phase two of direction (d),
        // todo/periodic_cap_cut_emission.md ). map_facets() aligns the two
        // facet lists index-wise, and each cap facet lies in exactly one cut
        // sector per side, so the k-th facet pair carries the sector-resolved
        // raw edge correspondence. Within one pair, edges match through their
        // ORIGINAL endpoints ( unique per facet ); the raw endpoints of the
        // matched partner are whatever the cut relinking chose.
        //   - STRAIGHT profile ( partner raw endpoints are exactly the node
        //     ties of the source raw endpoints ): tied directly.
        //   - PURE MIXED profile ( original <-> duplicate across the
        //     identification, both endpoints sharing the cut membership on
        //     each side ): tied DIRECTLY — the abstract current cancels in
        //     the Whitney circulation, the jump stays in the node hanging
        //     ( Grok EDGE-TIES combo III ).
        //   - HALF-CUT profile ( exactly one duplicated endpoint ): tied
        //     purely IFF both sides are half-cut and all four raw endpoints
        //     are periodically paired endpointwise ( then the identified H
        //     circulations are equal exactly; applied deferred, after the
        //     loop, so pure ties win under every encounter order — the corc
        //     cap-corner fix, 2026-08-25 ). Anything else is the genuinely
        //     affine class and stays untied.
        // Afterwards the edge containers are compacted to tied pairs only,
        // because every downstream consumer ( set_entity_dependencies, BFM
        // save, crosslink ) assumes matched pairs.

        bool
        PeriodicityFactory::match_edges(
            Cell< Facet * > & aSourceFacets,
            Cell< Facet * > & aTargetFacets,
            Cell< Node * >  & aSourceNodes,
            Cell< Edge * >  & aSourceEdges,
            Cell< Node * >  & aTargetNodes,
            Cell< Edge * >  & aTargetEdges )
        {
            if ( ! mMesh->edges_exist() ) return false ;

            // the folded originals of the two surfaces must still correspond
            BELFEM_ERROR( aSourceNodes.size() == aTargetNodes.size(),
                "number of nodes on source and target periodic surfaces do not match (%lu vs. %lu)",
                ( long unsigned int ) aSourceNodes.size(), ( long unsigned int ) aTargetNodes.size() );

            BELFEM_ERROR( aSourceFacets.size() == aTargetFacets.size(),
                "number of facets on source and target periodic surfaces do not match (%lu vs. %lu)",
                ( long unsigned int ) aSourceFacets.size(), ( long unsigned int ) aTargetFacets.size() );

            if ( aTargetEdges.size() == 0 && aSourceEdges.size() == 0 ) return false ;

            // per-facet edge collection, same two paths as collect_edges()
            auto tEdgesOfFacet = []( Facet * aFacet, Cell< Edge * > & aOut )
            {
                Element * tElement = aFacet->element() ;
                if ( tElement->has_edges() )
                {
                    aOut.set_size( tElement->number_of_edges(), nullptr );
                    for ( uint e = 0; e < tElement->number_of_edges(); ++e )
                    {
                        aOut( e ) = tElement->edge( e );
                    }
                }
                else if ( aFacet->has_master() && aFacet->master()->has_edges() )
                {
                    aFacet->master()->get_edges_of_facet( aFacet->index_on_master(), aOut );
                }
                else
                {
                    aOut.clear();
                }
            };

            auto tKeyOfIds = []( const luint aA, const luint aB ) -> luint
            {
                return aA < aB ? ( aA << 32 ) | aB : ( aB << 32 ) | aA ;
            };

            // a skipped trace twin must be Whitney-identical to the kept tie:
            // same original endpoint pair and NOT half-cut ( both endpoints
            // share the cut membership, so the abstract current cancels in
            // the circulation ). Anything else must fail loudly ( Codex
            // phase-two final audit, alias invariant )
            auto tIsPureTwin = [&tKeyOfIds]( Edge * aEdge, Edge * aKept ) -> bool
            {
                const bool tHalf =
                    ( aEdge->node( 0 ) != aEdge->node( 0 )->original() ) !=
                    ( aEdge->node( 1 ) != aEdge->node( 1 )->original() );
                if ( tHalf ) return false ;
                return tKeyOfIds(
                        ( luint ) aEdge->node( 0 )->original()->id(),
                        ( luint ) aEdge->node( 1 )->original()->id() )
                    == tKeyOfIds(
                        ( luint ) aKept->node( 0 )->original()->id(),
                        ( luint ) aKept->node( 1 )->original()->id() );
            };

            index_t tNumTied           = 0 ;
            index_t tNumMixed          = 0 ; // pure-mixed pairs, direct-tied
            index_t tNumHalfCut        = 0 ; // affine class, stays untied
            index_t tNumHalfCutTied    = 0 ; // pure half-cut pairs, tied
            index_t tNumAlias          = 0 ; // trace twins ( 1:2 realization )

            // deferred pure half-cut tie candidates: collected during the
            // loop, applied afterwards so that pure ties and aliases always
            // win regardless of encounter order ( pure-twin-wins )
            struct HalfCutCandidate { Edge * E ; Edge * F ; bool Swap ; };
            std::vector< HalfCutCandidate > tPureHalfCutCandidates ;
            index_t tNumFacetMismatch  = 0 ;

            // proposal bookkeeping: agreement check for edges shared by two
            // facets, and target coverage
            std::map< Edge *, Edge * > tProposals ;
            std::set< Edge * >         tProposedTargets ;

            Cell< Edge * > tEm ;
            Cell< Edge * > tEs ;

            for ( index_t k = 0; k < aSourceFacets.size(); ++k )
            {
                tEdgesOfFacet( aSourceFacets( k ), tEm );
                tEdgesOfFacet( aTargetFacets( k ), tEs );

                if ( tEm.size() != tEs.size() )
                {
                    ++tNumFacetMismatch ;
                    continue ;
                }

                for ( Edge * E : tEm )
                {
                    Node * tOA = E->node( 0 )->original() ;
                    Node * tOB = E->node( 1 )->original() ;

                    // originals on the surfaces are always tied
                    if ( ! tOA->is_periodic() || ! tOB->is_periodic() )
                    {
                        ++tNumFacetMismatch ;
                        continue ;
                    }

                    const luint tWant = tKeyOfIds(
                        ( luint ) tOA->periodic()->id(),
                        ( luint ) tOB->periodic()->id() );

                    Edge * F = nullptr ;
                    for ( Edge * tCandidate : tEs )
                    {
                        if ( tKeyOfIds(
                                ( luint ) tCandidate->node( 0 )->original()->id(),
                                ( luint ) tCandidate->node( 1 )->original()->id() ) == tWant )
                        {
                            F = tCandidate ;
                            break ;
                        }
                    }

                    // a null match most likely means a duplicate that never
                    // received set_original() — its original() returns itself
                    // and the key cannot match
                    if ( F == nullptr )
                    {
                        ++tNumFacetMismatch ;
                        continue ;
                    }

                    // an edge shared by two facets is proposed twice: the two
                    // sectors must agree on the partner
                    auto tHit = tProposals.find( E );
                    if ( tHit != tProposals.end() )
                    {
                        tProposedTargets.insert( F );

                        if ( tHit->second == F )
                        {
                            continue ;
                        }

                        // 1:2 realization along the cut trace. If this source
                        // edge already carries a tie, the newcomer must be a
                        // Whitney-identical pure twin ( benign alias ) or the
                        // half-cut twin ( circulation differs by I; belongs to
                        // the censused affine class ). If NO tie exists yet
                        // ( the first proposal was the half-cut twin ), fall
                        // through so the pure twin can still be tied.
                        if ( E->is_periodic() )
                        {
                            if ( tIsPureTwin( F, E->periodic() ) )
                            {
                                ++tNumAlias ;
                            }
                            else
                            {
                                const bool tHalfTwin =
                                    ( F->node( 0 ) != F->node( 0 )->original() ) !=
                                    ( F->node( 1 ) != F->node( 1 )->original() );
                                BELFEM_ERROR( tHalfTwin,
                                    "periodic trace-twin alias is neither a pure twin nor half-cut ( target edge %lu vs %lu )",
                                    ( long unsigned int ) F->id(),
                                    ( long unsigned int ) E->periodic()->id() );
                                ++tNumHalfCut ;
                            }
                            continue ;
                        }
                        // no tie yet: retry classification with the newcomer
                        tHit->second = F ;
                    }
                    else
                    {
                        tProposals[ E ] = F ;
                        tProposedTargets.insert( F );
                    }

                    // classify the raw profile ( Grok EDGE-TIES, combo III:
                    // straight node ties + DIRECT edge ties win ). For an edge
                    // whose endpoints are BOTH cut duplicates ( or both
                    // originals ) the abstract current cancels in the Whitney
                    // relation h = phi0 - phi1, so the direct period tie is
                    // exact — the jump stays in the node hanging. HALF-CUT
                    // edges ( exactly one duplicated endpoint on either
                    // side ) split: with all four raw endpoints paired
                    // endpointwise they tie purely too ( deferred, below );
                    // otherwise they are genuinely affine and stay untied,
                    // their closure running through node-source hanging.
                    Node * tA = E->node( 0 );
                    Node * tB = E->node( 1 );

                    const bool tHalfCutSource =
                        ( tA != tA->original() ) != ( tB != tB->original() );
                    const bool tHalfCutTarget =
                        ( F->node( 0 ) != F->node( 0 )->original() ) !=
                        ( F->node( 1 ) != F->node( 1 )->original() );

                    if ( ! tHalfCutSource && ! tHalfCutTarget )
                    {
                        // source-side trace twin: the target edge is already
                        // claimed by a Whitney-identical twin of this source
                        // edge ( the mirror image of the target-alias case ).
                        // Re-tying would overwrite the partner AND re-swap
                        // the direction alignment of the first tie. A pure
                        // twin is a benign alias; a half-cut twin belongs to
                        // the censused affine class ( note: E cannot reach
                        // here half-cut — the branch above filters those )
                        if ( F->is_periodic() )
                        {
                            BELFEM_ERROR( tIsPureTwin( E, F->periodic() ),
                                "periodic trace-twin alias is not Whitney-identical to the kept tie ( source edge %lu vs %lu )",
                                ( long unsigned int ) E->id(),
                                ( long unsigned int ) F->periodic()->id() );
                            ++tNumAlias ;
                            continue ;
                        }

                        // align the intrinsic direction through the ORIGINALS
                        // ( valid for straight and pure-mixed pairs alike ):
                        // F->node(0) must correspond to E->node(0)
                        if ( F->node( 0 )->original() == tOB->periodic() &&
                             F->node( 1 )->original() == tOA->periodic() )
                        {
                            BELFEM_ASSERT( F->number_of_nodes() < 4, "Node swapping is not implemented for higher-order edges" );

                            Node * P = F->node( 0 );
                            Node * Q = F->node( 1 );
                            F->insert_node( Q, 0 );
                            F->insert_node( P, 1 );
                        }
#if !defined( NDEBUG ) || defined( DEBUG )
                        else
                        {
                            BELFEM_ERROR(
                                F->node( 0 )->original() == tOA->periodic() &&
                                F->node( 1 )->original() == tOB->periodic(),
                                "Topology error in facet-mediated edge pair" );
                        }
#endif
                        E->set_periodic( F );
                        F->set_periodic( E );
                        ++tNumTied ;

                        // bookkeeping: how many ties are pure-mixed ( at
                        // least one side runs through duplicates )
                        const bool tDupSource =
                            tA != tA->original() || tB != tB->original() ;
                        const bool tDupTarget =
                            F->node( 0 ) != F->node( 0 )->original() ||
                            F->node( 1 ) != F->node( 1 )->original() ;
                        if ( tDupSource || tDupTarget )
                        {
                            ++tNumMixed ;
                        }
                    }
                    else
                    {
                        // half-cut profile: exactly one duplicated endpoint on
                        // a side. The affine caution ( abstract current
                        // surviving on the edge dof ) applies only when an
                        // endpoint carrier is NOT periodically resolved. When
                        // BOTH sides are half-cut and ALL FOUR raw endpoints
                        // are periodically paired — endpointwise, after
                        // alignment — the identified H circulations are equal
                        // exactly and the tie is pure: leaving such a pair
                        // untied gives each geometric cap an independent
                        // corner dof ( the corc cap-corner antisymmetry,
                        // 2026-08-25 ). Genuinely affine half-cuts stay
                        // untied as before.
                        bool tTiePure = false ;

                        bool tSwap = false ;

                        if ( tHalfCutSource && tHalfCutTarget &&
                             tA->is_periodic() && tB->is_periodic() &&
                             F->node( 0 )->is_periodic() &&
                             F->node( 1 )->is_periodic() )
                        {
                            // alignment decision through the originals ( same
                            // rule as the pure ties above ) — decided WITHOUT
                            // mutating F; mutation happens only if the
                            // deferred pass actually ties the pair
                            tSwap =
                                F->node( 0 )->original() == tOB->periodic() &&
                                F->node( 1 )->original() == tOA->periodic() ;

                            Node * tF0 = tSwap ? F->node( 1 ) : F->node( 0 );
                            Node * tF1 = tSwap ? F->node( 0 ) : F->node( 1 );

                            // endpointwise raw-carrier test: each endpoint's
                            // periodic partner must be the CORRESPONDING
                            // endpoint of F — anything else, e.g. a crossed
                            // (orig,dup)x(dup,orig) pair, is affine
                            tTiePure =
                                tA->periodic() == tF0 &&
                                tB->periodic() == tF1 ;
                        }

                        if ( tTiePure )
                        {
                            // DEFER the tie: pure ties and aliases of the
                            // main loop must win regardless of encounter
                            // order ( pure-twin-wins, order-independent ).
                            // Candidates are applied after the loop, only
                            // if neither edge has been claimed by then.
                            tPureHalfCutCandidates.push_back( { E, F, tSwap } );
                        }
                        else
                        {
                            // half-cut: the abstract current survives on the
                            // edge DOF ( affine relation ); stays untied —
                            // closure via node-source hanging where DOFs exist
                            ++tNumHalfCut ;
                        }
                    }
                }
            }

            // apply the deferred pure half-cut ties: only pairs that no tie
            // from the main loop has claimed in the meantime — pure ties
            // therefore always win, independent of encounter order
            for ( const HalfCutCandidate & tC : tPureHalfCutCandidates )
            {
                Edge * E = tC.E ;
                Edge * F = tC.F ;

                if ( E->is_periodic() || F->is_periodic() )
                {
                    // displaced by a main-loop tie: stays untied ( affine )
                    ++tNumHalfCut ;
                    continue ;
                }

                if ( tC.Swap )
                {
                    BELFEM_ASSERT( F->number_of_nodes() < 4,
                        "Node swapping is not implemented for higher-order edges" );

                    Node * P = F->node( 0 );
                    Node * Q = F->node( 1 );
                    F->insert_node( Q, 0 );
                    F->insert_node( P, 1 );
                }

                E->set_periodic( F );
                F->set_periodic( E );
                ++tNumTied ;
                ++tNumHalfCutTied ;
            }

            // coverage: every collected edge must appear in the proposals
            const index_t tNumSourceUncovered =
                aSourceEdges.size() - tProposals.size() ;
            const index_t tNumTargetUncovered =
                aTargetEdges.size() - tProposedTargets.size() ;

            gLog.message( InfoLevel::Detailed,
                "    periodic edges : %lu tied ( %lu through cut duplicates, %lu pure half-cut ), %lu half-cut untied, %lu trace-twin aliases, %lu facet mismatches",
                ( long unsigned int ) tNumTied,
                ( long unsigned int ) tNumMixed,
                ( long unsigned int ) tNumHalfCutTied,
                ( long unsigned int ) tNumHalfCut,
                ( long unsigned int ) tNumAlias,
                ( long unsigned int ) tNumFacetMismatch );

            // facet mismatches and uncovered edges are hard inconsistencies
            BELFEM_ERROR( tNumFacetMismatch == 0
                          && tNumSourceUncovered == 0
                          && tNumTargetUncovered == 0,
                "periodic facet-mediated edge matching inconsistent: %lu facet mismatches, %lu + %lu uncovered edges",
                ( long unsigned int ) tNumFacetMismatch,
                ( long unsigned int ) tNumSourceUncovered,
                ( long unsigned int ) tNumTargetUncovered );

            // compact the containers to tied pairs only ( downstream consumes
            // matched pairs: set_entity_dependencies, BFM save, crosslink )
            {
                Cell< Edge * > tTiedSource( tNumTied, nullptr );
                Cell< Edge * > tTiedTarget( tNumTied, nullptr );
                index_t tCount = 0 ;
                for ( Edge * E : aSourceEdges )
                {
                    if ( E->is_periodic() )
                    {
                        tTiedSource( tCount ) = E ;
                        tTiedTarget( tCount ) = E->periodic() ;
                        ++tCount ;
                    }
                }
                BELFEM_ERROR( tCount == tNumTied,
                    "periodic edge compaction mismatch ( %lu vs %lu )",
                    ( long unsigned int ) tCount,
                    ( long unsigned int ) tNumTied );
                aSourceEdges = tTiedSource ;
                aTargetEdges = tTiedTarget ;
            }

            return tNumTied > 0 ;
        }

        bool
        PeriodicityFactory::match_facets_and_faces(
                Cell< Node * >  & aSourceNodes,
                Cell< Facet * > & aSourceFacets,
                Cell< Face * >  & aSourceFaces,
                Cell< Node * >  & aTargetNodes,
                Cell< Facet * > & aTargetFacets,
                Cell< Face * >  & aTargetFaces )
        {

            BELFEM_ASSERT( aSourceFacets.size() == aTargetFacets.size(), "number of nodes on source and target facets do not match" );
            BELFEM_ASSERT( aSourceNodes.size() == aTargetNodes.size(), "number of nodes on source and target nodes do not match" );

            if ( aSourceFacets.size() == 0 || aSourceNodes.size() == 0 ) return false ;

            Map< key128_t, Facet * > tSources ;
            this->create_facet_map( aSourceNodes, aSourceFacets, tSources );

            Map< key128_t, Facet * > tTargets ;
            this->create_facet_map( aTargetNodes, aTargetFacets, tTargets );

            for ( auto tPair : tSources )
            {
                Facet * A = tPair.second ;
                Facet * B = tTargets( tPair.first ) ;

                A->set_periodic( B );
                B->set_periodic( A );
            }


            aSourceFaces.clear();
            aTargetFaces.clear();

            if ( mMesh->faces_exist() )
            {
                aSourceFaces.reserve( aSourceFacets.size() );
                aTargetFaces.reserve( aTargetFacets.size() );

                for ( Facet * A : aSourceFacets )
                {
                    Facet * B = A->periodic();

                    if ( A->master()->has_faces() && B->master()->has_faces() )
                    {
                        Face * a = A->master()->face( A->index_on_master() );
                        Face * b = B->master()->face( B->index_on_master() );

                        a->set_periodic( b );
                        b->set_periodic( a );

                        aSourceFaces.push( a );
                        aTargetFaces.push( b );
                    }

                }
                aSourceFaces.shrink_to_fit();
                aTargetFaces.shrink_to_fit();
            }

            return true ;
        }

        void
        PeriodicityFactory::compute_transformation_matrix(
            const Node * A,
            const Node * B,
            const Node * C,
                  real & D,
            Matrix< real > & T,
            Vector< real > & H )
        {
            // Build orthonormal basis from three non-collinear points:
            //   P = normalized(B - A)          in-plane direction 1
            //   N = normalized(cross(P, C-A))  surface normal
            //   Q = cross(N, P)                in-plane direction 2
            //
            // T = [P; Q; N] so that Y = T*X gives:
            //   Y(0) = dot(P, X)   in-plane coordinate 1
            //   Y(1) = dot(Q, X)   in-plane coordinate 2
            //   Y(2) = dot(N, X)   signed distance from origin along normal
            //
            // D = dot(A, N) is the plane offset: dot(X, N) = D for all
            // points X on the plane.
            // H = hesse normal form of the plane

            Vector< real > R( 3 );
            R( 0 ) = A->x();
            R( 1 ) = A->y();
            R( 2 ) = A->z();

            Vector< real > P( 3 );
            P( 0 ) = B->x();
            P( 1 ) = B->y();
            P( 2 ) = B->z();
            P-= R ;

            Vector< real > Q( 3 );
            Q( 0 ) = C->x();
            Q( 1 ) = C->y();
            Q( 2 ) = C->z();

            Q-= R;

            Vector< real > N( 3 );
            N = cross( P, Q );

            real tNormN = norm( N );
            BELFEM_ERROR( std::abs( tNormN ) > BELFEM_MESH_EPSILON, "could not determine surface normal" );


            N /= tNormN;

            D = dot( R, N );

            P/= norm( P );
            Q = cross( N, P );
            Q /= norm( Q );

            T.set_size( 3, 3 );

            T( 0, 0 ) = P( 0 );
            T( 1, 0 ) = Q( 0 );
            T( 2, 0 ) = N( 0 );

            T( 0, 1 ) = P( 1 );
            T( 1, 1 ) = Q( 1 );
            T( 2, 1 ) = N( 1 );

            T( 0, 2 ) = P( 2 );
            T( 1, 2 ) = Q( 2 );
            T( 2, 2 ) = N( 2 );

            H.set_size( 4 );
            H( 0 ) = N( 0 );
            H( 1 ) = N( 1 );
            H( 2 ) = N( 2 );
            H( 3 ) = D ;
        }

        void
        PeriodicityFactory::select_sidesets(
            const Matrix< real > & aTransform,
            const real aDistance,
             Cell< SideSet * > & aSideSets )
        {

            Vector< real > X( 3 );
            Vector< real > N( 3 );

            // extract the surface normal (row 2 of the transform matrix)
            N( 0 ) = aTransform( 2, 0 );
            N( 1 ) = aTransform( 2, 1 );
            N( 2 ) = aTransform( 2, 2 );

            // A sideset lies on the plane if ALL its nodes satisfy
            // dot(X, N) = D within tolerance. Collect matching sidesets.
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                Cell< Node * > & tNodes = tSideSet->nodes() ;

                index_t tCount = 0 ;
                for ( Node * tNode : tNodes )
                {
                    X( 0 ) = tNode->x();
                    X( 1 ) = tNode->y();
                    X( 2 ) = tNode->z();

                    if ( std::abs( dot( X, N ) - aDistance ) < BELFEM_MESH_EPSILON )
                    {
                        ++tCount ;
                    }
                    else
                    {
                        break ;
                    }
                }

                if ( tCount == tNodes.size() )
                {
                    aSideSets.push( tSideSet ) ;
                }
            }

        }

        void
        PeriodicityFactory::create_temporary_nodes(
                const Matrix< real >    & aTransform,
                const real aDistance,
                const Cell< SideSet * > & aSideSets,
                Cell< Node * > & aNodes )
        {
            index_t tCount = 0 ;
            Vector< real > X( 3 );
            Vector< real > M( 3 );

            for ( SideSet * tSideSet : aSideSets )
            {
                tCount += tSideSet->number_of_facets();
            }
            aNodes.set_size( tCount, nullptr );

            tCount = 0 ;
            for ( SideSet * tSideSet : aSideSets )
            {
                Cell< Facet * > & tFacets = tSideSet->facets() ;
                for ( Facet * tFacet : tFacets )
                {
                    M.fill( 0 );
                    for ( uint k=0; k<tFacet->number_of_corner_nodes(); ++k )
                    {
                        tFacet->node( k )->get_coords( X );
                        M += X ;
                    }
                    M/= tFacet->number_of_corner_nodes();

                    X = aTransform * M ;
                    X( 2 ) -= aDistance ; //<-- should be close to zero
                    if ( std::abs( X( 2 ) ) < BELFEM_MESH_EPSILON ) X( 2 ) = 0.0;

                    Node * tNode = new Node( tFacet->id(), X( 0 ), X( 1 ), X( 2 ) );
                    tNode->set_index( tCount );
                    aNodes( tCount++ ) = tNode ;
                }
            }
        }

        Node *
        PeriodicityFactory::find_closest_node(
            Node   * aNode,
            real     aX,
            real     aY,
            bool     aFlip,
            Node   * aBest,
            real   & aBestDist )
        {
            if ( aNode == nullptr )
            {
                return aBest ;
            }

            // update best candidate using the Euclidean distance
            real tDx = aX - aNode->x() ;
            real tDy = aY - aNode->y() ;
            real tDist = std::sqrt( tDx * tDx + tDy * tDy ) ;

            if ( tDist < aBestDist )
            {
                aBestDist = tDist ;
                aBest = aNode ;
            }

            if ( aNode->number_of_nodes() == 0 )
            {
                return aBest ;
            }

            // signed distance to the splitting plane (alternates x/y)
            real tDelta = aFlip ? ( aY - aNode->y() ) : ( aX - aNode->x() ) ;

            // search the near subtree first, then prune the far subtree
            // if the splitting plane is farther than the current best
            Node * tNear = tDelta < 0.0 ? aNode->node( 0 ) : aNode->node( 1 ) ;
            Node * tFar  = tDelta < 0.0 ? aNode->node( 1 ) : aNode->node( 0 ) ;

            bool tFlip = ! aFlip ;

            aBest = this->find_closest_node( tNear, aX, aY, tFlip, aBest, aBestDist ) ;

            if ( std::abs( tDelta ) < aBestDist )
            {
                aBest = this->find_closest_node( tFar, aX, aY, tFlip, aBest, aBestDist ) ;
            }

            return aBest ;
        }

        Node * PeriodicityFactory::create_kdtree( Cell< Node * > & aNodes )
        {
            if ( aNodes.size() > 0 )
            {
                return this->kdsort( aNodes, 0, aNodes.size(), false );
            }
            return nullptr;
        }

        Node *
        PeriodicityFactory::kdsort(
            Cell< Node * > & aNodes, index_t aStart,
            index_t aEnd, const bool aFlip )
        {
            if ( aStart >= aEnd )
            {
                return nullptr ;
            }

            // leaf node
            if ( aStart + 1 == aEnd )
            {
                return aNodes( aStart ) ;
            }

            // sort range by current dimension
            if ( aFlip )
            {
                sort( aNodes.begin() + aStart,
                      aNodes.begin() + aEnd,
                      []( Node * a, Node * b ) { return a->y() < b->y(); } ) ;
            }
            else
            {
                sort( aNodes.begin() + aStart,
                      aNodes.begin() + aEnd,
                      []( Node * a, Node * b ) { return a->x() < b->x(); } ) ;
            }

            // pick the median as the splitting node
            index_t tMid = ( aStart + aEnd ) / 2 ;
            Node * tNode = aNodes( tMid ) ;

            bool tFlip = ! aFlip ;

            // recurse into left [aStart, tMid) and right [tMid+1, aEnd)
            Node * tLeft  = this->kdsort( aNodes, aStart, tMid, tFlip ) ;
            Node * tRight = this->kdsort( aNodes, tMid + 1, aEnd, tFlip ) ;

            // attach children: index 0 = left, index 1 = right
            if ( tLeft != nullptr || tRight != nullptr )
            {
                tNode->allocate_node_container( 2 ) ;
                tNode->insert_node( tLeft, 0 ) ;
                tNode->insert_node( tRight, 1 ) ;
            }

            return tNode ;
        }

        void
        PeriodicityFactory::fix_face_slaves( Periodicity * aPeriodicity )
        {
            Cell< Face * > & tSlaveFaces = aPeriodicity->slave_faces();

            if ( tSlaveFaces.empty() ) return;

            // For each slave face, determine its orientation relative to
            // its volume element's facet ordering, then reassign it from
            // master to slave role with the correct orientation index.
            Cell< Node * > tNodes ;
            index_t tCount = 0 ;
            for ( Face * tFace : tSlaveFaces )
            {
                // idempotency guard: a periodic boundary face that was already
                // converted to slave-owned by a previous update() has a null
                // master. Skip it so a second pass can't dereference nullptr.
                if ( tFace->master() == nullptr ) continue ;

                tFace->flag_nodes() ;
                uint tOrientation = BELFEM_UINT_MAX ;

                // the volume element that owns this face
                Element * tElement = tFace->master() ;

                // find which facet of the volume element matches this face
                for ( uint f=0; f<tElement->number_of_faces(); ++f )
                {
                    tCount = 0 ;
                    tElement->get_corner_nodes_of_facet( f, tNodes ) ;
                    for ( Node * tNode : tNodes )
                    {
                        if ( tNode->is_flagged() )
                        {
                            ++tCount ;
                        }
                        else
                        {
                            break ;
                        }
                    }

                    if ( tCount == tFace->number_of_corner_nodes() )
                    {
                        // determine rotation: find which face node matches
                        // the first corner node of the element facet
                        id_t tID = tNodes.first()->original()->id();
                        for ( uint k=0; k<tCount; ++k )
                        {
                            if ( tID == tFace->node( k )->original()->id() )
                            {
                                tOrientation = k + 1 ;
                                tFace->set_master( nullptr, BELFEM_UINT_MAX ) ;
                                tFace->set_slave( tElement, f, tOrientation ) ;
                                break ;
                            }
                        }
                    }

                    if ( tOrientation != BELFEM_UINT_MAX )
                    {
                        break ;
                    }
                }

                BELFEM_ASSERT( tOrientation != BELFEM_UINT_MAX, "Could not find orientation" );

                tFace->unflag_nodes();
            }
        }

        void
        PeriodicityFactory::create_hesse( Cell< Node * > & aNodes , Vector< real > & aHesse )
        {
            BELFEM_ASSERT( aNodes.size() == 3, "Invalid number of nodes, need exactly 3." );

            Vector< real > tP( 3 );
            aNodes( 0 )->get_coords( tP );

            Vector< real > tQ( 3 );
            aNodes( 1 )->get_coords( tQ );

            Vector< real > tR( 3 );
            aNodes( 2 )->get_coords( tR );

            Vector< real > tM( 3, 0.0 );
            tM += tP ;
            tM += tQ ;
            tM += tR ;
            tM /= 3.0 ;

            tQ -= tP;
            tR -= tP;

            Vector< real > tN( cross( tQ, tR ) );

            real tNorm = norm( tN );
            BELFEM_ERROR( tNorm > BELFEM_MESH_EPSILON,
                "could not determine surface normal: plane points are collinear" );

            tN /= tNorm;

            aHesse.set_size( 4 );
            aHesse( 0 ) = tN( 0 );
            aHesse( 1 ) = tN( 1 );
            aHesse( 2 ) = tN( 2 );
            aHesse( 3 ) = dot( tN, tM );
        }

        void
        PeriodicityFactory::tag_periodic_sidesets()
        {
            BELFEM_ERROR( mMesh->is_finalized(), "PeriodicityFactory::tag_periodic_sidesets() must operate on a finalized mesh" );
            BELFEM_ERROR( mSourceHesse.length() == 4 && mTargetHesse.length() == 4,
                "planes must be set before tagging periodic sidesets" );

            Vector< real > tNm = { mSourceHesse( 0 ), mSourceHesse( 1 ), mSourceHesse( 2 ) };
            real tDm = mSourceHesse( 3 );

            Vector< real > tNs = { mTargetHesse( 0 ), mTargetHesse( 1 ), mTargetHesse( 2 ) };
            real tDs = mTargetHesse( 3 );

            Vector< real > tP( 3 );

            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                if ( tSideSet->nodes().size() == 0 ) continue ;

                bool tIsMaster = true ;
                bool tIsSlave  = true ;

                for ( Node * tNode : tSideSet->nodes() )
                {
                    tNode->get_coords( tP );

                    tIsMaster = tIsMaster && std::abs( dot( tNm, tP ) - tDm ) < BELFEM_MESH_EPSILON ;
                    tIsSlave  = tIsSlave  && std::abs( dot( tNs, tP ) - tDs ) < BELFEM_MESH_EPSILON ;

                    // most sidesets leave the planes at their first nodes
                    if ( ! ( tIsMaster || tIsSlave ) ) break ;
                }

                if ( tIsMaster || tIsSlave )
                {
                    tSideSet->set_domain_type( DomainType::Periodic );
                }

            }
        }

        void
        PeriodicityFactory::reset_bitsets( const index_t aNumNodes, const index_t aNumEdges, const index_t aNumFaces )
        {
            if ( mNodeBitset == nullptr )
            {
                mNodeBitset = new DynamicBitset( aNumNodes );
            }
            else  if ( mNodeBitset->size() == aNumNodes )
            {
                mNodeBitset->reset();
            }
            else
            {
                delete mNodeBitset ;
                mNodeBitset = new DynamicBitset( aNumNodes );
            }

            if ( mMesh->edges_exist() )
            {
                if ( mEdgeBitset == nullptr )
                {
                    mEdgeBitset = new DynamicBitset( aNumEdges );
                }
                else if ( mEdgeBitset->size() == aNumEdges )
                {
                    mEdgeBitset->reset();
                }
                else
                {
                    delete mEdgeBitset ;
                    mEdgeBitset = new DynamicBitset( aNumEdges );
                }
            }
            else if ( mEdgeBitset != nullptr )
            {
                delete mEdgeBitset ;
                mEdgeBitset = nullptr ;
            }

            if ( mMesh->faces_exist() )
            {
                if ( mFaceBitset == nullptr )
                {
                    mFaceBitset = new DynamicBitset( aNumFaces );
                }
                else if ( mFaceBitset->size() == aNumFaces )
                {
                    mFaceBitset->reset();
                }
                else
                {
                    delete mFaceBitset ;
                    mFaceBitset = new DynamicBitset( aNumFaces );
                }
            }
            else if ( mFaceBitset != nullptr )
            {
                delete mFaceBitset ;
                mFaceBitset = nullptr ;
            }
        }

        void
        PeriodicityFactory::crosslink( Periodicity * aPeriodicity )
        {
            if ( aPeriodicity->master_nodes().size() > 0 )
            {
                this->crosslink( aPeriodicity->master_nodes(), aPeriodicity->slave_nodes() );
                aPeriodicity->flag( EntityType::NODE );
            }
            if ( aPeriodicity->master_edges().size() > 0 )
            {
                this->crosslink( aPeriodicity->master_edges(), aPeriodicity->slave_edges() );
                aPeriodicity->flag( EntityType::EDGE );
            }
            if ( aPeriodicity->master_faces().size() > 0 )
            {
                this->crosslink( aPeriodicity->master_faces(), aPeriodicity->slave_faces() );
                aPeriodicity->flag( EntityType::FACE );
            }
            if ( aPeriodicity->master_facets().size() > 0 )
            {
                this->crosslink( aPeriodicity->master_facets(), aPeriodicity->slave_facets() );
                aPeriodicity->flag( EntityType::FACET );
            }
        }

    }
}