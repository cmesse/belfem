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

#include "assert.hpp"
#include "cl_Mesh_Periodicity.hpp"

#include "cl_Communicator.hpp"
#include "cl_Mesh.hpp"
#include "cl_Mesh_PeriodicityFactory.hpp"
namespace belfem
{
    namespace mesh
    {
        Periodicity::Periodicity( Mesh *aMesh ) :
            mMesh( aMesh )
        {

        }

        Periodicity::~Periodicity()
        {
            this->reset_nodes();
            this->reset_edges();
            this->reset_faces();
            this->reset_facets();
            this->clear_node_pair_backup();
        }

        void
        Periodicity::update()
        {
            PeriodicityFactory  tFactory( mMesh );

            tFactory.set_master_plane( mMasterPlane );
            tFactory.set_slave_plane( mSlavePlane );

            tFactory.update_periodicity( this );
        }

        void
        Periodicity::set_entity_dependencies()
        {
            BELFEM_ASSERT( gComm.rank() == 0, "Only rank 0 can set entity dependencies" );

            index_t n = mMasterNodes.size();
            for ( index_t k = 0; k < n; ++k )
            {
                Node * A = mMasterNodes( k );
                Node * B = mSlaveNodes( k );

                if ( B->is_hanging() ) B->reset_source_container();

                if ( A->is_hanging() )
                {
                    uint m = A->number_of_sources();
                    B->allocate_source_container( m );
                    for ( uint i = 0; i < m; ++i )
                    {
                        B->add_source( A->source( i ), A->weight( i ) );
                    }
                }
                else
                {
                    B->allocate_source_container( 1 );
                    B->add_source( A, 1.0 );
                }
            }

            if ( mMesh->edges_exist() )
            {
                n = mMasterEdges.size();
                for ( index_t e = 0 ; e<n; ++e )
                {
                    Edge * A = mMasterEdges( e );

                    // after match_edges() the edge lists are compacted to tied pairs,
                    // so mMasterEdges(e) <-> mSlaveEdges(e) is aligned. We still go
                    // through the partner pointer set by match_edges(), which is the
                    // authoritative link.
                    Edge * B = A->periodic();
                    BELFEM_ERROR( B != nullptr, "master periodic edge ( index %lu ) has no slave partner",
                                  ( long unsigned int ) A->index() );

                    if ( B->is_hanging() ) B->reset_source_container();

                    // orientation is original-aware: a cut-duplicate endpoint has no
                    // periodic() of its own, but its original is the periodic partner of
                    // the master edge's original ( for non-duplicates original()==self,
                    // so this is backward-compatible ).
                    BELFEM_ASSERT( A->node( 0 )->original()->periodic() == B->node( 0 )->original(), "invalid edge orientation at periodicity" );
                    BELFEM_ASSERT( A->node( 1 )->original()->periodic() == B->node( 1 )->original(), "invalid edge orientation at periodicity" );
                    if ( A->is_hanging() )
                    {
                        uint m = A->number_of_sources();
                        B->allocate_source_container( m );
                        for ( uint i = 0; i < m; ++i )
                        {
                            B->add_source( A->source( i ), A->weight( i ) );
                        }
                    }
                    else
                    {
                        B->allocate_source_container( 1 );
                        B->add_source( A, 1.0 );
                    }
                }
            }
            if ( mMesh->faces_exist() )
            {
                n = mMasterFaces.size();
                for ( index_t f = 0 ; f<n; ++f )
                {
                    Face * A = mMasterFaces( f );
                    Face * B = mSlaveFaces( f );

                    if ( B->is_hanging() ) B->reset_source_container();

                    // we need to set the weights to NAN because they depend on orientation
                    // and type of the dof. This must be handled specifically to the considered physics
                    if ( A->is_hanging() )
                    {
                        uint m = A->number_of_sources();
                        B->allocate_source_container( m );
                        for ( uint i = 0; i < m; ++i )
                        {
                            B->add_source( A->source( i ), BELFEM_QUIET_NAN );
                        }
                    }
                    else
                    {
                        B->allocate_source_container( 1 );
                        B->add_source( A, BELFEM_QUIET_NAN );
                    }
                }
            }

            n = mMasterFacets.size();
            for ( index_t f = 0 ; f<n; ++f )
            {
                Facet * A = mMasterFacets( f );
                Facet * B = mSlaveFacets( f );

                if ( B->is_hanging() ) B->reset_source_container();

                // we need to set the weights to NAN because they depend on orientation
                // and type of the dof. This must be handled specifically to the considered physics
                if ( A->is_hanging() )
                {
                    uint m = A->number_of_sources();
                    B->allocate_source_container( m );
                    for ( uint i = 0; i < m; ++i )
                    {
                        B->add_source( A->source( i ), BELFEM_QUIET_NAN );
                    }
                }
                else
                {
                    B->allocate_source_container( 1 );
                    B->add_source( A, BELFEM_QUIET_NAN );
                }
            }
        }

        void
        Periodicity::reset_nodes()
        {
            if ( this->is_flagged( EntityType::NODE  ) )
            {
                for ( Node * tNode : mMasterNodes )
                {
                    tNode->set_periodic( nullptr );
                }
                for ( Node * tNode : mSlaveNodes )
                {
                    tNode->set_periodic( nullptr );
                }

                this->unflag( EntityType::NODE );
            }

            mMasterNodes.clear();
            mSlaveNodes.clear();
        }

        void
        Periodicity::reset_edges()
        {
            if ( this->is_flagged( EntityType::EDGE ) )
            {
                for ( Edge * tEdge : mMasterEdges )
                {
                    tEdge->set_periodic( nullptr );
                }
                for ( Edge * tEdge : mSlaveEdges )
                {
                    tEdge->set_periodic( nullptr );
                }

                this->unflag( EntityType::EDGE );
            }

            mMasterEdges.clear();
            mSlaveEdges.clear();

        }

        void
        Periodicity::reset_faces()
        {
            if ( this->is_flagged( EntityType::FACE ) )
            {
                for ( Face * tFace : mMasterFaces )
                {
                    tFace->set_periodic( nullptr );
                }
                for ( Face * tFace : mSlaveFaces )
                {
                    tFace->set_periodic( nullptr );
                }

                this->unflag( EntityType::FACE );
            }

            mMasterFaces.clear();
            mSlaveFaces.clear();
        }

        void
        Periodicity::reset_facets()
        {
            if ( this->is_flagged( EntityType::FACET ) )
            {
                for ( Facet * tFacet : mMasterFacets )
                {
                    tFacet->set_periodic( nullptr );
                }
                for ( Facet * tFacet : mSlaveFacets )
                {
                    tFacet->set_periodic( nullptr );
                }

                this->unflag( EntityType::FACET );
            }

            mMasterFacets.clear();
            mSlaveFacets.clear();

        }
        void
        Periodicity::backup_node_pairs()
        {

            BELFEM_ASSERT( mMasterNodes.size() == mSlaveNodes.size(), "Size of node pairs does not match" );

            index_t n = mMasterNodes.size() ;

            mNodePairBackup.clear();
            mNodePairBackup.reserve( n );

            for ( index_t k=0; k<n; ++k )
            {
                this->add_node_pair_to_backup(mMasterNodes( k ),  mSlaveNodes ( k ) );
            }

        }

        void
        Periodicity::restore_node_pairs()
        {
            index_t n = mNodePairBackup.size();

            BELFEM_ASSERT( n > 0, "size of node pair backup is zero" );
            mMasterNodes.set_size( n, nullptr );
            mSlaveNodes.set_size( n, nullptr );
            for ( index_t k=0; k<n; ++k )
            {
                Node * A = mNodePairBackup(k).first;
                Node * B = mNodePairBackup(k).second;

                A->set_periodic( B );
                B->set_periodic( A );

                BELFEM_ASSERT( this->pair_orientation( A, B, mMaxTolerance ) == 1,
                    "node pair backup is corrupted ( nodes %lu and %lu )",
                    ( long unsigned int ) A->id(),
                    ( long unsigned int ) B->id() );

                mMasterNodes( k ) = A ;
                mSlaveNodes ( k ) = B ;
            }

            mNodePairsRestored = true ;
        }

        void
        Periodicity::clear_node_pair_backup()
        {
            mNodePairBackup.clear();
        }


        void
        Periodicity::add_node_pair_to_backup( Node * aNodeA, Node * aNodeB, const real aTolerance )
        {
            // needed for assert during backup
            mMaxTolerance = std::max( mMaxTolerance, aTolerance );

            switch ( this->pair_orientation( aNodeA, aNodeB, aTolerance ) )
            {
                case 1 :
                {
                    mNodePairBackup.push( std::make_pair( aNodeA, aNodeB ) );
                    break ;
                }
                case -1 :
                {
                    mNodePairBackup.push( std::make_pair( aNodeB, aNodeA ) );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Periodicity: nodes are not on the periodic surfaces" );
                }
            }
        }

    }
}