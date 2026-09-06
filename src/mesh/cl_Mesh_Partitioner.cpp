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

#include "cl_Mesh_Partitioner.hpp"
#include "random.hpp"

#include "cl_Timer.hpp"
#include "cl_Vector.hpp"
#include "cl_Logger.hpp"
#include "fn_Graph_METIS.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Partitioner::Partitioner( Mesh * aMesh,
                const uint aNumberOfPartitions,
                const bool aSetProcOwnerships,
                const bool aForceContiguousPartitions,
                const bool aResetVertexContainers ) :
            mMesh( aMesh ),
            mNumberOfPartitions( aNumberOfPartitions ),
            mForceContiguousPartitions( aForceContiguousPartitions ),
            mResetVertexContainers( aResetVertexContainers )
        {
#ifdef BELFEM_METIS

            // prepare the graph
            Graph tGraph ;
            this->create_graph( tGraph );

            // call metis
            this->run_metis( tGraph );

            if( aSetProcOwnerships )
            {
                // set element owners
                this->set_element_owners();

                // tell mesh how many partitions exist
                aMesh->set_number_of_partitions( mNumberOfPartitions );

                // fix an ownership bug
                this->fix_facet_related_ownerships( );

                // update nodes, edges, faces, control points and vertices
                aMesh->update_ownerships();
            }

#else
         message( InfoLevel::Minimal, "WARNING: The mesh partitioner was called, but the executable  is not" );
         message( InfoLevel::Minimal,"          linked against METIS. The mesh will not be partitioned!");
#endif
        }

//------------------------------------------------------------------------------b

        Partitioner::~Partitioner()
        {
        }

//------------------------------------------------------------------------------

        void
        Partitioner::create_graph( Graph & aGraph )
        {
            // start the timer
            Timer tTimer;

            // get pointer to mesh Cell
            Cell< Element * > & tElements = mMesh->elements();

            mMesh->unflag_all_nodes() ;

            index_t tCount = 0 ;
            for ( Element * tElement : tElements )
            {
                if ( tElement->is_flagged() )
                {
                    ++ tCount ;
                    tElement->reset_vertex_container() ;
                }
                else
                {
                    tElement->set_index( gNoIndex );
                }
            }

            DynamicBitset tBitset( mMesh->number_of_elements() ) ;
            aGraph.set_size( tCount, nullptr ) ;

            uint d = mMesh->number_of_dimensions() ;
            Cell< index_t > tIndices ;

            tCount = 0 ;
            for ( Element * tElement : tElements )
            {
                if ( ! tElement->is_flagged() ) continue;

                for ( uint k=0; k<tElement->number_of_corner_nodes(); ++k )
                {
                    tElement->node( k )->original()->flag() ;
                }

                tBitset.reset();
                tBitset.set( tElement->index() ) ;
                for ( uint i = 0 ; i < tElement->number_of_corner_nodes() ; ++i )
                {
                    uint n = tElement->node( i )->original()->number_of_duplicates() ;

                    for ( uint k = 0 ; k <= n ; ++k )
                    {
                        Node * tNode = k == n ? tElement->node( i )->original() : tElement->node( i )->original()->duplicate( k ) ;

                        for ( index_t e = 0 ; e<tNode->number_of_elements() ; ++e )
                        {
                            // get the neighbor
                            Element * tNeighbor = tNode->element( e ) ;

                            if ( ! tNeighbor->is_flagged() || tBitset.test( tNeighbor->index() )) continue ;

                            uint c = 0 ;
                            for ( uint j = 0 ; j < tNeighbor->number_of_corner_nodes() ; ++j )
                            {
                                if ( tNeighbor->node( j )->original()->is_flagged() ) ++c ;
                                if ( c == d )
                                {
                                    tBitset.set( tNeighbor->index() ) ;
                                    break ;
                                }
                            }
                        }
                    }
                }
                tBitset.reset( tElement->index() ) ;
                for ( uint k=0; k<tElement->number_of_corner_nodes(); ++k )
                {
                    tElement->node( k )->original()->unflag() ;
                }
                tBitset.where( tIndices ) ;

                tElement->init_vertex_container( tIndices.size() );
                for ( index_t e : tIndices )
                {
                    tElement->insert_vertex( tElements( e ) );
                }
                aGraph( tCount++ ) = tElement ;
            }

            tCount = 0 ;
            for ( graph::Vertex * tVertex : aGraph )
            {
                tVertex->set_index( tCount++ ) ;
            }

            message( InfoLevel::Verbose, "    Created Mesh Graph: %i ms\n", ( uint ) tTimer.stop() );

        }

//------------------------------------------------------------------------------

        void
        Partitioner::run_metis( Graph & aGraph )
        {
#ifdef BELFEM_METIS

            Timer tTimer;

            message( InfoLevel::Detailed, "    ... time for repartitioning         : %u ms",
                     ( unsigned int ) tTimer.stop() );

            graph::metis_partition( aGraph,
                mNumberOfPartitions,
                mForceContiguousPartitions,
                & mPartition );

            // tidy up element container
            Cell< Element * > & tElements = mMesh->elements();
            index_t tCount = 0 ;

            if ( mResetVertexContainers )
            {
                for ( Element * tElement : tElements )
                {
                    tElement->set_index( tCount++ );
                    if ( tElement->is_flagged() )
                    {
                        tElement->reset_vertex_container();
                    }
                }
            }
            else
            {
                for ( Element * tElement : tElements )
                {
                    tElement->set_index( tCount++ );
                }
            }

#endif
        }

//------------------------------------------------------------------------------

        void
        Partitioner::set_element_owners()
        {

            // get ref to elements
            Cell< Element * > & tElements = mMesh->elements();

            // reset counter
            index_t tCount = 0;

            // loop over all nodes
            for( Element * tElement : tElements )
            {
                if( tElement->is_flagged() )
                {
                    // set the owner
                    tElement->set_owner( mPartition( tCount++ ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Partitioner::fix_facet_related_ownerships( )
        {
            Cell< Facet * > & tFacets = mMesh->facets();
            for ( Facet * tFacet : tFacets )
            {
                tFacet->set_owner( gNoOwner );
            }

            index_t tSweep = 1 ;
            while ( tSweep != 0 )
            {
                tSweep = 0 ;
                for ( Facet * tFacet : tFacets )
                {
                    proc_t tO = tFacet->owner() ;
                    if ( tFacet->has_master() )
                    {
                        tO = std::min( tFacet->master()->owner(), tO );
                    }
                    if ( tFacet->has_slave() )
                    {
                        tO = std::min( tFacet->slave()->owner(), tO );
                    }

                    if ( tFacet->owner() != tO )
                    {
                        tFacet->set_owner( tO ) ;
                        ++tSweep ;
                    }
                    if ( tFacet->has_master() )
                    {
                        if ( tFacet->master()->owner() != tO )
                        {
                            tFacet->master()->set_owner( tO ) ;
                            ++tSweep ;
                        }
                    }
                    if ( tFacet->has_slave() )
                    {
                        if ( tFacet->slave()->owner() != tO )
                        {
                            tFacet->slave()->set_owner( tO ) ;
                            ++tSweep ;
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        const Vector < metis_t > &
        Partitioner::partition()
        {
            return mPartition;
        }

//------------------------------------------------------------------------------
    }
}
