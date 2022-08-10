//
// Created by Christian Messe on 27.10.19.
//
#include "cl_Mesh_Partitioner.hpp"

#include "commtools.hpp"
#include "random.hpp"

#include "cl_Timer.hpp"
#include "cl_Vector.hpp"
#include "cl_Logger.hpp"


namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Partitioner::Partitioner( Mesh * aMesh,
                const uint aNumberOfPartitions,
                bool aSetProcOwnerships ) :
            mMesh( aMesh ),
            mNumberOfPartitions( aNumberOfPartitions )
        {
#ifdef BELFEM_SUITESPARSE

            // prepare the graph
            this->create_graph();

            // call metis
            this->run_metis();

            if( aSetProcOwnerships )
            {
                // set element owners
                this->set_element_owners();

                // tell mesh how many partitions exist
                aMesh->set_number_of_partitions( mNumberOfPartitions );

                // fix a rare bug
                this->fix_facet_related_ownerships( );
                this->fix_cut_related_ownerships(  );


                // set node owners
                aMesh->set_node_owners();

                //this->fix_cut_related_ownerships();


                // set vertex owners
                aMesh->set_vertex_owners();

                //this->fix_facet_related_ownerships();


                // fix ghost sidesets
                this->fix_ghost_sideset_related_ownerships() ;

                // set owners of domain cuts
                aMesh->set_connector_owners();
            }

#elif
         message( 0, "WARNING: The mesh partitioner was called, but the executable  is not" );
         message( 0,"          linked against METIS. The mesh will not be partitioned!");
#endif
        }

//------------------------------------------------------------------------------b

        Partitioner::~Partitioner()
        {
            if(  mNumberOfElements > 0 )
            {
                free( mElementPointers );
            }
            if( mAdjacencySize > 0 )
            {
                free( mAdjacency );
            }
        }

//------------------------------------------------------------------------------

        void
        Partitioner::create_graph()
        {

            // start the timer
            Timer tTimer;

            // get pointer to mesh Cell
            Cell< Element * > & tElements = mMesh->elements();


            // number of flagged elements
            mNumberOfElements = 0 ;

            // siye of adjacency
            mAdjacencySize = 0;

            // counter for elements
            metis_t tCount = 0;

            // loop over all elements
            for( Element * tElement : tElements )
            {
                if( tElement->is_flagged() )
                {
                    // increment adjacency counter
                    uint tNumElements = tElement->number_of_elements();

                    // loop over all elements
                    for ( uint e = 0; e < tNumElements; ++e )
                    {
                        if ( tElement->element( e )->is_flagged())
                        {
                            mAdjacencySize += 1;
                        }
                    }

                    // also, update the element index
                    tElement->set_index( tCount++ );

                    // increment element counter
                    ++mNumberOfElements ;
                }
            }

            BELFEM_ERROR( mNumberOfElements > 0, "No elements were selected for partitioning");

            if( mNumberOfElements > 0 )
            {
                mElementPointers = ( metis_t * ) malloc(( mNumberOfElements + 1 ) * sizeof( metis_t ));
            }
            if( mAdjacencySize > 0 )
            {
                mAdjacency = ( metis_t * ) malloc( mAdjacencySize * sizeof( metis_t ));
            }

            // reset the counter
            tCount = 0;

            // loop over all elements
            for( Element * tElement : tElements )
            {
                if( tElement->is_flagged() )
                {
                    mElementPointers[ tElement->index() ] = tCount;

                    uint tNumConnectedElements = ( metis_t ) tElement->number_of_elements();

                    for ( uint k = 0; k < tNumConnectedElements; ++k )
                    {
                        if( tElement->element( k )->is_flagged() )
                        {
                            mAdjacency[ tCount++ ] = ( metis_t ) tElement->element( k )->index();
                        }
                    }
                }
            }

            mElementPointers[ mNumberOfElements ] = tCount;

            // reset element index
            tCount = 0 ;
            for( Element * tElement : tElements )
            {
                tElement->set_index( tCount++ );
            }

            message( 4, "    Created Mesh Graph: %i ms\n", ( uint ) tTimer.stop() );

        }

//------------------------------------------------------------------------------

        void
        Partitioner::run_metis()
        {
#ifdef BELFEM_SUITESPARSE

            // set options
            metis_t tOptions[ METIS_NOPTIONS ];
            METIS_SetDefaultOptions( tOptions );
            tOptions[ METIS_OPTION_COMPRESS ]  = 0; // Do not try to compress the graph.
            tOptions[ METIS_OPTION_CONTIG ]    = 1; // Forces contiguous partitions.
            tOptions[ METIS_OPTION_NUMBERING ] = 0; // Use zero based numbering
            tOptions[ METIS_OPTION_OBJTYPE ] = METIS_OBJTYPE_VOL;

            // set the verbose flag
            if( gLog.info_level() >= 5 )
            {
                tOptions[ METIS_OPTION_DBGLVL ] = METIS_DBG_INFO; // METIS_DBG_INFO;
            }

            tOptions[ METIS_OPTION_CTYPE]  = METIS_CTYPE_SHEM; // SHEM crashes or METIS_CTYPE_RM
            tOptions[ METIS_OPTION_RTYPE ] = METIS_RTYPE_SEP1SIDED;

            // the result
            metis_t tEdgeCut;

            // allocate the element data
            mPartition.set_size( mNumberOfElements );

            metis_t tNumberOfConstraints = 1;

            message( 3, "    Starting METIS ..." );

            // start the timer
            Timer tTimer;

            int tStatus =
                    METIS_PartGraphKway(
                            &mNumberOfElements,     // The number of vertices in the graph.
                            &tNumberOfConstraints,  // The number of balancing constraints. It should be at least 1.
                            mElementPointers,       // The adjacency structure of the graph
                            mAdjacency,             // The adjacency structure of the graph
                            NULL,                   // The weights of the vertices
                            NULL,                   // The size of the vertices for computing the total communication volume
                            NULL,                   // The weights of the edges
                            &mNumberOfPartitions,   // The number of parts to partition the graph.
                            NULL,                   // Target Partition weights
                            NULL,                   // the allowed load imbalance tolerance for each constraint
                            tOptions,               // Options for METIS
                            &tEdgeCut,
                            mPartition.data() );

            // interpret status from METIS
            string tMessage = metis_status( tStatus );

            BELFEM_ERROR( tStatus == METIS_OK,
                         "An error occured while calling METIS: %s",
                         tMessage.c_str() );

            message( 3, "    ... time for repartitioning         : %u ms",
                     ( unsigned int ) tTimer.stop() );

#endif
        }

//------------------------------------------------------------------------------

        string
        Partitioner::metis_status( const int aStatus )
        {
            switch ( aStatus )
            {
                case( METIS_OK ) :
                {
                    return "Success";
                }
                case( METIS_ERROR_INPUT ) :
                {
                    return "Invalid input";
                }
                case( METIS_ERROR_MEMORY ) :
                {
                    return "Memory Error";
                }
                default :
                {
                    return "Unknown Error";
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Partitioner::set_element_owners()
        {
            // start the timer
            Timer tTimer;

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
                else
                {
                    tElement->set_owner( mNumberOfPartitions ) ;
                }
            }

            message( 4, "    ... time for setting element owners : %u ms",
                     ( unsigned int ) tTimer.stop() );
        }

//------------------------------------------------------------------------------

        void
        Partitioner::fix_facet_related_ownerships( )
        {
            mMesh->unflag_all_nodes() ;

            Cell< Facet * > & tFacets = mMesh->facets() ;

            index_t tCount = 1 ;

            while( tCount > 0 )
            {
                tCount = 0 ;
                for ( Facet * tFacet: tFacets )
                {
                    proc_t tOwner = tFacet->master()->owner() ;

                    if( tFacet->has_slave() )
                    {
                        if( tFacet->slave()->owner() < tOwner )
                        {
                            tOwner = tFacet->slave()->owner() ;
                        }
                    }

                    tFacet->set_owner( tOwner );

                    if( tFacet->master()->is_flagged() && tFacet->master()->owner() > tOwner )
                    {
                       tFacet->master()->set_owner( tOwner );
                       ++tCount ;
                    }
                    if( tFacet->has_slave() )
                    {
                        if( tFacet->slave()->is_flagged() &&  tFacet->slave()->owner() > tOwner )
                        {
                            tFacet->slave()->set_owner( tOwner );
                            ++tCount ;
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Partitioner::fix_ghost_sideset_related_ownerships()
        {
            // note: this function is called after fix_facet_related_ownerships()
            // latter function sets the sideset ownership based on the master
            // and slave elements. What is left to be done is fix the remaining
            // node, face, and edge ownerships

            mMesh->unflag_all_nodes() ;
            mMesh->unflag_all_edges() ;

            for ( id_t tID : mMesh->ghost_sideset_ids() )
            {
                // grab ghost
                SideSet * tSideSet = mMesh->sideset( tID ) ;

                // get entity numbers
                uint tNumberOfNodes = mesh::number_of_nodes( tSideSet->element_type() );
                uint tNumberOfEdges = mesh::number_of_edges( tSideSet->element_type() );

                // loop over all facets
                for( Facet * tFacet : tSideSet->facets() )
                {
                    // the clone of the element
                    Element * tClone = tFacet->element() ;

                    // grab original facet
                    Element * tOriginal = tClone->element( 0 );

                    // should not be neccessary, but we do it anywase
                    tClone->set_owner( tOriginal->owner() );

                    for( uint k=0; k<tNumberOfNodes; ++k )
                    {
                        // grab node
                        Node * tNode = tClone->node( k );

                        if( ! tNode->is_flagged() )
                        {
                            tNode->set_owner( tOriginal->node( k )->owner() ) ;
                            tNode->flag() ;
                        }
                    }
                }

                // fix edge ownerships
                if( mMesh->edges_exist() )
                {
                    for( Facet * tFacet : tSideSet->facets() )
                    {
                        // the clone of the element
                        Element * tClone = tFacet->element() ;

                        // grab original facet
                        Element * tOriginal = tClone->element( 0 );

                        for( uint e=0; e<tNumberOfEdges; ++e )
                        {
                            // grab node
                            Edge * tEdge = tClone->edge( e );

                            if( ! tEdge->is_flagged() )
                            {
                                tEdge->set_owner( tOriginal->edge( e )->owner() ) ;
                                tEdge->flag() ;
                            }
                        }
                    }
                }

                if( mMesh->faces_exist() && mMesh->number_of_dimensions() == 3 )
                {
                    for( Facet * tFacet : tSideSet->facets() )
                    {

                        if( ! tFacet->is_flagged() )
                        {
                            tFacet->face()->set_owner( tFacet->owner() );
                        }
                    }
                }
            }

            mMesh->unflag_all_nodes() ;
            mMesh->unflag_all_edges() ;
        }

//------------------------------------------------------------------------------

        void
        Partitioner::fix_cut_related_ownerships()
        {
            /*Matrix< id_t > &  tNodeTable = mMesh->node_cut_table() ;

            index_t tNumNodes = tNodeTable.n_cols() ;

            index_t tCount = 1 ;

            while( tCount > 0 )
            {
                // reset the counter
                tCount = 0;

                // check element owners of nodes connected to cuts
                for ( index_t k = 0; k < tNumNodes; ++k )
                {
                    proc_t tOwner = gNoOwner;

                    for( uint j=0; j<2; ++j )
                    {
                        Node * tN = mMesh->node( tNodeTable( j, k ));
                        for ( uint e = 0; e < tN->number_of_elements(); ++e )
                        {
                            Element * tE = tN->element( e );

                            if ( tE->owner() < tOwner )
                            {
                                tOwner = tE->owner();
                            }
                        }
                    }

                    for( uint j=0; j<2; ++j )
                    {
                        Node * tN = mMesh->node( tNodeTable( j, k ));
                        for ( uint e = 0; e < tN->number_of_elements(); ++e )
                        {
                            Element * tE = tN->element( e );

                            if ( tOwner < tE->owner() )
                            {
                                tE->set_owner( tOwner );
                                ++tCount ;
                            }
                        }
                    }
                }
            } */
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