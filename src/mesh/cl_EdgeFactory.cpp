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
#include <cmath>

#include "commtools.hpp"
#include "cl_EdgeFactory.hpp"
#include "cl_Element.hpp"
#include "meshtools.hpp"
#include "fn_unique.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "meshtools.hpp"
#include "op_Graph_Vertex_ID.hpp"
#include "cl_Curve.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        EdgeFactory::EdgeFactory( Mesh & aMesh ) :
                mCommRank( comm_rank() ),
                mMesh( aMesh ),
                mNumberOfNodes( aMesh.number_of_nodes() )
        {
            // Validate that mesh size is compatible with key_t for edge hashing
            // Edge keys use: larger_index * N + smaller_index
            // Maximum safe N is sqrt(key_t_max)
            const key_t tMaxNodes = static_cast< key_t >( std::sqrt( static_cast< double >( std::numeric_limits< key_t >::max() ) ) );
            BELFEM_ERROR( mNumberOfNodes <= tMaxNodes,
                "Mesh has too many nodes (%llu) for edge key computation with current key_t type.\n"
                "Maximum allowed nodes: %llu\n"
                "Recommendation: Change typedef for key_t in src/core/typedefs.hpp to a larger unsigned integer type (e.g., __uint128_t)",
                ( long long unsigned int ) mNumberOfNodes,
                ( long long unsigned int ) tMaxNodes );
        }

//------------------------------------------------------------------------------

        EdgeFactory::EdgeFactory( Mesh * aMesh ) :
            mCommRank( comm_rank() ),
            mMesh( * aMesh ),
            mNumberOfNodes( aMesh->number_of_nodes() )
        {
            // Validate that mesh size is compatible with key_t for edge hashing
            // Edge keys use: larger_index * N + smaller_index
            // Maximum safe N is sqrt(key_t_max)
            const key_t tMaxNodes = static_cast< key_t >( std::sqrt( static_cast< double >( std::numeric_limits< key_t >::max() ) ) );
            BELFEM_ERROR( mNumberOfNodes <= tMaxNodes,
                "Mesh has too many nodes (%llu) for edge key computation with current key_t type.\n"
                "Maximum allowed nodes: %llu\n"
                "Recommendation: Change typedef for key_t in src/core/typedefs.hpp to a larger unsigned integer type (e.g., __uint128_t)",
                ( long long unsigned int ) mNumberOfNodes,
                ( long long unsigned int ) tMaxNodes );
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::
        create_edges(
                const Vector< id_t > aNedelecBlocks,
                const Vector< id_t > aNedelecSideSets,
                const bool aCreateEdgesOnAllSideSets )
        {
            // start timer
            Timer tTimer;

            if ( mCommRank == 0 )
            {
                mMesh.update_node_indices() ;

                message( InfoLevel::Default, "    Creating edges ...");

                // collect the elements from these blocks
                Cell< Element * > tElements ;

                this->collect_elements( aNedelecBlocks, aNedelecSideSets, aCreateEdgesOnAllSideSets );

                // create the edge IDs
                Vector< key_t > tEdgeKeys ;
                this->create_edge_keys( tEdgeKeys );

                this->create_edges_on_master( tEdgeKeys );

                this->set_edge_ids();

                this->link_elements_to_edges() ;

                //std::cout << "#warning: curves are not linked to edges" << std::endl ;
                //this->link_curves_to_edges() ;

                this->compute_edge_ownerships() ;

                mMesh.finalize_edges( mElements ) ;

                message( InfoLevel::Detailed, "    ... number of edges                 : %lu",
                         ( long unsigned int ) mMesh.number_of_edges() );

                message( InfoLevel::Detailed, "    ... time for creating edges         : %u ms\n",
                         ( unsigned int ) tTimer.stop() );
            }

        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::get_all_block_ids( Vector< id_t > & aBlockIDs )
        {

            // grab all ids from this block
            Cell< Block * > & tBlocks = mMesh.blocks();

            // ids for the blocks
            aBlockIDs.set_size( tBlocks.size() );

            // initialize counter
            uint tCount = 0 ;

            // collect the ids
            for ( Block * tBlock : tBlocks )
            {
                aBlockIDs( tCount++ ) = tBlock->id() ;
            }
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::collect_elements(
                const Vector< id_t > & aBlockIDs,
                const Vector< id_t > & aSideSetIDs,
                const bool aCreateEdgesOnAllSideSets )
        {

            mMesh.unflag_all_elements() ;
            mMesh.unflag_all_facets() ;

            // reset element order
            mElementOrder = 0 ;

            // count number of elements
            if( aBlockIDs.length() == 0 && aSideSetIDs.length() == 0 )
            {
                for( Block * tBlock : mMesh.blocks() )
                {
                    tBlock->flag_elements() ;

                    // get order
                    uint tOrder = interpolation_order_numeric( tBlock->element_type() );

                    mElementOrder = tOrder > mElementOrder ?
                                    tOrder : mElementOrder ;

                }
            }
            else
            {
                for ( id_t tID: aBlockIDs )
                {
                    mMesh.block( tID )->flag_elements();

                    // get order
                    uint tOrder = interpolation_order_numeric( mMesh.block( tID )->element_type());

                    mElementOrder = tOrder > mElementOrder ?
                                    tOrder : mElementOrder;

                    mMesh.block( tID )->set_faces_flag( true );
                }
            }

            if( aSideSetIDs.length() == 0 && aCreateEdgesOnAllSideSets )
            {
                for( SideSet * tSideSet : mMesh.sidesets() )
                {
                    Cell< Facet * > & tFacets = tSideSet->facets();
                    for( Facet * tFacet : tFacets )
                    {
                        if( tFacet->master()->is_flagged() )
                        {
                            tFacet->flag();

                            // get order
                            uint tOrder = interpolation_order_numeric( tSideSet->element_type() );

                            mElementOrder = tOrder > mElementOrder ?
                                            tOrder : mElementOrder ;
                        }
                    }
                }
            }
            else
            {
                for ( id_t tID: aSideSetIDs )
                {
                    Cell< Facet * > & tFacets = mMesh.sideset( tID )->facets();

                    for ( Facet * tFacet: tFacets )
                    {
                        tFacet->flag();
                    }

                    if ( tFacets.size() > 0 )
                    {
                        // get order
                        uint tOrder = interpolation_order_numeric( mMesh.sideset( tID )->element_type());

                        mElementOrder = tOrder > mElementOrder ?
                                        tOrder : mElementOrder;
                    }
                }
            }

            // initialize counter
            index_t tCount = 0 ;

            // count flagged elements
            Cell< Element * > & tAllElements = mMesh.elements() ;
            for( Element * tElement : tAllElements )
            {
                if( tElement->is_flagged() )
                {
                    ++tCount ;
                }
            }

            // count flagged facets
            Cell< Facet * > & tAllFacets = mMesh.facets() ;
            for( Facet * tFacet : tAllFacets )
            {
                if( tFacet->is_flagged() )
                {
                    ++tCount ;
                }
            }

            // allocate memory
            mElements.set_size( tCount, nullptr );

            // reset counter
            tCount = 0 ;

            for( Element * tElement : tAllElements )
            {
                if( tElement->is_flagged() )
                {
                    mElements( tCount++ ) = tElement;
                }
            }

            // add remaining facets
            for( Facet* tFacet : tAllFacets )
            {
                if( tFacet->is_flagged() )
                {
                    // increment the counter
                    mElements( tCount++ ) = tFacet->element();
                }
            }

            // fix node flipping on edges if this happens
            BELFEM_ERROR( mElementOrder < 3,
                         "Element order not supported. Need to add Node flipping routine." );
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::create_edge_keys( Vector< key_t >    & aKeys )
        {
            // count maximum number of edges
            index_t tCount = 0 ;

            // loop over all elements
            for ( Element * tElement : mElements )
            {
                tCount += tElement->number_of_edges();
            }

            // allocate memory
            aKeys.set_size( tCount );

            // reset counter
            tCount = 0 ;

            // work array for nodes
            Cell< Node * > tNodes ;

            // loop over all elements
            for ( Element * tElement : mElements )
            {
                uint tNumEdges = tElement->number_of_edges();
                for( uint e=0; e<tNumEdges; ++e )
                {
                    aKeys( tCount++ ) = this->edge_key( tElement, e, tNodes );
                }

            }

            // make edges unique
            unique( aKeys );
        }

//---------------------------------------------------------------------------

        key_t
        EdgeFactory::edge_key(
                Element           * aElement,
                const uint          aEdgeIndex,
                Cell< Node * >    & aNodes )
        {
            // grab nodes from edge
            aElement->get_nodes_of_edge( aEdgeIndex, aNodes );

            // grab indices from corner nodes
            key_t tA = aNodes( 0 )->index();
            key_t tB = aNodes( 1 )->index();

            if( tA > tB )
            {
                return tA * mNumberOfNodes + tB ;
            }
            else
            {
                return tB * mNumberOfNodes + tA ;
            }
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::create_edges_on_master(
                      const Vector< key_t >   & aKeys )
        {
            Cell< Edge * > & tEdges = mMesh.edges() ;

            BELFEM_ASSERT( tEdges.size() == 0, "Edges of mesh have already been created");

            tEdges.set_size( aKeys.length(), nullptr );

            // reset map
            mMap.clear() ;

            // edge counter
            index_t tCount = 0 ;

            // container for nodes
            Cell< Node * > tNodes ;


            for( key_t tKey : aKeys )
            {
                // create a new edge
                Edge * tEdge = new Edge() ;

                tEdge->set_index( tCount );

                // grab nodes from edge
                this->grab_nodes( tKey, tNodes );

                // get number of nodes from edge
                uint tNumNodes = tNodes.size() ;

                tEdge->allocate_node_container( tNumNodes );

                // link edge to nodes
                for( uint k=0; k< tNumNodes; ++k )
                {
                    tEdge->insert_node( tNodes( k ), k );
                }

                // add edge to map
                mMap[ tKey ] = tEdge ;

                // add edge to array
                tEdges( tCount++ ) = tEdge ;
            }
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::link_elements_to_edges()
        {
            Cell< Node * > tNodes ;

            for ( Element * tElement : mElements )
            {
                uint tNumEdges = tElement->number_of_edges();

                tElement->allocate_edge_container();
                for ( index_t k = 0; k < tNumEdges; ++k )
                {
                    // grab edge from map and insert
                    tElement->insert_edge( mMap( this->edge_key(
                            tElement, k, tNodes )), k );

                }
            }
        }

        void
        EdgeFactory::link_curves_to_edges()
        {
            for ( Curve * tCurve : mMesh.curves() )
            {
                Cell< Edge * > & tEdges = tCurve->edges();
                tEdges.set_size( tCurve->segments().size(), nullptr );

                index_t tCount = 0 ;

                for ( Segment * tSegment : tCurve->segments() )
                {
                    key_t tA = tSegment->node( 0 )->index();
                    key_t tB = tSegment->node( 1 )->index();

                    key_t tKey = tA > tB ? tA * mNumberOfNodes + tB
                        : tB * mNumberOfNodes + tA ;

                    // find the edge
                    Edge * tEdge = mMap( tKey ) ;

                    // add edge to container
                    tEdges( tCount++ ) = tEdge ;

                    // link edge to segment
                    tSegment->insert_edge(  tEdge );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::compute_edge_ownerships()
        {
            // grab container
            Cell< Edge * > & tEdges = mMesh.edges() ;

            if( comm_size() == 1 )
            {
                for( Edge * tEdge : tEdges )
                {
                    tEdge->set_owner( tEdge->node( 0 )->owner() );
                }
            }
            else
            {
                proc_t tOwnerA ;
                proc_t tOwnerB ;
                for( Edge * tEdge : tEdges )
                {
                    // grab owners from nodes
                    tOwnerA = tEdge->node( 0 )->owner() ;
                    tOwnerB = tEdge->node( 1 )->owner() ;

                    // set owner to smaller one
                    tEdge->set_owner( tOwnerA < tOwnerB ? tOwnerA : tOwnerB );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::grab_nodes( const key_t aKey, Cell< Node* > & aNodes )
        {
            // compute first and second node ids
            index_t tNodeA = aKey % mNumberOfNodes ;
            index_t tNodeB = ( aKey - tNodeA ) / mNumberOfNodes ;

            Cell< Node * > & tNodes = mMesh.nodes() ;

            if( mElementOrder == 1 )
            {
                // trival
                aNodes.set_size( 2, nullptr );
                aNodes( 0 ) = tNodes( tNodeA );
                aNodes( 1 ) = tNodes( tNodeB );
            }
            else
            {
                aNodes.set_size( 3, nullptr );
                aNodes( 0 ) = tNodes( tNodeA );
                aNodes( 1 ) = tNodes( tNodeB );

                Node * tNode = tNodes( tNodeB );
                Cell< Node * > tEdgeNodes( 3, nullptr );

                for ( uint e=0; e<tNode->number_of_elements(); ++e )
                {
                    Element * tElement = tNode->element( e );
                    for ( uint d=0; d<tElement->number_of_edges(); ++d )
                    {
                        tElement->get_nodes_of_edge( d, tEdgeNodes );
                        if (  ( tEdgeNodes( 0 )->index() == tNodeA && tEdgeNodes( 1 )->index() == tNodeB )
                            || ( tEdgeNodes( 1 )->index() == tNodeA && tEdgeNodes( 0 )->index() == tNodeB ) )
                        {
                            aNodes( 2 ) = tEdgeNodes( 2 );
                            return ;
                        }
                    }
                }
                // we should never end up here!
                BELFEM_ERROR( false, "failed to find the correct edge for %lu - %lu " ,
                    ( long unsigned int ) tNodeA , ( long unsigned int ) tNodeB);
            }
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::print()
        {
            if ( mCommRank == 0 )
            {
                Cell< Edge * >    & tEdges = mMesh.edges();

                std::cout << "EDGES:" << std::endl ;

                for( Edge * tEdge : tEdges )
                {
                    std::cout << tEdge->id() << " :" ;

                    for( uint k=0; k<tEdge->number_of_nodes(); ++k )
                    {
                        std::cout << " " << tEdge->node( k )->id() ;
                    }
                    std::cout << std::endl ;
                }
            }
        }


//------------------------------------------------------------------------------

        void
        EdgeFactory::set_edge_ids()
        {
            switch( mMesh.number_of_dimensions() )
            {
                case( 2 ) :
                {
                    for( mesh::Facet * tFacet : mMesh.facets() )
                    {
                        // // compute key
                        key_t tA = tFacet->element()->node( 0 )->index();
                        key_t tB =tFacet->element()->node( 1 )->index();
                        key_t tKey = tA > tB ? tA * mNumberOfNodes + tB
                                             :  tB * mNumberOfNodes + tA ;

                        // get facet if it has been created
                        if( mMap.key_exists( tKey ) )
                        {
                            mesh::Edge * tEdge = mMap[ tKey ];
                            tEdge->set_id( tFacet->id() );
                            tEdge->flag();
                        }
                    }
                    break ;
                }
                case( 3 ) :
                {
                    for( mesh::Element * tBoundaryEdge : mMesh.boundary_edges() )
                    {
                        // // compute key
                        key_t tA = tBoundaryEdge->node( 0 )->index();
                        key_t tB = tBoundaryEdge->node( 1 )->index();
                        key_t tKey = tA > tB ? tA * mNumberOfNodes + tB
                                             :  tB * mNumberOfNodes + tA ;

                        // get edge if it has been created
                        if( mMap.key_exists( tKey ) )
                        {
                            mesh::Edge * tEdge = mMap[ tKey ];
                            tEdge->set_id( tBoundaryEdge->id() );
                            tEdge->flag() ;
                        }
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid mesh dimension");
                }
            }

            // compute the maximum ID so far
            id_t tMaxID = 0 ;
            for( mesh::Element * tEdge : mMesh.boundary_edges() )
            {
                tMaxID = tMaxID < tEdge->id() ? tEdge->id() : tMaxID ;
            }
            for( mesh::Facet * tFacet : mMesh.facets() )
            {
                tMaxID = tMaxID < tFacet->id() ? tFacet->id() : tMaxID ;
            }
            for( mesh::Face * tFace : mMesh.faces() )
            {
                tMaxID = tMaxID < tFace->id() ? tFace->id() : tMaxID ;
            }
            for( mesh::Element * tElement : mMesh.elements() )
            {
                tMaxID = tMaxID < tElement->id() ? tElement->id() : tMaxID ;
            }

            Cell< Edge * >    & tEdges = mMesh.edges();

            for( Edge * tEdge : tEdges )
            {
                if( ! tEdge->is_flagged() )
                {
                    tEdge->set_id( ++tMaxID );
                }
            }

            // resort the array
            sort( tEdges, opVertexID );

        }

//------------------------------------------------------------------------------
    }
}
