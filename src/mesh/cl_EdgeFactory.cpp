//
// Created by christian on 6/9/21.
//
#include "commtools.hpp"
#include "cl_EdgeFactory.hpp"
#include "cl_Element.hpp"
#include "meshtools.hpp"
#include "fn_unique.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "meshtools.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        EdgeFactory::EdgeFactory( Mesh & aMesh ) :
                mRank( comm_rank() ),
                mMesh( aMesh ),
                mNumberOfNodes( aMesh.number_of_nodes() )
        {

        }

//------------------------------------------------------------------------------

        EdgeFactory::EdgeFactory( Mesh * aMesh ) :
            mRank( comm_rank() ),
            mMesh( * aMesh ),
            mNumberOfNodes( aMesh->number_of_nodes() )
        {

        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::
        create_edges(
                const Vector< id_t > aNedelecBlocks,
                const Vector< id_t > aNedelecSideSets )
        {
            // start timer
            Timer tTimer;

            if ( mRank == 0 )
            {
                mMesh.update_node_indices() ;

                Vector< id_t > tBlockIDs ;
                if( aNedelecBlocks.length() == 0 && aNedelecSideSets.length() == 0 )
                {
                    // select all blocks
                    this->get_all_block_ids( tBlockIDs );
                }
                else
                {
                    // select given blocks
                    tBlockIDs = aNedelecBlocks ;
                }

                message( 2, "Creating edges ...");

                // collect the elements from these blocks
                Cell< Element * > tElements ;

                this->collect_elements( tBlockIDs, aNedelecSideSets );

                // create the edge IDs
                Vector< luint > tEdgeKeys ;
                this->create_edge_keys( tEdgeKeys );

                this->create_edges_on_master( tEdgeKeys );

                this->set_edge_ids();

                this->link_elements_to_edges() ;

                this->compute_edge_ownerships() ;

                mMesh.finalize_edges( mElements ) ;

                message( 2, "    ... number of edges                 : %lu",
                         ( long unsigned int ) mMesh.number_of_edges() );

                message( 3, "    ... time for creating edges         : %u ms",
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
                const Vector< id_t > & aSideSetIDs )
        {
            mMesh.unflag_all_elements() ;
            mMesh.unflag_all_facets() ;

            // reset element order
            mElementOrder = 0 ;

            // count number of elements
            for( id_t tID : aBlockIDs )
            {
                mMesh.block( tID )->flag_elements() ;

                // get order
                uint tOrder = interpolation_order_numeric( mMesh.block( tID )->element_type() );

                mElementOrder = tOrder > mElementOrder ?
                                tOrder : mElementOrder ;
            }
            for( id_t tID : aSideSetIDs )
            {
                Cell< Facet * > & tFacets = mMesh.sideset( tID )->facets();

                for( Facet * tFacet : tFacets )
                {
                    tFacet->flag() ;
                }

                if( tFacets.size() > 0 )
                {
                    // get order
                    uint tOrder = interpolation_order_numeric( mMesh.sideset( tID )->element_type() );

                    mElementOrder = tOrder > mElementOrder ?
                                    tOrder : mElementOrder ;
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

            // add elements from blocks
            for( id_t tID : aBlockIDs )
            {
                // grab element container
                Cell< Element * > & tElements
                        = mMesh.block( tID )->elements();

                // add elements
                for ( Element * tElement : tElements )
                {
                    // increment the counter
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
                         "Element order not supported. Need to fix Node flipping routine." );
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::create_edge_keys( Vector< luint >    & aKeys )
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

        luint
        EdgeFactory::edge_key(
                Element           * aElement,
                const uint          aEdgeIndex,
                Cell< Node * >    & aNodes )
        {
            // grab nodes from edge
            aElement->get_nodes_of_edge( aEdgeIndex, aNodes );

            // grab indices from corner nodes
            luint tA = aNodes( 0 )->index();
            luint tB = aNodes( 1 )->index();

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
                      const Vector< luint >   & aKeys )
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


            for( luint tKey : aKeys )
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
                    tElement->insert_edge( mMap[ this->edge_key(
                            tElement, k, tNodes ) ], k );

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
        EdgeFactory::grab_nodes( const luint & aKey, Cell< Node* > & aNodes )
        {
            // compute first and second node ids
            id_t tNodeA = aKey % mNumberOfNodes ;
            id_t tNodeB = ( aKey - tNodeA ) / mNumberOfNodes ;

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
                // get first node
                Node * tNode = tNodes( tNodeA );

                Element * tElement = nullptr;

                // get number of elements for this node
                uint tNumElements = tNode->number_of_elements();

                // flag that tells if node has been found
                bool tFound = false;

                // find element with which code is shared
                for ( uint e = 0; e < tNumElements; ++e )
                {
                    // get pointer to element
                    tElement = tNode->element( e );

                    // get number of nodes
                    uint tNumNodes = tElement->number_of_nodes();

                    for ( uint k = 0; k < tNumNodes; ++k )
                    {
                        if ( tElement->node( k )->index() == tNodeB )
                        {
                            // set found flag and break
                            tFound = true;
                            break;
                        }
                    }

                    // check if element has been found
                    if ( tFound )
                    {
                        // break the loop
                        break;
                    }
                }

                BELFEM_ASSERT( tFound, "Something went wrong in edge generation" );

                // find correct edge
                uint tNumEdges = tElement->number_of_edges();
                for ( uint e = 0; e < tNumEdges; ++e )
                {
                    // grab nodes from edge
                    tElement->get_nodes_of_edge( e, aNodes );

                    // check if this is the correct edge
                    if (    aNodes( 0 )->index() == tNodeA
                         && aNodes( 1 )->index() == tNodeB )
                    {
                        // exit the function
                        return;
                    }
                    else if (    aNodes( 1 )->index() == tNodeA
                              && aNodes( 0 )->index() == tNodeB )
                    {
                        // flip nodes 0 and 1
                        aNodes( 0 ) = tNodes( tNodeA );
                        aNodes( 1 ) = tNodes( tNodeB );

                        BELFEM_ASSERT( aNodes.size() <= 3, "Flipping of third order edges is not implemented" );

                        // exit the function
                        return;;
                    }
                }

                // we should never end up here!
                BELFEM_ERROR( false, "failed to find the correct edge" );
            }
        }

//------------------------------------------------------------------------------

        void
        EdgeFactory::print()
        {
            if ( mRank == 0 )
            {
                Cell< Edge * >    & tEdges = mMesh.edges();
                Cell< Element * > & tElements = mMesh.elements();

                std::cout << "Edges:" << std::endl ;

                for( Edge * tEdge : tEdges )
                {
                    std::cout << tEdge->id() << " :" ;

                    for( uint k=0; k<tEdge->number_of_nodes(); ++k )
                    {
                        std::cout << " " << tEdge->node( k )->id() ;
                    }
                    std::cout << std::endl ;
                }

                std::cout << "Elements:" << std::endl ;

                for( Element * tElement : tElements )
                {
                    if( tElement->has_edges() )
                    {
                        std::cout << tElement->id() << " :" ;
                        for( uint k=0; k<tElement->number_of_edges(); ++k )
                        {
                            std::cout << " " << tElement->edge( k )->id() ;
                        }
                        std::cout << std::endl ;
                    }
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
                        luint tA = tFacet->element()->node( 0 )->index();
                        luint tB =tFacet->element()->node( 1 )->index();
                        luint tKey = tA > tB ? tA * mNumberOfNodes + tB
                                             :  tB * mNumberOfNodes + tA ;

                        // get facet if it has been created
                        if( mMap.key_exists( tKey ) )
                        {
                            mesh::Edge * tEdge = mMap[ tKey ];
                            tEdge->set_id( tFacet->id());
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
                        luint tA = tBoundaryEdge->node( 0 )->index();
                        luint tB = tBoundaryEdge->node( 1 )->index();
                        luint tKey = tA > tB ? tA * mNumberOfNodes + tB
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
        }

//------------------------------------------------------------------------------
    }
}
