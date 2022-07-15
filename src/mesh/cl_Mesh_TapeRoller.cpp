//
// Created by christian on 3/30/22.
//
#include "commtools.hpp"
#include "cl_Mesh_TapeRoller.hpp"
#include "cl_Element_Factory.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        TapeRoller::TapeRoller( Mesh * aMesh, const uint aNumberOfLayers ) :
                mMyRank( comm_rank()),
                mMesh( aMesh ),
                mNumberOfLayers( aNumberOfLayers )
        {

        }

//------------------------------------------------------------------------------

        TapeRoller::~TapeRoller()
        {
            for ( Layer * tLayer: mLayers )
            {
                delete tLayer;
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::add_sideset( const id_t aSideSetID )
        {
            if ( comm_rank() == 0 )
            {
                mSelectedSideSets.push( aSideSetID );
            }

        }

//------------------------------------------------------------------------------

        void
        TapeRoller::add_sidesets( const Vector< id_t > & aSideSetIDs )
        {
            if ( comm_rank() == 0 )
            {
                for ( id_t tID: aSideSetIDs )
                {
                    mSelectedSideSets.push( tID );
                }
            }
        }

//------------------------------------------------------------------------------

        id_t
        TapeRoller::run()
        {
            id_t aMaxBlockID = 0 ;

            if ( comm_rank() == 0 )
            {
                // make sidesets unique
                unique( mSelectedSideSets );

                mMesh->unflag_all_nodes();
                mMesh->unflag_all_edges();
                mMesh->unflag_all_faces();
                mMesh->unflag_all_elements();

                // determine max ids for mesh entities
                this->get_max_ids();
                this->create_layers();
                this->clone_nodes();
                this->clone_elements();
                this->clone_edges();
                this->clone_faces();

                // create the list with the sideset ids
                mGhostSideSetIDs.set_size( mNumberOfLayers );

                for ( uint l = 0; l < mNumberOfLayers; ++l )
                {
                    mGhostSideSetIDs( l ) = ++mMaxSideSetId;
                }

                aMaxBlockID = mMesh->create_ghost_sidesets( mGhostSideSetIDs, mElementIDs, mLayers );

            }

            return aMaxBlockID ;
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::get_max_ids()
        {
            // Node IDs
            mMaxNodeId = 0;
            Cell< Node * > & tNodes = mMesh->nodes();
            for ( Node * tNode: tNodes )
            {
                mMaxNodeId = tNode->id() > mMaxNodeId ? tNode->id() : mMaxNodeId;
            }

            // Edge IDs
            mMaxEdgeId = 0;
            if ( mMesh->edges_exist())
            {
                Cell< Edge * > & tEdges = mMesh->edges();
                for ( Edge * tEdge: tEdges )
                {
                    mMaxEdgeId = tEdge->id() > mMaxEdgeId ? tEdge->id() : mMaxEdgeId;
                }
            }

            // face IDs
            mMaxFaceId = 0;
            if ( mMesh->faces_exist( ))
            {
                Cell< Face * > & tFaces = mMesh->faces();
                for ( Face * tFace: tFaces )
                {
                    mMaxFaceId = tFace->id() > mMaxFaceId ? tFace->id() : mMaxFaceId;
                }
            }

            // element IDs
            mMaxElementID = 0;
            Cell< Element * > & tElements = mMesh->elements();
            for ( Element * tElement: tElements )
            {
                mMaxElementID = tElement->id() > mMaxElementID ? tElement->id() : mMaxElementID;
            }
            for( mesh::Facet * tConnector : mMesh->connectors() )
            {
                mMaxElementID = tConnector->id() > mMaxElementID ? tConnector->id() : mMaxElementID;
            }

            // compute maximum sideset id (including cuts)
            mMaxSideSetId = 0;

            for ( mesh::SideSet * tSideSet: mMesh->sidesets() )
            {
                mMaxSideSetId = tSideSet->id() > mMaxSideSetId ? tSideSet->id() : mMaxSideSetId;
            }
            for ( mesh::SideSet * tSideSet: mMesh->cuts() )
            {
                mMaxSideSetId = tSideSet->id() > mMaxSideSetId ? tSideSet->id() : mMaxSideSetId;
            }
        }

//-----------------------------------------------------------------------------

        void
        TapeRoller::create_layers()
        {
            mLayers.set_size( mNumberOfLayers, nullptr );
            for ( uint t = 0; t < mNumberOfLayers; ++t )
            {
                mLayers( t ) = new Layer;
            }
        }

//-----------------------------------------------------------------------------

        void
        TapeRoller::clone_nodes()
        {
            mMesh->unflag_all_nodes() ;

            for ( id_t tID: mSelectedSideSets )
            {
                mMesh->sideset( tID )->flag_all_nodes();
            }

            // count nodes
            index_t tCount = 0;

            Cell< Node * > & tNodes = mMesh->nodes();

            // reset the map
            mNodeMap.clear();

            // the node map is needed to assign the nodes correctly
            // with the new tapes
            for ( Node * tNode: tNodes )
            {
                if ( tNode->is_flagged() )
                {
                    mNodeMap[ tNode->id() ] = tCount++;
                }
            }

            // now we populate the nodes for each layer
            for ( uint l = 0; l < mNumberOfLayers; ++l )
            {
                // grab object
                Cell< Node * > & tLayerNodes = mLayers( l )->Nodes;

                // allocate the memory
                tLayerNodes.set_size( tCount, nullptr );

                // reset the counter
                tCount = 0;

                // loop over all nodes
                for ( Node * tNode: tNodes )
                {
                    if ( tNode->is_flagged())
                    {
                        // crate the new node
                        // create a new node
                        tLayerNodes( tCount++ ) = new Node(
                                ++mMaxNodeId,
                                tNode->x(),
                                tNode->y(),
                                tNode->z());
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        TapeRoller::clone_elements()
        {
            // count number of elements that must be cloned
            index_t tCount = 0;

            for ( id_t tID: mSelectedSideSets )
            {
                tCount += mMesh->sideset( tID )->number_of_facets();
            }

            // allocate id container
            mElementIDs.set_size( tCount );

            // reset the map
            mElementMap.clear();

            // reset the counter
            tCount = 0;

            // collect the element IDs and create the map
            for ( id_t tID: mSelectedSideSets )
            {
                // get the container with the facets
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets();
                {
                    // loop over all facets
                    for ( Facet * tFacet: tFacets )
                    {
                        mElementMap[ tFacet->id() ] = tCount;
                        mElementIDs( tCount++ ) = tFacet->id();
                    }
                }
            }

            // create the element factory
            ElementFactory tFactory;

            // next, we clone the elements
            for ( uint l = 0; l < mNumberOfLayers; ++l )
            {
                // grab node container
                Cell< Node * > & tLayerNodes = mLayers( l )->Nodes;

                // grab element container
                Cell< Element * > & tLayerElements = mLayers( l )->Elements;

                // allocate the memory
                tLayerElements.set_size( tCount, nullptr );

                // reset the counter
                tCount = 0;

                // loop over all sidesets
                for ( id_t tID: mSelectedSideSets )
                {
                    // get the container with the facets
                    Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets();

                    // get the type
                    ElementType tType = mMesh->sideset( tID )->element_type();

                    // number of nodes per element type
                    uint tNumNodes = number_of_nodes( tType );


                    // loop over all facets
                    for ( Facet * tFacet: tFacets )
                    {
                        // original element
                        Element * tElement = tFacet->element();

                        // create a new element
                        Element * tClone = tFactory.create_lagrange_element( tType, ++mMaxElementID );

                        // this is needed to get the ownerships right later on
                        tClone->allocate_element_container( 1 ) ;
                        tClone->insert_element( tElement );

                        // link the nodes
                        for ( uint k = 0; k < tNumNodes; ++k )
                        {
                            tClone->insert_node( tLayerNodes( mNodeMap( tElement->node( k )->id())), k );
                        }

                        // print
                        std::cout << "#elements " << tElement->id() << " - " << l << " : " << tClone->id() << " - " << tClone->node( 0 )->id() << " " << tClone->node( 1 )->id() << std::endl ;

                        tLayerElements( tCount++ ) = tClone;
                    } // end loop over all facets
                } // end loop over all sidesets
            } // end loop over all layers
        }

//--------------------------------------------------------------------------

        void
        TapeRoller::clone_edges()
        {
            if( mMesh->edges_exist() )
            {
                // flag edges that are to be cloned
                mMesh->unflag_all_edges() ;
                for ( id_t tID: mSelectedSideSets )
                {
                    // get the container with the facets
                    Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets();

                    // loop over all facets
                    for( Facet * tFacet : tFacets )
                    {
                        tFacet->element()->flag_edges() ;
                    }
                }

                // count number of edges
                index_t tCount = 0;

                // reset the map
                mEdgeMap.clear() ;

                // count cloned edges and create edge map
                Cell< Edge * > & tEdges = mMesh->edges() ;
                for( Edge * tEdge : tEdges )
                {
                    if( tEdge->is_flagged() )
                    {
                        mEdgeMap[ tEdge->id() ] = tCount++ ;
                    }
                }

                // clone the edges
                for( uint l=0; l<mNumberOfLayers; ++l )
                {
                    // get the layer containers
                    Cell< Node * >    & tLayerNodes    = mLayers( l )->Nodes ;
                    Cell< Edge * >    & tLayerEdges    = mLayers( l )->Edges ;
                    Cell< Element * > & tLayerElements = mLayers( l )->Elements ;

                    // allocate memory
                    tLayerEdges.set_size( tCount, nullptr );

                    // reset the counter
                    tCount = 0 ;

                    // clone edges
                    for( Edge * tEdge : tEdges )
                    {
                        if ( tEdge->is_flagged())
                        {
                            // crate a new edge
                            Edge * tClone = new Edge();
                            tClone->set_id( ++mMaxEdgeId );
                            tClone->allocate_node_container( tEdge->number_of_nodes());

                            // link nodes
                            for ( uint k = 0; k < tEdge->number_of_nodes(); ++k )
                            {
                                tClone->insert_node( tLayerNodes( mNodeMap( tEdge->node( k )->id())), k );
                            }

                            // add clone to list
                            tLayerEdges( tCount++ ) = tClone;
                        } // end loop over all edges

                    }

                    index_t e = 0 ;

                    // link new edges with elements
                    for( Element * tClone : tLayerElements )
                    {
                        // get original element
                        Element * tElement = mMesh->facet( mElementIDs( e++ ))->element() ;

                        tClone->allocate_edge_container();

                        // link element clone with new edges
                        for ( uint k = 0; k < tElement->number_of_edges(); ++k )
                        {
                            Edge * tEdge = tLayerEdges( mEdgeMap( tElement->edge( k )->id()));

                            tClone->insert_edge( tEdge, k );
                        }
                    } // end loop over layer elements

                } // end loop over all layers
            } // end edges exist
        }

//-----------------------------------------------------------------------------

        void
        TapeRoller::clone_faces()
        {
            if( mMesh->faces_exist() && mMesh->number_of_dimensions() == 3 )
            {
                // flag faces that are to be cloned
                mMesh->unflag_all_faces() ;
                for ( id_t tID: mSelectedSideSets )
                {
                    // get the container with the facets
                    Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets();

                    // loop over all facets
                    for( Facet * tFacet : tFacets )
                    {
                        tFacet->element()->flag_faces() ;
                    }
                }

                // count number of faces
                index_t tCount = 0;

                // grab container
                Cell< Face * > & tFaces = mMesh->faces() ;

                // clone the edges
                for( uint l=0; l<mNumberOfLayers; ++l )
                {
                    // get the layer containers
                    Cell< Face * > & tLayerFaces = mLayers( l )->Faces;
                    Cell< Element * > & tLayerElements = mLayers( l )->Elements ;

                    // allocate memory
                    tLayerFaces.set_size( tCount, nullptr );

                    // reset the counter
                    tCount = 0 ;

                    for ( Face * tFace: tFaces )
                    {

                        if ( tFace->is_flagged() )
                        {
                            Element * tParent = tLayerElements( mElementMap( tFace->master()->id() ) ) ;

                            Face * tClone = new Face( tParent );

                            tClone->set_id( ++mMaxFaceId );

                            // allocate memory
                            tParent->allocate_face_container();

                            // link face to element container
                            tParent->insert_face( tClone, 0 );

                            // add clone to list
                            tLayerFaces( tCount++ ) = tClone ;
                        }
                    } // end loop over all faces
                } // end loop over all layers
            } // end faces exist and model is 3D
        }

//-----------------------------------------------------------------------------
    }
}