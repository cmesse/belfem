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

        TapeRoller::TapeRoller( Mesh * aMesh, const uint aNumberOfLayers, const uint aElementOrder ) :
                mMyRank( comm_rank() ),
                mMesh( aMesh ),
                mNumberOfBlocks( aNumberOfLayers ),
                mElementOrder( aElementOrder ),
                mNumberOfGhostLayers( aNumberOfLayers * aElementOrder + 1 )
        {

        }

//------------------------------------------------------------------------------

        TapeRoller::~TapeRoller()
        {
            for ( Layer * tLayer: mGhostLayers )
            {
                delete tLayer;
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

        void
        TapeRoller::get_sidesets( Vector< id_t > & aSideSetIDs )
        {

            // we must convert the cell to a vector
            aSideSetIDs.set_size( mSelectedSideSets.size() );
            index_t tCount = 0 ;
            for( id_t s : mSelectedSideSets )
            {
                aSideSetIDs( tCount++ ) = s ;
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::add_master_blocks( const Vector< id_t > & aMasterBlockIDs )
        {
            if ( comm_rank() == 0 )
            {
                for ( id_t tID: aMasterBlockIDs )
                {
                    mMasterBlocks.push( tID );
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
                mGhostSideSetIDs.set_size( mNumberOfGhostLayers );

                for ( uint l = 0; l < mNumberOfGhostLayers; ++l )
                {
                    mGhostSideSetIDs( l ) = ++mMaxSideSetId;
                }

                // this command contains a finalize
                aMaxBlockID = mMesh->create_ghost_sidesets(
                        mGhostSideSetIDs,
                        mElementIDs,
                        mGhostLayers );

                // fix node connectivity, because the finalize command messes up with the linking
                this->fix_tape_to_node_connectivities( aMaxBlockID );
            }

            return aMaxBlockID ;
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::fix_tape_to_node_connectivities( const id_t aMaxBlockID )
        {
            // get number of elements per layer
            index_t tNumElements = mElementIDs.length() ;

            // loop over all laters
            for ( uint l = 0; l < mNumberOfGhostLayers; ++l )
            {
                uint tNumNodes = mMesh->facet(
                        mElementIDs( 0 ) )->element()->number_of_nodes() ;

                Cell< Node * > & tLayerNodes = mGhostLayers( l )->Nodes ;

                // loop over all elements
                for( index_t e=0; e<tNumElements; ++e )
                {
                    // grab the original element
                    Element * tElement = mMesh->facet(
                            mElementIDs( e ) )->element() ;

                    // grab the element from the block
                    Element * tClone = mMesh->ghost_facet(
                            mElementIDs( e ), l )->element() ;

                    // relink nodes
                    for( uint k=0; k<tNumNodes; ++k )
                    {
                        tClone->insert_node(
                                tLayerNodes( mNodeMap( tElement->node( k )->id())), k );
                    }

                    // this will be used later on to restore the element ownership
                    tClone->reset_element_container();
                    tClone->allocate_element_container( 1 );
                    tClone->insert_element( tElement );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::shift_nodes()
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // step 1: create backup of original master and slave setup
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        }

//------------------------------------------------------------------------------

        void
        TapeRoller::flip_element_orientation()
        {
            // unflag all elements on the mesh
            mMesh->unflag_all_elements();

            uint tNumSidesets = mSelectedSideSets.size();

            // count memory
            index_t tCount = 0;

            for ( uint s = 0; s < tNumSidesets; ++s )
            {

                // grab the facets
                Cell< Facet * > & tFacets
                        = mMesh->sideset( mSelectedSideSets( s ))->facets();

                index_t tNumNodes = mesh::number_of_nodes(
                        mMesh->sideset( mSelectedSideSets( s ) )->element_type());

                tCount += tNumNodes * tFacets.size() ;

            }

            // allocate backup container
            mOrientationBackup.set_size( tCount, nullptr );

            // reset counter
            tCount = 0;

            for ( uint s = 0; s < tNumSidesets; ++s )
            {
                // flag the master block
                mMesh->block( mMasterBlocks( s ))->flag_elements();

                // grab the facets
                Cell< Facet * > & tFacets
                        = mMesh->sideset( mSelectedSideSets( s ))->facets();

                // get the number of nodes per facet
                uint tNumNodes = mesh::number_of_nodes(
                        mMesh->sideset( mSelectedSideSets( s ) )->element_type());

                // loop over all facets
                for ( Facet * tFacet: tFacets )
                {
                    // backup nodes
                    for ( uint k = 0; k < tNumNodes; ++k )
                    {
                        mOrientationBackup( tCount++ ) = tFacet->element()->node( k );
                    }

                    // check if we need to swap
                    if ( tFacet->slave()->is_flagged() )
                    {
                        // backup element
                        Element * tMaster = tFacet->slave();
                        uint tMasterIndex = tFacet->slave_index();

                        // swap
                        tFacet->set_slave( tFacet->master(), tFacet->master_index() );
                        tFacet->set_master( tMaster, tMasterIndex );
                        tFacet->flag();
                    }
                    else
                    {
                        // relink nodes
                        tFacet->set_master( tFacet->master(), tFacet->master_index()  );
                        tFacet->unflag();
                    }
                }

                // flag the master block
                mMesh->block( mMasterBlocks( s ))->unflag_elements();
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::revert_element_orientation()
        {
            uint tNumSidesets = mSelectedSideSets.size() ;

            index_t tCount = 0 ;

            // grab the normal fields
            Vector< real > & tNx = mMesh->field_data("SurfaceNormalsx");
            Vector< real > & tNy = mMesh->field_data("SurfaceNormalsy");
            Vector< real > & tNz = mMesh->field_data("SurfaceNormalsz");

            // loop over all facets
            for ( uint s = 0; s < tNumSidesets; ++s )
            {
                // flag the master block
                mMesh->block( mMasterBlocks( s ))->flag_elements();

                // grab the facets
                Cell< Facet * > & tFacets
                        = mMesh->sideset( mSelectedSideSets( s ))->facets();

                // get the number of nodes per facet
                uint tNumNodes = mesh::number_of_nodes(
                        mMesh->sideset( mSelectedSideSets( s ))->element_type());

                Vector< index_t > tIndices( tNumNodes );

                for ( Facet * tFacet: tFacets )
                {
                    // check if we need to swap
                    if ( tFacet->slave()->is_flagged() )
                    {
                        // backup element
                        Element * tMaster = tFacet->slave();
                        uint tMasterIndex = tFacet->slave_index();

                        // swap
                        tFacet->set_slave( tFacet->master(), tFacet->master_index());
                        tFacet->set_master( tMaster, tMasterIndex, false );

                        // grab indices in reverse order
                        uint c = 0 ;
                        for( int k=tNumNodes-1; k>=0; k-- )
                        {
                            tIndices( c++ ) = tFacet->element()->node( k )->index() ;
                        }

                        // relink nodes
                        for ( uint k = 0; k < tNumNodes; ++k )
                        {
                            tFacet->element()->insert_node( mOrientationBackup( tCount++ ), k );
                        }

                        // fix normal data
                        for( uint k=0; k<tNumNodes; ++k )
                        {
                            index_t i = tFacet->element()->node( k )->index() ;
                            index_t j = tIndices( k );

                            tNx( i ) = tNx( j );
                            tNy( i ) = tNy( j );
                            tNz( i ) = tNz( j );
                        }

                        tFacet->unflag() ;
                    }
                    else
                    {
                        // grab indices in normal
                        for( uint k = 0; k < tNumNodes; ++k )
                        {
                            tIndices( k ) = tFacet->element()->node( k )->index() ;
                        }

                        // relink nodes
                        for ( uint k = 0; k < tNumNodes; ++k )
                        {
                            tFacet->element()->insert_node( mOrientationBackup( tCount++ ), k );
                        }

                        // fix normal data
                        for( uint k=0; k<tNumNodes; ++k )
                        {
                            index_t i = tFacet->element()->node( k )->index() ;
                            index_t j = tIndices( k );

                            tNx( i ) = tNx( j );
                            tNy( i ) = tNy( j );
                            tNz( i ) = tNz( j );
                        }

                    }
                }
            }
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
            mGhostLayers.set_size( mNumberOfGhostLayers, nullptr );
            for ( uint t = 0; t < mNumberOfGhostLayers; ++t )
            {
                mGhostLayers( t ) = new Layer;
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
            for ( uint l = 0; l < mNumberOfGhostLayers; ++l )
            {
                // grab object
                Cell< Node * > & tLayerNodes = mGhostLayers( l )->Nodes;

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
            for ( uint l = 0; l < mNumberOfGhostLayers; ++l )
            {
                // grab node container
                Cell< Node * > & tLayerNodes = mGhostLayers( l )->Nodes;

                // grab element container
                Cell< Element * > & tLayerElements = mGhostLayers( l )->Elements;

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
                for( uint l=0; l<mNumberOfGhostLayers; ++l )
                {
                    // get the layer containers
                    Cell< Node * >    & tLayerNodes    = mGhostLayers( l )->Nodes ;
                    Cell< Edge * >    & tLayerEdges    = mGhostLayers( l )->Edges ;
                    Cell< Element * > & tLayerElements = mGhostLayers( l )->Elements ;

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
                for( uint l=0; l<mNumberOfGhostLayers; ++l )
                {
                    // get the layer containers
                    Cell< Face * > & tLayerFaces = mGhostLayers( l )->Faces;
                    Cell< Element * > & tLayerElements = mGhostLayers( l )->Elements ;

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