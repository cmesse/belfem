//
// Created by christian on 3/30/22.
//
#include "commtools.hpp"
#include "cl_Mesh_TapeRoller.hpp"
#include "cl_Element_Factory.hpp"
#include "fn_sum.hpp"

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



                mMesh->unflag_everything();

                // determine max ids for mesh entities
                this->get_max_ids();
                this->create_layers();
                this->clone_nodes();
                this->clone_elements();
                this->clone_edges();
                this->clone_faces();

                // for visualization
                aMaxBlockID = this->create_tape_blocks() ;

                // create the list with the sideset ids
                mGhostSideSetIDs.set_size( mNumberOfGhostLayers );

                for ( uint l = 0; l < mNumberOfGhostLayers; ++l )
                {
                    mGhostSideSetIDs( l ) = ++mMaxSideSetId;
                }

                // this command contains a finalize
                mMesh->create_ghost_sidesets(
                        mGhostBlockIDs,
                        mGhostSideSetIDs,
                        mElementIDs,
                        mGhostLayers );

                // fix node connectivity, because the finalize command messes up with the linking
                this->fix_tape_to_node_connectivities();
            }

            return aMaxBlockID ;
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::fix_tape_to_node_connectivities()
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

        id_t
        TapeRoller::create_tape_blocks()
        {
            // compute the max block id
            id_t aMaxBlockID ;
            mMaxBlockID = 0 ;
            for( mesh::Block * tBlock : mMesh->blocks() )
            {
                mMaxBlockID = tBlock->id() > mMaxBlockID ?
                              tBlock->id() : mMaxBlockID ;
            }

            // save value for return
            aMaxBlockID = mMaxBlockID ;

            // number of elements per block
            index_t tNumElements = mElementIDs.length() ;

            // offset to bottom layer
            uint tLayerOffset = 0 ;

            // get element type
            ElementType tType = mGhostLayers( 0 )->Elements( 0 )->type() ;

            mGhostBlockIDs.set_size( mNumberOfBlocks );

            for( uint b=0; b<mNumberOfBlocks; ++b )
            {
                // create a new block
                Block * tBock = new mesh::Block( ++mMaxBlockID, tNumElements );
                mGhostBlockIDs( b ) = tBock->id() ;

                // populate the block
                switch( tType )
                {
                    case( ElementType::LINE2 ) :
                    {
                        this->create_tape_blocks_quad4( tBock, tLayerOffset );
                        break ;
                    }
                    case( ElementType::LINE3 ) :
                    {
                        this->create_tape_blocks_quad9( tBock, tLayerOffset );
                        break ;
                    }
                    case( ElementType::TRI3 ) :
                    {
                        this->create_tape_blocks_penta6( tBock, tLayerOffset );
                        break ;
                    }
                    case( ElementType::PENTA18 ) :
                    {
                        this->create_tape_blocks_penta18( tBock, tLayerOffset );
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid Element Type");
                    }
                }

                // add block to mesh
                mMesh->blocks().push( tBock );
            }

            // return the value
            return aMaxBlockID ;
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::create_tape_blocks_quad4(
                mesh::Block * aBlock, uint & aLayerOffset )
        {
            // Create an element factort
            ElementFactory tFactory ;

            // grab element container
            Cell< Element * > & tElements = aBlock->elements() ;

            // grab first layer
            Cell< Element * > & tLayer0 = mGhostLayers( aLayerOffset++ )->Elements ;

            // grab second layer
            Cell< Element * > & tLayer1 = mGhostLayers( aLayerOffset )->Elements ;

            // get number of elements
            index_t tNumElements = aBlock->number_of_elements() ;

            // loop over all elements
            for( index_t e=0; e<tNumElements; ++e )
            {
                // grab faces
                Element * tFacet0 = tLayer0( e );
                Element * tFacet1 = tLayer1( e );

                // create a new element
                Element * tElement = tFactory.create_lagrange_element( ElementType::QUAD4, ++mMaxElementID );

                // link nodes

                tElement->insert_node( tFacet0->node( 0 ), 0 );
                tElement->insert_node( tFacet0->node( 1 ), 1 );
                tElement->insert_node( tFacet1->node( 1 ), 2 );
                tElement->insert_node( tFacet1->node( 0 ), 3 );

                // add element to Container
                tElements( e ) = tElement ;
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::create_tape_blocks_quad9(
                mesh::Block * aBlock, uint & aLayerOffset )
        {
            // Create an element factort
            ElementFactory tFactory ;

            // grab element container
            Cell< Element * > & tElements = aBlock->elements() ;

            // grab first layer
            Cell< Element * > & tLayer0 = mGhostLayers( aLayerOffset++ )->Elements ;

            // grab second layer
            Cell< Element * > & tLayer1 = mGhostLayers( aLayerOffset++ )->Elements ;

            // grab third layer
            Cell< Element * > & tLayer2 = mGhostLayers( aLayerOffset )->Elements ;

            // get number of elements
            index_t tNumElements = aBlock->number_of_elements() ;

            // loop over all elements
            for( index_t e=0; e<tNumElements; ++e )
            {
                // grab faces
                Element * tFacet0 = tLayer0( e );
                Element * tFacet1 = tLayer1( e );
                Element * tFacet2 = tLayer2( e );

                // create a new element
                Element * tElement = tFactory.create_lagrange_element( ElementType::QUAD9, ++mMaxElementID );

                // link nodes
                tElement->insert_node( tFacet0->node( 0 ), 0 );
                tElement->insert_node( tFacet0->node( 1 ), 1 );
                tElement->insert_node( tFacet2->node( 1 ), 2 );
                tElement->insert_node( tFacet2->node( 0 ), 3 );
                tElement->insert_node( tFacet0->node( 2 ), 4 );
                tElement->insert_node( tFacet1->node( 1 ), 5 );
                tElement->insert_node( tFacet2->node( 2 ), 6 );
                tElement->insert_node( tFacet1->node( 0 ), 7 );
                tElement->insert_node( tFacet1->node( 2 ), 8 );

                // add element to Container
                tElements( e ) = tElement ;
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::create_tape_blocks_penta6(
                mesh::Block * aBlock, uint & aLayerOffset )
        {
            // Create an element factort
            ElementFactory tFactory ;

            // grab element container
            Cell< Element * > & tElements = aBlock->elements() ;

            // grab first layer
            Cell< Element * > & tLayer0 = mGhostLayers( aLayerOffset++ )->Elements ;

            // grab second layer
            Cell< Element * > & tLayer1 = mGhostLayers( aLayerOffset )->Elements ;

            // get number of elements
            index_t tNumElements = aBlock->number_of_elements() ;

            // loop over all elements
            for( index_t e=0; e<tNumElements; ++e )
            {
                // grab faces
                Element * tFacet0 = tLayer0( e );
                Element * tFacet1 = tLayer1( e );

                // create a new element
                Element * tElement = tFactory.create_lagrange_element( ElementType::PENTA6, ++mMaxElementID );

                // link nodes
                tElement->insert_node( tFacet0->node( 0 ), 0 );
                tElement->insert_node( tFacet0->node( 1 ), 1 );
                tElement->insert_node( tFacet0->node( 2 ), 2 );

                tElement->insert_node( tFacet1->node( 0 ), 3 );
                tElement->insert_node( tFacet1->node( 1 ), 4 );
                tElement->insert_node( tFacet1->node( 2 ), 5 );

                // add element to Container
                tElements( e ) = tElement ;
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::create_tape_blocks_penta18(
                mesh::Block * aBlock, uint & aLayerOffset )
        {
            // Create an element factort
            ElementFactory tFactory ;

            // grab element container
            Cell< Element * > & tElements = aBlock->elements() ;

            // grab first layer
            Cell< Element * > & tLayer0 = mGhostLayers( aLayerOffset++ )->Elements ;

            // grab second layer
            Cell< Element * > & tLayer1 = mGhostLayers( aLayerOffset++ )->Elements ;

            // grab third layer
            Cell< Element * > & tLayer2 = mGhostLayers( aLayerOffset )->Elements ;

            // get number of elements
            index_t tNumElements = aBlock->number_of_elements() ;

            // loop over all elements
            for( index_t e=0; e<tNumElements; ++e )
            {
                // grab faces
                Element * tFacet0 = tLayer0( e );
                Element * tFacet1 = tLayer1( e );
                Element * tFacet2 = tLayer2( e );

                // create a new element
                Element * tElement = tFactory.create_lagrange_element( ElementType::PENTA18, ++mMaxElementID );

                // link nodes
                tElement->insert_node( tFacet0->node( 0 ),  0 );
                tElement->insert_node( tFacet0->node( 1 ),  1 );
                tElement->insert_node( tFacet0->node( 2 ),  2 );
                tElement->insert_node( tFacet2->node( 0 ),  3 );
                tElement->insert_node( tFacet2->node( 1 ),  4 );
                tElement->insert_node( tFacet2->node( 2 ),  5 );
                tElement->insert_node( tFacet0->node( 3 ),   6 );
                tElement->insert_node( tFacet0->node( 4 ),   7 );
                tElement->insert_node( tFacet0->node( 5 ),   8 );
                tElement->insert_node( tFacet1->node( 0 ),   9 );
                tElement->insert_node( tFacet1->node( 1 ),  10 );
                tElement->insert_node( tFacet1->node( 2 ),  11 );
                tElement->insert_node( tFacet2->node( 3 ),  12 );
                tElement->insert_node( tFacet2->node( 4 ),  13 );
                tElement->insert_node( tFacet2->node( 5 ),  14 );
                tElement->insert_node( tFacet1->node( 3 ),  15 );
                tElement->insert_node( tFacet1->node( 4 ),  16 );
                tElement->insert_node( tFacet1->node( 5 ),  17 );

                // add element to Container
                tElements( e ) = tElement ;
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::shift_nodes(const Vector< real > & aLayerThicknesses )
        {
            // field data
            Vector< real > & tNx = mMesh->field_data( "SurfaceNormalsx" );
            Vector< real > & tNy = mMesh->field_data( "SurfaceNormalsy" );
            Vector< real > & tNz = mMesh->field_data( "SurfaceNormalsz" );

            // container for normals
            Vector< real > tX( 3 );
            Vector< real > tN( 3 );

            uint tNumLayers = mGhostLayers.size() ;

            Vector< real > tShift( tNumLayers, BELFEM_QUIET_NAN );

            uint tNumThick = aLayerThicknesses.length();

            // create shifts
            switch( mElementOrder )
            {
                case( 1 ) :
                {
                    tShift( 0 ) = -0.5 * sum( aLayerThicknesses );
                    for ( uint l = 1; l < tNumLayers; ++l )
                    {
                        tShift( l ) = tShift( l - 1 ) + aLayerThicknesses( l - 1 );
                    }
                    break;
                }
                case( 2 ) :
                {
                    uint tCount = 0;
                    tShift( 0 ) = -0.5 * sum( aLayerThicknesses );
                    for ( uint l = 0; l < tNumThick; ++l )
                    {
                        real tX0 = tShift( tCount );
                        tShift( ++tCount ) = tX0 + 0.5 * aLayerThicknesses( l );
                        tShift( ++tCount ) = tX0 +       aLayerThicknesses( l );
                    }
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "unsupported element order: %u", ( unsigned int ) mElementOrder );
                }
            }


            // loop over all layers
            for( uint l=0 ; l<tNumLayers; ++l )
            {
                index_t tIndex;

                Cell< Node * > & tNodes = mGhostLayers( l )->Nodes;
                index_t tNumNodes = tNodes.size();

                // loop over all nodes
                for ( index_t k = 0; k < tNumNodes; ++k )
                {
                    // get index from original node
                    tIndex = mMesh->node( mNodeIDs( k ) )->index();

                    // grab node on layer
                    Node * tNode = tNodes( k );

                    // grab node coords
                    tX( 0 ) = tNode->x();
                    tX( 1 ) = tNode->y();
                    tX( 2 ) = tNode->z();

                    // populate normal vector
                    tN( 0 ) = tNx( tIndex );
                    tN( 1 ) = tNy( tIndex );
                    tN( 2 ) = tNz( tIndex );

                    // shift node
                    tX += tShift( l ) * tN;

                    if( tNode->id() == 39574 || tNode->id() == 39575 )
                    {
                        std::cout << "check" << std::endl ;
                        tX.print("X");
                        tN.print("N");
                    }

                    // write node coordinates back
                    tNode->set_coords( tX );

                    // grab index from this node
                    tIndex = tNode->index();

                    // write normal into field
                    tNx( tIndex ) = tN( 0 );
                    tNy( tIndex ) = tN( 1 );
                    tNz( tIndex ) = tN( 2 );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::flip_element_orientation()
        {
            // unflag all elements on the mesh
            mMesh->unflag_all_elements();

            BELFEM_ASSERT( mSelectedSideSets.size() == mMasterBlocks.size(),
                           "Number of sidesets and master blocks does not match (%u vs %u)",
                           ( unsigned int ) mSelectedSideSets.size(),
                           ( unsigned int ) mMasterBlocks.size() );

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
                        tFacet->set_master( tMaster, tMasterIndex, true );

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

            // allocate container
            mNodeIDs.set_size( tCount );
            tCount = 0 ;
            for ( Node * tNode: tNodes )
            {
                if ( tNode->is_flagged() )
                {
                   mNodeIDs( tCount++ ) = tNode->id() ;
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

            // get the table
            Matrix< id_t > & tTable = mMesh->tape_facet_table() ;

            // allocate memory
            tTable.set_size( 2, tCount * mNumberOfGhostLayers );

            // offset for table
            index_t tOffset = 0 ;

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

                        tLayerElements( tCount++ ) = tClone;

                        // remember ids
                        tTable( 0, tOffset )   = tClone->id() ;
                        tTable( 1, tOffset++ ) = tElement->id() ;

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
                mMesh->unflag_all_edges();
                for ( id_t tID: mSelectedSideSets )
                {
                    // get the container with the facets
                    Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets();

                    // loop over all facets
                    for ( Facet * tFacet: tFacets )
                    {
                        tFacet->element()->flag_edges();
                    }
                }

                // count number of edges
                index_t tCount = 0;



                // count cloned edges and create edge map
                Cell< Edge * > & tEdges = mMesh->edges();

                for ( Edge * tEdge: tEdges )
                {
                    if ( tEdge->is_flagged())
                    {
                        ++tCount ;
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

                    // reset the map
                    Map< index_t, Edge * > tEdgeMap ;

                    // clone edges
                    for( Edge * tEdge : tEdges )
                    {
                        if ( tEdge->is_flagged() )
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

                            tEdgeMap[ tEdge->id() ] = tClone ;

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
                            Edge * tEdge = tEdgeMap( tElement->edge( k )->id());

                            tClone->insert_edge( tEdge, k );
                        }
                    } // end loop over layer elements

                } // end loop over all layers
            } // end edges exist
        }

//------------------------------------------------------------------------------

        void
        TapeRoller::compute_edge_signs_2d()
        {
            // the normal vector
            Vector< real > & tNx = mMesh->field_data( "SurfaceNormalsx" );
            Vector< real > & tNy = mMesh->field_data( "SurfaceNormalsy" );


            // finally, we set the tags also for the ghost elements
            index_t tCount = 0 ;

            // loop over all sidesets
            for( id_t tID : mSelectedSideSets )
            {
                // grab the facets on the sideset
                Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;

                // loop over all facets on this sideset
                for( Facet * tFacet : tFacets )
                {
                   // approximate normal of facet
                   real tnx =     tNx( tFacet->node( 0 )->index() )
                           + tNx( tFacet->node( 1 )->index() );

                   real tny =   tNy( tFacet->node( 0 )->index() )
                         + tNy( tFacet->node( 1 )->index() );

                   real tnorm = std::sqrt( tnx * tnx + tny * tny );
                   tnx /= tnorm ;
                   tny /= tnorm ;

                   // approximate direction vector of facet
                   real tdx = tFacet->node( 1 )->x() - tFacet->node( 0 )->x() ;
                   real tdy = tFacet->node( 1 )->y() - tFacet->node( 0 )->y() ;
                   tnorm = std::sqrt( tdx*tdx + tdy*tdy );
                   tdx /= tnorm ;
                   tdy /= tnorm ;

                   // compute expression d - ( n x z ) ;
                   // note that ( n x z ) = [ ny ; -nx ; 0 ]

                   tdx -= tny ;
                   tdy += tnx ;
                   tnorm = std::sqrt( tdx*tdx + tdy*tdy );

                   // we store the information on the tags of the element
                   // todo: this is wrong and must be fixed

                   uint tTag = tnorm > 1.0 ? 1 : 0 ;

                   tFacet->element()->set_physical_tag( tTag );

                   // also set the tag for the ghost layers
                   for( uint l=0 ; l<mNumberOfGhostLayers; ++l )
                   {
                       mMesh->sideset( mGhostSideSetIDs( l ) )
                       ->facet_by_index( tCount )->element()->set_physical_tag( tTag );
                   }

                   // increment the counter
                   ++tCount ;
                }
            }
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