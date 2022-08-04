//
// Created by christian on 7/12/21.
//
#include "cl_FaceFactory.hpp"

#include "assert.hpp"
#include "commtools.hpp"

#include "op_Graph_Vertex_Index.hpp"

#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Mesh.hpp"
#include "meshtools.hpp"
#include "fn_unique.hpp"
#include "cl_Face.hpp"
#include "op_Graph_Vertex_ID.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        FaceFactory::FaceFactory( Mesh & aMesh ):
                mRank( comm_rank() ),
                mMesh( aMesh ),
                mNumberOfNodes( aMesh.number_of_nodes() )
        {

        }

//------------------------------------------------------------------------------

        FaceFactory::FaceFactory( Mesh * aMesh ):
                mRank( comm_rank() ),
                mMesh( *aMesh ) ,
                mNumberOfNodes( aMesh->number_of_nodes() )
        {

        }

//------------------------------------------------------------------------------

        void
        FaceFactory::create_faces( const Vector< id_t > aNedelecBlocks, const Vector< id_t > aNedelecSideSets  )
        {
            if ( comm_rank() == 0 )
            {
                Vector< id_t > tNedelecBlocks;
                if( aNedelecBlocks.length() == 0 && aNedelecSideSets.length() == 0 )
                {
                    // select all blocks
                    this->get_all_block_ids( tNedelecBlocks );
                }
                else
                {
                    // select given blocks
                    tNedelecBlocks = aNedelecBlocks ;
                }

                if( mMesh.number_of_dimensions() == 2 )
                {
                    this->create_faces_2d( tNedelecBlocks );
                }
                else if( mMesh.number_of_dimensions() == 3 )
                {
                    BELFEM_ASSERT( mNumberOfNodes < std::pow( static_cast< real >( BELFEM_LUINT_MAX ), 1.0 / 3.0 ),
                                  "Too many nodes" );


                    // count faces and crate a temporary map
                    Map< luint, index_t > tFaceMap;
                    index_t tNumFaces = this->count_faces( tNedelecBlocks, aNedelecSideSets, tFaceMap );

                    Vector< id_t > tPrimaryOwner;
                    Vector< id_t > tSecondaryOwner;
                    Vector< index_t > tPrimaryIndex;
                    Vector< index_t > tSecondaryIndex;

                    this->find_face_owners(
                            tNedelecBlocks,
                            aNedelecSideSets,
                            tNumFaces,
                            tFaceMap,
                            tPrimaryOwner,
                            tPrimaryIndex,
                            tSecondaryOwner,
                            tSecondaryIndex );


                    // allocate containers of elements
                    this->allocate_face_containers( tNedelecBlocks );

                    this->create_faces_3d(
                            tPrimaryOwner,
                            tPrimaryIndex,
                            tSecondaryOwner,
                            tSecondaryIndex );

                    // compute the face IDs
                    this->set_face_ids( aNedelecSideSets, tFaceMap );
                }

                mMesh.finalize_faces();
            }
        }

//------------------------------------------------------------------------------

        index_t
        FaceFactory::count_faces(
                const Vector< id_t >  & aBlockIDs,
                const Vector< id_t >  & aSideSetIDs,
                Map< luint, index_t > & aFaceMap  )
        {

            index_t tCount = 0 ;
            for( id_t tID : aBlockIDs )
            {
                index_t tNumElems = mMesh.block( tID )->number_of_elements() ;

                // num faces
                switch( geometry_type( mMesh.block( tID )->element_type() ) )
                {
                    case( GeometryType::TET ) :
                    {
                        tCount += tNumElems * 4 ;
                        break ;
                    }
                    case( GeometryType::PENTA ) :
                    {
                        tCount += tNumElems * 5 ;
                        break ;
                    }
                    case( GeometryType::HEX ) :
                    {
                        tCount += tNumElems * 6 ;
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid Element Type");
                    }
                }
            }
            for( id_t tID : aSideSetIDs )
            {
                tCount += mMesh.sideset( tID )->number_of_facets() ;
            }

            // allocate vector
            Vector< luint > tFaceIDs( tCount );

            // reset counter
            tCount = 0 ;

            Cell< mesh::Node * > tAllNodes ;
            Cell< mesh::Node * > tTriNodes( 3, nullptr );
            Cell< mesh::Node * > tQuadNodes( 4, nullptr );

            // populate vector
            for( id_t tID : aBlockIDs )
            {
                // get block
                mesh::Block * tBlock = mMesh.block( tID );

                switch( geometry_type( tBlock->element_type() ) )
                {
                    case( GeometryType::TET ) :
                    {
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            for ( uint f = 0; f < 4; ++f )
                            {
                                tFaceIDs( tCount++ ) = this->face_key_tri(
                                        tElement,
                                        f,
                                        tAllNodes,
                                        tTriNodes );
                            }
                        }
                        break ;
                    }
                    case( GeometryType::PENTA ) :
                    {
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            for ( uint f = 0; f < 3; ++f )
                            {
                                tFaceIDs( tCount++ ) = this->face_key_tri(
                                        tElement,
                                        f,
                                        tAllNodes,
                                        tTriNodes );
                            }
                            for ( uint f = 3; f < 5; ++f )
                            {
                                tFaceIDs( tCount++ ) = this->face_key_quad(
                                        tElement,
                                        f,
                                        tAllNodes,
                                        tQuadNodes );
                            }
                        }
                        break ;
                    }
                    case( GeometryType::HEX ) :
                    {
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            for ( uint f = 0; f < 6; ++f )
                            {
                                tFaceIDs( tCount++ ) = this->face_key_quad(
                                        tElement,
                                        f,
                                        tAllNodes,
                                        tQuadNodes );
                            }
                        }
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid Element Type");
                    }
                }
           }

            for( id_t tID : aSideSetIDs )
            {
                // get facets on sideset
                Cell< Facet * > & tFacets = mMesh.sideset( tID )->facets() ;

                // loop over all facets
                for( Facet * tFacet : tFacets )
                {
                   // compute face ID
                   tFaceIDs( tCount++ ) = this->face_key_2d(
                           tFacet->element(),
                           tTriNodes,
                           tQuadNodes );
                }
            }

            // make sure that we have computed all faces
            BELFEM_ASSERT( tCount == tFaceIDs.length(), "Unknown Error" );

            // make ids unique
            unique( tFaceIDs );

            // reset the counter
            index_t aCount = 0 ;

            // create the map
            aFaceMap.clear() ;
            for( luint tID : tFaceIDs )
            {
                // write index into map
                aFaceMap[ tID ] = aCount++ ;
            }

            return aCount ;
        }

//------------------------------------------------------------------------------

        void
        FaceFactory::get_all_block_ids( Vector< id_t > & aBlockIDs )
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

        luint
        FaceFactory::face_key_2d(
                Element           * aElement,
                Cell< Node * >    & aTriNodes,
                Cell< Node * >    & aQuadNodes )
        {
            switch( geometry_type( aElement->type() ) )
            {
                case( GeometryType::TRI ) :
                {
                    // collect nodes
                    aTriNodes( 0 ) = aElement->node( 0 ) ;
                    aTriNodes( 1 ) = aElement->node( 1 ) ;
                    aTriNodes( 2 ) = aElement->node( 2 ) ;

                    // sort the nodes
                    sort( aTriNodes,opVertexIndex );

                    // return the ID
                    return (   aTriNodes( 2 )->index() * mNumberOfNodes
                             + aTriNodes( 1 )->index() ) * mNumberOfNodes
                             + aTriNodes( 0 )->index() ;

                }
                case( GeometryType::QUAD ) :
                {
                    aQuadNodes( 0 ) = aElement->node( 0 ) ;
                    aQuadNodes( 1 ) = aElement->node( 1 ) ;
                    aQuadNodes( 2 ) = aElement->node( 2 ) ;
                    aQuadNodes( 3 ) = aElement->node( 3 ) ;

                    sort( aQuadNodes,opVertexIndex  );

                    return (  aQuadNodes( 2 )->index() * mNumberOfNodes
                            + aQuadNodes( 1 )->index() ) * mNumberOfNodes
                            + aQuadNodes( 0 )->index() ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Element Type" );
                    return BELFEM_LUINT_MAX ;
                }
            }
        }

//------------------------------------------------------------------------------
        luint
        FaceFactory::face_key_tri(
                Element           * aElement,
                const uint          aFaceIndex,
                Cell< Node * >    & aAllNodes,
                Cell< Node * >    & aCornerNodes )
        {
            // grab face
            aElement->get_nodes_of_facet( aFaceIndex, aAllNodes );

            aCornerNodes( 0 ) = aAllNodes( 0 );
            aCornerNodes( 1 ) = aAllNodes( 1 );
            aCornerNodes( 2 ) = aAllNodes( 2 );

            // sort the nodes
            sort( aCornerNodes,opVertexIndex  );

            // return the ID
            return ( aCornerNodes( 2 )->index() * mNumberOfNodes
                      + aCornerNodes( 1 )->index() ) * mNumberOfNodes
                      + aCornerNodes( 0 )->index() ;
        }

//------------------------------------------------------------------------------

        luint
        FaceFactory::face_key_quad(
                Element           * aElement,
                const uint          aFaceIndex,
                Cell< Node * >    & aAllNodes,
                Cell< Node * >    & aCornerNodes )
        {
            // grab face
            aElement->get_nodes_of_facet( aFaceIndex, aAllNodes );

            aCornerNodes( 0 ) = aAllNodes( 0 );
            aCornerNodes( 1 ) = aAllNodes( 1 );
            aCornerNodes( 2 ) = aAllNodes( 2 );
            aCornerNodes( 3 ) = aAllNodes( 3 );

            // sort the nodes
            sort( aCornerNodes,opVertexIndex  );

            // return the ID
            return (   aCornerNodes( 3 )->index()   * mNumberOfNodes
                     + aCornerNodes( 2 )->index() ) * mNumberOfNodes
                     + aCornerNodes( 1 )->index() ;
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::find_face_owners(
                const Vector< id_t >        & aBlockIDs,
                const Vector< id_t >        & aSideSetIDs,
                const index_t               & aNumFaces,
                const Map< luint, index_t > & aFaceMap,
                Vector< id_t >              & aMasterOwner,
                Vector< index_t >           & aMasterIndex,
                Vector< id_t >              & aSlaveOwner,
                Vector< index_t >           & aSlaveIndex   )
        {
            // allocate memory
            aMasterOwner.set_size( aNumFaces, gNoID );
            aSlaveOwner.set_size( aNumFaces, gNoID );
            aMasterIndex.set_size( aNumFaces, gNoIndex );
            aSlaveIndex.set_size( aNumFaces, gNoIndex );

            Cell< mesh::Node * > tAllNodes ;
            Cell< mesh::Node * > tTriNodes( 3, nullptr );
            Cell< mesh::Node * > tQuadNodes( 4, nullptr );

            // loop over all blocks
            for( id_t tBlockID : aBlockIDs )
            {
                // get block
                mesh::Block * tBlock = mMesh.block( tBlockID );

                switch( geometry_type( tBlock->element_type() ) )
                {
                    case ( GeometryType::TET ) :
                    {
                        for ( mesh::Element * tElement : tBlock->elements() )
                        {
                            for ( uint f = 0; f < 4; ++f )
                            {
                                // get the index
                                index_t tIndex = aFaceMap( this->face_key_tri(
                                        tElement,
                                        f,
                                        tAllNodes,
                                        tTriNodes ));

                                // check if owner is present
                                if ( aMasterOwner( tIndex ) != gNoID )
                                {
                                    // shift owner
                                    aSlaveOwner( tIndex ) = aMasterOwner( tIndex );
                                    aSlaveIndex( tIndex ) = aMasterIndex( tIndex );
                                }

                                // remember stuff
                                aMasterOwner( tIndex ) = tElement->id();
                                aMasterIndex( tIndex ) = f;
                            }
                        }
                        break;
                    }
                    case ( GeometryType::PENTA ) :
                    {
                        for ( mesh::Element * tElement : tBlock->elements())
                        {
                            for ( uint f = 0; f < 3; ++f )
                            {
                                // get the index
                                index_t tIndex = aFaceMap( this->face_key_tri(
                                        tElement,
                                        f,
                                        tAllNodes,
                                        tTriNodes ));

                                // check if owner is present
                                if ( aMasterOwner( tIndex ) != gNoID )
                                {
                                    // shift owner
                                    aSlaveOwner( tIndex ) = aMasterOwner( tIndex );
                                    aSlaveIndex( tIndex ) = aMasterIndex( tIndex );
                                }

                                // remember stuff
                                aMasterOwner( tIndex ) = tElement->id();
                                aMasterIndex( tIndex ) = f;
                            }
                            for ( uint f = 3; f < 5; ++f )
                            {
                                // get the index
                                index_t tIndex = aFaceMap( this->face_key_quad(
                                        tElement,
                                        f,
                                        tAllNodes,
                                        tQuadNodes ));

                                // check if owner is present
                                if ( aMasterOwner( tIndex ) != gNoID )
                                {
                                    // shift owner
                                    aSlaveOwner( tIndex ) = aMasterOwner( tIndex );
                                    aSlaveIndex( tIndex ) = aMasterIndex( tIndex );
                                }

                                // remember stuff
                                aMasterOwner( tIndex ) = tElement->id();
                                aMasterIndex( tIndex ) = f;
                            }
                        }
                        break;
                    }
                    case ( GeometryType::HEX ) :
                    {
                        for ( mesh::Element * tElement : tBlock->elements())
                        {
                            for ( uint f = 0; f < 6; ++f )
                            {
                                // get the index
                                index_t tIndex = aFaceMap( this->face_key_quad(
                                        tElement,
                                        f,
                                        tAllNodes,
                                        tQuadNodes ));

                                // check if owner is present
                                if ( aMasterOwner( tIndex ) != gNoID )
                                {
                                    // shift owner
                                    aSlaveOwner( tIndex ) = aMasterOwner( tIndex );
                                    aSlaveIndex( tIndex ) = aMasterIndex( tIndex );
                                }

                                // remember stuff
                                aMasterOwner( tIndex ) = tElement->id();
                                aMasterIndex( tIndex ) = f;
                            }
                        }
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid Element Type" );
                    }
                }
            } // end loop over all blocks

            // loop over all sidesets
            for( id_t tSideSetID : aSideSetIDs )
            {
                // get facet
                Cell< Facet * > & tFacets = mMesh.sideset( tSideSetID )->facets() ;

                for( Facet * tFacet : tFacets )
                {
                    // get the index
                    index_t tIndex = aFaceMap(
                            this->face_key_2d(
                            tFacet->element(),
                            tTriNodes,
                            tQuadNodes ) );

                    // check if facet has been claimed
                    if( aMasterOwner( tIndex ) == gNoID )
                    {
                        aMasterOwner( tIndex ) = tFacet->id();

                        // we deliberateley don't write an index here
                    }
                }
            } // end loop over all sidesets

        }
//-----------------------------------------------------------------------

        void
        FaceFactory::allocate_face_containers( Vector< id_t > & aBlockIDs )
        {
            // loop over all blocks
            for( id_t tBlockID : aBlockIDs )
            {
                // get block
                mesh::Block * tBlock = mMesh.block( tBlockID );

                for( mesh::Element * tElement : tBlock->elements() )
                {
                    tElement->allocate_face_container() ;
                }
            }
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::create_faces_2d( const Vector< id_t > & aBlockIDs )
        {
            Cell< mesh::Face * > & tFaces = mMesh.faces() ;

            BELFEM_ASSERT( tFaces.size() == 0, "Faces of mesh have already been created");

            // count elements
            index_t tCount = 0 ;

            for( id_t b  : aBlockIDs )
            {
                tCount += mMesh.block( b )->number_of_elements() ;
            }

            // allocate container
            tFaces.set_size( tCount, nullptr );
            tCount = 0 ;

            for( id_t b  : aBlockIDs )
            {
                // get elements
                Cell< Element * > & tElements = mMesh.block( b )->elements() ;

                for( Element * tElement : tElements )
                {
                    // create new face
                    Face * tFace =  new Face( tElement ) ;

                    // in 2D, ids of elements and faces are identical
                    tFace->set_id( tElement->id() );

                    // set the index of the face
                    tFace->set_index( tCount );

                    // allocate memory
                    tElement->allocate_face_container();

                    // link face to element container
                    tElement->insert_face( tFace, 0 );

                    // add face to container
                    tFaces( tCount++ ) = tFace ;
                }
            }
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::create_faces_3d( Vector< id_t > & aMasterOwner,
                      Vector< index_t >             & aMasterIndex,
                      Vector< id_t >                & aSlaveOwner,
                      Vector< index_t >             & aSlaveIndex )
        {
            index_t tNumFaces = aMasterOwner.length();

            Cell< mesh::Face * > & tFaces = mMesh.faces() ;

            BELFEM_ASSERT( tFaces.size() == 0, "Faces of mesh have already been created");

            tFaces.set_size( tNumFaces, nullptr );

            // loop over all faces that are to be created
            for( index_t tIndex = 0; tIndex<tNumFaces; ++tIndex )
            {

                if( aSlaveOwner( tIndex ) < gNoID ) // face has primary and secondary elements
                {
                    // get primary element
                    Element * tMaster
                            = mMesh.element( aMasterOwner( tIndex ) );

                    // get secondary element
                    Element * tSlave
                            = mMesh.element( aSlaveOwner( tIndex ) );

                    // create new face
                    Face * tFace =  new Face(
                            tMaster,
                            aMasterIndex( tIndex ),
                            tSlave,
                            aSlaveIndex( tIndex ) ) ;

                    // link face to primary element
                    tMaster->insert_face( tFace, aMasterIndex( tIndex ) );

                    // link face to secondary element
                    tSlave->insert_face( tFace, aSlaveIndex( tIndex ) );

                    // add face to container
                    tFaces( tIndex ) = tFace ;
                }
                else if ( aMasterIndex( tIndex ) < gNoIndex ) // face is only connected to one element
                {
                    // get primary element
                    Element * tMaster
                            = mMesh.element( aMasterOwner( tIndex ) );

                    // create new face
                    Face * tFace =  new Face(
                            tMaster,
                            aMasterIndex( tIndex ),
                            nullptr,
                            gNoIndex) ;

                    // link face to primary element
                    tMaster->insert_face( tFace, aMasterIndex( tIndex ) );

                    // add face to container
                    tFaces( tIndex ) = tFace ;
                }
                else // this is a facet-only face
                {
                    // get Facet Element
                    Element * tElement = mMesh.facet( aMasterOwner( tIndex ) )->element() ;

                    // create a new face
                    Face * tFace = new Face( tElement, 0, nullptr, gNoIndex );

                    // link face to primary element
                    tElement->insert_face( tFace, aMasterIndex( tIndex ) );

                    // add face to container
                    tFaces( tIndex ) = tFace ;
                }
            }
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::set_face_ids(  const Vector< id_t > & aSideSets, const Map< luint, index_t > & aFaceMap )
        {

            Cell< mesh::Node * > tTriNodes( 3, nullptr );
            Cell< mesh::Node * > tQuadNodes( 4, nullptr );

            luint tKey ;

            Cell< mesh::Face * > & tFaces = mMesh.faces() ;

            // loop over all sidesets
            for ( id_t s : aSideSets )
            {
                // grab sideset
                mesh::SideSet * tSideSet = mMesh.sideset( s );

                // loop over all facets on this sideset
                for( mesh::Facet * tFacet : tSideSet->facets() )
                {
                    switch( geometry_type( tFacet->element()->type() ) )
                    {
                        case( GeometryType::TRI ) :
                        {
                            tTriNodes( 0 ) = tFacet->element()->node( 0 );
                            tTriNodes( 1 ) = tFacet->element()->node( 1 );
                            tTriNodes( 2 ) = tFacet->element()->node( 2 );
                            sort( tTriNodes, opVertexIndex  );

                            // compute key
                            tKey =  (   tTriNodes( 2 )->index() * mNumberOfNodes
                                      + tTriNodes( 1 )->index() ) * mNumberOfNodes
                                      + tTriNodes( 0 )->index() ;
                            break ;
                        }
                        case( GeometryType::QUAD ) :
                        {
                            tQuadNodes( 0 ) = tFacet->element()->node( 0 );
                            tQuadNodes( 1 ) = tFacet->element()->node( 1 );
                            tQuadNodes( 2 ) = tFacet->element()->node( 2 );
                            tQuadNodes( 3 ) = tFacet->element()->node( 3 );
                            sort( tTriNodes, opVertexIndex  );

                            // compute key
                            tKey =  (   tTriNodes( 3 )->index() * mNumberOfNodes
                                      + tTriNodes( 2 )->index() ) * mNumberOfNodes
                                      + tTriNodes( 1 )->index() ;
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Invalid geometry type");
                        }
                    }

                    // get index for face
                    index_t tIndex = aFaceMap( tKey ) ;

                    // set id of face to be identical to facet
                    tFaces( tIndex )->set_id( tFacet->id() );

                    // flag this face
                    tFaces( tIndex )->flag() ;
                }
            }

            // compute the maximum ID so far
            id_t tMaxID = 0 ;
            for( mesh::Element * tEdge : mMesh.boundary_edges() )
            {
                tMaxID = tMaxID < tEdge->id() ? tEdge->id() : tMaxID ;
            }
            for( mesh::Edge * tEdge : mMesh.edges() )
            {
                tMaxID = tMaxID < tEdge->id() ? tEdge->id() : tMaxID ;
            }
            for( mesh::Facet * tFacet : mMesh.facets() )
            {
                tMaxID = tMaxID < tFacet->id() ? tFacet->id() : tMaxID ;
            }
            for( mesh::Element * tElement : mMesh.elements() )
            {
                tMaxID = tMaxID < tElement->id() ? tElement->id() : tMaxID ;
            }

            for( mesh::Face * tFace : mMesh.faces() )
            {
                if( ! tFace->is_flagged() )
                {
                    tFace->set_id( ++tMaxID );
                }
            }

            // resort the array
            sort( tFaces, opVertexID );
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::print()
        {
            Cell< mesh::Node * > tNodes ;

            for( mesh::Face * tFace : mMesh.faces() )
            {
                if( tFace->slave() != nullptr )
                {
                    std::cout << "id : " << tFace->id() << " : " << tFace->master()->id() << " "
                     << tFace->slave()->id() << " :" ;
                }
                else
                {
                    std::cout << "id : " << tFace->id() << " : " << tFace->master()->id() << " "
                              << " :" ;
                }

                tFace->master()->get_nodes_of_facet( tFace->index_on_master(), tNodes );

                for( mesh::Node * tNode : tNodes )
                {
                    std::cout << " " << tNode->id() ;
                }

                std::cout << std::endl ;
            }
        }

//-----------------------------------------------------------------------
    }
}