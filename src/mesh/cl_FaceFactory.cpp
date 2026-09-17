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

#include "cl_FaceFactory.hpp"

#include "assert.hpp"
#include "commtools.hpp"

#include "cl_Timer.hpp"
#include "cl_Logger.hpp"

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
                mCommRank( comm_rank() ),
                mMesh( aMesh ),
                mNumberOfNodes( aMesh.number_of_nodes() )
        {
            // Validate that mesh size is compatible with key128_t for face hashing
            // Face keys use: (i * N + j) * N + k (3 nodes minimum)
            // Maximum safe N is cube_root(key128_t_max)
            const key128_t tMaxNodes = static_cast< key128_t >( std::pow( static_cast< double >( std::numeric_limits< key128_t >::max() ), 1.0 / 3.0 ) );
            BELFEM_ERROR( mNumberOfNodes <= tMaxNodes,
                "Mesh has too many nodes (%llu) for face key computation with current key128_t type.\n"
                "Maximum allowed nodes: %llu\n"
                "Recommendation: This should never happen with 128-bit keys. Please report this as a bug.",
                ( long long unsigned int ) mNumberOfNodes,
                ( long long unsigned int ) tMaxNodes );
        }

//------------------------------------------------------------------------------

        FaceFactory::FaceFactory( Mesh * aMesh ):
                mCommRank( comm_rank() ),
                mMesh( *aMesh ) ,
                mNumberOfNodes( aMesh->number_of_nodes() )
        {
            // Validate that mesh size is compatible with key128_t for face hashing
            // Face keys use: (i * N + j) * N + k (3 nodes minimum)
            // Maximum safe N is cube_root(key128_t_max)
            const key128_t tMaxNodes = static_cast< key128_t >( std::pow( static_cast< double >( std::numeric_limits< key128_t >::max() ), 1.0 / 3.0 ) );
            BELFEM_ERROR( mNumberOfNodes <= tMaxNodes,
                "Mesh has too many nodes (%llu) for face key computation with current key128_t type.\n"
                "Maximum allowed nodes: %llu\n"
                "Recommendation: This should never happen with 128-bit keys. Please report this as a bug.",
                ( long long unsigned int ) mNumberOfNodes,
                ( long long unsigned int ) tMaxNodes );
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
                    this->get_all_block_ids( tNedelecBlocks );
                }
                else
                {
                    tNedelecBlocks = aNedelecBlocks ;
                }

                if( mMesh.number_of_dimensions() == 2 )
                {
                    this->create_faces_2d( tNedelecBlocks );

                    mMesh.finalize_faces();
                }
                else if( mMesh.number_of_dimensions() == 3 )
                {
                    message( InfoLevel::Default, "    Creating faces ...");

                    Timer tTimer ;

                    Map< key128_t, index_t > tFaceMap;
                    index_t tNumFaces = this->count_faces( tNedelecBlocks, aNedelecSideSets, tFaceMap );

                    Vector< id_t > tMaster;
                    Vector< id_t > tSlave;
                    Vector< index_t > tIndexOnMaster;
                    Vector< index_t > tIndexOnSlave;

                    this->find_face_owners(
                            tNedelecBlocks,
                            aNedelecSideSets,
                            tNumFaces,
                            tFaceMap,
                            tMaster,
                            tIndexOnMaster,
                            tSlave,
                            tIndexOnSlave );

                    this->allocate_face_containers( tNedelecBlocks );

                    mMesh.unflag_all_elements() ;
                    for ( id_t tBlockID : tNedelecBlocks )
                    {
                        mMesh.block( tBlockID )->flag_elements() ;
                    }

                    this->create_faces_3d(
                            tMaster,
                            tIndexOnMaster,
                            tSlave,
                            tIndexOnSlave );

                    this->set_face_ids( aNedelecSideSets, tFaceMap );

                    mMesh.finalize_faces();

                    message( InfoLevel::Detailed, "    ... number of faces                 : %lu",
                             ( long unsigned int ) mMesh.number_of_faces() );

                    message( InfoLevel::Detailed, "    ... time for creating faces         : %u ms\n",
                             ( unsigned int ) tTimer.stop() );
                }

            }
        }

//------------------------------------------------------------------------------

        index_t
        FaceFactory::count_faces(
                const Vector< id_t >  & aBlockIDs,
                const Vector< id_t >  & aSideSetIDs,
                Map< key128_t, index_t > & aFaceMap  )
        {

            index_t tCount = 0 ;
            for( id_t tID : aBlockIDs )
            {
                index_t tNumElems = mMesh.block( tID )->number_of_elements() ;

                switch( geometry_type( mMesh.block( tID )->element_type() ) )
                {
                    case( GeometryType::TET ) :
                    {
                        tCount += tNumElems * 4 ;
                        break ;
                    }
                    case( GeometryType::PENTA ) :
                    case( GeometryType::PYRA ) :
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

            Cell< key128_t > tFaceKeys( tCount, 0 );

            tCount = 0 ;

            Cell< Node * > tWork ;

            for( id_t tID : aBlockIDs )
            {
                Block * tBlock = mMesh.block( tID );

                uint n = number_of_faces( tBlock->element_type() );
                for( Element * tElement : tBlock->elements() )
                {
                    for ( uint f = 0; f < n; ++f )
                    {
                        tFaceKeys( tCount++ ) = this->face_key_3d(
                                tElement,
                                f,
                                tWork );
                    }
                }
            }

            for( id_t tID : aSideSetIDs )
            {
                Cell< Facet * > & tFacets = mMesh.sideset( tID )->facets() ;

                for( Facet * tFacet : tFacets )
                {
                   tFaceKeys( tCount++ ) = this->face_key_2d(
                           tFacet->element(),
                           tWork );
                }
            }

            BELFEM_ASSERT( tCount == tFaceKeys.size(), "Unknown Error" );

            unique( tFaceKeys );

            index_t aCount = 0 ;

            aFaceMap.clear() ;
            for( key128_t tKey : tFaceKeys )
            {
                aFaceMap[ tKey ] = aCount++ ;
            }

            return aCount ;
        }

//------------------------------------------------------------------------------

        void
        FaceFactory::get_all_block_ids( Vector< id_t > & aBlockIDs )
        {

            Cell< Block * > & tBlocks = mMesh.blocks();

            aBlockIDs.set_size( tBlocks.size() );

            uint tCount = 0 ;

            for ( Block * tBlock : tBlocks )
            {
                aBlockIDs( tCount++ ) = tBlock->id() ;
            }
        }

//------------------------------------------------------------------------------

        key128_t
        FaceFactory::face_key_2d(
                Element           * aElement,
                Cell< Node * >    & aWork )
        {
            switch( geometry_type( aElement->type() ) )
            {
                case( GeometryType::TRI ) :
                {
                    // the scratch is shared with face_key_3d, whose corner getter
                    // sizes it to the facet at hand ( 3 or 4 ) -- size it here too
                    aWork.set_size( 3, nullptr );

                    aWork( 0 ) = aElement->node( 0 ) ;
                    aWork( 1 ) = aElement->node( 1 ) ;
                    aWork( 2 ) = aElement->node( 2 ) ;

                    sort( aWork,opVertexIndex );

                    return (   aWork( 2 )->index() * mNumberOfNodes
                             + aWork( 1 )->index() ) * mNumberOfNodes
                             + aWork( 0 )->index() ;

                }
                case( GeometryType::QUAD ) :
                {
                    aWork.set_size( 4, nullptr );

                    aWork( 0 ) = aElement->node( 0 ) ;
                    aWork( 1 ) = aElement->node( 1 ) ;
                    aWork( 2 ) = aElement->node( 2 ) ;
                    aWork( 3 ) = aElement->node( 3 ) ;

                    sort( aWork,opVertexIndex  );

                    return (  aWork( 2 )->index() * mNumberOfNodes
                            + aWork( 1 )->index() ) * mNumberOfNodes
                            + aWork( 0 )->index() ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Element Type" );
                    return BELFEM_LUINT_MAX ;
                }
            }
        }

//------------------------------------------------------------------------------
        key128_t
        FaceFactory::face_key_3d(
                Element           * aElement,
                const uint          aFaceIndex,
                Cell< Node * >    & aWork )
        {
            aElement->get_corner_nodes_of_facet( aFaceIndex, aWork );

            sort( aWork,opVertexIndex  );

            return ( aWork( 2 )->index() * mNumberOfNodes
                      + aWork( 1 )->index() ) * mNumberOfNodes
                      + aWork( 0 )->index() ;
        }

//------------------------------------------------------------------------------

        void
        FaceFactory::find_face_owners(
                const Vector< id_t >        & aBlockIDs,
                const Vector< id_t >        & aSideSetIDs,
                const index_t               & aNumFaces,
                const Map< key128_t, index_t > & aFaceMap,
                Vector< id_t >              & aMasterIDs,
                Vector< index_t >           & aMasterIndex,
                Vector< id_t >              & aSlaveIDs,
                Vector< index_t >           & aSlaveIndex   )
        {
            aMasterIDs.set_size( aNumFaces, gNoID );
            aSlaveIDs.set_size( aNumFaces, gNoID );
            aMasterIndex.set_size( aNumFaces, gNoIndex );
            aSlaveIndex.set_size( aNumFaces, gNoIndex );

            Cell< Node * > tWork ;

            for( id_t tBlockID : aBlockIDs )
            {
                Block * tBlock = mMesh.block( tBlockID );

                uint n = number_of_faces( tBlock->element_type() );

                for ( Element * tElement : tBlock->elements() )
                {
                    for ( uint f = 0; f < n; ++f )
                    {
                        index_t tIndex = aFaceMap( this->face_key_3d(
                                tElement,
                                f,
                                tWork ) );

                        if ( aMasterIDs( tIndex ) != gNoID )
                        {
                            if( tElement->id() < aMasterIDs( tIndex ) )
                            {
                                aSlaveIDs( tIndex ) = aMasterIDs( tIndex );
                                aSlaveIndex( tIndex ) = aMasterIndex( tIndex );

                                aMasterIDs( tIndex ) = tElement->id();
                                aMasterIndex( tIndex ) = f;
                            }
                            else
                            {
                                aSlaveIDs( tIndex ) = tElement->id();
                                aSlaveIndex( tIndex ) = f;
                            }
                        }
                        else
                        {
                            aMasterIDs( tIndex ) = tElement->id();
                            aMasterIndex( tIndex ) = f;
                        }
                    }
                }
            }

            for( id_t tSideSetID : aSideSetIDs )
            {
                Cell< Facet * > & tFacets = mMesh.sideset( tSideSetID )->facets() ;

                for( Facet * tFacet : tFacets )
                {
                    index_t tIndex = aFaceMap(
                            this->face_key_2d(
                            tFacet->element(),
                            tWork ) );

                    if( aMasterIDs( tIndex ) == gNoID )
                    {
                        aMasterIDs( tIndex ) = tFacet->id();

                        // we deliberately don't write an index here
                    }
                }
            } // end loop over all sidesets
        }
//-----------------------------------------------------------------------

        void
        FaceFactory::allocate_face_containers( Vector< id_t > & aBlockIDs )
        {
            for( id_t tBlockID : aBlockIDs )
            {
                Block * tBlock = mMesh.block( tBlockID );

                for( Element * tElement : tBlock->elements() )
                {
                    tElement->allocate_face_container() ;
                }

                tBlock->set_faces_flag( true );
            }
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::create_faces_2d( const Vector< id_t > & aBlockIDs )
        {
            Cell< Face * > & tFaces = mMesh.faces() ;

            BELFEM_ASSERT( tFaces.size() == 0, "Faces of mesh have already been created");

            index_t tCount = 0 ;

            for( id_t b  : aBlockIDs )
            {
                tCount += mMesh.block( b )->number_of_elements() ;
            }

            tFaces.set_size( tCount, nullptr );
            tCount = 0 ;

            for( id_t b  : aBlockIDs )
            {
                Cell< Element * > & tElements = mMesh.block( b )->elements() ;

                for( Element * tElement : tElements )
                {
                    Face * tFace =  new Face( tElement ) ;

                    // in 2D, ids of elements and faces are identical
                    tFace->set_id( tElement->id() );

                    tFace->set_index( tCount );

                    tElement->allocate_face_container();

                    tElement->insert_face( tFace, 0 );

                    tFaces( tCount++ ) = tFace ;
                }
            }
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::create_faces_3d( Vector< id_t > & aMasterIDs,
                      Vector< index_t >             & aMasterIndex,
                      Vector< id_t >                & aSlaveIDs,
                      Vector< index_t >             & aSlaveIndex )
        {
            index_t tNumFaces = aMasterIDs.length();

            Cell< Face * > & tFaces = mMesh.faces() ;

            BELFEM_ASSERT( tFaces.size() == 0, "Faces of mesh have already been created");

            tFaces.set_size( tNumFaces, nullptr );

            Cell< Node * > tNodes ;

            for( index_t tIndex = 0; tIndex<tNumFaces; ++tIndex )
            {

                if( aSlaveIDs( tIndex ) < gNoIndex ) // face has primary and secondary elements
                {
                    Element * tMaster
                            = mMesh.element( aMasterIDs( tIndex ) );

                    Element * tSlave
                            = mMesh.element( aSlaveIDs( tIndex ) );

                    BELFEM_ASSERT( tMaster->is_flagged(), "Master element %lu is not flagged",
                        ( long unsigned int ) tMaster->id()  );

                    // we can't link the face to a slave if the slave is not part
                    // of the selected groups
                    if ( ! tSlave->is_flagged() )
                    {
                        tSlave = nullptr ;
                        aSlaveIndex( tIndex ) = gNoIndex ;
                    }

                    Face * tFace =  new Face(
                            tMaster,
                            aMasterIndex( tIndex ),
                            tSlave,
                            aSlaveIndex( tIndex ) ) ;

                    tMaster->insert_face( tFace, aMasterIndex( tIndex ) );

                    if ( tSlave != nullptr )
                    {
                        tSlave->insert_face( tFace, aSlaveIndex( tIndex ) );
                    }

                    tFaces( tIndex ) = tFace ;

                }
                else if ( aMasterIndex( tIndex ) < gNoIndex ) // face is only connected to one element
                {
                    Element * tMaster
                            = mMesh.element( aMasterIDs( tIndex ) );

                    Face * tFace =  new Face(
                            tMaster,
                            aMasterIndex( tIndex ),
                            nullptr,
                            gNoIndex) ;

                    tMaster->insert_face( tFace, aMasterIndex( tIndex ) );

                    tFaces( tIndex ) = tFace ;
                }
                else
                {
                    Element * tElement = mMesh.facet( aMasterIDs( tIndex ) )->element() ;

                    Face * tFace = new Face( tElement, 0, nullptr, gNoIndex );

                    tElement->insert_face( tFace, aMasterIndex( tIndex ) );

                    tFaces( tIndex ) = tFace ;
                }
        }
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::set_face_ids(  const Vector< id_t > & aSideSets, const Map< key128_t, index_t > & aFaceMap )
        {

            // scratch for face_key_2d, which writes up to four corner nodes
            Cell< Node * > tWork( 4, nullptr );

            Cell< Face * > & tFaces = mMesh.faces() ;

            for ( id_t s : aSideSets )
            {
                SideSet * tSideSet = mMesh.sideset( s );

                for( Facet * tFacet : tSideSet->facets() )
                {
                    // The key MUST come from face_key_2d, the same function that
                    // built aFaceMap. This used to be a hand-copied duplicate of
                    // that logic and had drifted from it: the QUAD branch sorted
                    // tTriNodes instead of tQuadNodes, so the quad key was
                    // computed from UNSORTED nodes and could not match the sorted
                    // key in the map -- a wrong index or a bare "Key not found in
                    // map" for any quad facet. The local key was also a luint,
                    // truncating the 128-bit key that the map type exists to
                    // provide ( the key is O( numNodes^3 ), so it leaves 64 bits
                    // at ~2.1e6 nodes ).
                    const key128_t tKey =
                            this->face_key_2d( tFacet->element(), tWork );

                    index_t tIndex = aFaceMap( tKey ) ;

                    tFaces( tIndex )->set_id( tFacet->id() );

                    tFaces( tIndex )->flag() ;
                }
            }

            id_t tMaxID = 0 ;
            for( Element * tEdge : mMesh.boundary_edges() )
            {
                tMaxID = tMaxID < tEdge->id() ? tEdge->id() : tMaxID ;
            }
            for( Edge * tEdge : mMesh.edges() )
            {
                tMaxID = tMaxID < tEdge->id() ? tEdge->id() : tMaxID ;
            }
            for( Facet * tFacet : mMesh.facets() )
            {
                tMaxID = tMaxID < tFacet->id() ? tFacet->id() : tMaxID ;
            }
            for( Element * tElement : mMesh.elements() )
            {
                tMaxID = tMaxID < tElement->id() ? tElement->id() : tMaxID ;
            }

            for( Face * tFace : mMesh.faces() )
            {
                if( ! tFace->is_flagged() )
                {
                    tFace->set_id( ++tMaxID );
                }
            }

            sort( tFaces, opVertexID );
        }

//-----------------------------------------------------------------------

        void
        FaceFactory::print()
        {
            Cell< Node * > tNodes ;

            for( Face * tFace : mMesh.faces() )
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

                for( Node * tNode : tNodes )
                {
                    std::cout << " " << tNode->id() ;
                }

                std::cout << std::endl ;
            }
        }

//-----------------------------------------------------------------------
    }
}
