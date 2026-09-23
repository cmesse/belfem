/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for mesh::FaceFactory::create_faces() on a hand-built pair of
 * tetrahedra that share one face: the face key is built from the sorted
 * corner triple alone, so at order 2 no midside node can merge two facets of
 * one element; exactly the coincident facets are paired, the lower element id
 * is the master whichever element the factory visits first, and after the
 * connectivity pass every face carries its three edges and each neighbour
 * once. A HEX8 with one sideset facet pins that the face at a sideset facet
 * takes the facet's id, and a sideset facet on an unselected block is
 * refused by an error naming the facet ( the mesh main forces the throwing
 * reaction, so the refusal is a std::runtime_error in both builds ).
 *
 * Premises this file relies on: node indices are push positions ( the key is
 * O( N^3 ) in a 128-bit integer, far from overflow for 14 nodes ), and the
 * mesh is built with connectivity computation on, which is what makes
 * finalize_faces() run the face-edge and face-face passes once edges exist.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Face.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "meshtools.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    const real tTiny = 1e-12 ;

    // Two tetrahedra on one block. A = ( P0, P1, P2, P3 ), id 1, and
    // B = ( P1, P2, P3, P4 ), id 2, both positively oriented; they share
    // the triangle P1 P2 P3, which is facet 1 of A and facet 3 of B.
    //
    // Node push order is the index order the face key sees. Corners and
    // midsides are interleaved so that the two lowest corner indices of A's
    // facets 0 and 3 ( P0, P1 ) and the midside between them are the three
    // lowest indices of both facets: a key taken from the full facet node
    // set would merge the two facets and pair A with itself, a key taken
    // from the corners alone keeps them apart.
    //
    // At order 2 every midside is set from the element's own edge table, so
    // node k of an element is the midpoint of the edge the table names.
    //
    // aReverse pushes B before A ( ids kept ), so the factory meets the
    // higher id first and must displace it when A arrives.
    struct TetPair
    {
        Mesh    * mMesh ;
        Element * mA ;
        Element * mB ;
        uint      mNodesPerFace ;

        Node * mP0 ;
        Node * mP1 ;
        Node * mP2 ;
        Node * mP3 ;
        Node * mP4 ;

        TetPair( const uint aOrder, const bool aReverse = false ) :
            mNodesPerFace( aOrder == 2 ? 6 : 3 )
        {
            mMesh = new Mesh( 3, 0, true );

            Cell< Node * > & tNodes = mMesh->nodes();

            auto push = [ &tNodes ]( const real aX, const real aY, const real aZ ) -> Node *
            {
                Node * tNode = new Node( ( id_t ) tNodes.size() + 1, aX, aY, aZ );
                tNodes.push( tNode );
                return tNode ;
            };

            auto mid = [ &push ]( const Node * aP, const Node * aQ ) -> Node *
            {
                return push( 0.5 * ( aP->x() + aQ->x() ),
                             0.5 * ( aP->y() + aQ->y() ),
                             0.5 * ( aP->z() + aQ->z() ) );
            };

            mP0 = push( 0.0, 0.0, 0.0 );
            mP1 = push( 1.0, 0.0, 0.0 );

            Node * tM01 = nullptr ;
            if ( aOrder == 2 )
            {
                tM01 = mid( mP0, mP1 );
            }

            mP2 = push( 0.0, 1.0, 0.0 );
            mP3 = push( 0.0, 0.0, 1.0 );
            mP4 = push( 1.0, 1.0, 1.0 );

            ElementFactory tFactory ;

            const ElementType tType = ( aOrder == 2 ) ?
                    ElementType::TET10 : ElementType::TET4 ;

            mA = tFactory.create_element( tType, 1 );
            mB = tFactory.create_element( tType, 2 );

            mA->insert_node( mP0, 0 );
            mA->insert_node( mP1, 1 );
            mA->insert_node( mP2, 2 );
            mA->insert_node( mP3, 3 );

            mB->insert_node( mP1, 0 );
            mB->insert_node( mP2, 1 );
            mB->insert_node( mP3, 2 );
            mB->insert_node( mP4, 3 );

            if ( aOrder == 2 )
            {
                Node * tM12 = mid( mP1, mP2 );
                Node * tM02 = mid( mP0, mP2 );
                Node * tM03 = mid( mP0, mP3 );
                Node * tM13 = mid( mP1, mP3 );
                Node * tM23 = mid( mP2, mP3 );
                Node * tM14 = mid( mP1, mP4 );
                Node * tM24 = mid( mP2, mP4 );
                Node * tM34 = mid( mP3, mP4 );

                // TET10 edge table: 4 = ( 0, 1 ), 5 = ( 1, 2 ), 6 = ( 2, 0 ),
                // 7 = ( 0, 3 ), 8 = ( 1, 3 ), 9 = ( 2, 3 ), in local corners
                mA->insert_node( tM01, 4 );
                mA->insert_node( tM12, 5 );
                mA->insert_node( tM02, 6 );
                mA->insert_node( tM03, 7 );
                mA->insert_node( tM13, 8 );
                mA->insert_node( tM23, 9 );

                mB->insert_node( tM12, 4 );
                mB->insert_node( tM23, 5 );
                mB->insert_node( tM13, 6 );
                mB->insert_node( tM14, 7 );
                mB->insert_node( tM24, 8 );
                mB->insert_node( tM34, 9 );
            }

            Block * tBlock = new Block( 1, 2 );
            Cell< Element * > & tElements = mMesh->elements();

            Element * tFirst  = aReverse ? mB : mA ;
            Element * tSecond = aReverse ? mA : mB ;

            tFirst->set_block_id( 1 );
            tElements.push( tFirst );
            tBlock->insert_element( tFirst );

            tSecond->set_block_id( 1 );
            tElements.push( tSecond );
            tBlock->insert_element( tSecond );

            mMesh->blocks().push( tBlock );
        }

        ~TetPair()
        {
            delete mMesh ;   // owns nodes, blocks ( and their elements ), edges, faces
        }

        // owning fixture: a copy or a move would double-free
        TetPair( const TetPair & ) = delete ;
        TetPair & operator=( const TetPair & ) = delete ;
        TetPair( TetPair && ) = delete ;
        TetPair & operator=( TetPair && ) = delete ;

        // edges must exist before the faces, or the face connectivity
        // passes are skipped
        void
        build_faces()
        {
            mMesh->finalize();
            mMesh->create_edges();
            mMesh->create_faces();
        }

        // an edge of the shared triangle has both corners among P1, P2, P3
        bool
        is_shared_face_corner( const Node * aNode ) const
        {
            return aNode == mP1 || aNode == mP2 || aNode == mP3 ;
        }
    };

//------------------------------------------------------------------------------

    // the seven faces are the seven geometric faces: one two-owner face at
    // ( A, 1 ) / ( B, 3 ), six one-owner faces, every ( element, f ) slot
    // linked back to exactly the face that names it
    void
    expect_two_tets_share_one_face( TetPair & aPair )
    {
        Mesh    * tMesh = aPair.mMesh ;
        Element * tA    = aPair.mA ;
        Element * tB    = aPair.mB ;

        ASSERT_EQ( tMesh->number_of_faces(), ( index_t ) 7 );

        Cell< Face * > & tFaces = tMesh->faces();

        for ( uint k = 0; k < 7; ++k )
        {
            EXPECT_EQ( tFaces( k )->index(), ( index_t ) k );
        }

        // ids pairwise distinct and disjoint from the element and edge ids
        for ( uint k = 0; k < 7; ++k )
        {
            for ( uint l = k + 1; l < 7; ++l )
            {
                EXPECT_NE( tFaces( k )->id(), tFaces( l )->id() );
            }
            EXPECT_NE( tFaces( k )->id(), tA->id() );
            EXPECT_NE( tFaces( k )->id(), tB->id() );
            for ( Edge * tEdge : tMesh->edges() )
            {
                EXPECT_NE( tFaces( k )->id(), tEdge->id() );
            }
        }

        uint tNumShared = 0 ;

        for ( Face * tFace : tFaces )
        {
            ASSERT_NE( tFace->master(), nullptr );
            EXPECT_NE( tFace->master(), tFace->slave() );

            if ( tFace->number_of_elements() == 2 )
            {
                ++tNumShared ;
                EXPECT_EQ( tFace->master(), tA );
                EXPECT_EQ( tFace->index_on_master(), 1u );
                EXPECT_EQ( tFace->slave(), tB );
                EXPECT_EQ( tFace->index_on_slave(), 3u );
                EXPECT_EQ( tFace->orientation_on_slave(), 1u );

                // both facets are the same corner set
                Cell< Node * > tMasterCorners ;
                Cell< Node * > tSlaveCorners ;
                tA->get_corner_nodes_of_facet( 1, tMasterCorners );
                tB->get_corner_nodes_of_facet( 3, tSlaveCorners );
                ASSERT_EQ( tMasterCorners.size(), tSlaveCorners.size() );
                for ( uint i = 0; i < tMasterCorners.size(); ++i )
                {
                    bool tFound = false ;
                    for ( uint j = 0; j < tSlaveCorners.size(); ++j )
                    {
                        tFound = tFound || ( tMasterCorners( i ) == tSlaveCorners( j ) );
                    }
                    EXPECT_TRUE( tFound );
                }
            }
            else
            {
                EXPECT_EQ( tFace->number_of_elements(), 1u );
                EXPECT_EQ( tFace->slave(), nullptr );
            }

            // the face's nodes are, positionally, the master's facet nodes
            Cell< Node * > tNodes ;
            tFace->master()->get_nodes_of_facet( tFace->index_on_master(), tNodes );
            ASSERT_EQ( tFace->number_of_nodes(), aPair.mNodesPerFace );
            ASSERT_EQ( tNodes.size(), ( index_t ) aPair.mNodesPerFace );
            for ( uint k = 0; k < aPair.mNodesPerFace; ++k )
            {
                EXPECT_EQ( tFace->node( k ), tNodes( k ) );
            }

            // every node lies on the plane through the three corners
            const Node * tC0 = tFace->node( 0 );
            const Node * tC1 = tFace->node( 1 );
            const Node * tC2 = tFace->node( 2 );

            const real tUx = tC1->x() - tC0->x();
            const real tUy = tC1->y() - tC0->y();
            const real tUz = tC1->z() - tC0->z();
            const real tVx = tC2->x() - tC0->x();
            const real tVy = tC2->y() - tC0->y();
            const real tVz = tC2->z() - tC0->z();

            const real tNx = tUy * tVz - tUz * tVy ;
            const real tNy = tUz * tVx - tUx * tVz ;
            const real tNz = tUx * tVy - tUy * tVx ;

            EXPECT_GT( std::sqrt( tNx * tNx + tNy * tNy + tNz * tNz ), 0.1 );

            for ( uint k = 0; k < aPair.mNodesPerFace; ++k )
            {
                const Node * tNode = tFace->node( k );
                const real tDist =   tNx * ( tNode->x() - tC0->x() )
                                   + tNy * ( tNode->y() - tC0->y() )
                                   + tNz * ( tNode->z() - tC0->z() );
                EXPECT_NEAR( tDist, 0.0, tTiny );
            }
        }

        EXPECT_EQ( tNumShared, 1u );

        // back-links: every facet slot of both elements holds a mesh face
        // that names this element at this facet index
        Element * tElements[ 2 ] = { tA, tB };

        for ( Element * tElement : tElements )
        {
            ASSERT_EQ( tElement->number_of_faces(), 4u );

            for ( uint f = 0; f < 4; ++f )
            {
                Face * tFace = tElement->face( f );
                ASSERT_NE( tFace, nullptr );

                const bool tAsMaster = ( tFace->master() == tElement )
                                    && ( tFace->index_on_master() == f );
                const bool tAsSlave  = ( tFace->slave() == tElement )
                                    && ( tFace->index_on_slave() == f );
                EXPECT_TRUE( tAsMaster || tAsSlave )
                    << "element " << tElement->id() << " facet " << f ;

                bool tInMesh = false ;
                for ( Face * tOther : tFaces )
                {
                    tInMesh = tInMesh || ( tOther == tFace );
                }
                EXPECT_TRUE( tInMesh );
            }
        }
    }

//------------------------------------------------------------------------------

    // after the connectivity pass: three edges per face in the master's
    // facet order, an edge of the shared triangle meets three faces and
    // every other edge two, the shared face has six neighbours and each
    // outer face four, lists symmetric and without self or duplicates
    void
    expect_face_connectivity( TetPair & aPair )
    {
        Mesh * tMesh = aPair.mMesh ;

        ASSERT_EQ( tMesh->number_of_edges(), ( index_t ) 9 );
        ASSERT_EQ( tMesh->number_of_faces(), ( index_t ) 7 );

        Cell< Face * > & tFaces = tMesh->faces();

        for ( Face * tFace : tFaces )
        {
            ASSERT_EQ( tFace->number_of_edges(), 3u );

            Cell< Edge * > tEdges ;
            tFace->master()->get_edges_of_facet( tFace->index_on_master(), tEdges );
            ASSERT_EQ( tEdges.size(), ( index_t ) 3 );
            for ( uint e = 0; e < 3; ++e )
            {
                EXPECT_EQ( tFace->edge( e ), tEdges( e ) );
            }
        }

        uint tNumThree = 0 ;
        uint tNumTwo   = 0 ;

        for ( Edge * tEdge : tMesh->edges() )
        {
            const bool tOnShared = aPair.is_shared_face_corner( tEdge->node( 0 ) )
                                && aPair.is_shared_face_corner( tEdge->node( 1 ) );

            const uint tExpected = tOnShared ? 3 : 2 ;
            EXPECT_EQ( tEdge->number_of_faces(), tExpected )
                << "edge " << tEdge->id();

            if ( tEdge->number_of_faces() == 3 ) ++tNumThree ;
            if ( tEdge->number_of_faces() == 2 ) ++tNumTwo ;
        }
        EXPECT_EQ( tNumThree, 3u );
        EXPECT_EQ( tNumTwo, 6u );

        for ( Face * tFace : tFaces )
        {
            const uint tExpected = ( tFace->number_of_elements() == 2 ) ? 6 : 4 ;
            EXPECT_EQ( tFace->number_of_faces(), tExpected )
                << "face " << tFace->id();

            for ( uint n = 0; n < tFace->number_of_faces(); ++n )
            {
                Face * tOther = tFace->face( n );
                ASSERT_NE( tOther, nullptr );
                EXPECT_NE( tOther, tFace );

                for ( uint m = n + 1; m < tFace->number_of_faces(); ++m )
                {
                    EXPECT_NE( tOther, tFace->face( m ) );
                }

                bool tBack = false ;
                for ( uint m = 0; m < tOther->number_of_faces(); ++m )
                {
                    tBack = tBack || ( tOther->face( m ) == tFace );
                }
                EXPECT_TRUE( tBack )
                    << "face " << tFace->id() << " lists " << tOther->id()
                    << " but not the reverse" ;
            }
        }
    }
}

//------------------------------------------------------------------------------

TEST( FaceFactory, TwoTet10sShareExactlyOneFace )
{
    TetPair tPair( 2 );
    tPair.build_faces();
    expect_two_tets_share_one_face( tPair );
}

//------------------------------------------------------------------------------

TEST( FaceFactory, FacesShareEdgesAndNeighboursAsExpected )
{
    TetPair tPair( 2 );
    tPair.build_faces();
    expect_face_connectivity( tPair );
}

//------------------------------------------------------------------------------

TEST( FaceFactory, LowerIdLaterElementBecomesMaster )
{
    TetPair tPair( 2, true );
    tPair.build_faces();
    expect_two_tets_share_one_face( tPair );
}

//------------------------------------------------------------------------------

TEST( FaceFactory, TwoTet4sShareExactlyOneFace )
{
    TetPair tPair( 1 );
    tPair.build_faces();
    expect_two_tets_share_one_face( tPair );

    // a fatal failure inside a helper returns from the helper only
    if ( ::testing::Test::HasFatalFailure() ) return ;

    expect_face_connectivity( tPair );
}

//------------------------------------------------------------------------------

// one HEX8 with a QUAD4 sideset facet on its bottom: the face at that
// facet carries the facet's id, the other five get fresh distinct ids
TEST( FaceFactory, SidesetFacetIdIsCarriedByTheFace )
{
    const id_t tFacetID = 500 ;

    // owner: a fatal assert leaves the body without reaching the delete
    std::unique_ptr< Mesh > tOwner( new Mesh( 3, 0, true ) );
    Mesh * tMesh = tOwner.get();

    Cell< Node * > & tNodes = tMesh->nodes();

    // unit cube, gmsh corner order: bottom quad ccw, then top
    tNodes.push( new Node( 1, 0.0, 0.0, 0.0 ) );
    tNodes.push( new Node( 2, 1.0, 0.0, 0.0 ) );
    tNodes.push( new Node( 3, 1.0, 1.0, 0.0 ) );
    tNodes.push( new Node( 4, 0.0, 1.0, 0.0 ) );
    tNodes.push( new Node( 5, 0.0, 0.0, 1.0 ) );
    tNodes.push( new Node( 6, 1.0, 0.0, 1.0 ) );
    tNodes.push( new Node( 7, 1.0, 1.0, 1.0 ) );
    tNodes.push( new Node( 8, 0.0, 1.0, 1.0 ) );

    ElementFactory tFactory ;

    Element * tHex = tFactory.create_element( ElementType::HEX8, 1 );
    for ( uint k = 0; k < 8; ++k )
    {
        tHex->insert_node( tNodes( k ), k );
    }
    tHex->set_block_id( 1 );

    Block * tBlock = new Block( 1, 1 );
    tMesh->elements().push( tHex );
    tBlock->insert_element( tHex );
    tMesh->blocks().push( tBlock );

    // the facet owns the wrapper, the sideset owns the facet
    const uint tBottom = bottom_facet_index( ElementType::HEX8 );

    Element * tSurf = tFactory.create_element( ElementType::QUAD4, tFacetID );
    Facet * tFacet = new Facet( tSurf );
    tFacet->set_master( tHex, tBottom );
    tFacet->set_sideset_id( 1 );

    SideSet * tSideSet = new SideSet( 1, 1 );
    tSideSet->insert_facet( tFacet );
    tMesh->sidesets().push( tSideSet );

    tMesh->finalize();

    // a non-empty sideset list needs the block list passed alongside
    Vector< id_t > tBlocks( 1, 1 );
    Vector< id_t > tSideSets( 1, 1 );
    tMesh->create_faces( false, tBlocks, tSideSets );

    ASSERT_EQ( tMesh->number_of_faces(), ( index_t ) 6 );

    Cell< Face * > & tFaces = tMesh->faces();

    uint tNumAtFacet = 0 ;

    for ( uint k = 0; k < 6; ++k )
    {
        Face * tFace = tFaces( k );

        EXPECT_EQ( tFace->index(), ( index_t ) k );
        EXPECT_EQ( tFace->master(), tHex );
        EXPECT_EQ( tFace->slave(), nullptr );

        if ( tFace->index_on_master() == tBottom )
        {
            ++tNumAtFacet ;
            EXPECT_EQ( tFace->id(), tFacetID );
        }
        else
        {
            EXPECT_NE( tFace->id(), tFacetID );
        }

        for ( uint l = k + 1; l < 6; ++l )
        {
            EXPECT_NE( tFace->id(), tFaces( l )->id() );
        }
    }

    EXPECT_EQ( tNumAtFacet, 1u );
}

//------------------------------------------------------------------------------

// a sideset facet whose element lies in no selected block has no volume
// facet to become a face of, so create_faces refuses it by id
TEST( FaceFactory, SidesetFacetOutsideTheSelectedBlocks )
{
    const id_t tFacetID = 500 ;

    // owner: the mesh is deleted on the way out of the throw, together
    // with any face created before the refusal
    std::unique_ptr< Mesh > tOwner( new Mesh( 3, 0, true ) );
    Mesh * tMesh = tOwner.get();

    Cell< Node * > & tNodes = tMesh->nodes();

    // two disjoint unit cubes, gmsh corner order: bottom quad ccw, then top
    for ( uint k = 0; k < 2; ++k )
    {
        const real tX = 2.0 * k ;
        const id_t tID = 8 * k ;
        tNodes.push( new Node( tID + 1, tX + 0.0, 0.0, 0.0 ) );
        tNodes.push( new Node( tID + 2, tX + 1.0, 0.0, 0.0 ) );
        tNodes.push( new Node( tID + 3, tX + 1.0, 1.0, 0.0 ) );
        tNodes.push( new Node( tID + 4, tX + 0.0, 1.0, 0.0 ) );
        tNodes.push( new Node( tID + 5, tX + 0.0, 0.0, 1.0 ) );
        tNodes.push( new Node( tID + 6, tX + 1.0, 0.0, 1.0 ) );
        tNodes.push( new Node( tID + 7, tX + 1.0, 1.0, 1.0 ) );
        tNodes.push( new Node( tID + 8, tX + 0.0, 1.0, 1.0 ) );
    }

    ElementFactory tFactory ;

    Element * tHexes[ 2 ] ;
    for ( uint k = 0; k < 2; ++k )
    {
        tHexes[ k ] = tFactory.create_element( ElementType::HEX8, k + 1 );
        for ( uint i = 0; i < 8; ++i )
        {
            tHexes[ k ]->insert_node( tNodes( 8 * k + i ), i );
        }
        tHexes[ k ]->set_block_id( k + 1 );

        Block * tBlock = new Block( k + 1, 1 );
        tMesh->elements().push( tHexes[ k ] );
        tBlock->insert_element( tHexes[ k ] );
        tMesh->blocks().push( tBlock );
    }

    // the facet sits on the top of the second cube, whose block is not selected
    const uint tTop = 5 ;

    Element * tSurf = tFactory.create_element( ElementType::QUAD4, tFacetID );
    Facet * tFacet = new Facet( tSurf );
    tFacet->set_master( tHexes[ 1 ], tTop );
    tFacet->set_sideset_id( 1 );

    SideSet * tSideSet = new SideSet( 1, 1 );
    tSideSet->insert_facet( tFacet );
    tMesh->sidesets().push( tSideSet );

    tMesh->finalize();

    Vector< id_t > tBlocks( 1, 1 );
    Vector< id_t > tSideSets( 1, 1 );

    bool tThrown = false ;
    try
    {
        tMesh->create_faces( false, tBlocks, tSideSets );
    }
    catch( const std::runtime_error & aError )
    {
        tThrown = true ;
        const std::string tMessage( aError.what() );
        EXPECT_NE( tMessage.find( "Sideset facet 500" ), std::string::npos )
            << "message was: " << tMessage ;
        EXPECT_NE( tMessage.find( "owned by no selected block" ), std::string::npos )
            << "message was: " << tMessage ;
    }
    EXPECT_TRUE( tThrown ) << "expected create_faces to refuse the facet" ;

    // the refusal comes before the surface element is given a face
    EXPECT_FALSE( tSurf->has_faces() );
    EXPECT_FALSE( tHexes[ 1 ]->has_faces() );
}
