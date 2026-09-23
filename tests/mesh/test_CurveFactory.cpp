/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for mesh::CurveFactory::sort_end_nodes() on hand-built
 * polylines: a corner-node set whose degrees pass the branch guard but that
 * is not one connected chain or one closed loop is refused by an error naming
 * the curve ( a branch or an isolated node still dies at the older degree
 * guard ), a single open chain sorts to the walk order from its smaller-id tip
 * whatever order it was handed in, and a single closed loop is accepted and
 * walked. On a 2-D mesh, from_2d_sidesets() builds one segment per sideset
 * facet, also where the sideset turns a corner inside one element, and gives
 * every quadratic segment the midpoint of that facet. On 3-D meshes,
 * intersect() pairs only the edges that a facet of each sideset shares and
 * refuses midpoints that disagree, a closed curve comes out as one head-to-
 * tail loop in 2-D and 3-D, thin-shell side curves take their midpoints from
 * their own facet and refuse a closed rim, and mixed-order, non-line and
 * empty 2-D sidesets are refused.
 *
 * The mesh main forces the throwing reaction, so every BELFEM_ERROR on the
 * refusal paths is observable as a std::runtime_error in both builds. A
 * polyline is never sorted twice after a refusal: the sorter leaves node
 * flags dirty when it throws.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "cl_Segment.hpp"
#include "cl_Curve.hpp"
#include "cl_CurveFactory.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    // runs aCall, requires a std::runtime_error whose what() carries every
    // literal part of the BELFEM_ERROR text that is listed
    template< typename F >
    void
    expect_refusal( F && aCall, const Cell< std::string > & aSubstrings )
    {
        bool tThrown = false ;
        try
        {
            aCall();
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true ;
            for ( const std::string & tSubstring : aSubstrings )
            {
                EXPECT_NE( std::string( aError.what() ).find( tSubstring ), std::string::npos )
                    << "message was: " << aError.what();
            }
        }
        EXPECT_TRUE( tThrown ) << "expected a refusal from sort_end_nodes" ;
    }

    // hand-built polyline: an empty mesh for the factory, a curve of LINE2
    // segments and the test's own nodes, none of which are inserted into the
    // mesh. mNodes is the presentation order and stays fixed; sort_end_nodes
    // permutes the working copy mEndNodes. A Segment owns its element and the
    // Curve owns its segments, so the destructor deletes the nodes and the
    // curve. Scope-bound, so a failing expectation cannot skip the deletes
    struct Polyline
    {
        Mesh         * mMesh ;
        CurveFactory * mFactory ;
        Curve        * mCurve ;

        Cell< Node * > mNodes ;
        Cell< Node * > mEndNodes ;

        Matrix< index_t > mAdjacency ;

        Polyline( const id_t aCurveID,
                  const Cell< id_t > & aNodeIDs,
                  const Cell< std::pair< id_t, id_t > > & aSegments ) :
            mMesh( new Mesh( 3, 0, true ) ),
            mFactory( nullptr ),
            mCurve( new Curve( aCurveID, ElementType::LINE2 ) )
        {
            mFactory = new CurveFactory( mMesh );

            mNodes.set_size( aNodeIDs.size(), nullptr );

            index_t tCount = 0 ;
            for ( id_t tID : aNodeIDs )
            {
                // the sorter reads no coordinate; the zig-zag only keeps the
                // node objects distinguishable under a debugger
                mNodes( tCount ) = new Node( tID,
                    ( real ) tCount,
                    ( tCount % 2 == 0 ) ? 0.0 : 1.0,
                    0.0 );
                ++tCount ;
            }

            ElementFactory tFactory ;

            Cell< Segment * > & tSegments = mCurve->segments() ;
            tSegments.set_size( aSegments.size(), nullptr );

            tCount = 0 ;
            for ( auto tPair : aSegments )
            {
                Element * tElement = tFactory.create_element( ElementType::LINE2, tCount + 1 );
                tElement->insert_node( this->node( tPair.first ), 0 );
                tElement->insert_node( this->node( tPair.second ), 1 );

                Segment * tSegment = new Segment( tElement );
                tSegment->set_index( tCount );

                tSegments( tCount++ ) = tSegment ;
            }
        }

        ~Polyline()
        {
            delete mCurve ;

            for ( Node * tNode : mNodes )
            {
                delete tNode ;
            }

            delete mFactory ;
            delete mMesh ;
        }

        Polyline( const Polyline & ) = delete ;
        Polyline( Polyline && ) = delete ;
        Polyline & operator=( const Polyline & ) = delete ;
        Polyline & operator=( Polyline && ) = delete ;

        Node *
        node( const id_t aID )
        {
            for ( Node * tNode : mNodes )
            {
                if ( tNode->id() == aID )
                {
                    return tNode ;
                }
            }
            BELFEM_ERROR( false, "no node %lu in this polyline", ( luint ) aID );
            return nullptr ;
        }

        // restores the sorter's preconditions: the working cell is the
        // presentation order again, every node carries its slot of that cell
        // as index() and is unflagged, and every segment carries its slot of
        // the segment cell
        void
        present()
        {
            mEndNodes = mNodes ;

            index_t tCount = 0 ;
            for ( Node * tNode : mEndNodes )
            {
                tNode->set_index( tCount++ );
                tNode->unflag();
            }

            tCount = 0 ;
            for ( Segment * tSegment : mCurve->segments() )
            {
                tSegment->set_index( tCount++ );
            }
        }

        void
        sort_end_nodes()
        {
            this->present();
            mFactory->sort_end_nodes( mCurve, mAdjacency, mEndNodes );
        }
    };

    // how many segments of aCurve carry exactly the node pair { aOne, aTwo }
    index_t
    count_segments( Curve & aCurve, const id_t aOne, const id_t aTwo )
    {
        index_t tCount = 0 ;

        for ( Segment * tSegment : aCurve.segments() )
        {
            const id_t tFirst = tSegment->node( 0 )->id() ;
            const id_t tSecond = tSegment->node( 1 )->id() ;

            if ( ( tFirst == aOne && tSecond == aTwo ) ||
                 ( tFirst == aTwo && tSecond == aOne ) )
            {
                ++tCount ;
            }
        }

        return tCount ;
    }

    // the chain 40-10-50-20-60-30 walked from its smaller-id tip
    const Cell< id_t > tChain = { 30, 60, 20, 50, 10, 40 } ;

    // row k of the returned adjacency holds the positions ( in the segment
    // cell ) of the segments incident to the node now at index k, gNoIndex
    // in the second slot of a tip; the next sorter reads these rows
    void
    check_adjacency( Polyline & aLine )
    {
        const index_t tN = aLine.mEndNodes.size() ;
        ASSERT_EQ( aLine.mAdjacency.n_rows(), tN );
        ASSERT_EQ( aLine.mAdjacency.n_cols(), ( index_t ) 2 );

        Cell< Segment * > & tSegments = aLine.mCurve->segments() ;

        for ( index_t k=0; k<tN; ++k )
        {
            const id_t tID = aLine.mEndNodes( k )->id() ;

            for ( index_t j=0; j<2; ++j )
            {
                const index_t tPos = aLine.mAdjacency( k, j );

                if ( tPos == gNoIndex )
                {
                    // only a tip has an empty slot
                    EXPECT_TRUE( k == 0 || k == tN - 1 )
                        << "empty adjacency slot at interior index " << k ;
                    continue ;
                }

                ASSERT_LT( tPos, tSegments.size() );
                Segment * tSegment = tSegments( tPos );
                EXPECT_TRUE( tSegment->node( 0 )->id() == tID ||
                             tSegment->node( 1 )->id() == tID )
                    << "adjacency row " << k << " names segment " << tPos
                    << " which does not touch node " << tID ;
            }

            EXPECT_NE( aLine.mAdjacency( k, 0 ), aLine.mAdjacency( k, 1 ) );
        }
    }

    void
    check_chain( Polyline & aLine )
    {
        EXPECT_FALSE( aLine.mCurve->is_closed() );

        ASSERT_EQ( aLine.mEndNodes.size(), tChain.size() );

        for ( index_t k=0; k<tChain.size(); ++k )
        {
            EXPECT_EQ( aLine.mEndNodes( k )->id(), tChain( k ) );
            EXPECT_EQ( aLine.mEndNodes( k )->index(), k );
        }

        for ( index_t k=1; k<tChain.size(); ++k )
        {
            EXPECT_EQ( count_segments( *aLine.mCurve, tChain( k-1 ), tChain( k ) ),
                       ( index_t ) 1 );
        }

        check_adjacency( aLine );
    }
}

//------------------------------------------------------------------------------

TEST( CurveFactory, TwoChainsAreRefusedNamingTheCurve )
{
    // 1-5-3 and 2-6-4: four tips and two interior nodes pass the degree
    // guard, but the walk only covers the component of the smaller-id tip
    Polyline tLine( 7,
        { 1, 5, 3, 2, 6, 4 },
        { { 1, 5 }, { 5, 3 }, { 2, 6 }, { 6, 4 } } );

    expect_refusal( [ & ]() { tLine.sort_end_nodes() ; },
        { "Curve 7", "not one connected chain" } );
}

//------------------------------------------------------------------------------

TEST( CurveFactory, LoopAndChainAreRefused )
{
    // 1-5-3 open beside the closed 2-6-4-2: two tips and four interior nodes
    Polyline tLine( 7,
        { 1, 5, 3, 2, 6, 4 },
        { { 1, 5 }, { 5, 3 }, { 2, 6 }, { 6, 4 }, { 4, 2 } } );

    expect_refusal( [ & ]() { tLine.sort_end_nodes() ; },
        { "Curve 7", "not one connected chain" } );
}

//------------------------------------------------------------------------------

TEST( CurveFactory, TwoLoopsAreRefused )
{
    // two closed triangles: no tip at all, so the curve reads as closed
    Polyline tLine( 7,
        { 1, 2, 3, 4, 5, 6 },
        { { 1, 2 }, { 2, 3 }, { 3, 1 }, { 4, 5 }, { 5, 6 }, { 6, 4 } } );

    expect_refusal( [ & ]() { tLine.sort_end_nodes() ; },
        { "Curve 7", "not one connected chain" } );
}

//------------------------------------------------------------------------------

TEST( CurveFactory, SingleChainIsMonotone )
{
    // the smallest id sits in the interior, so the expected order cannot be
    // mistaken for a sort by node id
    Polyline tLine( 3,
        { 20, 60, 40, 30, 10, 50 },
        { { 40, 10 }, { 10, 50 }, { 50, 20 }, { 20, 60 }, { 60, 30 } } );

    tLine.sort_end_nodes();

    check_chain( tLine );
}

//------------------------------------------------------------------------------

TEST( CurveFactory, ReversalGivesTheIdenticalChain )
{
    // the same chain with the node cell, the segment cell and every segment's
    // node pair reversed: the start is the smaller-id tip, not an input slot,
    // so the result is the chain itself and not its reverse
    Polyline tLine( 3,
        { 50, 10, 30, 40, 60, 20 },
        { { 30, 60 }, { 60, 20 }, { 20, 50 }, { 50, 10 }, { 10, 40 } } );

    tLine.sort_end_nodes();

    check_chain( tLine );
}

//------------------------------------------------------------------------------

TEST( CurveFactory, SingleLoopIsAccepted )
{
    Polyline tLine( 3,
        { 11, 12, 13, 14 },
        { { 11, 12 }, { 12, 13 }, { 13, 14 }, { 14, 11 } } );

    tLine.sort_end_nodes();

    EXPECT_TRUE( tLine.mCurve->is_closed() );

    // the walk starts at the smallest id and takes the larger-id neighbour
    // first, so the loop reads 11-14-13-12 and consecutive nodes share a
    // segment, the last one closing back to the first
    const Cell< id_t > tLoop = { 11, 14, 13, 12 } ;

    ASSERT_EQ( tLine.mEndNodes.size(), tLoop.size() );

    for ( index_t k=0; k<tLoop.size(); ++k )
    {
        EXPECT_EQ( tLine.mEndNodes( k )->id(), tLoop( k ) );
        EXPECT_EQ( tLine.mEndNodes( k )->index(), k );
        EXPECT_EQ( count_segments( *tLine.mCurve, tLoop( k ),
                                   tLoop( ( k + 1 ) % tLoop.size() ) ),
                   ( index_t ) 1 );
    }

    check_adjacency( tLine );
}

//------------------------------------------------------------------------------

// a 2-D curve carries no second sideset; its order follows the one it has,
// and the 3-D rule ( both sidesets quadratic ) is unchanged
TEST( CurveFactory, TwoDCurveOrderFollowsItsOnlySideset )
{
    ElementFactory tFactory ;

    SideSet tLinear( 1, 1 );
    tLinear.insert_facet( new Facet( tFactory.create_element( ElementType::LINE2, 1 ) ) );

    SideSet tQuadratic( 2, 1 );
    tQuadratic.insert_facet( new Facet( tFactory.create_element( ElementType::LINE3, 2 ) ) );

    ASSERT_EQ( tQuadratic.element_type(), ElementType::LINE3 );

    Curve tQuadCurve( 1, &tQuadratic, nullptr );
    EXPECT_EQ( tQuadCurve.element_type(), ElementType::LINE3 );

    Curve tLinCurve( 2, &tLinear, nullptr );
    EXPECT_EQ( tLinCurve.element_type(), ElementType::LINE2 );

    Curve tMixed( 3, &tQuadratic, &tLinear );
    EXPECT_EQ( tMixed.element_type(), ElementType::LINE2 );

    Curve tBoth( 4, &tQuadratic, &tQuadratic );
    EXPECT_EQ( tBoth.element_type(), ElementType::LINE3 );
}

//------------------------------------------------------------------------------

// two QUAD8 side by side with a LINE3 sideset on their bottom edges. The
// factory orders its pairs by a key over the end nodes' slots in the cell
// that unique() sorts by address, so the facets are inserted in descending
// key order: a midpoint gathered in facet order would land on the wrong
// segment, and a midpoint written into the Segment instead of its line
// element is never seen by node( 2 )
TEST( CurveFactory, TwoDQuadraticSegmentsTakeTheirFacetsMidpoints )
{
    // owner: a fatal assert leaves the body without reaching the delete
    std::unique_ptr< Mesh > tOwner( new Mesh( 2, 0, true ) );
    Mesh * tMesh = tOwner.get();

    Cell< Node * > & tNodes = tMesh->nodes();

    auto push = [ &tNodes ]( const real aX, const real aY ) -> Node *
    {
        Node * tNode = new Node( ( id_t ) tNodes.size() + 1, aX, aY, 0.0 );
        tNodes.push( tNode );
        return tNode ;
    };

    // corners: bottom row P0 P1 P2, top row P5 P4 P3
    Node * tP0 = push( 0.0, 0.0 );
    Node * tP1 = push( 1.0, 0.0 );
    Node * tP2 = push( 2.0, 0.0 );
    Node * tP3 = push( 2.0, 1.0 );
    Node * tP4 = push( 1.0, 1.0 );
    Node * tP5 = push( 0.0, 1.0 );

    Node * tM01 = push( 0.5, 0.0 );
    Node * tM12 = push( 1.5, 0.0 );
    Node * tM14 = push( 1.0, 0.5 );
    Node * tM45 = push( 0.5, 1.0 );
    Node * tM50 = push( 0.0, 0.5 );
    Node * tM23 = push( 2.0, 0.5 );
    Node * tM34 = push( 1.5, 1.0 );

    ElementFactory tFactory ;

    // QUAD8: corners ccw, then the midsides of edges 0-1, 1-2, 2-3, 3-0
    Element * tA = tFactory.create_element( ElementType::QUAD8, 1 );
    Node * tNodesA[ 8 ] = { tP0, tP1, tP4, tP5, tM01, tM14, tM45, tM50 };

    Element * tB = tFactory.create_element( ElementType::QUAD8, 2 );
    Node * tNodesB[ 8 ] = { tP1, tP2, tP3, tP4, tM12, tM23, tM34, tM14 };

    for ( uint k = 0; k < 8; ++k )
    {
        tA->insert_node( tNodesA[ k ], k );
        tB->insert_node( tNodesB[ k ], k );
    }

    Block * tBlock = new Block( 1, 2 );
    for ( Element * tElement : { tA, tB } )
    {
        tElement->set_block_id( 1 );
        tMesh->elements().push( tElement );
        tBlock->insert_element( tElement );
    }
    tMesh->blocks().push( tBlock );

    // slot of each end node in the address-sorted end-node cell
    Cell< Node * > tSorted = { tP0, tP1, tP2 };
    std::sort( tSorted.vector_data().begin(), tSorted.vector_data().end() );

    auto slot = [ &tSorted ]( const Node * aNode ) -> belfem::key_t
    {
        for ( index_t k = 0; k < tSorted.size(); ++k )
        {
            if ( tSorted( k ) == aNode )
            {
                return k ;
            }
        }
        return gNoIndex ;
    };

    auto key = [ &slot ]( const Node * aOne, const Node * aTwo ) -> belfem::key_t
    {
        const belfem::key_t tOne = slot( aOne );
        const belfem::key_t tTwo = slot( aTwo );
        return tOne > tTwo ? tOne * 3 + tTwo : tTwo * 3 + tOne ;
    };

    // facet 0 of either QUAD8 is its bottom edge: corners 0, 1, midside 4
    Facet * tFacetA = new Facet( tFactory.create_element( ElementType::LINE3, 101 ) );
    tFacetA->set_master( tA, 0 );
    tFacetA->set_sideset_id( 1 );

    Facet * tFacetB = new Facet( tFactory.create_element( ElementType::LINE3, 102 ) );
    tFacetB->set_master( tB, 0 );
    tFacetB->set_sideset_id( 1 );

    const bool tAFirst = key( tP0, tP1 ) > key( tP1, tP2 );

    SideSet * tSideSet = new SideSet( 1, 2 );
    tSideSet->insert_facet( tAFirst ? tFacetA : tFacetB );
    tSideSet->insert_facet( tAFirst ? tFacetB : tFacetA );
    tMesh->sidesets().push( tSideSet );

    tMesh->finalize();

    CurveFactory tCurveFactory( tMesh );

    Vector< id_t > tSideSets( 1, 1 );
    Curve * tCurve = tCurveFactory.from_2d_sidesets( tSideSets, 7 );
    tMesh->curves().push( tCurve );

    ASSERT_EQ( tCurve->element_type(), ElementType::LINE3 );
    ASSERT_EQ( tCurve->segments().size(), ( index_t ) 2 );

    for ( Segment * tSegment : tCurve->segments() )
    {
        const Node * tFirst  = tSegment->node( 0 );
        const Node * tSecond = tSegment->node( 1 );

        const bool tIsA = ( tFirst == tP0 && tSecond == tP1 ) || ( tFirst == tP1 && tSecond == tP0 );
        const bool tIsB = ( tFirst == tP1 && tSecond == tP2 ) || ( tFirst == tP2 && tSecond == tP1 );

        ASSERT_TRUE( tIsA || tIsB ) << "segment " << tFirst->id() << "-" << tSecond->id()
                                    << " joins no facet" ;

        ASSERT_NE( tSegment->node( 2 ), nullptr )
            << "segment " << tFirst->id() << "-" << tSecond->id() << " has no midpoint" ;

        EXPECT_EQ( tSegment->node( 2 ), tIsA ? tM01 : tM12 )
            << "segment " << tFirst->id() << "-" << tSecond->id()
            << " carries midpoint " << tSegment->node( 2 )->id() ;
    }
}

//------------------------------------------------------------------------------

namespace
{
    // one QUAD4 or QUAD8 with a sideset on its edges 0 and 1, so the single
    // element holds all three end nodes of the curve. Builds the curve and
    // checks that it has one segment per facet, each joining that facet's
    // corners and, if quadratic, carrying its midside
    void
    check_corner_curve( const bool aQuadratic )
    {
        std::unique_ptr< Mesh > tOwner( new Mesh( 2, 0, true ) );
        Mesh * tMesh = tOwner.get();

        Cell< Node * > & tNodes = tMesh->nodes();

        auto push = [ &tNodes ]( const real aX, const real aY ) -> Node *
        {
            Node * tNode = new Node( ( id_t ) tNodes.size() + 1, aX, aY, 0.0 );
            tNodes.push( tNode );
            return tNode ;
        };

        Node * tP0 = push( 0.0, 0.0 );
        Node * tP1 = push( 1.0, 0.0 );
        Node * tP2 = push( 1.0, 1.0 );
        Node * tP3 = push( 0.0, 1.0 );

        ElementFactory tFactory ;

        Element * tQuad = tFactory.create_element(
            aQuadratic ? ElementType::QUAD8 : ElementType::QUAD4, 1 );

        Node * tCorners[ 4 ] = { tP0, tP1, tP2, tP3 };
        for ( uint k = 0; k < 4; ++k )
        {
            tQuad->insert_node( tCorners[ k ], k );
        }

        Node * tM01 = nullptr ;
        Node * tM12 = nullptr ;
        if ( aQuadratic )
        {
            tM01 = push( 0.5, 0.0 );
            tM12 = push( 1.0, 0.5 );
            Node * tM23 = push( 0.5, 1.0 );
            Node * tM30 = push( 0.0, 0.5 );

            Node * tMidsides[ 4 ] = { tM01, tM12, tM23, tM30 };
            for ( uint k = 0; k < 4; ++k )
            {
                tQuad->insert_node( tMidsides[ k ], k + 4 );
            }
        }

        Block * tBlock = new Block( 1, 1 );
        tQuad->set_block_id( 1 );
        tMesh->elements().push( tQuad );
        tBlock->insert_element( tQuad );
        tMesh->blocks().push( tBlock );

        const ElementType tLine = aQuadratic ? ElementType::LINE3 : ElementType::LINE2 ;

        SideSet * tSideSet = new SideSet( 1, 2 );
        for ( uint k = 0; k < 2; ++k )
        {
            Facet * tFacet = new Facet( tFactory.create_element( tLine, 101 + k ) );
            tFacet->set_master( tQuad, k );
            tFacet->set_sideset_id( 1 );
            tSideSet->insert_facet( tFacet );
        }
        tMesh->sidesets().push( tSideSet );

        tMesh->finalize();

        CurveFactory tCurveFactory( tMesh );

        Vector< id_t > tSideSets( 1, 1 );
        Curve * tCurve = tCurveFactory.from_2d_sidesets( tSideSets, 7 );
        tMesh->curves().push( tCurve );

        ASSERT_EQ( tCurve->element_type(), tLine );
        ASSERT_EQ( tCurve->segments().size(), tSideSet->number_of_facets() );
        EXPECT_FALSE( tCurve->is_closed() );

        for ( Segment * tSegment : tCurve->segments() )
        {
            const Node * tFirst  = tSegment->node( 0 );
            const Node * tSecond = tSegment->node( 1 );

            const bool tIs01 = ( tFirst == tP0 && tSecond == tP1 ) || ( tFirst == tP1 && tSecond == tP0 );
            const bool tIs12 = ( tFirst == tP1 && tSecond == tP2 ) || ( tFirst == tP2 && tSecond == tP1 );

            ASSERT_TRUE( tIs01 || tIs12 ) << "segment " << tFirst->id() << "-" << tSecond->id()
                                          << " joins no facet" ;

            if ( aQuadratic )
            {
                ASSERT_NE( tSegment->node( 2 ), nullptr )
                    << "segment " << tFirst->id() << "-" << tSecond->id() << " has no midpoint" ;

                EXPECT_EQ( tSegment->node( 2 ), tIs01 ? tM01 : tM12 )
                    << "segment " << tFirst->id() << "-" << tSecond->id()
                    << " carries midpoint " << tSegment->node( 2 )->id() ;
            }
        }
    }
}

// a sideset that turns a corner inside one element: the element holds three
// end nodes, and the curve must not gain the pair of nodes 0 and 2, which
// share the element but no facet
TEST( CurveFactory, TwoDCornerInOneLinearElementGivesOneSegmentPerFacet )
{
    check_corner_curve( false );
}

TEST( CurveFactory, TwoDCornerInOneQuadraticElementGivesOneSegmentPerFacet )
{
    check_corner_curve( true );
}

//------------------------------------------------------------------------------

namespace
{
    // a hand-built mesh: nodes numbered 1, 2, ... in push order, all elements
    // in block 1, and sidesets of facets linked to their master elements.
    // The mesh deletes everything pushed into it
    struct Fixture
    {
        std::unique_ptr< Mesh > mOwner ;
        Mesh * mMesh ;
        ElementFactory mFactory ;
        id_t mFacetID = 1000 ;

        Fixture( const uint aDimension, const bool aConnectivities ) :
            mOwner( new Mesh( aDimension, 0, aConnectivities ) ),
            mMesh( mOwner.get() )
        {
        }

        Node *
        node( const real aX, const real aY, const real aZ )
        {
            Node * tNode = new Node( ( id_t ) mMesh->nodes().size() + 1, aX, aY, aZ );
            mMesh->nodes().push( tNode );
            return tNode ;
        }

        Element *
        element( const ElementType aType, const Cell< Node * > & aNodes )
        {
            Element * tElement = mFactory.create_element( aType, ( id_t ) mMesh->elements().size() + 1 );

            uint tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tElement->insert_node( tNode, tCount++ );
            }

            tElement->set_block_id( 1 );
            mMesh->elements().push( tElement );
            return tElement ;
        }

        // a HEX8 or HEX20 on the corners in HEX8 order ( bottom 0-3, top
        // 4-7 ); a HEX20 gets twelve new edge midpoints, in slots 8 ... 19
        Element *
        hex( const Cell< Node * > & aCorners, const bool aQuadratic )
        {
            if ( ! aQuadratic )
            {
                return this->element( ElementType::HEX8, aCorners );
            }

            const uint tEdges[ 12 ][ 2 ] = {
                { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
                { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
                { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 } };

            Cell< Node * > tNodes( 20, nullptr );
            for ( uint k = 0; k < 8; ++k )
            {
                tNodes( k ) = aCorners( k );
            }
            for ( uint k = 0; k < 12; ++k )
            {
                const Node * tA = aCorners( tEdges[ k ][ 0 ] );
                const Node * tB = aCorners( tEdges[ k ][ 1 ] );
                tNodes( 8 + k ) = this->node(
                    0.5 * ( tA->x() + tB->x() ),
                    0.5 * ( tA->y() + tB->y() ),
                    0.5 * ( tA->z() + tB->z() ) );
            }

            return this->element( ElementType::HEX20, tNodes );
        }

        // a sideset of facet elements of type aType, facet k on facet
        // aFacets( k ).second of element aFacets( k ).first. Before
        // finalize() only
        void
        sideset( const id_t aID, const ElementType aType,
                 const Cell< std::pair< Element *, uint > > & aFacets )
        {
            SideSet * tSideSet = new SideSet( aID, aFacets.size() );

            for ( const std::pair< Element *, uint > & tMaster : aFacets )
            {
                Facet * tFacet = new Facet( mFactory.create_element( aType, ++mFacetID ) );
                tFacet->set_master( tMaster.first, tMaster.second );
                tFacet->set_sideset_id( aID );
                tSideSet->insert_facet( tFacet );
            }

            mMesh->sidesets().push( tSideSet );
        }

        void
        finalize()
        {
            Block * tBlock = new Block( 1, mMesh->elements().size() );
            for ( Element * tElement : mMesh->elements() )
            {
                tBlock->insert_element( tElement );
            }
            mMesh->blocks().push( tBlock );
            mMesh->finalize();
        }
    };

    // the unit cube: bottom 0-3 counter-clockwise at z = 0, top 4-7 above
    Cell< Node * >
    unit_cube( Fixture & aFixture )
    {
        return {
            aFixture.node( 0.0, 0.0, 0.0 ), aFixture.node( 1.0, 0.0, 0.0 ),
            aFixture.node( 1.0, 1.0, 0.0 ), aFixture.node( 0.0, 1.0, 0.0 ),
            aFixture.node( 0.0, 0.0, 1.0 ), aFixture.node( 1.0, 0.0, 1.0 ),
            aFixture.node( 1.0, 1.0, 1.0 ), aFixture.node( 0.0, 1.0, 1.0 ) };
    }

    // aCurve runs through aNodes in this order ( corner, midpoint, corner,
    // ... for LINE3; the start repeated at the end of a closed curve ), its
    // segments run head to tail along it, and every edge has unit length
    void
    check_curve( const Curve * aCurve, const Cell< Node * > & aNodes, const bool aClosed )
    {
        const bool tQuadratic = aCurve->element_type() == ElementType::LINE3 ;
        const index_t tStep = tQuadratic ? 2 : 1 ;

        EXPECT_EQ( aCurve->is_closed(), aClosed );

        ASSERT_EQ( aCurve->nodes().size(), aNodes.size() );
        ASSERT_EQ( aCurve->segments().size(), ( aNodes.size() - 1 ) / tStep );

        for ( index_t k = 0; k < aNodes.size(); ++k )
        {
            EXPECT_EQ( aCurve->nodes()( k ), aNodes( k ) )
                << "node " << k << " is " << aCurve->nodes()( k )->id()
                << ", expected " << aNodes( k )->id() ;
        }

        for ( index_t k = 0; k < aCurve->segments().size(); ++k )
        {
            const Segment * tSegment = aCurve->segments()( k );

            EXPECT_EQ( tSegment->node( 0 ), aNodes( k * tStep ) ) << "segment " << k ;
            EXPECT_EQ( tSegment->node( 1 ), aNodes( ( k + 1 ) * tStep ) ) << "segment " << k ;

            if ( tQuadratic )
            {
                EXPECT_EQ( tSegment->node( 2 ), aNodes( k * tStep + 1 ) ) << "segment " << k ;
            }
        }

        ASSERT_EQ( aCurve->arclength().length(), aNodes.size() );
        for ( index_t k = 0; k < aNodes.size(); ++k )
        {
            EXPECT_NEAR( aCurve->arclength()( k ), ( real ) k / ( real ) tStep, 1e-12 )
                << "arclength at node " << k ;
        }
    }

    // one HEX8 or HEX20: A is the top face, B the faces aFacesB. Builds
    // the intersection curve and returns it, owned by the mesh
    Curve *
    hex_curve( Fixture & aFixture, Element * aHex, const bool aQuadratic, const Cell< uint > & aFacesB )
    {
        const ElementType tFacet = aQuadratic ? ElementType::QUAD8 : ElementType::QUAD4 ;

        aFixture.sideset( 1, tFacet, { { aHex, 5 } } );

        Cell< std::pair< Element *, uint > > tB ;
        for ( uint tFace : aFacesB )
        {
            tB.push( std::make_pair( aHex, tFace ) );
        }
        aFixture.sideset( 2, tFacet, tB );

        aFixture.finalize();

        CurveFactory tFactory( aFixture.mMesh );
        Curve * tCurve = tFactory.intersect( 1, 2, 7 );
        aFixture.mMesh->curves().push( tCurve );
        return tCurve ;
    }

    // the top face meets the front and the left face: the curve turns the
    // corner at node 4, and the top face diagonal 5-7, whose nodes share
    // the element but no facet edge of A and B, is not a pair
    void
    check_hex_corner( const bool aQuadratic )
    {
        Fixture tFixture( 3, true );
        Cell< Node * > tP = unit_cube( tFixture );
        Element * tHex = tFixture.hex( tP, aQuadratic );

        const Curve * tCurve = hex_curve( tFixture, tHex, aQuadratic, { 0, 3 } );

        ASSERT_EQ( tCurve->element_type(), aQuadratic ? ElementType::LINE3 : ElementType::LINE2 );

        // the walk starts at the smaller-id tip, node 5
        if ( aQuadratic )
        {
            check_curve( tCurve, { tP( 5 ), tHex->node( 16 ), tP( 4 ), tHex->node( 19 ), tP( 7 ) }, false );
        }
        else
        {
            check_curve( tCurve, { tP( 5 ), tP( 4 ), tP( 7 ) }, false );
        }
    }

    // the top face meets all four side faces: one closed loop around the
    // rim, walked from the smallest id toward its larger-id neighbour
    void
    check_hex_rim( const bool aQuadratic )
    {
        Fixture tFixture( 3, true );
        Cell< Node * > tP = unit_cube( tFixture );
        Element * tHex = tFixture.hex( tP, aQuadratic );

        const Curve * tCurve = hex_curve( tFixture, tHex, aQuadratic, { 0, 1, 2, 3 } );

        if ( aQuadratic )
        {
            check_curve( tCurve,
                { tP( 4 ), tHex->node( 19 ), tP( 7 ), tHex->node( 18 ), tP( 6 ),
                  tHex->node( 17 ), tP( 5 ), tHex->node( 16 ), tP( 4 ) }, true );
        }
        else
        {
            check_curve( tCurve, { tP( 4 ), tP( 7 ), tP( 6 ), tP( 5 ), tP( 4 ) }, true );
        }
    }

    // one QUAD4 or QUAD8 on the unit square with a sideset on all four
    // edges: one closed loop
    void
    check_square_loop( const bool aQuadratic )
    {
        Fixture tFixture( 2, true );

        Node * tP0 = tFixture.node( 0.0, 0.0, 0.0 );
        Node * tP1 = tFixture.node( 1.0, 0.0, 0.0 );
        Node * tP2 = tFixture.node( 1.0, 1.0, 0.0 );
        Node * tP3 = tFixture.node( 0.0, 1.0, 0.0 );

        Cell< Node * > tNodes = { tP0, tP1, tP2, tP3 };

        Node * tM01 = nullptr ;
        Node * tM12 = nullptr ;
        Node * tM23 = nullptr ;
        Node * tM30 = nullptr ;
        if ( aQuadratic )
        {
            tM01 = tFixture.node( 0.5, 0.0, 0.0 );
            tM12 = tFixture.node( 1.0, 0.5, 0.0 );
            tM23 = tFixture.node( 0.5, 1.0, 0.0 );
            tM30 = tFixture.node( 0.0, 0.5, 0.0 );
            tNodes = { tP0, tP1, tP2, tP3, tM01, tM12, tM23, tM30 };
        }

        Element * tQuad = tFixture.element( aQuadratic ? ElementType::QUAD8 : ElementType::QUAD4, tNodes );

        tFixture.sideset( 1, aQuadratic ? ElementType::LINE3 : ElementType::LINE2,
            { { tQuad, 0 }, { tQuad, 1 }, { tQuad, 2 }, { tQuad, 3 } } );

        tFixture.finalize();

        CurveFactory tFactory( tFixture.mMesh );
        Vector< id_t > tSideSets( 1, 1 );
        Curve * tCurve = tFactory.from_2d_sidesets( tSideSets, 7 );
        tFixture.mMesh->curves().push( tCurve );

        if ( aQuadratic )
        {
            check_curve( tCurve, { tP0, tM30, tP3, tM23, tP2, tM12, tP1, tM01, tP0 }, true );
        }
        else
        {
            check_curve( tCurve, { tP0, tP3, tP2, tP1, tP0 }, true );
        }
    }
}

//------------------------------------------------------------------------------

TEST( CurveFactory, ThreeDCornerInOneLinearHexGivesTwoSegments )
{
    check_hex_corner( false );
}

TEST( CurveFactory, ThreeDCornerInOneQuadraticHexGivesTwoMidpoints )
{
    check_hex_corner( true );
}

//------------------------------------------------------------------------------

// the triangle of A has all three corners on B, but its edge 0-2 lies on no
// facet of B: the curve is 0-1-2, not a triangle
TEST( CurveFactory, ThreeDEarFacetOfOneTetGivesTwoSegments )
{
    Fixture tFixture( 3, true );

    Node * tP0 = tFixture.node( 0.0, 0.0, 0.0 );
    Node * tP1 = tFixture.node( 1.0, 0.0, 0.0 );
    Node * tP2 = tFixture.node( 1.0, 1.0, 0.0 );
    Node * tP3 = tFixture.node( 0.0, 0.0, 1.0 );

    Element * tTet = tFixture.element( ElementType::TET4, { tP0, tP1, tP2, tP3 } );

    // TET4 facets: 0 = 0-1-3, 1 = 1-2-3, 3 = 0-2-1
    tFixture.sideset( 1, ElementType::TRI3, { { tTet, 3 } } );
    tFixture.sideset( 2, ElementType::TRI3, { { tTet, 0 }, { tTet, 1 } } );
    tFixture.finalize();

    CurveFactory tFactory( tFixture.mMesh );
    Curve * tCurve = tFactory.intersect( 1, 2, 7 );
    tFixture.mMesh->curves().push( tCurve );

    check_curve( tCurve, { tP0, tP1, tP2 }, false );
}

//------------------------------------------------------------------------------

TEST( CurveFactory, ThreeDTopRimOfOneLinearHexIsOneClosedLoop )
{
    check_hex_rim( false );
}

TEST( CurveFactory, ThreeDTopRimOfOneQuadraticHexIsOneClosedLoop )
{
    check_hex_rim( true );
}

//------------------------------------------------------------------------------

// two HEX20 side by side that share the corners of a face but not its edge
// midpoints. B holds that face from both sides, so its two facets join
// nodes 4 and 5 through different midpoints: refused whichever facet of B
// comes first, including the order in which the last one agrees with A
TEST( CurveFactory, ThreeDSharedEdgeMidpointsMustAgreeInEitherOrder )
{
    Fixture tFixture( 3, true );
    Cell< Node * > tP = unit_cube( tFixture );
    Element * tFront = tFixture.hex( tP, true );

    // the neighbour at y < 0; its facet 2 ( 2-3-7-6 ) is facet 0 of tFront
    Cell< Node * > tQ = {
        tFixture.node( 0.0, -1.0, 0.0 ), tFixture.node( 1.0, -1.0, 0.0 ), tP( 1 ), tP( 0 ),
        tFixture.node( 0.0, -1.0, 1.0 ), tFixture.node( 1.0, -1.0, 1.0 ), tP( 5 ), tP( 4 ) };
    Element * tBack = tFixture.hex( tQ, true );

    tFixture.sideset( 1, ElementType::QUAD8, { { tFront, 5 } } );
    tFixture.sideset( 2, ElementType::QUAD8, { { tFront, 0 }, { tBack, 2 } } );
    tFixture.sideset( 3, ElementType::QUAD8, { { tBack, 2 }, { tFront, 0 } } );
    tFixture.finalize();

    CurveFactory tFactory( tFixture.mMesh );

    expect_refusal( [ & ]() { delete tFactory.intersect( 1, 2, 7 ) ; },
        { "different midpoints" } );

    expect_refusal( [ & ]() { delete tFactory.intersect( 1, 3, 8 ) ; },
        { "different midpoints" } );
}

//------------------------------------------------------------------------------

TEST( CurveFactory, ThreeDDisjointSidesetsAreRefused )
{
    Fixture tFixture( 3, true );
    Cell< Node * > tP = unit_cube( tFixture );
    Element * tHex = tFixture.hex( tP, false );

    tFixture.sideset( 1, ElementType::QUAD4, { { tHex, 5 } } );
    tFixture.sideset( 2, ElementType::QUAD4, { { tHex, 4 } } );
    tFixture.finalize();

    CurveFactory tFactory( tFixture.mMesh );

    expect_refusal( [ & ]() { delete tFactory.intersect( 1, 2, 7 ) ; },
        { "Curve 7", "share no facet edge" } );
}

//------------------------------------------------------------------------------

TEST( CurveFactory, TwoDClosedLoopLinear )
{
    check_square_loop( false );
}

TEST( CurveFactory, TwoDClosedLoopQuadratic )
{
    check_square_loop( true );
}

//------------------------------------------------------------------------------

namespace
{
    // one QUAD8 on the unit square with a LINE3 sideset 2 on its edge 1;
    // the test adds its own sidesets after finalize()
    struct QuadEight
    {
        Fixture mFixture ;
        Cell< Node * > mP ;
        Element * mQuad ;

        QuadEight() :
            mFixture( 2, true )
        {
            mP = { mFixture.node( 0.0, 0.0, 0.0 ), mFixture.node( 1.0, 0.0, 0.0 ),
                   mFixture.node( 1.0, 1.0, 0.0 ), mFixture.node( 0.0, 1.0, 0.0 ),
                   mFixture.node( 0.5, 0.0, 0.0 ), mFixture.node( 1.0, 0.5, 0.0 ),
                   mFixture.node( 0.5, 1.0, 0.0 ), mFixture.node( 0.0, 0.5, 0.0 ) };

            mQuad = mFixture.element( ElementType::QUAD8, mP );
            mFixture.sideset( 2, ElementType::LINE3, { { mQuad, 1 } } );
            mFixture.finalize();
        }

        // registers a sideset of one facet of type aType carrying aNodes,
        // with no master element
        void
        add_sideset( const id_t aID, const ElementType aType, const Cell< Node * > & aNodes )
        {
            SideSet * tSideSet = new SideSet( aID, 1 );
            Facet * tFacet = new Facet( mFixture.mFactory.create_element( aType, ++mFixture.mFacetID ) );

            uint tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tFacet->element()->insert_node( tNode, tCount++ );
            }

            tFacet->set_sideset_id( aID );
            tSideSet->insert_facet( tFacet );
            mFixture.mMesh->add_sideset( tSideSet );
        }

        void
        refuses( const Cell< id_t > & aSideSets, const Cell< std::string > & aSubstrings )
        {
            Vector< id_t > tSideSets( aSideSets.size() );
            for ( index_t k = 0; k < aSideSets.size(); ++k )
            {
                tSideSets( k ) = aSideSets( k );
            }

            CurveFactory tFactory( mFixture.mMesh );
            expect_refusal( [ & ]() { delete tFactory.from_2d_sidesets( tSideSets, 7 ) ; },
                aSubstrings );
        }
    };
}

// a LINE2 sideset before a LINE3 one would make a linear curve and drop the
// midpoint; the reverse order is refused as well
TEST( CurveFactory, TwoDMixedOrderSidesetsAreRefused )
{
    QuadEight tQuad ;
    tQuad.add_sideset( 1, ElementType::LINE2, { tQuad.mP( 0 ), tQuad.mP( 1 ) } );

    tQuad.refuses( { 1, 2 }, { "Curve 7", "same element type" } );
    tQuad.refuses( { 2, 1 }, { "Curve 7", "same element type" } );
}

// the facets carry no nodes: the refusal comes before any node is read
TEST( CurveFactory, TwoDNonLineSidesetsAreRefused )
{
    QuadEight tQuad ;
    tQuad.add_sideset( 3, ElementType::LINE4, {} );
    tQuad.add_sideset( 4, ElementType::QUAD4, {} );

    tQuad.refuses( { 3 }, { "Curve 7", "line4", "LINE2 or LINE3" } );
    tQuad.refuses( { 4 }, { "Curve 7", "quad4", "LINE2 or LINE3" } );
}

TEST( CurveFactory, TwoDEmptyCurvesAreRefused )
{
    QuadEight tQuad ;
    tQuad.mFixture.mMesh->add_sideset( new SideSet( 5, 0 ) );

    tQuad.refuses( { 5 }, { "Curve 7", "no facets" } );
    tQuad.refuses( {}, { "Curve 7", "no sideset" } );

    // next to a nonempty sideset, in either order
    tQuad.refuses( { 5, 2 }, { "Curve 7", "sideset 5", "no facets" } );
    tQuad.refuses( { 2, 5 }, { "Curve 7", "sideset 5", "no facets" } );
}

//------------------------------------------------------------------------------

// a mesh without node-to-element links: the midpoints must come from the
// facets. The right face is a domain boundary, so its edge 5-6 is not part
// of the side curve, which runs 5-4-7-6
TEST( CurveFactory, ThinShellQuadraticSideCurveTakesItsFacetsMidpoints )
{
    Fixture tFixture( 3, false );
    Cell< Node * > tP = unit_cube( tFixture );
    Element * tHex = tFixture.hex( tP, true );

    tFixture.sideset( 1, ElementType::QUAD8, { { tHex, 5 } } );
    tFixture.sideset( 2, ElementType::QUAD8, { { tHex, 1 } } );
    tFixture.finalize();

    ASSERT_EQ( tP( 5 )->number_of_elements(), ( uint ) 0 )
        << "the fixture must not link nodes to elements" ;

    CurveFactory tFactory( tFixture.mMesh );
    Vector< id_t > tShells( 1, 1 );
    Vector< id_t > tBoundaries( 1, 2 );
    Cell< Curve * > tCurves = tFactory.thin_shell_side_curves( tShells, tBoundaries );

    for ( Curve * tCurve : tCurves )
    {
        tFixture.mMesh->curves().push( tCurve );
    }

    ASSERT_EQ( tCurves.size(), ( index_t ) 1 );
    ASSERT_EQ( tCurves( 0 )->element_type(), ElementType::LINE3 );

    check_curve( tCurves( 0 ),
        { tP( 5 ), tHex->node( 16 ), tP( 4 ), tHex->node( 19 ), tP( 7 ),
          tHex->node( 18 ), tP( 6 ) }, false );
}

//------------------------------------------------------------------------------

// no domain boundary: the side curve is the whole closed rim, which is
// refused, and the mesh nodes get their indices and vertex containers back
TEST( CurveFactory, ThinShellClosedSideCurveIsRefused )
{
    Fixture tFixture( 3, true );
    Cell< Node * > tP = unit_cube( tFixture );
    Element * tHex = tFixture.hex( tP, false );

    tFixture.sideset( 1, ElementType::QUAD4, { { tHex, 5 } } );
    tFixture.finalize();

    CurveFactory tFactory( tFixture.mMesh );
    Vector< id_t > tShells( 1, 1 );
    Vector< id_t > tBoundaries ;

    expect_refusal( [ & ]()
        {
            Cell< Curve * > tCurves = tFactory.thin_shell_side_curves( tShells, tBoundaries );
            for ( Curve * tCurve : tCurves )
            {
                delete tCurve ;
            }
        },
        { "is closed", "not supported" } );

    index_t tCount = 0 ;
    for ( Node * tNode : tFixture.mMesh->nodes() )
    {
        EXPECT_EQ( tNode->index(), tCount++ );
    }

    for ( index_t k = 4; k < 8; ++k )
    {
        EXPECT_EQ( tP( k )->number_of_vertices(), ( uint ) 0 ) << "node " << tP( k )->id() ;
    }
}
