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
 * walked.
 *
 * The mesh main forces the throwing reaction, so every BELFEM_ERROR on the
 * refusal paths is observable as a std::runtime_error in both builds. A
 * polyline is never sorted twice after a refusal: the sorter leaves node
 * flags dirty when it throws.
 */

#include <gtest/gtest.h>
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
