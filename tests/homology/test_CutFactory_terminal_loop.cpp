/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * Unit tests for CutFactory::close_terminal_loop() on hand-built open
 * curves of two segments of unequal length ( 1 and 3 ). The return pass
 * must reverse every segment with the corners swapped and a LINE3 midside
 * kept at node 2, take the thin-shell duplicate of every indexed node and
 * keep a node at gNoIndex, append the arclength of the return pass in
 * reversed segment order, rebuild the node list in loop order without the
 * closing node, and set the closed flag. A Curve::reverse() after the
 * close must leave every midside at node 2. A node without a duplicate
 * and a mis-sized arclength are refused before anything is allocated:
 * the homology main forces the throwing reaction, so the BELFEM_ERROR is
 * observable as a std::runtime_error in both builds, and the curve and
 * the id counter are untouched afterwards. The curve owns its segments
 * and a segment its element; the fixture owns the nodes.
 */

#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Curve.hpp"
#include "cl_Segment.hpp"
#include "cl_CutFactory.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    // an open two-segment curve along x. For LINE3 the nodes are
    // 1, m12, 2, m23, 3 at x = 0, 0.5, 1, 2.5, 4 ( ids 1..5, indices
    // 0..4 ); for LINE2 they are 1, 2, 3 at x = 0, 1, 4 ( ids 1..3,
    // indices 0..2 ). The duplicates carry ids 11..15 at the same
    // coordinates. The arclength holds one value per node, as
    // CurveFactory::compute_coordinates writes it
    struct Loop
    {
        const ElementType mType ;
        Cell< Node * > mNodes ;
        Cell< Node * > mDuplicates ;
        std::unique_ptr< Curve > mCurve ;
        id_t mMaxID = 100 ;

        explicit Loop( const ElementType aType, const bool aIndexed = true ) :
            mType( aType )
        {
            const bool tQuadratic = aType == ElementType::LINE3 ;

            Cell< real > tX = tQuadratic ?
                Cell< real >( { 0.0, 0.5, 1.0, 2.5, 4.0 } ) :
                Cell< real >( { 0.0, 1.0, 4.0 } );

            for ( uint k = 0; k < tX.size(); ++k )
            {
                Node * tNode = new Node( k + 1, tX( k ), 0.0, 0.0 );
                tNode->set_index( aIndexed ? k : gNoIndex );
                mNodes.push( tNode );

                if ( aIndexed )
                {
                    Node * tDup = new Node( k + 11, tX( k ), 0.0, 0.0 );
                    tDup->set_index( gNoIndex );
                    mDuplicates.push( tDup );
                }
            }

            mCurve.reset( new Curve( 7, aType ) );
            mCurve->set_closed_flag( false );

            ElementFactory tFactory ;
            const uint tStride = tQuadratic ? 2 : 1 ;

            for ( uint s = 0; s < 2; ++s )
            {
                Element * tElement = tFactory.create_element( aType, s + 1 );
                tElement->insert_node( mNodes( tStride * s ), 0 );
                tElement->insert_node( mNodes( tStride * s + tStride ), 1 );
                if ( tQuadratic )
                {
                    tElement->insert_node( mNodes( tStride * s + 1 ), 2 );
                }
                Segment * tSegment = new Segment( tElement );
                tSegment->set_index( s );
                mCurve->segments().push( tSegment );
            }

            for ( Node * tNode : mNodes )
            {
                mCurve->nodes().push( tNode );
            }

            Vector< real > & tS = mCurve->arclength();
            tS.set_size( tX.size() );
            for ( uint k = 0; k < tX.size(); ++k )
            {
                tS( k ) = tX( k );
            }
        }

        ~Loop()
        {
            mCurve.reset();
            for ( Node * tNode : mNodes )
            {
                delete tNode ;
            }
            for ( Node * tNode : mDuplicates )
            {
                delete tNode ;
            }
        }

        Loop( const Loop & ) = delete ;
        Loop & operator=( const Loop & ) = delete ;

        void
        close()
        {
            CutFactory::close_terminal_loop( mCurve.get(), mDuplicates, mMaxID );
        }

        // the node the return pass must carry for original k
        Node *
        back( const uint k ) const
        {
            return mDuplicates.size() > 0 ? mDuplicates( k ) : mNodes( k );
        }

        // true if aNode is one of the midside nodes, original or duplicate
        bool
        is_midside( const Node * aNode ) const
        {
            return aNode == mNodes( 1 ) || aNode == mNodes( 3 )
                || aNode == this->back( 1 ) || aNode == this->back( 3 );
        }
    };

    void
    expect_arclength( const Vector< real > & aS, const Cell< real > & aExpected )
    {
        ASSERT_EQ( aS.length(), aExpected.size() );
        for ( index_t k = 0; k < aExpected.size(); ++k )
        {
            EXPECT_DOUBLE_EQ( aS( k ), aExpected( k ) ) << "entry " << k ;
        }
    }
}

//------------------------------------------------------------------------------

TEST( CutFactoryTerminalLoop, Line3KeepsMidsideAtNodeTwo )
{
    Loop tLoop( ElementType::LINE3 );
    tLoop.close();

    const Cell< Segment * > & tSegments = tLoop.mCurve->segments();
    ASSERT_EQ( tSegments.size(), ( index_t ) 4 );

    // the originals are untouched
    EXPECT_EQ( tSegments( 0 )->node( 0 ), tLoop.mNodes( 0 ) );
    EXPECT_EQ( tSegments( 0 )->node( 1 ), tLoop.mNodes( 2 ) );
    EXPECT_EQ( tSegments( 0 )->node( 2 ), tLoop.mNodes( 1 ) );

    // the return pass runs 3 -> 2 -> 1 on the duplicates, midsides at node 2
    EXPECT_EQ( tSegments( 2 )->node( 0 ), tLoop.back( 4 ) );
    EXPECT_EQ( tSegments( 2 )->node( 1 ), tLoop.back( 2 ) );
    EXPECT_EQ( tSegments( 2 )->node( 2 ), tLoop.back( 3 ) );
    EXPECT_EQ( tSegments( 3 )->node( 0 ), tLoop.back( 2 ) );
    EXPECT_EQ( tSegments( 3 )->node( 1 ), tLoop.back( 0 ) );
    EXPECT_EQ( tSegments( 3 )->node( 2 ), tLoop.back( 1 ) );

    // node list in loop order, corner then midside, the closing node omitted
    const Cell< Node * > & tNodes = tLoop.mCurve->nodes();
    ASSERT_EQ( tNodes.size(), ( index_t ) 8 );
    EXPECT_EQ( tNodes( 0 ), tLoop.mNodes( 0 ) );
    EXPECT_EQ( tNodes( 1 ), tLoop.mNodes( 1 ) );
    EXPECT_EQ( tNodes( 2 ), tLoop.mNodes( 2 ) );
    EXPECT_EQ( tNodes( 3 ), tLoop.mNodes( 3 ) );
    EXPECT_EQ( tNodes( 4 ), tLoop.back( 4 ) );
    EXPECT_EQ( tNodes( 5 ), tLoop.back( 3 ) );
    EXPECT_EQ( tNodes( 6 ), tLoop.back( 2 ) );
    EXPECT_EQ( tNodes( 7 ), tLoop.back( 1 ) );

    // return pass: the 3-long segment first, then the 1-long one
    expect_arclength( tLoop.mCurve->arclength(),
        { 0.0, 0.5, 1.0, 2.5, 4.0, 5.5, 7.0, 7.5, 8.0 } );

    EXPECT_TRUE( tLoop.mCurve->is_closed() );
    EXPECT_EQ( tLoop.mMaxID, ( id_t ) 102 );
}

//------------------------------------------------------------------------------

TEST( CutFactoryTerminalLoop, Line3WithoutIndexReusesOriginals )
{
    Loop tLoop( ElementType::LINE3, false );
    tLoop.close();

    const Cell< Segment * > & tSegments = tLoop.mCurve->segments();
    ASSERT_EQ( tSegments.size(), ( index_t ) 4 );

    EXPECT_EQ( tSegments( 2 )->node( 0 ), tLoop.mNodes( 4 ) );
    EXPECT_EQ( tSegments( 2 )->node( 1 ), tLoop.mNodes( 2 ) );
    EXPECT_EQ( tSegments( 2 )->node( 2 ), tLoop.mNodes( 3 ) );
    EXPECT_EQ( tSegments( 3 )->node( 0 ), tLoop.mNodes( 2 ) );
    EXPECT_EQ( tSegments( 3 )->node( 1 ), tLoop.mNodes( 0 ) );
    EXPECT_EQ( tSegments( 3 )->node( 2 ), tLoop.mNodes( 1 ) );

    EXPECT_TRUE( tLoop.mCurve->is_closed() );
}

//------------------------------------------------------------------------------

TEST( CutFactoryTerminalLoop, Line3SurvivesReverse )
{
    Loop tLoop( ElementType::LINE3 );
    tLoop.close();
    tLoop.mCurve->reverse();

    const Cell< Segment * > & tSegments = tLoop.mCurve->segments();
    ASSERT_EQ( tSegments.size(), ( index_t ) 4 );

    for ( const Segment * tSegment : tSegments )
    {
        EXPECT_FALSE( tLoop.is_midside( tSegment->node( 0 ) ) );
        EXPECT_FALSE( tLoop.is_midside( tSegment->node( 1 ) ) );
        EXPECT_TRUE( tLoop.is_midside( tSegment->node( 2 ) ) );
    }
}

//------------------------------------------------------------------------------

TEST( CutFactoryTerminalLoop, Line2ReversesCorners )
{
    Loop tLoop( ElementType::LINE2 );
    tLoop.close();

    const Cell< Segment * > & tSegments = tLoop.mCurve->segments();
    ASSERT_EQ( tSegments.size(), ( index_t ) 4 );

    EXPECT_EQ( tSegments( 2 )->node( 0 ), tLoop.back( 2 ) );
    EXPECT_EQ( tSegments( 2 )->node( 1 ), tLoop.back( 1 ) );
    EXPECT_EQ( tSegments( 3 )->node( 0 ), tLoop.back( 1 ) );
    EXPECT_EQ( tSegments( 3 )->node( 1 ), tLoop.back( 0 ) );

    const Cell< Node * > & tNodes = tLoop.mCurve->nodes();
    ASSERT_EQ( tNodes.size(), ( index_t ) 4 );
    EXPECT_EQ( tNodes( 0 ), tLoop.mNodes( 0 ) );
    EXPECT_EQ( tNodes( 1 ), tLoop.mNodes( 1 ) );
    EXPECT_EQ( tNodes( 2 ), tLoop.back( 2 ) );
    EXPECT_EQ( tNodes( 3 ), tLoop.back( 1 ) );

    expect_arclength( tLoop.mCurve->arclength(), { 0.0, 1.0, 4.0, 7.0, 8.0 } );

    EXPECT_TRUE( tLoop.mCurve->is_closed() );
    EXPECT_EQ( tLoop.mMaxID, ( id_t ) 102 );
}

//------------------------------------------------------------------------------

TEST( CutFactoryTerminalLoop, EmptyCurveIsLeftAlone )
{
    Loop tLoop( ElementType::LINE2 );
    Curve tEmpty( 8, ElementType::LINE2 );
    tEmpty.set_closed_flag( false );

    CutFactory::close_terminal_loop( &tEmpty, tLoop.mDuplicates, tLoop.mMaxID );

    EXPECT_EQ( tEmpty.segments().size(), ( index_t ) 0 );
    EXPECT_FALSE( tEmpty.is_closed() );
    EXPECT_EQ( tLoop.mMaxID, ( id_t ) 100 );
}

//------------------------------------------------------------------------------

namespace
{
    // the helper must refuse before it allocates: two segments, the id
    // counter and the open flag survive
    void
    expect_untouched( const Loop & aLoop )
    {
        EXPECT_EQ( aLoop.mCurve->segments().size(), ( index_t ) 2 );
        EXPECT_EQ( aLoop.mMaxID, ( id_t ) 100 );
        EXPECT_FALSE( aLoop.mCurve->is_closed() );
    }
}

TEST( CutFactoryTerminalLoop, RefusesIndexWithoutDuplicate )
{
    Loop tLoop( ElementType::LINE2 );
    tLoop.mNodes( 2 )->set_index( 7 );

    EXPECT_THROW( tLoop.close(), std::runtime_error );
    expect_untouched( tLoop );
}

//------------------------------------------------------------------------------

TEST( CutFactoryTerminalLoop, RefusesNullDuplicate )
{
    Loop tLoop( ElementType::LINE3 );
    Node * tSaved = tLoop.mDuplicates( 3 );
    tLoop.mDuplicates( 3 ) = nullptr ;

    EXPECT_THROW( tLoop.close(), std::runtime_error );
    expect_untouched( tLoop );

    tLoop.mDuplicates( 3 ) = tSaved ;
}

//------------------------------------------------------------------------------

TEST( CutFactoryTerminalLoop, RefusesMisSizedArclength )
{
    Loop tLoop( ElementType::LINE3 );
    tLoop.mCurve->arclength().set_size( 2 );

    EXPECT_THROW( tLoop.close(), std::runtime_error );
    expect_untouched( tLoop );
}
