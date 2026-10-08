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
 * Unit tests for CutFactory::refuse_closed_side_curve_sharing_nodes() on two
 * hand-built Protoshells. The mesh holds two disjoint HEX8 cubes; the side
 * curve of shell A is the closed rim of the top face of the first cube, built
 * by CurveFactory::thin_shell_side_curves() without domain boundaries. Shell
 * B gets, per test, the closed rim of the second cube, a hand-built side
 * curve, or nothing; terminal curves are hand-built node lists, which is all
 * the helper reads. As in the cut factory, the side curves are pushed onto
 * the mesh before the check, behind a curve the mesh already owns. A closed
 * rim that shares a node with a terminal curve of either shell, or with a
 * side curve of the other shell, is refused: the mesh is back to its one
 * earlier curve, the side curves of both shells are gone and the node flags
 * are clear. Isolated closed rims pass, stay intact, sit in the mesh exactly
 * once, and the node flags are clear.
 *
 * The homology main forces the throwing reaction, so the BELFEM_ERROR of the
 * refusal is observable as a std::runtime_error in both builds.
 */

#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "cl_Curve.hpp"
#include "cl_CurveFactory.hpp"
#include "cl_Protoshell.hpp"
#include "cl_CutFactory.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    // two unit cubes as HEX8s in block 1, the first at the origin and the
    // second shifted by 3 in x, with their top faces as sidesets 1 and 2.
    // Corners 0-7 belong to the first cube, 8-15 to the second, each bottom
    // 0-3 and top 4-7 in HEX8 order. Shell A ( "rimtape" ) lies on sideset
    // 1 and has its closed rim 5-8-7-6-5 as side curve; shell B
    // ( "othertape" ) lies on sideset 2 and starts without curves. The mesh
    // owns nodes, elements, the open curve 1-4 ( id 99 ) and every side
    // curve pushed into it; the fixture owns the terminal curves
    struct Rims
    {
        std::unique_ptr< Mesh > mMesh ;
        Cell< Node * > mP ;
        Protoshell mA ;
        Protoshell mB ;
        Cell< Protoshell * > mShells ;
        Cell< Curve * > mTerminals ;
        std::unique_ptr< CurveFactory > mCurveFactory ;

        Rims() :
            mMesh( new Mesh( 3, 0, true ) ),
            mA( 1 ),
            mB( 2 )
        {
            for ( index_t k = 0; k < 16; ++k )
            {
                const index_t c = k % 8 ;
                const real tX = ( ( c == 1 || c == 2 || c == 5 || c == 6 ) ? 1.0 : 0.0 ) + ( k < 8 ? 0.0 : 3.0 ) ;
                const real tY = ( c == 2 || c == 3 || c == 6 || c == 7 ) ? 1.0 : 0.0 ;
                const real tZ = c < 4 ? 0.0 : 1.0 ;
                Node * tNode = new Node( k + 1, tX, tY, tZ );
                mMesh->nodes().push( tNode );
                mP.push( tNode );
            }

            ElementFactory tFactory ;
            Block * tBlock = new Block( 1, 2 );

            for ( index_t h = 0; h < 2; ++h )
            {
                Element * tHex = tFactory.create_element( ElementType::HEX8, h + 1 );
                for ( uint k = 0; k < 8; ++k )
                {
                    tHex->insert_node( mP( 8 * h + k ), k );
                }
                tHex->set_block_id( 1 );
                mMesh->elements().push( tHex );
                tBlock->insert_element( tHex );

                // HEX8 facet 5 is the top face
                SideSet * tSideSet = new SideSet( h + 1, 1 );
                Facet * tFacet = new Facet( tFactory.create_element( ElementType::QUAD4, 1001 + h ) );
                tFacet->set_master( tHex, 5 );
                tFacet->set_sideset_id( h + 1 );
                tSideSet->insert_facet( tFacet );
                mMesh->sidesets().push( tSideSet );
            }

            mMesh->blocks().push( tBlock );
            mMesh->finalize();

            mA.label() = "rimtape" ;
            mA.sidesets().set_size( 1, 1 );
            mB.label() = "othertape" ;
            mB.sidesets().set_size( 1, 2 );
            mShells = { &mA, &mB };

            mMesh->curves().push( this->curve( 99, { 0, 3 } ) );

            mCurveFactory.reset( new CurveFactory( mMesh.get() ) );
            this->build_rim( mA );
        }

        ~Rims()
        {
            for ( Curve * tTerminal : mTerminals )
            {
                delete tTerminal ;
            }
        }

        Rims( const Rims & ) = delete ;
        Rims( Rims && ) = delete ;
        Rims & operator=( const Rims & ) = delete ;
        Rims & operator=( Rims && ) = delete ;

        // the side curves of aShell from its own sideset, no boundaries
        void
        build_rim( Protoshell & aShell )
        {
            Vector< id_t > tBoundaries ;
            aShell.side_curves() = mCurveFactory->thin_shell_side_curves( aShell.sidesets(), tBoundaries );
        }

        // a hand-built open curve through the corners aCorners
        Curve *
        curve( const id_t aID, const Cell< index_t > & aCorners )
        {
            Curve * tCurve = new Curve( aID, ElementType::LINE2 );
            tCurve->set_closed_flag( false );
            for ( index_t k : aCorners )
            {
                tCurve->nodes().push( mP( k ) );
            }
            return tCurve ;
        }

        void
        add_terminal( Protoshell & aShell, const id_t aID, const Cell< index_t > & aCorners )
        {
            Curve * tCurve = this->curve( aID, aCorners );
            mTerminals.push( tCurve );
            aShell.terminal_curves().push( tCurve );
        }

        void
        add_side_curve( Protoshell & aShell, const id_t aID, const Cell< index_t > & aCorners )
        {
            aShell.side_curves().push( this->curve( aID, aCorners ) );
        }

        // the message of the refusal, empty when the helper returned. The
        // side curves are pushed onto the mesh first, as in the cut factory
        std::string
        refusal()
        {
            const index_t tNumMeshCurves = mMesh->curves().size() ;

            for ( Protoshell * tShell : mShells )
            {
                for ( Curve * tCurve : tShell->side_curves() )
                {
                    mMesh->curves().push( tCurve );
                }
            }

            std::string tMessage ;
            try
            {
                CutFactory::refuse_closed_side_curve_sharing_nodes( mMesh.get(), mShells, tNumMeshCurves );
            }
            catch ( const std::runtime_error & aError )
            {
                tMessage = aError.what();
            }
            return tMessage ;
        }

        void
        expect_refused_and_clean()
        {
            EXPECT_TRUE( mA.side_curves().empty() );
            EXPECT_TRUE( mB.side_curves().empty() );
            ASSERT_EQ( mMesh->curves().size(), ( index_t ) 1 );
            EXPECT_EQ( mMesh->curves()( 0 )->id(), ( id_t ) 99 );
            this->expect_flags_clear();
        }

        // every surviving side curve is in the mesh exactly once
        void
        expect_owned_once()
        {
            index_t tNumSideCurves = 0 ;
            for ( Protoshell * tShell : mShells )
            {
                for ( Curve * tSideCurve : tShell->side_curves() )
                {
                    ++tNumSideCurves ;
                    index_t tCount = 0 ;
                    for ( Curve * tCurve : mMesh->curves() )
                    {
                        if ( tCurve == tSideCurve ) ++tCount ;
                    }
                    EXPECT_EQ( tCount, ( index_t ) 1 );
                }
            }
            EXPECT_EQ( mMesh->curves().size(), tNumSideCurves + 1 );
        }

        void
        expect_flags_clear()
        {
            for ( Node * tNode : mMesh->nodes() )
            {
                EXPECT_FALSE( tNode->is_flagged() ) << "node " << tNode->id() ;
            }
        }
    };

    bool
    contains( const std::string & aText, const std::string & aPart )
    {
        return aText.find( aPart ) != std::string::npos ;
    }

    void
    expect_closed_rim( const Curve * aCurve )
    {
        EXPECT_TRUE( aCurve->is_closed() );
        EXPECT_EQ( aCurve->nodes().size(), ( index_t ) 5 );
        EXPECT_EQ( aCurve->segments().size(), ( index_t ) 4 );
    }
}

//------------------------------------------------------------------------------

TEST( CutFactoryClosedRim, FixtureRimIsOneClosedCurve )
{
    Rims tRims ;

    ASSERT_EQ( tRims.mA.side_curves().size(), ( index_t ) 1 );
    EXPECT_TRUE( tRims.mA.side_curves()( 0 )->is_closed() );

    tRims.refusal();
}

//------------------------------------------------------------------------------

// the terminal of shell A runs along the rim edge 5-6
TEST( CutFactoryClosedRim, TerminalOnTheRimIsRefused )
{
    Rims tRims ;
    tRims.add_terminal( tRims.mA, 9, { 4, 5 } );
    const std::string tMessage = tRims.refusal();

    EXPECT_TRUE( contains( tMessage, "thin shell rimtape" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "with terminal curve 9 of thin shell rimtape" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "side curve" ) ) << tMessage ;
    tRims.expect_refused_and_clean();
}

//------------------------------------------------------------------------------

// the terminal 1-5 of shell A rises from the bottom and touches the rim at
// corner 5 only
TEST( CutFactoryClosedRim, TerminalTouchingOneCornerIsRefused )
{
    Rims tRims ;
    tRims.add_terminal( tRims.mA, 9, { 0, 4 } );
    const std::string tMessage = tRims.refusal();

    EXPECT_TRUE( contains( tMessage, "thin shell rimtape shares node 5 with terminal curve 9 of thin shell rimtape" ) ) << tMessage ;
    tRims.expect_refused_and_clean();
}

//------------------------------------------------------------------------------

// the terminal 1-2 of shell A lies on the bottom face, away from the rim
TEST( CutFactoryClosedRim, TerminalOffTheRimPasses )
{
    Rims tRims ;
    tRims.add_terminal( tRims.mA, 9, { 0, 1 } );
    const std::string tMessage = tRims.refusal();

    EXPECT_TRUE( tMessage.empty() ) << tMessage ;
    ASSERT_EQ( tRims.mA.side_curves().size(), ( index_t ) 1 );
    expect_closed_rim( tRims.mA.side_curves()( 0 ) );
    tRims.expect_owned_once();
    tRims.expect_flags_clear();
}

//------------------------------------------------------------------------------

// a terminal of the other shell, 6-7, lies on the rim of shell A, which
// walks 5-8-7-6, so node 7 is the first shared one
TEST( CutFactoryClosedRim, OtherShellsTerminalOnTheRimIsRefused )
{
    Rims tRims ;
    tRims.build_rim( tRims.mB );
    tRims.add_terminal( tRims.mA, 9, { 0, 1 } );
    tRims.add_terminal( tRims.mB, 10, { 5, 6 } );
    const std::string tMessage = tRims.refusal();

    EXPECT_TRUE( contains( tMessage, "of thin shell rimtape shares node 7 with terminal curve 10 of thin shell othertape" ) ) << tMessage ;
    tRims.expect_refused_and_clean();
}

//------------------------------------------------------------------------------

// an open side curve of the other shell, 8-9, ends on corner 8 of the rim
// of shell A
TEST( CutFactoryClosedRim, OtherShellsSideCurveOnTheRimIsRefused )
{
    Rims tRims ;
    tRims.add_side_curve( tRims.mB, 19, { 7, 8 } );
    const std::string tMessage = tRims.refusal();

    EXPECT_TRUE( contains( tMessage, "of thin shell rimtape shares node 8 with a side curve of thin shell othertape" ) ) << tMessage ;
    tRims.expect_refused_and_clean();
}

//------------------------------------------------------------------------------

// two closed rims that share no node with each other or with a terminal
// curve both pass and stay in the mesh
TEST( CutFactoryClosedRim, DisjointClosedRimsPass )
{
    Rims tRims ;
    tRims.build_rim( tRims.mB );
    tRims.add_terminal( tRims.mA, 9, { 0, 1 } );
    tRims.add_terminal( tRims.mB, 10, { 8, 9 } );
    const std::string tMessage = tRims.refusal();

    EXPECT_TRUE( tMessage.empty() ) << tMessage ;
    ASSERT_EQ( tRims.mA.side_curves().size(), ( index_t ) 1 );
    ASSERT_EQ( tRims.mB.side_curves().size(), ( index_t ) 1 );
    expect_closed_rim( tRims.mA.side_curves()( 0 ) );
    expect_closed_rim( tRims.mB.side_curves()( 0 ) );
    EXPECT_NE( tRims.mA.side_curves()( 0 ), tRims.mB.side_curves()( 0 ) );
    tRims.expect_owned_once();
    tRims.expect_flags_clear();
}
