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
 * A deck curve that lies on a thin-shell sideset becomes a terminal curve of
 * that shell, and a terminal curve must be open. MaxwellFactory refuses a
 * closed one while it reads the thin-shell data, on the 2-D branch and on
 * both 3-D branches ( the tape is sideset A of the curve, or sideset B ).
 *
 * The association is private and runs inside create_magnetic_kernel(), so
 * each case builds a factory from a deck and a hand-written Gmsh 2.2 mesh in
 * the test temp directory. The decks carry no homology section on purpose:
 * right after the thin-shell data and the topology scan, create_cuts()
 * refuses a deck without one. That refusal is the stop sign of an accepted
 * curve: it shows that reading the thin-shell data returned without a
 * refusal, and it keeps the case from building cuts and a kernel. The
 * terminal list itself is private, so an accepted case cannot see the
 * curve being stored; the 3-D cases see the sideset swap that precedes it.
 * Each case works in its own fresh directory, removed with its contents.
 */

#include <gtest/gtest.h>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_Curve.hpp"
#include "cl_SideSet.hpp"
#include "cl_MaxwellFactory.hpp"

using namespace belfem;

namespace
{
    // a fresh directory under the test temp directory, unique to this
    // process and case, removed with everything in it on scope exit
    class ScratchDir
    {
        std::string mPath ;

    public:

        explicit ScratchDir( const std::string & aName )
        {
            std::string tTemplate = ::testing::TempDir() + aName + "_XXXXXX" ;
            if ( mkdtemp( tTemplate.data() ) == nullptr )
            {
                throw std::runtime_error( "could not create a directory from " + tTemplate );
            }
            mPath = tTemplate + "/" ;
        }

        ScratchDir( const ScratchDir & ) = delete ;
        ScratchDir( ScratchDir && ) = delete ;
        ScratchDir & operator=( const ScratchDir & ) = delete ;
        ScratchDir & operator=( ScratchDir && ) = delete ;

        ~ScratchDir()
        {
            std::error_code tIgnored ;
            std::filesystem::remove_all( mPath, tIgnored );
        }

        // writes aBody to aName in this directory and returns its path
        std::string
        write( const std::string & aName, const std::string & aBody ) const
        {
            const std::string tPath = mPath + aName ;
            std::ofstream tFile( tPath );
            tFile << aBody ;
            if ( ! tFile )
            {
                throw std::runtime_error( "could not write " + tPath );
            }
            return tPath ;
        }
    };

    const char * gClosedRefusal = "is closed, but a terminal curve must be open" ;

    const char * gHomologyStop = "No homology section found in input file" ;

    const id_t gTape  = 5 ;
    const id_t gOther = 6 ;

    // one .msh 2.2 element line; both tags are the group id, the reader
    // keys blocks and sidesets by the second ( geometry ) tag
    std::string
    element_line( id_t & aID, const uint aGmshType, const id_t aTag,
                  const Cell< id_t > & aNodes )
    {
        std::string tLine = std::to_string( aID++ ) + " "
                          + std::to_string( aGmshType ) + " 2 "
                          + std::to_string( aTag ) + " " + std::to_string( aTag );
        for ( id_t tNode : aNodes )
        {
            tLine += " " + std::to_string( tNode );
        }
        return tLine + "\n" ;
    }

    std::string
    msh_text( const std::string & aNodes, const index_t aNumNodes,
              const std::string & aElements, const index_t aNumElements )
    {
        return "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n"
             + std::to_string( aNumNodes ) + "\n" + aNodes
             + "$EndNodes\n$Elements\n"
             + std::to_string( aNumElements ) + "\n" + aElements
             + "$EndElements\n" ;
    }

    // 2-D: the square [0,3]^2 as a 3x3 grid of TRI3 ( block 1 ), with the
    // tape as LINE2 elements of sideset 5 along grid edges. Node ( i, j )
    // has the id 1 + i + 4 j
    std::string
    mesh_2d( const Cell< Cell< id_t > > & aTape )
    {
        std::string tNodes ;
        for ( index_t j = 0 ; j < 4 ; ++j )
        {
            for ( index_t i = 0 ; i < 4 ; ++i )
            {
                tNodes += std::to_string( 1 + i + 4 * j ) + " "
                        + std::to_string( i ) + " " + std::to_string( j ) + " 0\n" ;
            }
        }

        id_t tID = 1 ;
        std::string tElements ;
        for ( id_t j = 0 ; j < 3 ; ++j )
        {
            for ( id_t i = 0 ; i < 3 ; ++i )
            {
                const id_t n00 = 1 + i + 4 * j ;
                const id_t n10 = n00 + 1 ;
                const id_t n01 = n00 + 4 ;
                const id_t n11 = n01 + 1 ;

                // counterclockwise
                tElements += element_line( tID, 2, 1, { n00, n10, n11 } );
                tElements += element_line( tID, 2, 1, { n00, n11, n01 } );
            }
        }
        for ( const Cell< id_t > & tLine : aTape )
        {
            tElements += element_line( tID, 1, gTape, tLine );
        }

        return msh_text( tNodes, 16, tElements, tID - 1 );
    }

    // corner v = i + 2 j + 4 k of the unit cube, coordinate aAxis
    uint
    corner_coordinate( const id_t aCorner, const uint aAxis )
    {
        return ( aCorner >> aAxis ) & 1u ;
    }

    // 3-D: the unit cube as the six Kuhn tetrahedra around the diagonal
    // from corner 0 to corner 7 ( block 1 ). The tape is the face z = 0
    // ( sideset 5 ); the faces listed in aOtherFaces as ( axis, value ) form
    // sideset 6. Corner v has the node id v + 1
    std::string
    mesh_3d( const Cell< std::pair< uint, uint > > & aOtherFaces )
    {
        std::string tNodes ;
        for ( id_t v = 0 ; v < 8 ; ++v )
        {
            tNodes += std::to_string( v + 1 ) + " "
                    + std::to_string( corner_coordinate( v, 0 ) ) + " "
                    + std::to_string( corner_coordinate( v, 1 ) ) + " "
                    + std::to_string( corner_coordinate( v, 2 ) ) + "\n" ;
        }

        const uint tPermutations[ 6 ][ 3 ] = {
            { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 },
            { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };

        Cell< Cell< id_t > > tTets ;
        for ( const auto & tAxes : tPermutations )
        {
            Cell< id_t > tTet = { 0, 1u << tAxes[ 0 ],
                                  ( 1u << tAxes[ 0 ] ) | ( 1u << tAxes[ 1 ] ), 7 };

            // positive volume: det( v1 - v0, v2 - v0, v3 - v0 ) > 0
            int e[ 3 ][ 3 ];
            for ( uint r = 0 ; r < 3 ; ++r )
            {
                for ( uint c = 0 ; c < 3 ; ++c )
                {
                    e[ r ][ c ] = ( int ) corner_coordinate( tTet( r + 1 ), c ) ;
                }
            }
            const int tDet = e[ 0 ][ 0 ] * ( e[ 1 ][ 1 ] * e[ 2 ][ 2 ] - e[ 1 ][ 2 ] * e[ 2 ][ 1 ] )
                           - e[ 0 ][ 1 ] * ( e[ 1 ][ 0 ] * e[ 2 ][ 2 ] - e[ 1 ][ 2 ] * e[ 2 ][ 0 ] )
                           + e[ 0 ][ 2 ] * ( e[ 1 ][ 0 ] * e[ 2 ][ 1 ] - e[ 1 ][ 1 ] * e[ 2 ][ 0 ] );
            if ( tDet < 0 )
            {
                std::swap( tTet( 2 ), tTet( 3 ) );
            }
            tTets.push( tTet );
        }

        id_t tID = 1 ;
        std::string tElements ;
        for ( const Cell< id_t > & tTet : tTets )
        {
            tElements += element_line( tID, 4, 1,
                { tTet( 0 ) + 1, tTet( 1 ) + 1, tTet( 2 ) + 1, tTet( 3 ) + 1 } );
        }

        // every tetrahedron face that lies in the cube face ( aAxis, aValue )
        auto tFaceTriangles = [ & ]( const uint aAxis, const uint aValue, const id_t aTag )
        {
            for ( const Cell< id_t > & tTet : tTets )
            {
                for ( uint tOmit = 0 ; tOmit < 4 ; ++tOmit )
                {
                    Cell< id_t > tTriangle ;
                    for ( uint k = 0 ; k < 4 ; ++k )
                    {
                        if ( k != tOmit && corner_coordinate( tTet( k ), aAxis ) == aValue )
                        {
                            tTriangle.push( tTet( k ) + 1 );
                        }
                    }
                    if ( tTriangle.size() == 3 )
                    {
                        tElements += element_line( tID, 2, aTag, tTriangle );
                    }
                }
            }
        };

        tFaceTriangles( 2, 0, gTape );
        for ( const std::pair< uint, uint > & tFace : aOtherFaces )
        {
            tFaceTriangles( tFace.first, tFace.second, gOther );
        }

        return msh_text( tNodes, 8, tElements, tID - 1 );
    }

    // a deck with one thin shell on sideset 5, one air block and the single
    // deck curve 1 : aCurve. Nothing downstream of the association is set up
    std::string
    deck( const std::string & aMeshPath, const std::string & aCurve )
    {
        return "mesh\n{\n    file : " + aMeshPath + " ;\n    unit : m ;\n}\n\n"
               "layers : tape\n{\n    copper : 20 mum ;\n}\n\n"
               "topology\n{\n"
               "    thinshell : tape\n    {\n        sidesets : 5 ;\n    }\n"
               "    air\n    {\n        blocks : 1 ;\n    }\n"
               "    curves\n    {\n        1 : " + aCurve + " ;\n    }\n"
               "}\n\n"
               "boundary conditions\n{\n}\n" ;
    }

    // what create_magnetic_kernel() threw on this mesh and curve; the
    // sideset A and B ids of deck curve 1 afterwards land in aSideSetA and
    // aSideSetB. Only std::runtime_error is caught, so any other exception
    // fails the test
    std::string
    refusal( const std::string & aName,
             const std::string & aMesh,
             const std::string & aCurve,
             id_t              & aSideSetA,
             id_t              & aSideSetB )
    {
        // a fresh directory also holds no processed mesh ( .bfm ) the
        // factory could reuse instead of reading the .msh
        ScratchDir tDir( aName );
        const std::string tMesh = tDir.write( aName + ".msh", aMesh );
        const std::string tDeck = tDir.write( aName + ".conf", deck( tMesh, aCurve ) );

        fem::MaxwellFactory tFactory( tDeck );

        std::string tMessage ;
        try
        {
            tFactory.create_magnetic_kernel();
        }
        catch ( const std::runtime_error & aError )
        {
            tMessage = aError.what();
        }

        aSideSetA = gNoID ;
        aSideSetB = gNoID ;
        if ( tFactory.mesh()->curves().size() != 1 )
        {
            ADD_FAILURE() << "expected one mesh curve, found "
                          << tFactory.mesh()->curves().size() ;
            return tMessage ;
        }

        mesh::Curve * tCurve = tFactory.mesh()->curves()( 0 );
        aSideSetA = tCurve->sideset_a() == nullptr ? gNoID : tCurve->sideset_a()->id() ;
        aSideSetB = tCurve->sideset_b() == nullptr ? gNoID : tCurve->sideset_b()->id() ;

        return tMessage ;
    }

    bool
    contains( const std::string & aText, const std::string & aPart )
    {
        return aText.find( aPart ) != std::string::npos ;
    }

    // the loop around the centre cell [1,2]^2, and two of its edges
    const Cell< Cell< id_t > > gClosedTape = { { 6, 7 }, { 7, 11 }, { 11, 10 }, { 10, 6 } };
    const Cell< Cell< id_t > > gOpenTape   = { { 6, 7 }, { 7, 11 } };

    // x = 0 alone meets the face z = 0 in one edge; the four side faces
    // together meet it in its closed boundary
    const Cell< std::pair< uint, uint > > gOneSide  = { { 0, 0 } };
    const Cell< std::pair< uint, uint > > gAllSides = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 } };

    const std::string gTapeFirst = "5 @ 6" ;
    const std::string gTapeLast  = "6 @ 5" ;
}

//------------------------------------------------------------------------------
// 2-D
//------------------------------------------------------------------------------

TEST( ThinShellTerminalCurves, ClosedCurveIsRefused2D )
{
    id_t tA, tB ;
    const std::string tMessage = refusal( "tstc_2d_closed", mesh_2d( gClosedTape ), "5", tA, tB );

    EXPECT_TRUE( contains( tMessage, gClosedRefusal ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "curve 1 on thin shell tape" ) ) << tMessage ;
}

//------------------------------------------------------------------------------

TEST( ThinShellTerminalCurves, OpenCurveIsAccepted2D )
{
    id_t tA, tB ;
    const std::string tMessage = refusal( "tstc_2d_open", mesh_2d( gOpenTape ), "5", tA, tB );

    EXPECT_TRUE( contains( tMessage, gHomologyStop ) ) << tMessage ;
    EXPECT_EQ( tA, gTape );
}

//------------------------------------------------------------------------------
// 3-D, the tape is sideset A of the curve
//------------------------------------------------------------------------------

TEST( ThinShellTerminalCurves, ClosedCurveIsRefused3DTapeFirst )
{
    id_t tA, tB ;
    const std::string tMessage = refusal( "tstc_3d_closed_a", mesh_3d( gAllSides ), gTapeFirst, tA, tB );

    EXPECT_TRUE( contains( tMessage, gClosedRefusal ) ) << tMessage ;

    // refused before the swap that makes the shell sideset B
    EXPECT_EQ( tA, gTape );
    EXPECT_EQ( tB, gOther );
}

//------------------------------------------------------------------------------

TEST( ThinShellTerminalCurves, OpenCurveIsAccepted3DTapeFirst )
{
    id_t tA, tB ;
    const std::string tMessage = refusal( "tstc_3d_open_a", mesh_3d( gOneSide ), gTapeFirst, tA, tB );

    EXPECT_TRUE( contains( tMessage, gHomologyStop ) ) << tMessage ;

    // the association swapped the sidesets: the shell is now sideset B
    EXPECT_EQ( tA, gOther );
    EXPECT_EQ( tB, gTape );
}

//------------------------------------------------------------------------------
// 3-D, the tape is sideset B of the curve
//------------------------------------------------------------------------------

TEST( ThinShellTerminalCurves, ClosedCurveIsRefused3DTapeLast )
{
    id_t tA, tB ;
    const std::string tMessage = refusal( "tstc_3d_closed_b", mesh_3d( gAllSides ), gTapeLast, tA, tB );

    EXPECT_TRUE( contains( tMessage, gClosedRefusal ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "curve 1 on thin shell tape" ) ) << tMessage ;

    // this branch never swaps: the curve keeps the deck's order
    EXPECT_EQ( tA, gOther );
    EXPECT_EQ( tB, gTape );
}

//------------------------------------------------------------------------------

TEST( ThinShellTerminalCurves, OpenCurveIsAccepted3DTapeLast )
{
    id_t tA, tB ;
    const std::string tMessage = refusal( "tstc_3d_open_b", mesh_3d( gOneSide ), gTapeLast, tA, tB );

    EXPECT_TRUE( contains( tMessage, gHomologyStop ) ) << tMessage ;
    EXPECT_EQ( tA, gOther );
    EXPECT_EQ( tB, gTape );
}
