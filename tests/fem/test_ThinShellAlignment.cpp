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
 * The orientation of a 3-D thin-shell sideset between two conductor blocks,
 * from MaxwellFactory::fix_facet_masters() through the cut and the thin-shell
 * factory. Each case is a hand-written Gmsh 2.2 mesh of a box in unit cells,
 * two cells deep in z, every cell split into the six Kuhn tetrahedra, with
 * air around a conductor column that runs through the box along z. The shell
 * is sideset 10 in the plane y = 2. Every conductor|air face and every outer
 * face of the box is in a sideset: the cut through the air must end on
 * sidesets, or it is pruned to nothing.
 *
 * Conductor|conductor shell facets have no domain-type order, so the
 * same-type walk orients them from the conductor|air walls that meet the
 * shell at a T-junction, x = const. A wall above the shell and a wall below
 * it give the shell opposite windings, and the walk takes the one it pops
 * first, which follows the sideset order. The wall sidesets are numbered so
 * that the walk reaches the shell from the wall above at one end and from
 * the wall below at the other ( mixed ), or from the walls above at both
 * ends ( uniform ).
 *
 * Every case runs create_magnetic_kernel() in its own fresh directory, which
 * is also the working directory for the run, so that the rho database of the
 * conductor material does not land in the build tree. A refused case returns
 * the message; an accepted case returns the log and the block ids on both
 * sides of every shell facet after the kernel was built.
 */

#include <gtest/gtest.h>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_Element.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "cl_MaxwellFactory.hpp"
#include "en_DomainType.hpp"

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

        const std::string &
        path() const
        {
            return mPath ;
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

    // the working directory for one scope, restored on exit
    class WorkingDir
    {
        std::filesystem::path mPrevious ;

    public:

        explicit WorkingDir( const std::string & aPath ) :
            mPrevious( std::filesystem::current_path() )
        {
            std::filesystem::current_path( aPath );
        }

        WorkingDir( const WorkingDir & ) = delete ;
        WorkingDir( WorkingDir && ) = delete ;
        WorkingDir & operator=( const WorkingDir & ) = delete ;
        WorkingDir & operator=( WorkingDir && ) = delete ;

        ~WorkingDir()
        {
            std::error_code tIgnored ;
            std::filesystem::current_path( mPrevious, tIgnored );
        }
    };

    const id_t gAir   = 1 ;
    const id_t gShell = 10 ;

    // the number of cells along z: two, so that the air has a node off
    // every sideset for the gauge pin
    const uint gNz = 2 ;

    // the faces of the cells in the plane aAxis = aPlane whose cell indices
    // along the other two axes, in increasing axis order, lie in
    // [ aLo0, aHi0 ) x [ aLo1, aHi1 ); all of them form sideset aTag
    struct Face
    {
        id_t  mTag ;
        uint  mAxis ;
        uint  mPlane ;
        uint  mLo0 ;
        uint  mHi0 ;
        uint  mLo1 ;
        uint  mHi1 ;
    };

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

    // the box [ 0, aNx ] x [ 0, aNy ] x [ 0, gNz ] in unit cells; cell
    // ( i, j, k ) belongs to block aBlock( i, j ). The faces z = 0 and z = gNz
    // are the boundary sidesets 21 and 22, the four lateral faces sideset
    // 23, all written after aFaces
    std::string
    grid_mesh( const uint aNx, const uint aNy,
               const std::function< id_t( uint, uint ) > & aBlock,
               const Cell< Face > & aFaces )
    {
        const uint tN[ 3 ] = { aNx, aNy, gNz };

        auto tNodeID = [ & ]( const uint aI, const uint aJ, const uint aK ) -> id_t
        {
            return 1 + aI + ( tN[ 0 ] + 1 ) * ( aJ + ( tN[ 1 ] + 1 ) * aK );
        };

        std::string tNodes ;
        index_t tNumNodes = 0 ;
        for ( uint k = 0; k <= tN[ 2 ]; ++k )
        {
            for ( uint j = 0; j <= tN[ 1 ]; ++j )
            {
                for ( uint i = 0; i <= tN[ 0 ]; ++i )
                {
                    tNodes += std::to_string( tNodeID( i, j, k ) ) + " "
                            + std::to_string( i ) + " " + std::to_string( j ) + " "
                            + std::to_string( k ) + "\n" ;
                    ++tNumNodes ;
                }
            }
        }

        // the six Kuhn tetrahedra of the unit cube as corner bit sets
        // v = x + 2 y + 4 z, each with positive volume
        Cell< Cell< uint > > tKuhn ;
        {
            const uint tPermutations[ 6 ][ 3 ] = {
                { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 },
                { 1, 2, 0 }, { 2, 0, gNz }, { 2, 1, 0 } };

            for ( const auto & tAxes : tPermutations )
            {
                Cell< uint > tTet = { 0u, 1u << tAxes[ 0 ],
                                      ( 1u << tAxes[ 0 ] ) | ( 1u << tAxes[ 1 ] ), 7u };

                int e[ 3 ][ 3 ];
                for ( uint r = 0; r < 3; ++r )
                {
                    for ( uint c = 0; c < 3; ++c )
                    {
                        e[ r ][ c ] = ( int ) ( ( tTet( r + 1 ) >> c ) & 1u );
                    }
                }
                const int tDet = e[ 0 ][ 0 ] * ( e[ 1 ][ 1 ] * e[ 2 ][ 2 ] - e[ 1 ][ 2 ] * e[ 2 ][ 1 ] )
                               - e[ 0 ][ 1 ] * ( e[ 1 ][ 0 ] * e[ 2 ][ 2 ] - e[ 1 ][ 2 ] * e[ 2 ][ 0 ] )
                               + e[ 0 ][ 2 ] * ( e[ 1 ][ 0 ] * e[ 2 ][ 1 ] - e[ 1 ][ 1 ] * e[ 2 ][ 0 ] );
                if ( tDet < 0 )
                {
                    std::swap( tTet( 2 ), tTet( 3 ) );
                }
                tKuhn.push( tTet );
            }
        }

        auto tCorner = [ & ]( const uint aCell[ 3 ], const uint aV ) -> id_t
        {
            return tNodeID( aCell[ 0 ] + ( aV & 1u ),
                            aCell[ 1 ] + ( ( aV >> 1 ) & 1u ),
                            aCell[ 2 ] + ( ( aV >> 2 ) & 1u ) );
        };

        id_t tID = 1 ;
        std::string tElements ;

        for ( uint k = 0; k < tN[ 2 ]; ++k )
        {
            for ( uint j = 0; j < tN[ 1 ]; ++j )
            {
                for ( uint i = 0; i < tN[ 0 ]; ++i )
                {
                    const uint tCell[ 3 ] = { i, j, k };
                    for ( const Cell< uint > & tTet : tKuhn )
                    {
                        tElements += element_line( tID, 4, aBlock( i, j ),
                            { tCorner( tCell, tTet( 0 ) ), tCorner( tCell, tTet( 1 ) ),
                              tCorner( tCell, tTet( 2 ) ), tCorner( tCell, tTet( 3 ) ) } );
                    }
                }
            }
        }

        // the triangles of the cells next to the plane that lie in it
        auto tFaceTriangles = [ & ]( const Face & aFace )
        {
            const uint tB = aFace.mAxis == 0 ? 1 : 0 ;
            const uint tC = aFace.mAxis == 2 ? 1 : 2 ;
            const uint tLayer = aFace.mPlane < tN[ aFace.mAxis ] ? aFace.mPlane : aFace.mPlane - 1 ;

            for ( uint b = aFace.mLo0; b < aFace.mHi0; ++b )
            {
                for ( uint c = aFace.mLo1; c < aFace.mHi1; ++c )
                {
                    uint tCell[ 3 ];
                    tCell[ aFace.mAxis ] = tLayer ;
                    tCell[ tB ] = b ;
                    tCell[ tC ] = c ;

                    for ( const Cell< uint > & tTet : tKuhn )
                    {
                        for ( uint tOmit = 0; tOmit < 4; ++tOmit )
                        {
                            Cell< id_t > tTriangle ;
                            for ( uint k = 0; k < 4; ++k )
                            {
                                const uint tCoord = tCell[ aFace.mAxis ]
                                    + ( ( tTet( k ) >> aFace.mAxis ) & 1u );
                                if ( k != tOmit && tCoord == aFace.mPlane )
                                {
                                    tTriangle.push( tCorner( tCell, tTet( k ) ) );
                                }
                            }
                            if ( tTriangle.size() == 3 )
                            {
                                tElements += element_line( tID, 2, aFace.mTag, tTriangle );
                            }
                        }
                    }
                }
            }
        };

        for ( const Face & tFace : aFaces )
        {
            tFaceTriangles( tFace );
        }
        tFaceTriangles( { 21, 2, 0, 0, tN[ 0 ], 0, tN[ 1 ] } );
        tFaceTriangles( { 22, 2, gNz, 0, tN[ 0 ], 0, tN[ 1 ] } );
        tFaceTriangles( { 23, 0, 0, 0, tN[ 1 ], 0, gNz } );
        tFaceTriangles( { 23, 0, tN[ 0 ], 0, tN[ 1 ], 0, gNz } );
        tFaceTriangles( { 23, 1, 0, 0, tN[ 0 ], 0, gNz } );
        tFaceTriangles( { 23, 1, tN[ 1 ], 0, tN[ 0 ], 0, gNz } );

        return "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n"
             + std::to_string( tNumNodes ) + "\n" + tNodes
             + "$EndNodes\n$Elements\n"
             + std::to_string( tID - 1 ) + "\n" + tElements
             + "$EndElements\n" ;
    }

    // one thin shell on sideset 10, air in block 1, and copper in the
    // conductor blocks aConductors
    std::string
    deck( const std::string & aMeshPath, const std::string & aConductors )
    {
        return "mesh\n{\n    file : " + aMeshPath + " ;\n    unit : m ;\n}\n\n"
               "solver\n{\n    linear\n    {\n        library : mumps ;\n    }\n}\n\n"
               "materials\n{\n    copper\n    {\n        builtin : copper ;\n        RRR : 50 ;\n    }\n}\n\n"
               "layers : tape\n{\n    copper : 20 mum ;\n}\n\n"
               "homology\n{\n    algorithm : generalized pellikka ;\n}\n\n"
               "topology\n{\n"
               "    thinshell : tape\n    {\n        sidesets : 10 ;\n    }\n"
               "    air\n    {\n        blocks : 1 ;\n    }\n"
               "    conductor\n    {\n        blocks : " + aConductors + " ;\n"
               "        material : copper ;\n    }\n"
               "}\n\n"
               "boundary conditions\n{\n}\n" ;
    }

    struct Outcome
    {
        // what create_magnetic_kernel() threw, empty if it returned
        std::string mMessage ;

        // what the run wrote to stdout
        std::string mLog ;

        // the block ids on the master and on the slave side of every shell
        // facet
        Cell< id_t > mMasterBlocks ;
        Cell< id_t > mSlaveBlocks ;
    };

    Outcome
    run( const std::string & aName, const std::string & aMesh, const std::string & aConductors )
    {
        // a fresh directory also holds no processed mesh ( .bfm ) the
        // factory could reuse instead of reading the .msh
        ScratchDir tDir( aName );
        const std::string tMesh = tDir.write( aName + ".msh", aMesh );
        const std::string tDeck = tDir.write( aName + ".conf", deck( tMesh, aConductors ) );

        WorkingDir tWorkingDir( tDir.path() );

        Outcome tOutcome ;

        fem::MaxwellFactory tFactory( tDeck );

        ::testing::internal::CaptureStdout();
        try
        {
            tFactory.create_magnetic_kernel();
        }
        catch ( const std::runtime_error & aError )
        {
            tOutcome.mMessage = aError.what();
        }
        tOutcome.mLog = ::testing::internal::GetCapturedStdout();

        // the thin-shell factory moves the facets of the shell sidesets, in
        // their order, into one geometry-only sideset per shell; a refused
        // run stops before that and leaves them in the shell sideset
        for ( mesh::SideSet * tSideSet : tFactory.mesh()->sidesets() )
        {
            if ( tSideSet->id() != gShell && tSideSet->domain_type() != DomainType::GeometryOnly ) continue ;

            for ( mesh::Facet * tFacet : tSideSet->facets() )
            {
                tOutcome.mMasterBlocks.push( tFacet->master()->block_id() );
                tOutcome.mSlaveBlocks.push( tFacet->has_slave() ? tFacet->slave()->block_id() : gNoID );
            }
        }

        return tOutcome ;
    }

    bool
    contains( const std::string & aText, const std::string & aPart )
    {
        return aText.find( aPart ) != std::string::npos ;
    }

    // the shell facets of aBlockA | aBlockB at positions [ aFirst, aEnd )
    // all have the same master block, which is one of the two
    id_t
    common_master( const Outcome & aOutcome, const index_t aFirst, const index_t aEnd,
                   const id_t aBlockA, const id_t aBlockB )
    {
        const id_t tMaster = aOutcome.mMasterBlocks( aFirst );
        const id_t tSlave = tMaster == aBlockA ? aBlockB : aBlockA ;

        EXPECT_TRUE( tMaster == aBlockA || tMaster == aBlockB ) << "master block " << tMaster ;

        for ( index_t f = aFirst; f < aEnd; ++f )
        {
            EXPECT_EQ( aOutcome.mMasterBlocks( f ), tMaster ) << "shell facet " << f ;
            EXPECT_EQ( aOutcome.mSlaveBlocks( f ), tSlave ) << "shell facet " << f ;
        }
        return tMaster ;
    }

    // conductor A ( block 2 ) above the shell, B ( block 3 ) below it, both
    // spanning x in [ 1, 3 ); the walls at x = 1 and x = 3 are numbered so
    // that the walk reaches the left end of the shell first from the wall
    // above it when aAboveLeft is true, and the right end when aAboveRight is
    std::string
    one_patch_mesh( const bool aAboveLeft, const bool aAboveRight )
    {
        auto tBlock = []( const uint aI, const uint aJ ) -> id_t
        {
            if ( aI < 1 || aI >= 3 ) return gAir ;
            if ( aJ == 2 ) return 2 ;
            if ( aJ == 1 ) return 3 ;
            return gAir ;
        };

        // a wall at x = aPlane, above ( y in [ 2, 3 ) ) or below ( [ 1, 2 ) )
        auto tWall = []( const id_t aTag, const uint aPlane, const bool aAbove ) -> Face
        {
            return { aTag, 0, aPlane, aAbove ? 2u : 1u, aAbove ? 3u : 2u, 0, gNz };
        };

        return grid_mesh( 5, 5, tBlock, {
            { gShell, 1, 2, 1, 3, 0, gNz },
            tWall( 11, 1, aAboveLeft ),
            tWall( 12, 3, aAboveRight ),
            tWall( 13, 1, ! aAboveLeft ),
            tWall( 14, 3, ! aAboveRight ),
            { 15, 1, 3, 1, 3, 0, gNz },
            { 16, 1, 1, 1, 3, 0, gNz } } );
    }

    // two shell pieces in the plane y = 2, for x in [ 1, 2 ) and [ 3, 4 ),
    // which share no edge; the conductor column spans x in [ 1, 4 ). The
    // walk reaches the left piece from the wall above it and the right piece
    // from the wall below it, so the pieces come out anti-parallel
    std::string
    two_patch_mesh( const std::function< id_t( uint, uint ) > & aBlock )
    {
        return grid_mesh( 6, 5, aBlock, {
            { gShell, 1, 2, 1, 2, 0, gNz },
            { gShell, 1, 2, 3, 4, 0, gNz },
            { 11, 0, 1, 2, 3, 0, gNz },
            { 12, 0, 4, 1, 2, 0, gNz },
            { 13, 0, 1, 1, 2, 0, gNz },
            { 14, 0, 4, 2, 3, 0, gNz },
            { 15, 1, 3, 1, 4, 0, gNz },
            { 16, 1, 1, 1, 4, 0, gNz } } );
    }
}

//------------------------------------------------------------------------------

// the walk gives the left and the right half of the shell opposite
// windings; the alignment flips one half, and the run completes with one
// master block on the whole sideset
TEST( ThinShellAlignment, MixedWindingIsAligned )
{
    const Outcome tOutcome = run( "tsa_mixed", one_patch_mesh( true, false ), "2, 3" );

    EXPECT_TRUE( tOutcome.mMessage.empty() ) << tOutcome.mMessage ;
    EXPECT_TRUE( contains( tOutcome.mLog, "thin shell tape: aligned " ) ) << tOutcome.mLog ;
    EXPECT_TRUE( contains( tOutcome.mLog, " of 8 facets of sideset 10 to one orientation" ) ) << tOutcome.mLog ;

    ASSERT_EQ( tOutcome.mMasterBlocks.size(), ( index_t ) 8 );
    common_master( tOutcome, 0, 8, 2, 3 );
}

//------------------------------------------------------------------------------

// one block on both sides of the shell is refused, whatever the winding
TEST( ThinShellAlignment, OneBlockOnBothSidesIsRefused )
{
    auto tBlock = []( const uint aI, const uint aJ ) -> id_t
    {
        return ( aI >= 1 && aI < 3 && ( aJ == 1 || aJ == 2 ) ) ? 2 : gAir ;
    };

    const std::string tMesh = grid_mesh( 5, 5, tBlock, {
        { gShell, 1, 2, 1, 3, 0, gNz },
        { 11, 0, 1, 1, 3, 0, gNz },
        { 12, 0, 3, 1, 3, 0, gNz },
        { 15, 1, 3, 1, 3, 0, gNz },
        { 16, 1, 1, 1, 3, 0, gNz } } );

    const Outcome tOutcome = run( "tsa_one_block", tMesh, "2" );

    EXPECT_TRUE( contains( tOutcome.mMessage, "thin shell tape: block 2 lies on both sides of sideset 10" ) )
        << tOutcome.mMessage ;
    EXPECT_TRUE( contains( tOutcome.mMessage, "a thin shell must separate two different volume blocks" ) )
        << tOutcome.mMessage ;
    EXPECT_TRUE( contains( tOutcome.mMessage, "Split the volume block, not the sideset" ) ) << tOutcome.mMessage ;
}

//------------------------------------------------------------------------------

// the control for MixedWindingIsAligned: the walk reaches both halves from
// the walls above, so nothing is aligned
TEST( ThinShellAlignment, UniformWindingIsLeftAlone )
{
    const Outcome tOutcome = run( "tsa_uniform", one_patch_mesh( true, true ), "2, 3" );

    EXPECT_TRUE( tOutcome.mMessage.empty() ) << tOutcome.mMessage ;
    EXPECT_FALSE( contains( tOutcome.mLog, "aligned" ) ) << tOutcome.mLog ;
    EXPECT_FALSE( contains( tOutcome.mLog, "patches" ) ) << tOutcome.mLog ;

    ASSERT_EQ( tOutcome.mMasterBlocks.size(), ( index_t ) 8 );
    common_master( tOutcome, 0, 8, 2, 3 );
}

//------------------------------------------------------------------------------

// two pieces with the same two blocks on their sides, anti-parallel
// after the walk, get the same master block
TEST( ThinShellAlignment, PatchesWithOneBlockPairAreAligned )
{
    auto tBlock = []( const uint aI, const uint aJ ) -> id_t
    {
        if ( aI < 1 || aI >= 4 ) return gAir ;
        if ( aJ == 2 ) return 2 ;
        if ( aJ == 1 ) return 3 ;
        return gAir ;
    };

    const Outcome tOutcome = run( "tsa_same_pair", two_patch_mesh( tBlock ), "2, 3" );

    EXPECT_TRUE( tOutcome.mMessage.empty() ) << tOutcome.mMessage ;
    EXPECT_TRUE( contains( tOutcome.mLog,
        "sideset 10 falls apart into 2 patches; each took the majority winding of its own facets "
        "( the seed's on a tie ); 1 were then aligned to the master block of a patch with the same "
        "two blocks ( 1 flipped ), 0 kept their winding because no other patch has their block pair, "
        "0 kept their winding because their facets do not share one block pair" ) ) << tOutcome.mLog ;

    ASSERT_EQ( tOutcome.mMasterBlocks.size(), ( index_t ) 8 );
    common_master( tOutcome, 0, 8, 2, 3 );
}

//------------------------------------------------------------------------------

// two pieces with different block pairs keep the winding the walk gave
// them, anti-parallel, and are reported as unrelated
TEST( ThinShellAlignment, PatchesWithDistinctBlockPairsAreUnrelated )
{
    auto tBlock = []( const uint aI, const uint aJ ) -> id_t
    {
        if ( aI < 1 || aI >= 4 ) return gAir ;
        const bool tLeft = aI < 3 ;
        if ( aJ == 2 ) return tLeft ? 2 : 4 ;
        if ( aJ == 1 ) return tLeft ? 3 : 5 ;
        return gAir ;
    };

    const Outcome tOutcome = run( "tsa_distinct_pairs", two_patch_mesh( tBlock ), "2:5" );

    EXPECT_TRUE( tOutcome.mMessage.empty() ) << tOutcome.mMessage ;
    EXPECT_TRUE( contains( tOutcome.mLog,
        "sideset 10 falls apart into 2 patches; each took the majority winding of its own facets "
        "( the seed's on a tie ); 0 were then aligned to the master block of a patch with the same "
        "two blocks ( 0 flipped ), 2 kept their winding because no other patch has their block pair, "
        "0 kept their winding because their facets do not share one block pair" ) ) << tOutcome.mLog ;
    EXPECT_FALSE( contains( tOutcome.mLog, "facets of sideset 10 to one orientation" ) ) << tOutcome.mLog ;

    ASSERT_EQ( tOutcome.mMasterBlocks.size(), ( index_t ) 8 );

    // the left piece is the first four facets, the right piece the last four
    const id_t tLeft = common_master( tOutcome, 0, 4, 2, 3 );
    const id_t tRight = common_master( tOutcome, 4, 8, 4, 5 );

    // anti-parallel: the master is above the shell on one piece and below
    // it on the other
    EXPECT_NE( tLeft == 2, tRight == 4 );
}

//------------------------------------------------------------------------------

// a shell on the domain boundary has no element on its outer side, and the
// cut cannot open it
TEST( ThinShellAlignment, ShellOnTheDomainBoundaryIsRefused )
{
    auto tBlock = []( const uint aI, const uint aJ ) -> id_t
    {
        return ( aI >= 1 && aI < 3 && aJ >= 1 && aJ < 3 ) ? 2 : gAir ;
    };

    const std::string tMesh = grid_mesh( 5, 5, tBlock, {
        { gShell, 1, 0, 1, 3, 0, gNz },
        { 11, 0, 1, 1, 3, 0, gNz },
        { 12, 0, 3, 1, 3, 0, gNz },
        { 15, 1, 3, 1, 3, 0, gNz },
        { 16, 1, 1, 1, 3, 0, gNz } } );

    const Outcome tOutcome = run( "tsa_boundary", tMesh, "2" );

    EXPECT_TRUE( contains( tOutcome.mMessage, "of sideset 10 needs an element on both sides" ) )
        << tOutcome.mMessage ;
}

//------------------------------------------------------------------------------

// one connected shell, x in [ 1, 3 ), between a conductor above it on the
// left half and a different conductor below it on the right half, with air
// on the other side of each half. Its cross-type facets disagree on the
// side of the conductor, the vote keeps one winding, and air ends up on
// both sides of the sideset, which the block-set check refuses
TEST( ThinShellAlignment, SplitVoteIsRefusedByTheBlockSets )
{
    auto tBlock = []( const uint aI, const uint aJ ) -> id_t
    {
        if ( aI == 1 && aJ == 2 ) return 2 ;
        if ( aI == 2 && aJ == 1 ) return 3 ;
        return gAir ;
    };

    const std::string tMesh = grid_mesh( 5, 5, tBlock, {
        { gShell, 1, 2, 1, 3, 0, gNz },
        { 11, 0, 1, 2, 3, 0, gNz },
        { 12, 0, 2, 2, 3, 0, gNz },
        { 13, 1, 3, 1, 2, 0, gNz },
        { 14, 0, 2, 1, 2, 0, gNz },
        { 15, 0, 3, 1, 2, 0, gNz },
        { 16, 1, 1, 2, 3, 0, gNz } } );

    const Outcome tOutcome = run( "tsa_split_vote", tMesh, "2, 3" );

    EXPECT_TRUE( contains( tOutcome.mLog,
        "thin shell tape: the cross-type facets of sideset 10 disagreed on which side the higher "
        "domain type lies; the vote kept the majority winding" ) ) << tOutcome.mLog ;
    EXPECT_TRUE( contains( tOutcome.mMessage, "thin shell tape: block 1 lies on both sides of sideset 10" ) )
        << tOutcome.mMessage ;
}

//------------------------------------------------------------------------------

// one connected shell between conductor 2 above it and conductors 3 and 4
// below it, wound uniformly by the walk from the walls above: its facets do
// not share one block pair, which is valid. Nothing is aligned, a single
// patch prints no patch line, and the master and slave block sets stay
// disjoint
TEST( ThinShellAlignment, MixedPairPatchIsAccepted )
{
    auto tBlock = []( const uint aI, const uint aJ ) -> id_t
    {
        if ( aI < 1 || aI >= 3 ) return gAir ;
        if ( aJ == 2 ) return 2 ;
        if ( aJ == 1 ) return aI == 1 ? 3 : 4 ;
        return gAir ;
    };

    const std::string tMesh = grid_mesh( 5, 5, tBlock, {
        { gShell, 1, 2, 1, 3, 0, gNz },
        { 11, 0, 1, 2, 3, 0, gNz },
        { 12, 0, 3, 2, 3, 0, gNz },
        { 13, 0, 1, 1, 2, 0, gNz },
        { 14, 0, 3, 1, 2, 0, gNz },
        { 15, 1, 3, 1, 3, 0, gNz },
        { 16, 1, 1, 1, 3, 0, gNz } } );

    const Outcome tOutcome = run( "tsa_mixed_pair", tMesh, "2:4" );

    EXPECT_TRUE( tOutcome.mMessage.empty() ) << tOutcome.mMessage ;
    EXPECT_FALSE( contains( tOutcome.mLog, "aligned" ) ) << tOutcome.mLog ;
    EXPECT_FALSE( contains( tOutcome.mLog, "patches" ) ) << tOutcome.mLog ;

    ASSERT_EQ( tOutcome.mMasterBlocks.size(), ( index_t ) 8 );

    // block 2 is on the same side of every facet: all masters or all slaves
    const bool tAboveIsMaster = tOutcome.mMasterBlocks( 0 ) == 2 ;
    for ( index_t f = 0; f < 8; ++f )
    {
        const id_t tAbove = tAboveIsMaster ? tOutcome.mMasterBlocks( f ) : tOutcome.mSlaveBlocks( f );
        const id_t tBelow = tAboveIsMaster ? tOutcome.mSlaveBlocks( f ) : tOutcome.mMasterBlocks( f );
        EXPECT_EQ( tAbove, ( id_t ) 2 ) << "shell facet " << f ;
        EXPECT_TRUE( tBelow == 3 || tBelow == 4 ) << "shell facet " << f << ", block " << tBelow ;
    }
}

//------------------------------------------------------------------------------

// two shell pieces in the plane y = 2, for x in [ 1, 3 ) and [ 4, 5 ), with
// conductor 2 above both. Below the left piece lie conductors 4 and 3, below
// the right piece conductor 3, so the left piece is a mixed-pair patch and
// the right piece the only patch with the pair 2|3. The walk reaches both
// from the walls above; the multi-patch line counts one patch of each kind
TEST( ThinShellAlignment, MixedPairPieceIsCountedSeparately )
{
    auto tBlock = []( const uint aI, const uint aJ ) -> id_t
    {
        if ( aI < 1 || aI >= 5 ) return gAir ;
        if ( aJ == 2 ) return 2 ;
        if ( aJ == 1 ) return aI == 1 ? 4 : 3 ;
        return gAir ;
    };

    const std::string tMesh = grid_mesh( 7, 5, tBlock, {
        { gShell, 1, 2, 1, 3, 0, gNz },
        { gShell, 1, 2, 4, 5, 0, gNz },
        { 11, 0, 1, 2, 3, 0, gNz },
        { 12, 0, 5, 2, 3, 0, gNz },
        { 13, 0, 1, 1, 2, 0, gNz },
        { 14, 0, 5, 1, 2, 0, gNz },
        { 15, 1, 3, 1, 5, 0, gNz },
        { 16, 1, 1, 1, 5, 0, gNz } } );

    const Outcome tOutcome = run( "tsa_mixed_piece", tMesh, "2:4" );

    EXPECT_TRUE( tOutcome.mMessage.empty() ) << tOutcome.mMessage ;
    EXPECT_TRUE( contains( tOutcome.mLog,
        "sideset 10 falls apart into 2 patches; each took the majority winding of its own facets "
        "( the seed's on a tie ); 0 were then aligned to the master block of a patch with the same "
        "two blocks ( 0 flipped ), 1 kept their winding because no other patch has their block pair, "
        "1 kept their winding because their facets do not share one block pair" ) ) << tOutcome.mLog ;

    ASSERT_EQ( tOutcome.mMasterBlocks.size(), ( index_t ) 12 );

    // block 2 is on the same side of every facet of both pieces
    const bool tAboveIsMaster = tOutcome.mMasterBlocks( 0 ) == 2 ;
    for ( index_t f = 0; f < 12; ++f )
    {
        const id_t tAbove = tAboveIsMaster ? tOutcome.mMasterBlocks( f ) : tOutcome.mSlaveBlocks( f );
        EXPECT_EQ( tAbove, ( id_t ) 2 ) << "shell facet " << f ;
    }
}
