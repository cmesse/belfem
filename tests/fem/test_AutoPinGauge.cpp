/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_MaxwellFactory.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Dof.hpp"

using namespace belfem;

namespace
{
    // runs aCall, requires a std::runtime_error whose what() carries the
    // literal part of the BELFEM_ERROR text
    template< typename F >
    void
    expect_refusal( F && aCall, const char * aSubstring )
    {
        bool tThrown = false ;
        try
        {
            aCall();
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true ;
            EXPECT_NE( std::string( aError.what() ).find( aSubstring ), std::string::npos )
                << "message was: " << aError.what();
        }
        EXPECT_TRUE( tThrown ) << "expected a refusal carrying \"" << aSubstring << "\"";
    }
}

namespace
{
    // the contract under test: with no bearing in the deck, the factory fixes
    // one automatic gauge pin per connected phi component, EXCEPT in a
    // component whose phi is already fixed by a symmetry plane or a dirichlet
    // sideset; a deck bearing disables every automatic pin

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

    const uint gNx = 5 ;
    const uint gNy = 2 ;

    //! the strip [ 0, 5 ] x [ 0, 2 ] in unit cells, two TRI3 per cell: air
    //! in the columns x < 2 ( block 1 ) and x > 3 ( block 3 ), copper in
    //! the column 2 <= x <= 3 ( block 2 ) across the full height, so the air
    //! falls into two phi components. The edge x = 0 is sideset 24, the edge
    //! x = 5 sideset 25; the edges y = 0 and y = 2 are untagged. Node
    //! ( i, j ) has id 1 + i + 6 j; node 11 ( 4, 1 ) also carries a point
    //! element with tag 11, the vertex a deck bearing can name
    std::string
    strip_mesh()
    {
        auto tNodeID = []( const uint aI, const uint aJ ) -> id_t
        {
            return 1 + aI + ( gNx + 1 ) * aJ ;
        };

        std::string tNodes ;
        for ( uint j = 0; j <= gNy; ++j )
        {
            for ( uint i = 0; i <= gNx; ++i )
            {
                tNodes += std::to_string( tNodeID( i, j ) ) + " "
                        + std::to_string( i ) + " " + std::to_string( j ) + " 0\n" ;
            }
        }

        // one .msh 2.2 element line; both tags are the group id
        id_t tID = 1 ;
        std::string tElements ;
        auto tLine = [ & ]( const uint aGmshType, const id_t aTag, const Cell< id_t > & aNodes )
        {
            tElements += std::to_string( tID++ ) + " " + std::to_string( aGmshType ) + " 2 "
                       + std::to_string( aTag ) + " " + std::to_string( aTag );
            for ( id_t tNode : aNodes )
            {
                tElements += " " + std::to_string( tNode );
            }
            tElements += "\n" ;
        };

        for ( uint j = 0; j < gNy; ++j )
        {
            for ( uint i = 0; i < gNx; ++i )
            {
                const id_t tBlock = i < 2 ? 1 : ( i == 2 ? 2 : 3 );
                tLine( 2, tBlock, { tNodeID( i, j ), tNodeID( i + 1, j ), tNodeID( i + 1, j + 1 ) } );
                tLine( 2, tBlock, { tNodeID( i, j ), tNodeID( i + 1, j + 1 ), tNodeID( i, j + 1 ) } );
            }
        }
        for ( uint j = 0; j < gNy; ++j )
        {
            tLine( 1, 24, { tNodeID( 0, j ), tNodeID( 0, j + 1 ) } );
            tLine( 1, 25, { tNodeID( gNx, j ), tNodeID( gNx, j + 1 ) } );
        }
        tLine( 15, 11, { tNodeID( 4, 1 ) } );

        return "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n"
             + std::to_string( ( gNx + 1 ) * ( gNy + 1 ) ) + "\n" + tNodes
             + "$EndNodes\n$Elements\n"
             + std::to_string( tID - 1 ) + "\n" + tElements
             + "$EndElements\n" ;
    }

    //! aTopology goes into the topology section, aConditions into the
    //! boundary conditions section
    std::string
    deck( const std::string & aMeshPath, const std::string & aTopology, const std::string & aConditions )
    {
        return "mesh\n{\n    file : " + aMeshPath + " ;\n    unit : m ;\n}\n\n"
               "solver\n{\n    linear\n    {\n        library : mumps ;\n    }\n}\n\n"
               "materials\n{\n    copper\n    {\n        builtin : copper ;\n        RRR : 50 ;\n    }\n}\n\n"
               "homology\n{\n    algorithm : generalized pellikka ;\n}\n\n"
               "topology\n{\n"
               "    air\n    {\n        blocks : 1, 3 ;\n    }\n"
               "    conductor\n    {\n        blocks : 2 ;\n        material : copper ;\n    }\n"
               + aTopology +
               "}\n\n"
               "boundary conditions\n{\n" + aConditions + "}\n" ;
    }

    struct Outcome
    {
        // what the run wrote to stdout
        std::string mLog ;

        // the pin of the left ( x < 2.5 ) and the right air component, and
        // whether its phi dof is fixed
        id_t mLeftPin  = gNoID ;
        id_t mRightPin = gNoID ;
        bool mLeftFixed  = false ;
        bool mRightFixed = false ;
        uint mNumPins = 0 ;

        // whether node 11 ( x = 4, y = 1 ) is fixed
        bool mNode11Fixed = false ;

        // the kernel, kept alive so a test can initialize its field and
        // impose the conditions of a later time
        std::shared_ptr< fem::Kernel > mKernel ;
    };

    // builds the kernel from the deck in aDir ( the mesh is written there
    // once; a second call in the same directory reuses the .bfm the first
    // one wrote ) and reads the pin state through the dof manager: the
    // factory has disconnected the dofs from the mesh when it returns
    Outcome
    run( const ScratchDir & aDir, const std::string & aName,
         const std::string & aTopology, const std::string & aConditions )
    {
        const std::string tMesh = aDir.path() + "strip.msh" ;
        if ( ! std::filesystem::exists( tMesh ) )
        {
            aDir.write( "strip.msh", strip_mesh() );
        }
        const std::string tDeck = aDir.write( aName + ".conf", deck( tMesh, aTopology, aConditions ) );

        WorkingDir tWorkingDir( aDir.path() );

        Outcome tOutcome ;

        // the constructor reads ( or reuses ) the mesh and says so
        ::testing::internal::CaptureStdout();
        fem::MaxwellFactory tFactory( tDeck );

        std::shared_ptr< fem::Kernel > tKernel ;
        try
        {
            tKernel = tFactory.create_magnetic_kernel();
        }
        catch ( const std::runtime_error & aError )
        {
            tOutcome.mLog = ::testing::internal::GetCapturedStdout();
            throw ;
        }
        tOutcome.mLog = ::testing::internal::GetCapturedStdout();

        fem::DofManager * tField = tKernel->dofmgr();

        auto tIsFixed = [ & ]( mesh::Node * aNode ) -> bool
        {
            return tField->dof( tField->calculate_dof_id( aNode, 0 ) )->is_fixed();
        };

        for ( mesh::Node * tNode : tFactory.mesh()->autopins() )
        {
            ++tOutcome.mNumPins ;
            if ( tNode->x() < 2.5 )
            {
                EXPECT_EQ( tOutcome.mLeftPin, gNoID ) << "two pins in the left component" ;
                tOutcome.mLeftPin   = tNode->id();
                tOutcome.mLeftFixed = tIsFixed( tNode );
            }
            else
            {
                EXPECT_EQ( tOutcome.mRightPin, gNoID ) << "two pins in the right component" ;
                tOutcome.mRightPin   = tNode->id();
                tOutcome.mRightFixed = tIsFixed( tNode );
            }
        }

        tOutcome.mNode11Fixed = tIsFixed( tFactory.mesh()->node( 11 ) );
        tOutcome.mKernel = tKernel ;

        return tOutcome ;
    }

    // the phi dof of node aNodeID of a finished run
    fem::Dof *
    node_dof( const Outcome & aOutcome, const id_t aNodeID )
    {
        fem::DofManager * tField = aOutcome.mKernel->dofmgr();
        return tField->dof( tField->calculate_dof_id( tField->mesh()->node( aNodeID ), 0 ) );
    }

    bool
    contains( const std::string & aText, const std::string & aPart )
    {
        return aText.find( aPart ) != std::string::npos ;
    }

    const std::string gSymmetry  = "    air symmetry\n    {\n        sidesets : 24 ;\n    }\n" ;
    const std::string gDirichlet = "    dirichlet\n    {\n        sidesets  : 24 ;\n        type      : constant ;\n        amplitude : 0 ;\n    }\n" ;
    const std::string gBearing   = "    bearing\n    {\n        nodes : 11 ;\n    }\n" ;

    // a gauge is a bearing with a value: point 11, constant 0, constant 3,
    // a ramp, and the untyped form the factory must refuse
    const std::string gGauge        = "    gauge\n    {\n        nodes     : 11 ;\n        type      : constant ;\n        amplitude : 0 ;\n    }\n" ;
    const std::string gGaugeThree   = "    gauge\n    {\n        nodes     : 11 ;\n        type      : constant ;\n        amplitude : 3 ;\n    }\n" ;
    const std::string gGaugeRamp    = "    gauge\n    {\n        nodes     : 11 ;\n        type      : ramp ;\n        amplitude : 2 ;\n        period    : 1 s ;\n        offset    : 0 s ;\n    }\n" ;
    const std::string gGaugeUntyped = "    gauge\n    {\n        nodes : 11 ;\n    }\n" ;
    const std::string gDirichletUntyped = "    dirichlet\n    {\n        sidesets : 24 ;\n    }\n" ;

    // writes the deck and constructs the factory ( the refusal of a deck
    // comes from the constructor, before any kernel exists )
    void
    construct( const ScratchDir & aDir, const std::string & aName, const std::string & aConditions )
    {
        if ( ! std::filesystem::exists( aDir.path() + "strip.msh" ) )
        {
            aDir.write( "strip.msh", strip_mesh() );
        }
        const std::string tDeck = aDir.write( aName + ".conf", deck( aDir.path() + "strip.msh", "", aConditions ) );
        WorkingDir tWorkingDir( aDir.path() );
        ::testing::internal::CaptureStdout();
        try
        {
            fem::MaxwellFactory tFactory( tDeck );
        }
        catch ( ... )
        {
            ::testing::internal::GetCapturedStdout();
            throw ;
        }
        ::testing::internal::GetCapturedStdout();
    }
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, both_components_pinned_without_a_fixing_plane )
{
    ScratchDir tDir( "autopin_none" );
    Outcome tOutcome = run( tDir, "none", "", "" );

    EXPECT_EQ( tOutcome.mNumPins, 2u );
    ASSERT_NE( tOutcome.mLeftPin,  gNoID );
    ASSERT_NE( tOutcome.mRightPin, gNoID );
    EXPECT_TRUE( tOutcome.mLeftFixed )  << "left pin " << tOutcome.mLeftPin ;
    EXPECT_TRUE( tOutcome.mRightFixed ) << "right pin " << tOutcome.mRightPin ;
    EXPECT_FALSE( contains( tOutcome.mLog, "no automatic gauge pin" ) ) << tOutcome.mLog ;
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, symmetry_plane_takes_the_pin_of_its_component )
{
    ScratchDir tDir( "autopin_symmetry" );
    Outcome tOutcome = run( tDir, "symmetry", gSymmetry, "" );

    EXPECT_EQ( tOutcome.mNumPins, 2u );
    ASSERT_NE( tOutcome.mLeftPin,  gNoID );
    ASSERT_NE( tOutcome.mRightPin, gNoID );
    EXPECT_FALSE( tOutcome.mLeftFixed ) << "left pin " << tOutcome.mLeftPin << " is fixed although sideset 24 fixes phi there" ;
    EXPECT_TRUE( tOutcome.mRightFixed ) << "right pin " << tOutcome.mRightPin ;
    EXPECT_TRUE( contains( tOutcome.mLog, "no automatic gauge pin for the phi component of node "
                                          + std::to_string( tOutcome.mLeftPin ) ) ) << tOutcome.mLog ;
    EXPECT_TRUE( contains( tOutcome.mLog, "( sideset 24 )" ) ) << tOutcome.mLog ;
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, dirichlet_sideset_takes_the_pin_of_its_component )
{
    ScratchDir tDir( "autopin_dirichlet" );
    Outcome tOutcome = run( tDir, "dirichlet", "", gDirichlet );

    EXPECT_EQ( tOutcome.mNumPins, 2u );
    ASSERT_NE( tOutcome.mLeftPin,  gNoID );
    ASSERT_NE( tOutcome.mRightPin, gNoID );
    EXPECT_FALSE( tOutcome.mLeftFixed ) << "left pin " << tOutcome.mLeftPin ;
    EXPECT_TRUE( tOutcome.mRightFixed ) << "right pin " << tOutcome.mRightPin ;
    EXPECT_TRUE( contains( tOutcome.mLog, "( sideset 24 )" ) ) << tOutcome.mLog ;
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, cached_pins_follow_the_deck_not_the_bfm )
{
    // the first run writes the .bfm with both pins; the second reuses it
    // under a deck that fixes phi on sideset 24 — the stored pin of the left
    // component must then stay free
    ScratchDir tDir( "autopin_cache" );
    Outcome tFirst = run( tDir, "none", "", "" );
    ASSERT_TRUE( tFirst.mLeftFixed && tFirst.mRightFixed );

    ASSERT_TRUE( std::filesystem::exists( tDir.path() + "strip.bfm" ) ) << "the first run wrote no .bfm" ;

    Outcome tSecond = run( tDir, "dirichlet", "", gDirichlet );
    EXPECT_TRUE( contains( tSecond.mLog, "Reusing processed mesh" ) ) << tSecond.mLog ;
    EXPECT_EQ( tSecond.mLeftPin, tFirst.mLeftPin );
    EXPECT_FALSE( tSecond.mLeftFixed ) << "left pin " << tSecond.mLeftPin ;
    EXPECT_TRUE( tSecond.mRightFixed );

    // and back: the same cache under the plain deck pins both again
    Outcome tThird = run( tDir, "none", "", "" );
    EXPECT_TRUE( tThird.mLeftFixed && tThird.mRightFixed );
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, a_deck_bearing_disables_every_automatic_pin )
{
    ScratchDir tDir( "autopin_bearing" );
    Outcome tOutcome = run( tDir, "bearing", "", gBearing );

    EXPECT_EQ( tOutcome.mNumPins, 2u );
    EXPECT_TRUE( tOutcome.mNode11Fixed );
    if ( tOutcome.mLeftPin != 11 )
    {
        EXPECT_FALSE( tOutcome.mLeftFixed ) << "left pin " << tOutcome.mLeftPin ;
    }
    if ( tOutcome.mRightPin != 11 )
    {
        EXPECT_FALSE( tOutcome.mRightFixed ) << "right pin " << tOutcome.mRightPin ;
    }
    EXPECT_FALSE( contains( tOutcome.mLog, "automatic gauge pin" ) ) << tOutcome.mLog ;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

TEST( AutoPinGauge, gauge_point_is_a_contact_of_its_component )
{
    ScratchDir tDir( "autopin_gauge" );
    Outcome tOutcome = run( tDir, "gauge", "", gGauge );

    EXPECT_EQ( tOutcome.mNumPins, 2u );
    ASSERT_NE( tOutcome.mLeftPin,  gNoID );
    ASSERT_NE( tOutcome.mRightPin, gNoID );
    EXPECT_TRUE( tOutcome.mNode11Fixed ) << "the gauge point is fixed at setup" ;
    EXPECT_TRUE( tOutcome.mLeftFixed )   << "left pin " << tOutcome.mLeftPin ;
    EXPECT_FALSE( tOutcome.mRightFixed ) << "right pin " << tOutcome.mRightPin << " is fixed although point 11 fixes phi in that component" ;
    EXPECT_TRUE( contains( tOutcome.mLog, "( gauge point 11 )" ) ) << tOutcome.mLog ;
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, gauge_value_reaches_the_dof_after_initialization )
{
    // the factory fixes the point to 0; the per-timestep imposition after
    // the field is initialized is a value update, not a first fix
    ScratchDir tDir( "autopin_gauge_value" );
    Outcome tOutcome = run( tDir, "gauge3", "", gGaugeThree );
    ASSERT_TRUE( tOutcome.mNode11Fixed );

    fem::Dof * tDof = node_dof( tOutcome, 11 );
    EXPECT_EQ( tDof->value(), 0.0 );

    ::testing::internal::CaptureStdout();
    tOutcome.mKernel->dofmgr()->initialize();
    tOutcome.mKernel->compute_boundary_conditions( 0.5 );
    ::testing::internal::GetCapturedStdout();

    EXPECT_TRUE( tDof->is_fixed() );
    EXPECT_EQ( tDof->value(), 3.0 );
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, gauge_ramp_updates_the_value_every_timestep )
{
    ScratchDir tDir( "autopin_gauge_ramp" );
    Outcome tOutcome = run( tDir, "ramp", "", gGaugeRamp );
    ASSERT_TRUE( tOutcome.mNode11Fixed );

    fem::Dof * tDof = node_dof( tOutcome, 11 );

    ::testing::internal::CaptureStdout();
    tOutcome.mKernel->dofmgr()->initialize();
    tOutcome.mKernel->compute_boundary_conditions( 0.25 );
    const real tFirst = tDof->value();
    tOutcome.mKernel->compute_boundary_conditions( 0.75 );
    const real tSecond = tDof->value();
    ::testing::internal::GetCapturedStdout();

    // amplitude 2 over a 1 s period from offset 0: 2 * t
    EXPECT_TRUE( tDof->is_fixed() );
    EXPECT_DOUBLE_EQ( tFirst,  0.5 );
    EXPECT_DOUBLE_EQ( tSecond, 1.5 );
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, cached_pins_honor_a_gauge_added_to_the_deck )
{
    ScratchDir tDir( "autopin_gauge_cache" );
    Outcome tFirst = run( tDir, "none", "", "" );
    ASSERT_TRUE( tFirst.mLeftFixed && tFirst.mRightFixed );

    Outcome tSecond = run( tDir, "gauge", "", gGauge );
    EXPECT_TRUE( contains( tSecond.mLog, "Reusing processed mesh" ) ) << tSecond.mLog ;
    EXPECT_TRUE( tSecond.mNode11Fixed );
    EXPECT_TRUE( tSecond.mLeftFixed );
    EXPECT_FALSE( tSecond.mRightFixed ) << "right pin " << tSecond.mRightPin ;
}

//------------------------------------------------------------------------------

TEST( AutoPinGauge, a_condition_without_a_type_is_refused_at_deck_read )
{
    ScratchDir tDir( "autopin_untyped" );
    expect_refusal( [ & ]() { construct( tDir, "gauge_untyped", gGaugeUntyped ); },
                    "'type' is required" );
    expect_refusal( [ & ]() { construct( tDir, "dirichlet_untyped", gDirichletUntyped ); },
                    "'type' is required" );
}

//------------------------------------------------------------------------------
