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
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "globals.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Facet.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_ThinShell.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_MeshChecker.hpp"
#include "meshtools.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_DofMgr_BlockData.hpp"
#include "cl_Map.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "en_Maxwell_Formulations.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"
#include "cl_Matrix.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_TimestepMatrices.hpp"
#include "support/cl_TS_TestStack.hpp"
#include <memory>
#include "cl_Protoshell.hpp"
#include "cl_ThinShellFactory.hpp"
#include "cl_Material.hpp"
#include "cl_IWG_MaxwellThermal.hpp"

using namespace belfem;

namespace
{
    //! two thin-shell stacks side by side, each the TS_TestGhostStack
    //! construction: nine nodes, two PENTA6TS layer elements in two layer
    //! blocks, a ghost facet between them, a tape ( mid-surface ) facet,
    //! and a ThinShell record. Stack s uses blocks 2s+1, 2s+2, ghost sideset
    //! 10(s+1), tape sideset 10(s+1)+1, ghost facet id 100(s+1), tape facet
    //! id 100(s+1)+1; node, element and edge ids continue across the stacks.
    //! The IWG selects both shells or only the second one. Serial only: the
    //! fixture never enters Kernel::partition_mesh, so the parallel linker
    //! is not exercised here
    class TwoShellStacks
    {
        //! h_ghost falls back to the nodal "T" field when gTbulk is NaN
        struct TbulkGuard
        {
            const real mSaved ;
            TbulkGuard() : mSaved( gTbulk ) { gTbulk = 77.0 ; }
            ~TbulkGuard() { gTbulk = mSaved ; }
        };

        TbulkGuard mTbulkGuard ;

        public:

            Mesh                  * mMesh   = nullptr ;
            fem::KernelParameters * mParams = nullptr ;
            fem::Kernel           * mKernel = nullptr ;
            fem::DofManager       * mField  = nullptr ;   // owned by the kernel
            fem::IWG_Maxwell      * mIWG    = nullptr ;   // owned by the kernel

            mesh::Facet   * mTapeFacet[ 2 ] = { nullptr, nullptr };
            mesh::Element * mLayer[ 2 ][ 2 ] = { { nullptr, nullptr }, { nullptr, nullptr } };

//------------------------------------------------------------------------------

            explicit TwoShellStacks( const bool aSelectAllShells )
            {
                const real tThickness = 0.1 ;

                // third argument false: the MeshChecker below must run
                // before edges exist
                mMesh = new Mesh( 3, 0, false );

                mesh::ElementFactory tFactory ;

                // triangle footprint, deliberately not equilateral
                const real tX[ 3 ] = { 0.0, 1.3, 0.9 };
                const real tY[ 3 ] = { 0.0, 0.1, 1.1 };

                Cell< mesh::Node * >    & tNodes    = mMesh->nodes();
                Cell< mesh::Element * > & tElements = mMesh->elements();

                id_t tNodeID    = 1 ;
                id_t tElementID = 1 ;

                for ( uint s = 0; s < 2; ++s )
                {
                    const real tOffset = 3.0 * s ;
                    const uint tBase   = tNodes.size();

                    for ( uint r = 0; r < 3; ++r )
                    {
                        for ( uint k = 0; k < 3; ++k )
                        {
                            tNodes.push( new mesh::Node( tNodeID++,
                                    tX[ k ] + tOffset, tY[ k ], tThickness * r ) );
                        }
                    }

                    const id_t tGhostSetID = 10 * ( s + 1 );
                    const id_t tTapeSetID  = tGhostSetID + 1 ;
                    const id_t tGhostFacetID = 100 * ( s + 1 );
                    const id_t tTapeFacetID  = tGhostFacetID + 1 ;

                    mesh::Block * tBlocks[ 2 ] = { nullptr, nullptr };
                    for ( uint l = 0; l < 2; ++l )
                    {
                        tBlocks[ l ] = new mesh::Block( 2 * s + l + 1, 1 );
                        tBlocks[ l ]->set_thickness( tThickness );
                        tBlocks[ l ]->set_domain_type( DomainType::ThinShell );

                        mesh::Element * tPrism =
                                tFactory.create_element( ElementType::PENTA6TS, tElementID++ );
                        for ( uint k = 0; k < 6; ++k )
                        {
                            tPrism->insert_node( tNodes( tBase + 3 * l + k ), k );
                        }
                        tPrism->set_block_id( tBlocks[ l ]->id() );
                        tElements.push( tPrism );
                        tBlocks[ l ]->insert_element( tPrism );
                        mMesh->blocks().push( tBlocks[ l ] );
                        mLayer[ s ][ l ] = tPrism ;
                    }

                    // ghost facet on the middle node row, master = lower
                    // prism's top face, slave = upper prism's bottom face
                    mesh::Element * tSurf = tFactory.create_element( ElementType::TRI3, tGhostFacetID );
                    for ( uint k = 0; k < 3; ++k )
                    {
                        tSurf->insert_node( tNodes( tBase + 3 + k ), k );
                    }
                    mesh::Facet * tGhost = new mesh::Facet( tSurf );
                    tGhost->set_master( mLayer[ s ][ 0 ], mesh::top_facet_index( ElementType::PENTA6TS ) );
                    tGhost->set_slave( mLayer[ s ][ 1 ], mesh::bottom_facet_index( ElementType::PENTA6TS ), 1 );
                    tGhost->set_sideset_id( tGhostSetID );

                    mesh::SideSet * tGhostSet = new mesh::SideSet( tGhostSetID, 1 );
                    tGhostSet->set_domain_type( DomainType::Ghost );
                    tGhostSet->insert_facet( tGhost );
                    mMesh->sidesets().push( tGhostSet );

                    // tape ( mid-surface ) facet: what link_thin_shell_facets
                    // hangs on every layer element of this shell
                    mesh::Element * tMidSurf = tFactory.create_element( ElementType::TRI3, tTapeFacetID );
                    for ( uint k = 0; k < 3; ++k )
                    {
                        tMidSurf->insert_node( tNodes( tBase + 3 + k ), k );
                    }
                    mesh::Facet * tMid = new mesh::Facet( tMidSurf );
                    tMid->set_sideset_id( tTapeSetID );
                    mesh::SideSet * tTapeSet = new mesh::SideSet( tTapeSetID, 1 );
                    tTapeSet->set_domain_type( DomainType::ThinShell );
                    tTapeSet->insert_facet( tMid );
                    mMesh->sidesets().push( tTapeSet );
                    mTapeFacet[ s ] = tMid ;

                    mesh::ThinShell * tShell = new mesh::ThinShell( tTapeSet, tGhostSet );
                    tShell->blocks().push( tBlocks[ 0 ] );
                    tShell->blocks().push( tBlocks[ 1 ] );
                    {
                        Vector< real > tThicknesses( 2, tThickness );
                        tShell->set_thicknesses( tThicknesses );
                    }
                    mMesh->thin_shells().push( tShell );
                }

                mMesh->finalize();

                // reorientation check runs BEFORE edges exist; it sets the
                // mesh checker flag so the Kernel ctor skips its own run
                {
                    MeshChecker tCheck( mMesh );
                }

                // per-element edge objects in canonical get_nodes_of_edge() order
                {
                    id_t tEdgeID = 1 ;
                    Cell< mesh::Node * > tEdgeNodes ;
                    for ( uint s = 0; s < 2; ++s )
                    {
                        for ( uint l = 0; l < 2; ++l )
                        {
                            mesh::Element * tElement = mLayer[ s ][ l ];
                            const uint tNumEdges = tElement->number_of_edges();
                            tElement->allocate_edge_container();
                            for ( uint e = 0; e < tNumEdges; ++e )
                            {
                                tElement->get_nodes_of_edge( e, tEdgeNodes );
                                mesh::Edge * tEdge = new mesh::Edge();
                                tEdge->set_id( tEdgeID++ );
                                tEdge->allocate_node_container( 2 );
                                tEdge->insert_node( tEdgeNodes( 0 ), 0 );
                                tEdge->insert_node( tEdgeNodes( 1 ), 1 );
                                tEdge->set_owner( 0 );
                                tElement->insert_edge( tEdge, e );
                                mMesh->edges().push( tEdge );
                            }
                        }
                    }
                    mMesh->finalize_edges();
                }

                mParams = new fem::KernelParameters( mMesh );
                mKernel = new fem::Kernel( mParams );

                mIWG = new fem::IWG_Maxwell(
                        fem::maxwell::Formulation::HPhi,
                        ModelDimensionality::ThreeD,
                        false, false );

                if ( aSelectAllShells )
                {
                    Vector< id_t > tBlockIDs( 4 );
                    tBlockIDs( 0 ) = 1 ; tBlockIDs( 1 ) = 2 ; tBlockIDs( 2 ) = 3 ; tBlockIDs( 3 ) = 4 ;
                    Cell< DomainType > tBlockTypes( 4, DomainType::ThinShell );
                    mIWG->set_blocks( tBlockIDs, tBlockTypes );

                    Vector< id_t > tSideSetIDs( 2 );
                    tSideSetIDs( 0 ) = 10 ; tSideSetIDs( 1 ) = 20 ;
                    Cell< DomainType > tSideSetTypes( 2, DomainType::Ghost );
                    mIWG->set_sidesets( tSideSetIDs, tSideSetTypes );
                }
                else
                {
                    Vector< id_t > tBlockIDs( 2 );
                    tBlockIDs( 0 ) = 3 ; tBlockIDs( 1 ) = 4 ;
                    Cell< DomainType > tBlockTypes( 2, DomainType::ThinShell );
                    mIWG->set_blocks( tBlockIDs, tBlockTypes );

                    Vector< id_t > tSideSetIDs( 1 );
                    tSideSetIDs( 0 ) = 20 ;
                    Cell< DomainType > tSideSetTypes( 1, DomainType::Ghost );
                    mIWG->set_sidesets( tSideSetIDs, tSideSetTypes );
                }

                mKernel->add_equation( mIWG );
            }

//------------------------------------------------------------------------------

            //! the step under test: create_field runs BlockData::create_blocks,
            //! which collects the thin-shell facet ids and links them
            void
            create_field()
            {
                mField = mKernel->create_field( mIWG );
                mField->create_fields( mField->iwg() );
                mMesh->field_data( "element_rho" ).fill( 1.0e-8 );
                mField->set_solver( gDefaultSolver );
                mField->initialize();
            }

//------------------------------------------------------------------------------

            ~TwoShellStacks()
            {
                // kernel owns the dof manager and the equation; the mesh
                // owns nodes, blocks, elements, sidesets, facets, edges and
                // the thin-shell records
                delete mKernel ;
                delete mParams ;
                delete mMesh ;
            }

            // owning fixture: a copy or a move would double-free
            TwoShellStacks( const TwoShellStacks & ) = delete ;
            TwoShellStacks & operator=( const TwoShellStacks & ) = delete ;
            TwoShellStacks( TwoShellStacks && ) = delete ;
            TwoShellStacks & operator=( TwoShellStacks && ) = delete ;

//------------------------------------------------------------------------------

            const mesh::Facet *
            linked_facet( const id_t aBlockID, const id_t aElementID )
            {
                return mField->block( aBlockID )->element( aElementID )->facet();
            }
    };
}

//------------------------------------------------------------------------------

// the collector seam, called without a kernel: only the shells that own a
// selected block contribute their facet ids, in shell order; the count and
// the fill must follow one predicate, or the vector is overrun
TEST( ThinShellSelection, CollectorReturnsSelectedShellsOnly )
{
    TwoShellStacks tStacks( false );

    Map< id_t, index_t > tSecondOnly ;
    tSecondOnly[ 3 ] = 0 ;
    tSecondOnly[ 4 ] = 1 ;

    Vector< id_t > tIDs ;
    fem::dofmgr::BlockData::selected_thin_shell_facet_ids( tStacks.mMesh, tSecondOnly, tIDs );
    ASSERT_EQ( tIDs.length(), 1u );
    EXPECT_EQ( tIDs( 0 ), tStacks.mTapeFacet[ 1 ]->id() );

    Map< id_t, index_t > tBoth ;
    tBoth[ 1 ] = 0 ;
    tBoth[ 2 ] = 1 ;
    tBoth[ 3 ] = 2 ;
    tBoth[ 4 ] = 3 ;

    fem::dofmgr::BlockData::selected_thin_shell_facet_ids( tStacks.mMesh, tBoth, tIDs );
    ASSERT_EQ( tIDs.length(), 2u );
    EXPECT_EQ( tIDs( 0 ), tStacks.mTapeFacet[ 0 ]->id() );
    EXPECT_EQ( tIDs( 1 ), tStacks.mTapeFacet[ 1 ]->id() );

    Map< id_t, index_t > tNone ;
    fem::dofmgr::BlockData::selected_thin_shell_facet_ids( tStacks.mMesh, tNone, tIDs );
    EXPECT_EQ( tIDs.length(), 0u );

    // one selected block is enough to select its shell
    Map< id_t, index_t > tOneBlockOfSecond ;
    tOneBlockOfSecond[ 3 ] = 0 ;

    fem::dofmgr::BlockData::selected_thin_shell_facet_ids( tStacks.mMesh, tOneBlockOfSecond, tIDs );
    ASSERT_EQ( tIDs.length(), 1u );
    EXPECT_EQ( tIDs( 0 ), tStacks.mTapeFacet[ 1 ]->id() );
}

//------------------------------------------------------------------------------

// an IWG that selects the layer blocks of the second shell only: the first
// shell must be skipped by the facet-id collector and by the serial linker,
// and the second shell's elements must be linked to the second shell's tape
// facet
TEST( ThinShellSelection, UnselectedShellIsSkipped )
{
    TwoShellStacks tStacks( false );

    ASSERT_NO_THROW( tStacks.create_field() );

    EXPECT_FALSE( tStacks.mField->block_exists( 1 ) );
    EXPECT_FALSE( tStacks.mField->block_exists( 2 ) );
    ASSERT_TRUE( tStacks.mField->block_exists( 3 ) );
    ASSERT_TRUE( tStacks.mField->block_exists( 4 ) );

    EXPECT_EQ( tStacks.linked_facet( 3, 3 ), tStacks.mTapeFacet[ 1 ] );
    EXPECT_EQ( tStacks.linked_facet( 4, 4 ), tStacks.mTapeFacet[ 1 ] );
}

//------------------------------------------------------------------------------

// control: with both shells selected every layer element is linked to its
// own shell's tape facet
TEST( ThinShellSelection, AllShellsSelectedLinkOwnFacet )
{
    TwoShellStacks tStacks( true );

    ASSERT_NO_THROW( tStacks.create_field() );

    EXPECT_EQ( tStacks.linked_facet( 1, 1 ), tStacks.mTapeFacet[ 0 ] );
    EXPECT_EQ( tStacks.linked_facet( 2, 2 ), tStacks.mTapeFacet[ 0 ] );
    EXPECT_EQ( tStacks.linked_facet( 3, 3 ), tStacks.mTapeFacet[ 1 ] );
    EXPECT_EQ( tStacks.linked_facet( 4, 4 ), tStacks.mTapeFacet[ 1 ] );
}

//------------------------------------------------------------------------------
// the buffer layer: a thin-shell layer whose material has no resistivity is
// retagged DomainType::Buffer and takes the scalar-phi role of air. The
// cases below assert that role on a linked group, not a universal rule
//------------------------------------------------------------------------------

namespace
{
    //! the element matrices the magnetic IWG assembles for tri B
    //! ( block 2 ) when block 2 carries aTypeB
    void
    tri_b_matrices( const DomainType aTypeB,
                    Matrix< real > & aK, Matrix< real > & aM, Vector< real > & aF )
    {
        fem::test::TS_TestConductorShell2D tShell( DomainType::Air, aTypeB );

        // the Maxwell factory copies the mesh block's domain type onto the
        // fem group after the field is built; this fixture has no factory
        fem::Block * tBlock = tShell.field()->block( 2 );
        tBlock->set_domain_type( tShell.mesh()->block( 2 )->domain_type() );
        ASSERT_EQ( tBlock->domain_type(), aTypeB );

        fem::IWG_Maxwell * tIWG = tShell.iwg();
        tIWG->link_to_group( tBlock );

        const uint n = tIWG->number_of_dofs_per_element( tBlock );
        Matrix< real > tJ( n, n, 0.0 );
        Vector< real > tRHS( n, 0.0 );
        tIWG->compute_jacobian_and_rhs( tBlock->element( tShell.tri_b()->id() ), tJ, tRHS );

        aK = tIWG->matrices()->K();
        aM = tIWG->matrices()->M();
        aF = tIWG->matrices()->f();
    }
}

TEST( BufferDomain, BufferBlocksTakeThePhiDispatch )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    Matrix< real > tKair, tMair, tKbuf, tMbuf ;
    Vector< real > tFair, tFbuf ;

    tri_b_matrices( DomainType::Air,    tKair, tMair, tFair );
    tri_b_matrices( DomainType::Buffer, tKbuf, tMbuf, tFbuf );

    // the same kernel on the same element: bitwise the same operator
    ASSERT_EQ( tKbuf.n_rows(), tKair.n_rows() );
    ASSERT_EQ( tKbuf.n_cols(), tKair.n_cols() );
    ASSERT_GT( tKair.n_rows(), 0u );

    real tMaxK = 0.0 ;
    for ( uint j = 0; j < tKair.n_cols(); ++j )
    {
        for ( uint i = 0; i < tKair.n_rows(); ++i )
        {
            EXPECT_EQ( tKbuf( i, j ), tKair( i, j ) ) << "K( " << i << ", " << j << " )" ;
            tMaxK = std::max( tMaxK, std::abs( tKair( i, j ) ) );
        }
    }

    real tMaxM = 0.0 ;
    ASSERT_EQ( tMbuf.n_rows(), tMair.n_rows() );
    ASSERT_EQ( tMbuf.n_cols(), tMair.n_cols() );
    for ( uint j = 0; j < tMair.n_cols(); ++j )
    {
        for ( uint i = 0; i < tMair.n_rows(); ++i )
        {
            EXPECT_EQ( tMbuf( i, j ), tMair( i, j ) ) << "M( " << i << ", " << j << " )" ;
            tMaxM = std::max( tMaxM, std::abs( tMair( i, j ) ) );
        }
    }
    EXPECT_GT( std::max( tMaxK, tMaxM ), 0.0 )
        << "all-zero operators would compare equal to anything" ;

    ASSERT_EQ( tFbuf.length(), tFair.length() );
    for ( uint i = 0; i < tFair.length(); ++i )
    {
        EXPECT_EQ( tFbuf( i ), tFair( i ) ) << "f( " << i << " )" ;
    }
}

TEST( BufferDomain, ThermalAcceptsBuffer )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    fem::test::TS_TestConductorShell2D tShell( DomainType::Air, DomainType::Buffer );
    fem::Block * tBlock = tShell.field()->block( 2 );
    tBlock->set_domain_type( tShell.mesh()->block( 2 )->domain_type() );
    ASSERT_EQ( tBlock->domain_type(), DomainType::Buffer );

    // link_to_group refuses, in every build, a group it has no kernel for
    fem::IWG_MaxwellThermal tThermal( ModelDimensionality::TwoD );
    EXPECT_NO_THROW( tThermal.link_to_group( tBlock ) );
}

TEST( BufferDomain, BufferSymmetryLinksAntiSymmetryIsRefused )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    fem::test::TS_TestConductorShell2D tShell( DomainType::Air, DomainType::Buffer );

    // the tape sideset is a LINE2 group; retag it to exercise the
    // magnetic IWG's dispatch for boundary sidesets of a buffer
    fem::SideSet * tSideSet = tShell.field()->sideset( 10 );
    ASSERT_NE( tSideSet, nullptr );
    const DomainType tSaved = tSideSet->domain_type() ;

    tSideSet->set_domain_type( DomainType::BufferSymmetry );
    EXPECT_NO_THROW( tShell.iwg()->link_to_group( tSideSet ) );

    // antisymmetry must be disabled for the phi side, as for air: the
    // refusal must be the deliberate one, not the catch-all for an unknown
    // domain type
    for ( const DomainType tType : { DomainType::BufferAntiSymmetry, DomainType::AirAntiSymmetry } )
    {
        tSideSet->set_domain_type( tType );
        bool tThrown = false ;
        try
        {
            tShell.iwg()->link_to_group( tSideSet );
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true ;
            EXPECT_NE( std::string( aError.what() ).find( "Anti Symmetry" ), std::string::npos )
                << "message was: " << aError.what() ;
        }
        EXPECT_TRUE( tThrown ) << to_string( tType ) ;
    }

    tSideSet->set_domain_type( tSaved );
}

//------------------------------------------------------------------------------
// the thin-shell factory retags a layer block as DomainType::Buffer if and
// only if the layer's material declares no resistivity, and gives each edge
// of that block one source per edge node
//------------------------------------------------------------------------------

namespace
{
    //! a material that declares rho or not, and nothing else
    class FlagMaterial : public Material
    {
        public:
            explicit FlagMaterial( const bool aHaveRho ) :
                Material( MaterialType::NonMetal )
            {
                this->set_have( MaterialProperty::rho, aHaveRho );
            }
    };

    /**
     * two TRI3 on either side of a straight tape LINE2 along y = 0, the
     * tape facet in sideset 10. No cut duplicates: every tape node is its
     * own slave. The layer stack is built by ThinShellFactory::create()
     * from aMaterials, one layer per material
     *
     *        4
     *       / \      B = ( 1, 2, 4 ) above, A = ( 2, 1, 3 ) below
     *      1---2
     *       \ /
     *        3
     */
    class TapeStack
    {
        public:

            std::unique_ptr< Mesh >  mMesh ;
            FlagMaterial             mConductor ;
            FlagMaterial             mStabilizer ;
            FlagMaterial             mInsulator ;
            Map< string, Material * > mMaterials ;
            mesh::ThinShell        * mShell = nullptr ;   // owned by the mesh

            explicit
            TapeStack( const Cell< string > & aMaterials ) :
                mConductor( true ),
                mStabilizer( true ),
                mInsulator( false )
            {
                mMaterials[ "cu" ]     = & mConductor ;
                mMaterials[ "ag" ]     = & mStabilizer ;
                mMaterials[ "buffer" ] = & mInsulator ;
                mMaterials[ "kapton" ] = & mInsulator ;

                mMesh.reset( new Mesh( 2, 0, true ) );

                mesh::ElementFactory tFactory ;

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                tNodes.push( new mesh::Node( 1, 0.0,  0.0 ) );
                tNodes.push( new mesh::Node( 2, 2.0,  0.0 ) );
                tNodes.push( new mesh::Node( 3, 1.0, -1.0 ) );
                tNodes.push( new mesh::Node( 4, 1.0,  1.0 ) );

                mesh::Block * tBlockA = new mesh::Block( 1, 1 );
                tBlockA->set_domain_type( DomainType::Air );
                mesh::Element * tA = tFactory.create_element( ElementType::TRI3, 1 );
                tA->insert_node( tNodes( 1 ), 0 );
                tA->insert_node( tNodes( 0 ), 1 );
                tA->insert_node( tNodes( 2 ), 2 );
                tA->set_block_id( 1 );
                mMesh->elements().push( tA );
                tBlockA->insert_element( tA );
                mMesh->blocks().push( tBlockA );

                mesh::Block * tBlockB = new mesh::Block( 2, 1 );
                tBlockB->set_domain_type( DomainType::Air );
                mesh::Element * tB = tFactory.create_element( ElementType::TRI3, 2 );
                tB->insert_node( tNodes( 0 ), 0 );
                tB->insert_node( tNodes( 1 ), 1 );
                tB->insert_node( tNodes( 3 ), 2 );
                tB->set_block_id( 2 );
                mMesh->elements().push( tB );
                tBlockB->insert_element( tB );
                mMesh->blocks().push( tBlockB );

                mesh::Element * tLine = tFactory.create_element( ElementType::LINE2, 100 );
                tLine->insert_node( tNodes( 0 ), 0 );
                tLine->insert_node( tNodes( 1 ), 1 );
                mesh::Facet * tFacet = new mesh::Facet( tLine );
                tFacet->set_master( tB, 0 );
                tFacet->set_slave( tA, 0 );
                tFacet->set_sideset_id( 10 );

                mesh::SideSet * tTape = new mesh::SideSet( 10, 1 );
                tTape->insert_facet( tFacet );
                mMesh->sidesets().push( tTape );

                mMesh->finalize();

                Cell< mesh::Node * > tMasterNodes = { tNodes( 0 ), tNodes( 1 ) };
                Cell< mesh::Node * > tSlaveNodes  = { tNodes( 0 ), tNodes( 1 ) };

                Protoshell tProto( 1 );
                tProto.sidesets() = { 10 };
                tProto.materials() = aMaterials ;
                tProto.thicknesses().set_size( aMaterials.size(), 0.01 );

                mesh::ThinShellFactory tShellFactory( mMesh.get(), tMasterNodes, tSlaveNodes, & mMaterials, true );
                mShell = tShellFactory.create( & tProto );
            }
    };
}

TEST( ThinShellSelection, RhoLessLayerIsTaggedBuffer )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // the rule is the missing resistivity, not the material's name
    TapeStack tStack( { "cu", "kapton", "cu" } );
    ASSERT_NE( tStack.mShell, nullptr );
    ASSERT_EQ( tStack.mShell->blocks().size(), 3u );

    const DomainType tExpected[ 3 ] = { DomainType::ThinShell, DomainType::Buffer, DomainType::ThinShell };
    for ( uint b = 0; b < 3; ++b )
    {
        EXPECT_EQ( tStack.mShell->blocks()( b )->domain_type(), tExpected[ b ] )
            << "layer " << b << ": " << to_string( tStack.mShell->blocks()( b )->domain_type() ) ;
    }
}

TEST( ThinShellSelection, BufferEdgesCarryOneSourcePerNode )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    TapeStack tStack( { "cu", "buffer", "cu" } );
    mesh::Block * tBuffer = tStack.mShell->blocks()( 1 );
    ASSERT_EQ( tBuffer->domain_type(), DomainType::Buffer );
    ASSERT_GT( tBuffer->number_of_elements(), 0u );

    uint tNumEdges = 0 ;
    for ( mesh::Element * tElement : tBuffer->elements() )
    {
        for ( uint e = 0; e < tElement->number_of_edges(); ++e )
        {
            mesh::Edge * tEdge = tElement->edge( e );
            ++tNumEdges ;

            ASSERT_EQ( tEdge->number_of_sources(), tEdge->number_of_nodes() )
                << "edge " << tEdge->id() ;
            for ( uint k = 0; k < tEdge->number_of_nodes(); ++k )
            {
                EXPECT_EQ( tEdge->source( k ), tEdge->node( k ) )
                    << "edge " << tEdge->id() << ", source " << k ;
            }
        }
    }
    EXPECT_GT( tNumEdges, 0u );
}

TEST( ThinShellSelection, NoGhostFacetAcrossABuffer )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // controls: two and three resistive layers do get their ghost facets
    {
        TapeStack tStack( { "cu", "ag" } );
        EXPECT_NE( tStack.mShell->ghost_id(), gNoID );
        EXPECT_GT( tStack.mShell->ghost_facets().size(), 0u );
    }
    {
        TapeStack tStack( { "cu", "ag", "cu" } );
        EXPECT_NE( tStack.mShell->ghost_id(), gNoID );
        EXPECT_GT( tStack.mShell->ghost_facets().size(), 0u );
    }

    // with a buffer between them, neither interface is a resistive pair
    {
        TapeStack tStack( { "cu", "buffer", "ag" } );
        EXPECT_EQ( tStack.mShell->ghost_id(), gNoID );
    }
}
