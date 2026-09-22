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
