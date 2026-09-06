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

#ifndef BELFEM_CL_TS_TESTSTACK_HPP
#define BELFEM_CL_TS_TESTSTACK_HPP

#include <cmath>

#include "typedefs.hpp"
#include "globals.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Facet.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_ThinShell.hpp"
#include "cl_Element_Factory.hpp"
#include "meshtools.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_MeshChecker.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "en_Maxwell_Formulations.hpp"
#include "en_SolverEnums.hpp"

namespace belfem
{
    namespace fem
    {
        namespace test
        {
//------------------------------------------------------------------------------

            /**
             * Code-built N-layer 2-D thin-shell stack for edge-function tests.
             *
             * Geometry: N QUAD4TS layers stacked through-thickness, sharing the
             * interface node rows (the greg2 identical-material configuration:
             * shared edges, strong continuity, no ghost facets). The tape runs
             * along the unit vector (cos a, sin a); the stacking direction is
             * its +90-degree normal. aAngle != 0 exercises the non-axis-aligned
             * paths (the historical bug class assumed axis alignment).
             *
             * Node ordering per layer follows the ThinShellFactory convention
             * (cl_ThinShellFactory.cpp:1513): element node 3 sits above facet
             * node 0, element node 2 above facet node 1.
             *
             * Edge objects are created per interface ROW (row r shared by layer
             * r-1 top slot and layer r bottom slot). By default the edge node
             * order matches the canonical get_nodes_of_edge() order, so all
             * edge directions are positive; aFlipRow reverses one row's edge
             * node order to exercise the direction logic.
             *
             * The fem::Elements are built with the aura constructor against a
             * minimal real Kernel/DofManagerBase/Block chain, so the designed
             * thickness data path ( parent()->parent()->mesh()->block()->
             * thickness() ) is the one under test — no seams.
             */
            class TS_TestStack2D
            {
                Mesh                   * mMesh ;
                KernelParameters       * mParams ;
                Kernel                 * mKernel ;
                DofManagerBase         * mDofBase ;
                Block                  * mBlock ;

                Cell< fem::Element * >   mElements ;   // one per layer, owned
                Cell< mesh::Facet * >    mFacets ;     // bottom facet per layer, owned
                Cell< mesh::Edge * >     mEdges ;      // one per interface row, owned
                                                       // (fixture-owned: the MeshChecker
                                                       // inside the Kernel ctor rejects
                                                       // meshes that already carry edges)

                const uint mNumLayers ;
                const real mLength ;
                const real mThickness ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                TS_TestStack2D(
                        const uint aNumLayers,
                        const real aLength     = 2.0,
                        const real aThickness  = 0.1,
                        const real aAngle      = 0.0,
                        const int  aFlipRow    = -1 );

                ~TS_TestStack2D();

                uint
                num_layers() const { return mNumLayers ; }

                real
                length() const { return mLength ; }

                real
                thickness() const { return mThickness ; }

                fem::Element *
                element( const uint aLayer ) { return mElements( aLayer ); }

                Mesh *
                mesh() { return mMesh ; }
            };

//------------------------------------------------------------------------------

            inline
            TS_TestStack2D::TS_TestStack2D(
                    const uint aNumLayers,
                    const real aLength,
                    const real aThickness,
                    const real aAngle,
                    const int  aFlipRow ) :
                    mNumLayers( aNumLayers ),
                    mLength( aLength ),
                    mThickness( aThickness )
            {
                // tape direction and stacking normal
                const real tCa = std::cos( aAngle );
                const real tSa = std::sin( aAngle );

                mMesh = new Mesh( 2, 0, false );

                mesh::ElementFactory tFactory ;

                // node rows 0 .. N (row r = interface r), 2 nodes per row
                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                id_t tID = 1 ;
                for ( uint r = 0; r <= aNumLayers; ++r )
                {
                    const real tOffX = -tSa * aThickness * r ;
                    const real tOffY =  tCa * aThickness * r ;
                    tNodes.push( new mesh::Node( tID++, tOffX, tOffY ) );
                    tNodes.push( new mesh::Node( tID++,
                            tOffX + tCa * aLength,
                            tOffY + tSa * aLength ) );
                }

                // block with the layer elements
                mesh::Block * tBlock = new mesh::Block( 1, aNumLayers );
                tBlock->set_thickness( aThickness );

                Cell< mesh::Element * > & tElements = mMesh->elements();

                for ( uint l = 0; l < aNumLayers; ++l )
                {
                    mesh::Element * tElement =
                            tFactory.create_element( ElementType::QUAD4TS, l + 1 );

                    // bottom row l, top row l+1
                    tElement->insert_node( tNodes( 2 * l ),     0 );
                    tElement->insert_node( tNodes( 2 * l + 1 ), 1 );
                    tElement->insert_node( tNodes( 2 * l + 3 ), 2 );
                    tElement->insert_node( tNodes( 2 * l + 2 ), 3 );

                    tElement->set_block_id( 1 );
                    tElements.push( tElement );
                    tBlock->insert_element( tElement );
                }

                mMesh->blocks().push( tBlock );

                mMesh->finalize();

                // minimal real kernel chain for the thickness data path
                mParams  = new KernelParameters( mMesh );
                mKernel  = new Kernel( mParams );
                mDofBase = new DofManagerBase( DofManagerType::UNDEFINED, mKernel );
                mBlock   = new Block( mDofBase, ElementType::QUAD4TS );

                // one shared Edge object per interface row, AFTER the kernel
                // ( the MeshChecker refuses pre-existing edges ). Canonical
                // node order of the QUAD4TS bottom edge is { node0, node1 }
                // and of the top edge { node3, node2 } — both run along
                // +tape, so the row edge is { left, right } unless flipped.
                mEdges.set_size( aNumLayers + 1, nullptr );
                for ( uint r = 0; r <= aNumLayers; ++r )
                {
                    mesh::Edge * tEdge = new mesh::Edge();
                    tEdge->set_id( r + 1 );
                    tEdge->allocate_node_container( 2 );
                    mesh::Node * tLeft  = tNodes( 2 * r );
                    mesh::Node * tRight = tNodes( 2 * r + 1 );
                    if ( ( int ) r == aFlipRow )
                    {
                        tEdge->insert_node( tRight, 0 );
                        tEdge->insert_node( tLeft, 1 );
                    }
                    else
                    {
                        tEdge->insert_node( tLeft, 0 );
                        tEdge->insert_node( tRight, 1 );
                    }
                    mEdges( r ) = tEdge ;
                }

                for ( uint l = 0; l < aNumLayers; ++l )
                {
                    mesh::Element * tElement = mMesh->elements()( l );
                    tElement->allocate_edge_container();
                    tElement->insert_edge( mEdges( l ),     0 );   // bottom
                    tElement->insert_edge( mEdges( l + 1 ), 1 );   // top
                }

                // bottom facet (LINE2) per layer + fem element
                mFacets.set_size( aNumLayers, nullptr );
                mElements.set_size( aNumLayers, nullptr );
                for ( uint l = 0; l < aNumLayers; ++l )
                {
                    mesh::Element * tLine =
                            tFactory.create_element( ElementType::LINE2,
                                                     1000 + l + 1 );
                    tLine->insert_node( mMesh->elements()( l )->node( 0 ), 0 );
                    tLine->insert_node( mMesh->elements()( l )->node( 1 ), 1 );

                    mFacets( l ) = new mesh::Facet( tLine );

                    fem::Element * tFem =
                            new fem::Element( mBlock, mMesh->elements()( l ) );
                    tFem->set_facet( mFacets( l ) );

                    mElements( l ) = tFem ;
                }
            }

//------------------------------------------------------------------------------

            inline
            TS_TestStack2D::~TS_TestStack2D()
            {
                for ( fem::Element * tElement : mElements )
                {
                    delete tElement ;
                }
                for ( mesh::Facet * tFacet : mFacets )
                {
                    delete tFacet ;
                }
                for ( mesh::Edge * tEdge : mEdges )
                {
                    delete tEdge ;
                }
                delete mBlock ;
                delete mDofBase ;
                delete mKernel ;
                delete mParams ;
                delete mMesh ;
            }

//------------------------------------------------------------------------------

            /**
             * Single 3-D thin-shell prism ( PENTA6TS or HEX8TS ) with its
             * mid-surface facet, canonical edges, and the same minimal real
             * kernel chain. Edges are built from get_nodes_of_edge() in
             * canonical order, so all directions are positive by
             * construction; aFlipEdge reverses one edge.
             *
             * The default mid-surface is a flat GENERAL quad in the z = 0
             * plane ( NOT a parallelogram — the HEX8TS nablas vary per
             * point on it ). aCorners ( 2 x nFacetNodes, x row 0, y row 1,
             * z = 0 plane, counterclockwise from +z ) overrides the corner
             * coordinates, e.g. with a rotated rectangle when a test needs
             * a constant, orthogonal in-plane frame.
             */
            class TS_TestPrism
            {
                Mesh                   * mMesh ;
                KernelParameters       * mParams ;
                Kernel                 * mKernel ;
                DofManagerBase         * mDofBase ;
                Block                  * mBlock ;

                fem::Element           * mElement = nullptr ;
                mesh::Facet            * mFacet   = nullptr ;
                Cell< mesh::Edge * >     mEdges ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                TS_TestPrism(
                        const ElementType      aType,
                        const real             aThickness = 0.1,
                        const int              aFlipEdge  = -1,
                        const Matrix< real > * aCorners   = nullptr );

                ~TS_TestPrism();

                fem::Element *
                element() { return mElement ; }

                uint
                num_edges() const { return mEdges.size(); }

                mesh::Edge *
                edge( const uint aIndex ) { return mEdges( aIndex ); }
            };

//------------------------------------------------------------------------------

            inline
            TS_TestPrism::TS_TestPrism(
                    const ElementType      aType,
                    const real             aThickness,
                    const int              aFlipEdge,
                    const Matrix< real > * aCorners )
            {
                BELFEM_ERROR( aType == ElementType::PENTA6TS
                           || aType == ElementType::HEX8TS,
                    "TS_TestPrism: unsupported element type" );

                const bool tIsPenta = ( aType == ElementType::PENTA6TS );
                const uint tNumFacetNodes = tIsPenta ? 3 : 4 ;

                mMesh = new Mesh( 3, 0, false );

                mesh::ElementFactory tFactory ;

                // mid-surface corner coordinates ( z = 0 plane ),
                // counterclockwise as seen from +z; the default quad is
                // deliberately NOT a parallelogram
                real tX[ 4 ] = { 0.0, 1.3, 0.9, 0.0 };
                real tY[ 4 ] = { 0.0, 0.1, 1.1, 0.8 };

                if ( aCorners != nullptr )
                {
                    for ( uint k = 0; k < tNumFacetNodes; ++k )
                    {
                        tX[ k ] = ( *aCorners )( 0, k );
                        tY[ k ] = ( *aCorners )( 1, k );
                    }
                }

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                id_t tID = 1 ;
                for ( uint tLayer = 0; tLayer < 2; ++tLayer )
                {
                    const real tZ = tLayer == 0 ? 0.0 : aThickness ;
                    for ( uint k = 0; k < tNumFacetNodes; ++k )
                    {
                        tNodes.push( new mesh::Node( tID++,
                                tX[ k ], tY[ k ], tZ ) );
                    }
                }

                mesh::Block * tBlock = new mesh::Block( 1, 1 );
                tBlock->set_thickness( aThickness );

                mesh::Element * tElement = tFactory.create_element( aType, 1 );
                for ( uint k = 0; k < 2 * tNumFacetNodes; ++k )
                {
                    tElement->insert_node( tNodes( k ), k );
                }
                tElement->set_block_id( 1 );
                mMesh->elements().push( tElement );
                tBlock->insert_element( tElement );
                mMesh->blocks().push( tBlock );

                mMesh->finalize();

                mParams  = new KernelParameters( mMesh );
                mKernel  = new Kernel( mParams );
                mDofBase = new DofManagerBase( DofManagerType::UNDEFINED, mKernel );
                mBlock   = new Block( mDofBase, aType );

                // canonical edges from the element's own edge tables
                const uint tNumEdges = tElement->number_of_edges();
                mEdges.set_size( tNumEdges, nullptr );
                Cell< mesh::Node * > tEdgeNodes ;
                tElement->allocate_edge_container();
                for ( uint e = 0; e < tNumEdges; ++e )
                {
                    tElement->get_nodes_of_edge( e, tEdgeNodes );
                    mesh::Edge * tEdge = new mesh::Edge();
                    tEdge->set_id( e + 1 );
                    tEdge->allocate_node_container( 2 );
                    if ( ( int ) e == aFlipEdge )
                    {
                        tEdge->insert_node( tEdgeNodes( 1 ), 0 );
                        tEdge->insert_node( tEdgeNodes( 0 ), 1 );
                    }
                    else
                    {
                        tEdge->insert_node( tEdgeNodes( 0 ), 0 );
                        tEdge->insert_node( tEdgeNodes( 1 ), 1 );
                    }
                    tEdge->set_id( e + 1 );
                    mEdges( e ) = tEdge ;
                    tElement->insert_edge( tEdge, e );
                }

                // mid-surface facet element
                mesh::Element * tSurf = tFactory.create_element(
                        tIsPenta ? ElementType::TRI3 : ElementType::QUAD4, 1001 );
                for ( uint k = 0; k < tNumFacetNodes; ++k )
                {
                    tSurf->insert_node( tNodes( k ), k );
                }
                mFacet = new mesh::Facet( tSurf );

                mElement = new fem::Element( mBlock, tElement );
                mElement->set_facet( mFacet );
            }

//------------------------------------------------------------------------------

            inline
            TS_TestPrism::~TS_TestPrism()
            {
                delete mElement ;
                delete mFacet ;
                for ( mesh::Edge * tEdge : mEdges )
                {
                    delete tEdge ;
                }
                delete mBlock ;
                delete mDofBase ;
                delete mKernel ;
                delete mParams ;
                delete mMesh ;
            }

//------------------------------------------------------------------------------

            /**
             * HEX8TB side-connector wall with its designed data paths:
             * block 1 = one HEX8TS layer element ( thickness ), block 2 =
             * the wall ( block thickness carries the WIDTH ), and the
             * READ-ONLY recovery facet with id = wall id + 1 whose master
             * is the layer element ( cl_EF_HEX8TB.cpp:60-73 ).
             * Wall cuboid: length L along x, width along y, layer thickness
             * along z; nodes 0-3 bottom, 4-7 top.
             */
            class TS_TestWall
            {
                Mesh                   * mMesh ;
                KernelParameters       * mParams ;
                Kernel                 * mKernel ;
                DofManagerBase         * mDofBase ;
                Block                  * mBlock ;

                fem::Element           * mElement = nullptr ;
                Cell< mesh::Edge * >     mEdges ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                TS_TestWall(
                        const real aLength    = 2.0,
                        const real aWidth     = 0.05,
                        const real aThickness = 0.1,
                        const int  aFlipEdge  = -1 );

                ~TS_TestWall();

                fem::Element *
                element() { return mElement ; }
            };

//------------------------------------------------------------------------------

            inline
            TS_TestWall::TS_TestWall(
                    const real aLength,
                    const real aWidth,
                    const real aThickness,
                    const int  aFlipEdge )
            {
                mMesh = new Mesh( 3, 0, false );
                mesh::ElementFactory tFactory ;

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                id_t tID = 1 ;

                // layer element ( HEX8TS ) sits next to the wall in -y
                const real tYL0 = -3.0 * aWidth ;
                for ( uint tLayer = 0; tLayer < 2; ++tLayer )
                {
                    const real tZ = tLayer == 0 ? 0.0 : aThickness ;
                    tNodes.push( new mesh::Node( tID++, 0.0,     tYL0,   tZ ) );
                    tNodes.push( new mesh::Node( tID++, aLength, tYL0,   tZ ) );
                    tNodes.push( new mesh::Node( tID++, aLength, 0.0,    tZ ) );
                    tNodes.push( new mesh::Node( tID++, 0.0,     0.0,    tZ ) );
                }
                // wall nodes ( 0..3 bottom z=0, 4..7 top ), width in +y
                for ( uint tLayer = 0; tLayer < 2; ++tLayer )
                {
                    const real tZ = tLayer == 0 ? 0.0 : aThickness ;
                    tNodes.push( new mesh::Node( tID++, 0.0,     0.0,    tZ ) );
                    tNodes.push( new mesh::Node( tID++, aLength, 0.0,    tZ ) );
                    tNodes.push( new mesh::Node( tID++, aLength, aWidth, tZ ) );
                    tNodes.push( new mesh::Node( tID++, 0.0,     aWidth, tZ ) );
                }

                // block 1: the layer ( thickness )
                mesh::Block * tLayerBlock = new mesh::Block( 1, 1 );
                tLayerBlock->set_thickness( aThickness );
                mesh::Element * tLayer =
                        tFactory.create_element( ElementType::HEX8TS, 1 );
                for ( uint k = 0; k < 8; ++k )
                {
                    tLayer->insert_node( tNodes( k ), k );
                }
                tLayer->set_block_id( 1 );
                mMesh->elements().push( tLayer );
                tLayerBlock->insert_element( tLayer );
                mMesh->blocks().push( tLayerBlock );

                // block 2: the wall ( block thickness = WIDTH )
                mesh::Block * tWallBlock = new mesh::Block( 2, 1 );
                tWallBlock->set_thickness( aWidth );
                mesh::Element * tWall =
                        tFactory.create_element( ElementType::HEX8TB, 2 );
                for ( uint k = 0; k < 8; ++k )
                {
                    tWall->insert_node( tNodes( 8 + k ), k );
                }
                tWall->set_block_id( 2 );
                mMesh->elements().push( tWall );
                tWallBlock->insert_element( tWall );
                mMesh->blocks().push( tWallBlock );

                // recovery facet: id = wall id + 1, master = layer element
                mesh::Element * tSurf =
                        tFactory.create_element( ElementType::QUAD4, 3 );
                mesh::Facet * tFacet = new mesh::Facet( tSurf );
                tFacet->set_master( tLayer, 0 );
                mMesh->facets().push( tFacet );

                mMesh->finalize();

                mParams  = new KernelParameters( mMesh );
                mKernel  = new Kernel( mParams );
                mDofBase = new DofManagerBase( DofManagerType::UNDEFINED, mKernel );
                mBlock   = new Block( mDofBase, ElementType::HEX8TB );

                // the four longitudinal edges, canonical order
                const uint tNumEdges = tWall->number_of_edges();
                mEdges.set_size( tNumEdges, nullptr );
                Cell< mesh::Node * > tEdgeNodes ;
                tWall->allocate_edge_container();
                for ( uint e = 0; e < tNumEdges; ++e )
                {
                    tWall->get_nodes_of_edge( e, tEdgeNodes );
                    mesh::Edge * tEdge = new mesh::Edge();
                    tEdge->set_id( e + 1 );
                    tEdge->allocate_node_container( 2 );
                    if ( ( int ) e == aFlipEdge )
                    {
                        tEdge->insert_node( tEdgeNodes( 1 ), 0 );
                        tEdge->insert_node( tEdgeNodes( 0 ), 1 );
                    }
                    else
                    {
                        tEdge->insert_node( tEdgeNodes( 0 ), 0 );
                        tEdge->insert_node( tEdgeNodes( 1 ), 1 );
                    }
                    mEdges( e ) = tEdge ;
                    tWall->insert_edge( tEdge, e );
                }

                mElement = new fem::Element( mBlock, tWall );
            }

//------------------------------------------------------------------------------

            inline
            TS_TestWall::~TS_TestWall()
            {
                delete mElement ;
                for ( mesh::Edge * tEdge : mEdges )
                {
                    delete tEdge ;
                }
                delete mBlock ;
                delete mDofBase ;
                delete mKernel ;
                delete mParams ;
                delete mMesh ;
            }

//------------------------------------------------------------------------------

            /**
             * Two stacked PENTA6TS layer blocks with a ghost facet between
             * them, over the full production Kernel/DofManager/IWG_Maxwell
             * chain — the Calculator-level fixture for the h_ghost contract
             * Mirrors the post-ThinShellFactory mesh state:
             * shared middle node layer, DUPLICATED interface edges ( master
             * top and slave bottom are distinct objects, as after
             * cl_ThinShellFactory.cpp create() with hasDuplicates ), ghost
             * sideset ( DomainType::Ghost ), tape sideset carrying the
             * mid-surface reference facet the layer edge functions read
             * their in-plane frame from, and a mesh::ThinShell container
             * linking the layer blocks to it ( the ThinShell ctor hides
             * both sidesets from mesh output — harmless here ).
             *
             * aFlipEdge ( 0..2 ): reverses the node order of local edge
             * aFlipEdge in ALL edge rows of BOTH prisms. This is the only
             * legal flip: ThinShellFactory extrudes every row from the same
             * mid-surface edges, so direction patterns are identical per
             * local index across rows and elements ( the ghost element ctor
             * asserts the master/slave direction bitsets match,
             * cl_FEM_Element.cpp, and the kernel algebra additionally needs
             * co-located duplicates to agree ).
             *
             * The ctor saves gTbulk and sets 77 K so h_ghost skips its
             * temperature branch; the dtor restores the saved value.
             */
            class TS_TestGhostStack
            {
                //! saves gTbulk on construction, pins 77 K, restores on
                //! destruction. A member object (declared first) so the
                //! restore also runs when the enclosing ctor throws
                //! mid-construction — later tests must not inherit 77 K
                struct TbulkGuard
                {
                    const real mSaved ;
                    TbulkGuard() : mSaved( gTbulk ) { gTbulk = 77.0 ; }
                    ~TbulkGuard() { gTbulk = mSaved ; }
                };

                TbulkGuard               mTbulkGuard ;

                Mesh                   * mMesh ;
                KernelParameters       * mParams ;
                Kernel                 * mKernel ;
                DofManager             * mField ;      // owned by the kernel
                IWG_Maxwell            * mIWG ;        // owned by the kernel
                SideSet                * mGhostGroup ; // owned by the field

                mesh::Element          * mMaster = nullptr ;
                mesh::Element          * mSlave  = nullptr ;

                const real               mThickness ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                TS_TestGhostStack(
                        const real aThickness = 0.1,
                        const int  aFlipEdge  = -1 );

                // owning raw pointers: neither copyable nor movable
                TS_TestGhostStack( const TS_TestGhostStack & ) = delete ;
                TS_TestGhostStack &
                operator=( const TS_TestGhostStack & ) = delete ;

                ~TS_TestGhostStack();

                Calculator *
                calculator() { return mGhostGroup->calculator() ; }

                fem::Element *
                ghost_element() { return mGhostGroup->elements()( 0 ) ; }

                mesh::Element *
                master() { return mMaster ; }

                mesh::Element *
                slave() { return mSlave ; }

                IWG_Maxwell *
                iwg() { return mIWG ; }

                Mesh *
                mesh() { return mMesh ; }

                /**
                 * rewrite element_rho for both layer elements; h_ghost
                 * re-reads the field on every call, so no re-link is needed
                 */
                void
                set_rho( const real aRho )
                {
                    Vector< real > & tRho = mMesh->field_data( "element_rho" );
                    tRho( mMaster->index() ) = aRho ;
                    tRho( mSlave->index() )  = aRho ;
                }
            };

//------------------------------------------------------------------------------

            inline
            TS_TestGhostStack::TS_TestGhostStack(
                    const real aThickness,
                    const int  aFlipEdge ) :
                    mThickness( aThickness )
            {
                // mTbulkGuard has already pinned gTbulk = 77 K ( h_ghost
                // falls back to the nodal "T" field when gTbulk is NaN )

                // third argument false: the MeshChecker below must run
                // before edges exist, and connectivity computation inside
                // finalize() would create them
                mMesh = new Mesh( 3, 0, false );

                mesh::ElementFactory tFactory ;

                // triangle footprint, deliberately not equilateral
                const real tX[ 3 ] = { 0.0, 1.3, 0.9 };
                const real tY[ 3 ] = { 0.0, 0.1, 1.1 };

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                id_t tID = 1 ;
                for ( uint r = 0; r < 3; ++r )
                {
                    for ( uint k = 0; k < 3; ++k )
                    {
                        tNodes.push( new mesh::Node( tID++,
                                tX[ k ], tY[ k ], aThickness * r ) );
                    }
                }

                Cell< mesh::Element * > & tElements = mMesh->elements();

                // lower prism, block 1: node rows 0-1
                mesh::Block * tBlock1 = new mesh::Block( 1, 1 );
                tBlock1->set_thickness( aThickness );
                tBlock1->set_domain_type( DomainType::ThinShell );
                mMaster = tFactory.create_element( ElementType::PENTA6TS, 1 );
                for ( uint k = 0; k < 6; ++k )
                {
                    mMaster->insert_node( tNodes( k ), k );
                }
                mMaster->set_block_id( 1 );
                tElements.push( mMaster );
                tBlock1->insert_element( mMaster );
                mMesh->blocks().push( tBlock1 );

                // upper prism, block 2: node rows 1-2
                mesh::Block * tBlock2 = new mesh::Block( 2, 1 );
                tBlock2->set_thickness( aThickness );
                tBlock2->set_domain_type( DomainType::ThinShell );
                mSlave = tFactory.create_element( ElementType::PENTA6TS, 2 );
                for ( uint k = 0; k < 6; ++k )
                {
                    mSlave->insert_node( tNodes( k + 3 ), k );
                }
                mSlave->set_block_id( 2 );
                tElements.push( mSlave );
                tBlock2->insert_element( mSlave );
                mMesh->blocks().push( tBlock2 );

                // ghost facet: TRI3 on the middle node row, master = lower
                // prism's top face, slave = upper prism's bottom face with
                // orientation 1 ( cl_ThinShellFactory.cpp convention )
                mesh::Element * tSurf =
                        tFactory.create_element( ElementType::TRI3, 100 );
                for ( uint k = 0; k < 3; ++k )
                {
                    tSurf->insert_node( tNodes( k + 3 ), k );
                }
                mesh::Facet * tGhost = new mesh::Facet( tSurf );
                tGhost->set_master( mMaster,
                        mesh::top_facet_index( ElementType::PENTA6TS ) );
                tGhost->set_slave( mSlave,
                        mesh::bottom_facet_index( ElementType::PENTA6TS ), 1 );
                tGhost->set_sideset_id( 10 );

                mesh::SideSet * tGhostSet = new mesh::SideSet( 10, 1 );
                tGhostSet->set_domain_type( DomainType::Ghost );
                tGhostSet->insert_facet( tGhost );
                mMesh->sidesets().push( tGhostSet );

                // tape ( mid-surface ) sideset: carries the reference facet
                // the layer edge functions read their in-plane frame from;
                // deliberately NOT in the IWG sideset list — only
                // BlockData::link_thin_shell_facets consumes it
                mesh::Element * tMidSurf =
                        tFactory.create_element( ElementType::TRI3, 101 );
                for ( uint k = 0; k < 3; ++k )
                {
                    tMidSurf->insert_node( tNodes( k + 3 ), k );
                }
                mesh::Facet * tMid = new mesh::Facet( tMidSurf );
                tMid->set_sideset_id( 11 );
                mesh::SideSet * tTapeSet = new mesh::SideSet( 11, 1 );
                tTapeSet->set_domain_type( DomainType::ThinShell );
                tTapeSet->insert_facet( tMid );
                mMesh->sidesets().push( tTapeSet );

                // thin-shell container linking layer blocks to the
                // reference facets ( mesh owns and deletes it )
                mesh::ThinShell * tShell =
                        new mesh::ThinShell( tTapeSet, tGhostSet );
                tShell->blocks().push( tBlock1 );
                tShell->blocks().push( tBlock2 );
                {
                    Vector< real > tThicknesses( 2, aThickness );
                    tShell->set_thicknesses( tThicknesses );
                }
                mMesh->thin_shells().push( tShell );

                mMesh->finalize();

                // reorientation check must run BEFORE edges exist; it sets
                // the mesh checker flag so the Kernel ctor skips its own run
                {
                    MeshChecker tCheck( mMesh );
                }

                // per-element edge objects in canonical get_nodes_of_edge()
                // order; the interface edges are deliberately DUPLICATED
                // ( master ids 1-6, slave ids 7-12 ). aFlipEdge reverses
                // local edge aFlipEdge in every row of both elements
                {
                    id_t tEdgeID = 1 ;
                    Cell< mesh::Node * > tEdgeNodes ;
                    mesh::Element * tPrisms[ 2 ] = { mMaster, mSlave };
                    for ( uint p = 0; p < 2; ++p )
                    {
                        mesh::Element * tElement = tPrisms[ p ];
                        const uint tNumEdges = tElement->number_of_edges();
                        tElement->allocate_edge_container();
                        for ( uint e = 0; e < tNumEdges; ++e )
                        {
                            tElement->get_nodes_of_edge( e, tEdgeNodes );
                            mesh::Edge * tEdge = new mesh::Edge();
                            tEdge->set_id( tEdgeID++ );
                            tEdge->allocate_node_container( 2 );
                            if ( ( int )( e % 3 ) == aFlipEdge )
                            {
                                tEdge->insert_node( tEdgeNodes( 1 ), 0 );
                                tEdge->insert_node( tEdgeNodes( 0 ), 1 );
                            }
                            else
                            {
                                tEdge->insert_node( tEdgeNodes( 0 ), 0 );
                                tEdge->insert_node( tEdgeNodes( 1 ), 1 );
                            }
                            tEdge->set_owner( 0 );
                            tElement->insert_edge( tEdge, e );
                            mMesh->edges().push( tEdge );
                        }
                    }
                    mMesh->finalize_edges();
                }

                mParams = new KernelParameters( mMesh );
                mKernel = new Kernel( mParams );

                mIWG = new IWG_Maxwell(
                        maxwell::Formulation::HPhi,
                        ModelDimensionality::ThreeD,
                        false, false );

                Vector< id_t > tBlockIDs( 2 );
                tBlockIDs( 0 ) = 1 ;
                tBlockIDs( 1 ) = 2 ;
                Cell< DomainType > tBlockTypes( 2, DomainType::ThinShell );
                mIWG->set_blocks( tBlockIDs, tBlockTypes );

                Vector< id_t > tSideSetIDs( 1 );
                tSideSetIDs( 0 ) = 10 ;
                Cell< DomainType > tSideSetTypes( 1, DomainType::Ghost );
                mIWG->set_sidesets( tSideSetIDs, tSideSetTypes );

                mKernel->add_equation( mIWG );

                mField = mKernel->create_field( mIWG );
                mField->create_fields( mField->iwg() );
                this->set_rho( 1.0e-8 );

                mField->set_solver( gDefaultSolver );
                mField->initialize();

                mGhostGroup = mField->sideset( 10 );

                // gate: this call used to be skipped because on a
                // kernel without a Controller it read an uninitialized
                // mController ( heap-lottery jump, battery-order dependent ).
                // The read is guarded now — the link must succeed on a
                // controller-less kernel and must leave have_thermal()
                // false ( asserted by the consuming test )
                mIWG->link_to_group( mGhostGroup );
            }

//------------------------------------------------------------------------------

            inline
            TS_TestGhostStack::~TS_TestGhostStack()
            {
                // kernel owns the dof manager and the equation; the mesh
                // owns nodes, blocks, elements, sidesets, facets, edges,
                // and the thin-shell record ( test_DofSeeding order )
                delete mKernel ;
                delete mParams ;
                delete mMesh ;

                // gTbulk restored by mTbulkGuard's destructor
            }


//------------------------------------------------------------------------------

            /**
             * One PENTA6TS layer between two TET4 volumes, over the full
             * production Kernel/DofManager/IWG_Maxwell chain — the fixture
             * for the thin-shell normal-field recovery ( compute_hn and
             * compute_h_trace ) next to an h-conductor.
             *
             * Geometry: triangle T in the z = 0 plane ( CCW from +z ) is the
             * tape facet. Block 1 = tet A below T ( apex at z = -aHeight ),
             * block 2 = tet B above T ( apex at z = +aHeight ), block 3 = the
             * layer prism extruded from T to z = aThickness ( it overlaps tet
             * B geometrically, exactly as a ThinShellFactory layer overlaps
             * the volume mesh it was extruded into ). The tape sideset ( id
             * 10, DomainType::ThinShell, GeometryOnly in the dof manager )
             * holds one facet with master = tet A ( local face 3, outward
             * normal +z ) and slave = tet B ( local face 3 ), set explicitly:
             * every master/slave kind combination is reachable without the
             * factory's domain-type priority rule. aTypeA / aTypeB are
             * DomainType::Conductor or DomainType::Air.
             *
             * Edge objects exist on every Conductor tet and on the prism; an
             * Air tet has none, as in production. aFlipEdge ( 0..5 )
             * reverses the node order of that local edge on both tets ( where
             * they carry edges ) so the dof sign convention is exercised.
             *
             * No dof values are solved: the tests seed the edge_h and phi
             * mesh fields directly, which is the live storage compute_hn
             * reads ( Calculator::q() reads the same fields ).
             */
            class TS_TestConductorShell
            {
                Mesh                   * mMesh ;
                KernelParameters       * mParams ;
                Kernel                 * mKernel ;
                DofManager             * mField ;   // owned by the kernel
                IWG_Maxwell            * mIWG ;     // owned by the kernel

                mesh::Element          * mTetA  = nullptr ;
                mesh::Element          * mTetB  = nullptr ;
                mesh::Element          * mPrism = nullptr ;
                mesh::Facet            * mFacet = nullptr ;

                const real               mHeight ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                TS_TestConductorShell(
                        const DomainType aTypeA,
                        const DomainType aTypeB,
                        const int        aFlipEdge  = -1,
                        const real       aThickness = 0.1,
                        const real       aHeight    = 0.5 );

                TS_TestConductorShell( const TS_TestConductorShell & ) = delete ;
                TS_TestConductorShell &
                operator=( const TS_TestConductorShell & ) = delete ;

                ~TS_TestConductorShell();

                Mesh *
                mesh() { return mMesh ; }

                DofManager *
                field() { return mField ; }

                mesh::Element *
                tet_a() { return mTetA ; }

                mesh::Element *
                tet_b() { return mTetB ; }

                mesh::Facet *
                facet() { return mFacet ; }

                //! the layer block's calculator, linked to the prism
                Calculator *
                layer_calculator()
                {
                    Calculator * aCalc = mField->block( 3 )->calculator() ;
                    aCalc->link( mField->block( 3 )->element( mPrism->id() ) );
                    return aCalc ;
                }

                //! apex node of tet A ( belongs to tet A only )
                mesh::Node *
                apex_a() { return mTetA->node( 3 ) ; }

                //! apex node of tet B ( belongs to tet B only )
                mesh::Node *
                apex_b() { return mTetB->node( 3 ) ; }

                //! centroid of the tape facet
                void
                facet_centroid( real aX[ 3 ] ) const ;

                /**
                 * edge dofs of h = a + b x r on every edge of aElement:
                 * the line integral in the edge's own node order,
                 * a . ( q - p ) + b . ( p x q ) ( exact for a straight edge )
                 */
                void
                seed_edge_field( mesh::Element * aElement,
                                 const real aA[ 3 ],
                                 const real aB[ 3 ] );

                //! phi = -a . r on every node of aElement ( -grad phi = a )
                void
                seed_phi_uniform( mesh::Element * aElement, const real aA[ 3 ] );

                //! phi = aValue on one node
                void
                set_phi( mesh::Node * aNode, const real aValue );
            };

//------------------------------------------------------------------------------

            inline
            TS_TestConductorShell::TS_TestConductorShell(
                    const DomainType aTypeA,
                    const DomainType aTypeB,
                    const int        aFlipEdge,
                    const real       aThickness,
                    const real       aHeight ) :
                    mHeight( aHeight )
            {
                // no connectivity computation inside finalize(): facets and
                // edges are built by hand
                mMesh = new Mesh( 3, 0, false );

                mesh::ElementFactory tFactory ;

                // triangle footprint, deliberately not equilateral
                const real tX[ 3 ] = { 0.0, 1.3, 0.9 };
                const real tY[ 3 ] = { 0.0, 0.1, 1.1 };

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                id_t tID = 1 ;

                // nodes 0-2: the tape facet at z = 0
                for ( uint k = 0; k < 3; ++k )
                {
                    tNodes.push( new mesh::Node( tID++, tX[ k ], tY[ k ], 0.0 ) );
                }

                // node 3: apex of tet A, node 4: apex of tet B
                tNodes.push( new mesh::Node( tID++, 0.7, 0.4, -aHeight ) );
                tNodes.push( new mesh::Node( tID++, 0.7, 0.4,  aHeight ) );

                // nodes 5-7: top of the layer prism
                for ( uint k = 0; k < 3; ++k )
                {
                    tNodes.push( new mesh::Node( tID++, tX[ k ], tY[ k ], aThickness ) );
                }

                Cell< mesh::Element * > & tElements = mMesh->elements();

                // tet A below the facet: node order ( n0, n2, n1, apex ) gives
                // a positive volume with the apex below; its local face 3 is
                // ( mNodes[0], mNodes[2], mNodes[1] ) = ( n0, n1, n2 ), the
                // facet CCW from +z, outward normal +z
                mesh::Block * tBlockA = new mesh::Block( 1, 1 );
                tBlockA->set_domain_type( aTypeA );
                mTetA = tFactory.create_element( ElementType::TET4, 1 );
                mTetA->insert_node( tNodes( 0 ), 0 );
                mTetA->insert_node( tNodes( 2 ), 1 );
                mTetA->insert_node( tNodes( 1 ), 2 );
                mTetA->insert_node( tNodes( 3 ), 3 );
                mTetA->set_block_id( 1 );
                tElements.push( mTetA );
                tBlockA->insert_element( mTetA );
                mMesh->blocks().push( tBlockA );

                // tet B above the facet: ( n0, n1, n2, apex ) is positive with
                // the apex above; its local face 3 is ( n0, n2, n1 ), outward
                // normal -z
                mesh::Block * tBlockB = new mesh::Block( 2, 1 );
                tBlockB->set_domain_type( aTypeB );
                mTetB = tFactory.create_element( ElementType::TET4, 2 );
                mTetB->insert_node( tNodes( 0 ), 0 );
                mTetB->insert_node( tNodes( 1 ), 1 );
                mTetB->insert_node( tNodes( 2 ), 2 );
                mTetB->insert_node( tNodes( 4 ), 3 );
                mTetB->set_block_id( 2 );
                tElements.push( mTetB );
                tBlockB->insert_element( mTetB );
                mMesh->blocks().push( tBlockB );

                // layer prism, block 3: bottom = facet nodes, top = nodes 5-7
                mesh::Block * tBlockL = new mesh::Block( 3, 1 );
                tBlockL->set_thickness( aThickness );
                tBlockL->set_domain_type( DomainType::ThinShell );
                mPrism = tFactory.create_element( ElementType::PENTA6TS, 3 );
                for ( uint k = 0; k < 3; ++k )
                {
                    mPrism->insert_node( tNodes( k ),     k );
                    mPrism->insert_node( tNodes( k + 5 ), k + 3 );
                }
                mPrism->set_block_id( 3 );
                tElements.push( mPrism );
                tBlockL->insert_element( mPrism );
                mMesh->blocks().push( tBlockL );

                // tape facet: TRI3 on the z = 0 nodes; set_master relinks its
                // nodes to the master's face order, the slave orientation is
                // computed from the shared nodes
                mesh::Element * tSurf =
                        tFactory.create_element( ElementType::TRI3, 100 );
                for ( uint k = 0; k < 3; ++k )
                {
                    tSurf->insert_node( tNodes( k ), k );
                }
                mFacet = new mesh::Facet( tSurf );
                mFacet->set_master( mTetA, 3 );
                mFacet->set_slave( mTetB, 3 );
                mFacet->compute_orientation();
                mFacet->set_sideset_id( 10 );

                mesh::SideSet * tTapeSet = new mesh::SideSet( 10, 1 );
                tTapeSet->set_domain_type( DomainType::ThinShell );
                tTapeSet->insert_facet( mFacet );
                mMesh->sidesets().push( tTapeSet );

                // thin-shell container: BlockData::link_thin_shell_facets
                // hands the facet to the layer element through it ( the ctor
                // hides the sideset from mesh output — harmless here )
                mesh::ThinShell * tShell = new mesh::ThinShell( tTapeSet, nullptr );
                tShell->blocks().push( tBlockL );
                {
                    Vector< real > tThicknesses( 1, aThickness );
                    tShell->set_thicknesses( tThicknesses );
                }
                mMesh->thin_shells().push( tShell );

                mMesh->finalize();

                // reorientation check must run BEFORE edges exist; it sets
                // the mesh checker flag so the Kernel ctor skips its own run
                {
                    MeshChecker tCheck( mMesh );
                }

                // per-element edge objects in canonical get_nodes_of_edge()
                // order on the conductor tets and the prism; an air tet
                // carries none. aFlipEdge reverses that local edge on the tets
                {
                    id_t tEdgeID = 1 ;
                    Cell< mesh::Node * > tEdgeNodes ;

                    mesh::Element * tWithEdges[ 3 ] = {
                        aTypeA == DomainType::Conductor ? mTetA : nullptr,
                        aTypeB == DomainType::Conductor ? mTetB : nullptr,
                        mPrism };

                    for ( mesh::Element * tElement : tWithEdges )
                    {
                        if ( tElement == nullptr ) continue ;

                        const bool tIsTet = ( tElement != mPrism );
                        const uint tNumEdges = tElement->number_of_edges();
                        tElement->allocate_edge_container();
                        for ( uint e = 0; e < tNumEdges; ++e )
                        {
                            tElement->get_nodes_of_edge( e, tEdgeNodes );
                            mesh::Edge * tEdge = new mesh::Edge();
                            tEdge->set_id( tEdgeID++ );
                            tEdge->allocate_node_container( 2 );
                            if ( tIsTet && ( int ) e == aFlipEdge )
                            {
                                tEdge->insert_node( tEdgeNodes( 1 ), 0 );
                                tEdge->insert_node( tEdgeNodes( 0 ), 1 );
                            }
                            else
                            {
                                tEdge->insert_node( tEdgeNodes( 0 ), 0 );
                                tEdge->insert_node( tEdgeNodes( 1 ), 1 );
                            }
                            tEdge->set_owner( 0 );
                            tElement->insert_edge( tEdge, e );
                            mMesh->edges().push( tEdge );
                        }
                    }
                    mMesh->finalize_edges();
                }

                mParams = new KernelParameters( mMesh );
                mKernel = new Kernel( mParams );

                mIWG = new IWG_Maxwell(
                        maxwell::Formulation::HPhi,
                        ModelDimensionality::ThreeD,
                        false, false );

                Vector< id_t > tBlockIDs( 3 );
                tBlockIDs( 0 ) = 1 ;
                tBlockIDs( 1 ) = 2 ;
                tBlockIDs( 2 ) = 3 ;
                Cell< DomainType > tBlockTypes( 3, DomainType::ThinShell );
                tBlockTypes( 0 ) = aTypeA ;
                tBlockTypes( 1 ) = aTypeB ;
                mIWG->set_blocks( tBlockIDs, tBlockTypes );

                Vector< id_t > tSideSetIDs( 1 );
                tSideSetIDs( 0 ) = 10 ;
                Cell< DomainType > tSideSetTypes( 1, DomainType::ThinShell );
                mIWG->set_sidesets( tSideSetIDs, tSideSetTypes );

                mKernel->add_equation( mIWG );

                mField = mKernel->create_field( mIWG );
                mField->create_fields( mField->iwg() );

                mField->set_solver( gDefaultSolver );
                mField->initialize();
            }

//------------------------------------------------------------------------------

            inline
            TS_TestConductorShell::~TS_TestConductorShell()
            {
                // kernel owns the dof manager and the equation; the mesh
                // owns nodes, blocks, elements, sidesets, facets, edges,
                // and the thin-shell record
                delete mKernel ;
                delete mParams ;
                delete mMesh ;
            }

//------------------------------------------------------------------------------

            inline void
            TS_TestConductorShell::facet_centroid( real aX[ 3 ] ) const
            {
                aX[ 0 ] = 0.0 ; aX[ 1 ] = 0.0 ; aX[ 2 ] = 0.0 ;
                for ( uint k = 0; k < 3; ++k )
                {
                    aX[ 0 ] += mFacet->node( k )->x() / 3.0 ;
                    aX[ 1 ] += mFacet->node( k )->y() / 3.0 ;
                    aX[ 2 ] += mFacet->node( k )->z() / 3.0 ;
                }
            }

//------------------------------------------------------------------------------

            inline void
            TS_TestConductorShell::seed_edge_field(
                    mesh::Element * aElement,
                    const real aA[ 3 ],
                    const real aB[ 3 ] )
            {
                Vector< real > & tField = mMesh->field_data( "edge_h" );

                for ( uint e = 0; e < aElement->number_of_edges(); ++e )
                {
                    mesh::Edge * tEdge = aElement->edge( e );
                    const real p[ 3 ] = { tEdge->node( 0 )->x(), tEdge->node( 0 )->y(), tEdge->node( 0 )->z() };
                    const real q[ 3 ] = { tEdge->node( 1 )->x(), tEdge->node( 1 )->y(), tEdge->node( 1 )->z() };

                    // p x q
                    const real c[ 3 ] = {
                        p[ 1 ] * q[ 2 ] - p[ 2 ] * q[ 1 ],
                        p[ 2 ] * q[ 0 ] - p[ 0 ] * q[ 2 ],
                        p[ 0 ] * q[ 1 ] - p[ 1 ] * q[ 0 ] };

                    real tDof = 0.0 ;
                    for ( uint i = 0; i < 3; ++i )
                    {
                        tDof += aA[ i ] * ( q[ i ] - p[ i ] ) + aB[ i ] * c[ i ];
                    }
                    tField( tEdge->index() ) = tDof ;
                }
            }

//------------------------------------------------------------------------------

            inline void
            TS_TestConductorShell::seed_phi_uniform(
                    mesh::Element * aElement,
                    const real aA[ 3 ] )
            {
                Vector< real > & tPhi = mMesh->field_data( "phi" );
                for ( uint k = 0; k < aElement->number_of_nodes(); ++k )
                {
                    mesh::Node * tNode = aElement->node( k );
                    tPhi( tNode->index() ) = -( aA[ 0 ] * tNode->x()
                                              + aA[ 1 ] * tNode->y()
                                              + aA[ 2 ] * tNode->z() );
                }
            }

//------------------------------------------------------------------------------

            inline void
            TS_TestConductorShell::set_phi( mesh::Node * aNode, const real aValue )
            {
                mMesh->field_data( "phi" )( aNode->index() ) = aValue ;
            }

//------------------------------------------------------------------------------

            /**
             * 2-D twin of TS_TestConductorShell: one QUAD4TS layer between
             * two TRI3 volumes. The tape facet is the LINE2 from ( 0, 0 ) to
             * ( aLength, 0 ); tri A below it ( apex at y = -aHeight, node
             * order ( n1, n0, apex ) for a CCW triangle, local facet 0 =
             * ( n1, n0 ) with outward normal +y ), tri B above ( ( n0, n1,
             * apex ), facet 0 = ( n0, n1 ) ), the layer from y = 0 to
             * y = aThickness. Same block ids ( 1, 2, 3 ), sideset id 10, and
             * the same seeding helpers, with h = a + b ( -y, x ) so that
             * the normal trace a_y + b x varies along the facet: the 2-D
             * default line rule has no point at the midpoint, which is
             * what the weighted mean in compute_h_trace is for.
             */
            class TS_TestConductorShell2D
            {
                Mesh                   * mMesh ;
                KernelParameters       * mParams ;
                Kernel                 * mKernel ;
                DofManager             * mField ;
                IWG_Maxwell            * mIWG ;

                mesh::Element          * mTriA  = nullptr ;
                mesh::Element          * mTriB  = nullptr ;
                mesh::Element          * mQuad  = nullptr ;
                mesh::Facet            * mFacet = nullptr ;

                const real               mLength ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                TS_TestConductorShell2D(
                        const DomainType aTypeA,
                        const DomainType aTypeB,
                        const int        aFlipEdge  = -1,
                        const real       aLength    = 2.0,
                        const real       aThickness = 0.1,
                        const real       aHeight    = 0.7 );

                TS_TestConductorShell2D( const TS_TestConductorShell2D & ) = delete ;
                TS_TestConductorShell2D &
                operator=( const TS_TestConductorShell2D & ) = delete ;

                ~TS_TestConductorShell2D();

                Mesh *
                mesh() { return mMesh ; }

                real
                length() const { return mLength ; }

                mesh::Element *
                tri_a() { return mTriA ; }

                mesh::Element *
                tri_b() { return mTriB ; }

                mesh::Node *
                apex_a() { return mTriA->node( 2 ) ; }

                mesh::Node *
                apex_b() { return mTriB->node( 2 ) ; }

                Calculator *
                layer_calculator()
                {
                    Calculator * aCalc = mField->block( 3 )->calculator() ;
                    aCalc->link( mField->block( 3 )->element( mQuad->id() ) );
                    return aCalc ;
                }

                /**
                 * edge dofs of h = a + b ( -y, x ) on every edge of aElement:
                 * a . ( q - p ) + b ( p x q ) in the edge's own node order
                 */
                void
                seed_edge_field( mesh::Element * aElement,
                                 const real aA[ 2 ],
                                 const real aB );

                //! phi = -a . r on every node of aElement
                void
                seed_phi_uniform( mesh::Element * aElement, const real aA[ 2 ] );

                void
                set_phi( mesh::Node * aNode, const real aValue )
                {
                    mMesh->field_data( "phi" )( aNode->index() ) = aValue ;
                }
            };

//------------------------------------------------------------------------------

            inline
            TS_TestConductorShell2D::TS_TestConductorShell2D(
                    const DomainType aTypeA,
                    const DomainType aTypeB,
                    const int        aFlipEdge,
                    const real       aLength,
                    const real       aThickness,
                    const real       aHeight ) :
                    mLength( aLength )
            {
                mMesh = new Mesh( 2, 0, false );

                mesh::ElementFactory tFactory ;

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                id_t tID = 1 ;

                // nodes 0-1: the tape facet on y = 0
                tNodes.push( new mesh::Node( tID++, 0.0,     0.0 ) );
                tNodes.push( new mesh::Node( tID++, aLength, 0.0 ) );

                // node 2: apex of tri A, node 3: apex of tri B
                tNodes.push( new mesh::Node( tID++, 0.5 * aLength, -aHeight ) );
                tNodes.push( new mesh::Node( tID++, 0.5 * aLength,  aHeight ) );

                // nodes 4-5: top of the layer
                tNodes.push( new mesh::Node( tID++, 0.0,     aThickness ) );
                tNodes.push( new mesh::Node( tID++, aLength, aThickness ) );

                Cell< mesh::Element * > & tElements = mMesh->elements();

                // tri A below: ( n1, n0, apex ) is CCW, facet 0 = ( n1, n0 )
                mesh::Block * tBlockA = new mesh::Block( 1, 1 );
                tBlockA->set_domain_type( aTypeA );
                mTriA = tFactory.create_element( ElementType::TRI3, 1 );
                mTriA->insert_node( tNodes( 1 ), 0 );
                mTriA->insert_node( tNodes( 0 ), 1 );
                mTriA->insert_node( tNodes( 2 ), 2 );
                mTriA->set_block_id( 1 );
                tElements.push( mTriA );
                tBlockA->insert_element( mTriA );
                mMesh->blocks().push( tBlockA );

                // tri B above: ( n0, n1, apex ) is CCW, facet 0 = ( n0, n1 )
                mesh::Block * tBlockB = new mesh::Block( 2, 1 );
                tBlockB->set_domain_type( aTypeB );
                mTriB = tFactory.create_element( ElementType::TRI3, 2 );
                mTriB->insert_node( tNodes( 0 ), 0 );
                mTriB->insert_node( tNodes( 1 ), 1 );
                mTriB->insert_node( tNodes( 3 ), 2 );
                mTriB->set_block_id( 2 );
                tElements.push( mTriB );
                tBlockB->insert_element( mTriB );
                mMesh->blocks().push( tBlockB );

                // layer quad, block 3: ( bottom0, bottom1, top1, top0 ) as in
                // TS_TestStack2D
                mesh::Block * tBlockL = new mesh::Block( 3, 1 );
                tBlockL->set_thickness( aThickness );
                tBlockL->set_domain_type( DomainType::ThinShell );
                mQuad = tFactory.create_element( ElementType::QUAD4TS, 3 );
                mQuad->insert_node( tNodes( 0 ), 0 );
                mQuad->insert_node( tNodes( 1 ), 1 );
                mQuad->insert_node( tNodes( 5 ), 2 );
                mQuad->insert_node( tNodes( 4 ), 3 );
                mQuad->set_block_id( 3 );
                tElements.push( mQuad );
                tBlockL->insert_element( mQuad );
                mMesh->blocks().push( tBlockL );

                // tape facet
                mesh::Element * tLine =
                        tFactory.create_element( ElementType::LINE2, 100 );
                tLine->insert_node( tNodes( 0 ), 0 );
                tLine->insert_node( tNodes( 1 ), 1 );
                mFacet = new mesh::Facet( tLine );
                mFacet->set_master( mTriA, 0 );
                mFacet->set_slave( mTriB, 0 );
                mFacet->compute_orientation();
                mFacet->set_sideset_id( 10 );

                mesh::SideSet * tTapeSet = new mesh::SideSet( 10, 1 );
                tTapeSet->set_domain_type( DomainType::ThinShell );
                tTapeSet->insert_facet( mFacet );
                mMesh->sidesets().push( tTapeSet );

                mesh::ThinShell * tShell = new mesh::ThinShell( tTapeSet, nullptr );
                tShell->blocks().push( tBlockL );
                {
                    Vector< real > tThicknesses( 1, aThickness );
                    tShell->set_thicknesses( tThicknesses );
                }
                mMesh->thin_shells().push( tShell );

                mMesh->finalize();

                {
                    MeshChecker tCheck( mMesh );
                }

                {
                    id_t tEdgeID = 1 ;
                    Cell< mesh::Node * > tEdgeNodes ;

                    mesh::Element * tWithEdges[ 3 ] = {
                        aTypeA == DomainType::Conductor ? mTriA : nullptr,
                        aTypeB == DomainType::Conductor ? mTriB : nullptr,
                        mQuad };

                    for ( mesh::Element * tElement : tWithEdges )
                    {
                        if ( tElement == nullptr ) continue ;

                        const bool tIsTri = ( tElement != mQuad );
                        const uint tNumEdges = tElement->number_of_edges();
                        tElement->allocate_edge_container();
                        for ( uint e = 0; e < tNumEdges; ++e )
                        {
                            tElement->get_nodes_of_edge( e, tEdgeNodes );
                            mesh::Edge * tEdge = new mesh::Edge();
                            tEdge->set_id( tEdgeID++ );
                            tEdge->allocate_node_container( 2 );
                            if ( tIsTri && ( int ) e == aFlipEdge )
                            {
                                tEdge->insert_node( tEdgeNodes( 1 ), 0 );
                                tEdge->insert_node( tEdgeNodes( 0 ), 1 );
                            }
                            else
                            {
                                tEdge->insert_node( tEdgeNodes( 0 ), 0 );
                                tEdge->insert_node( tEdgeNodes( 1 ), 1 );
                            }
                            tEdge->set_owner( 0 );
                            tElement->insert_edge( tEdge, e );
                            mMesh->edges().push( tEdge );
                        }
                    }
                    mMesh->finalize_edges();
                }

                mParams = new KernelParameters( mMesh );
                mKernel = new Kernel( mParams );

                mIWG = new IWG_Maxwell(
                        maxwell::Formulation::HPhi,
                        ModelDimensionality::TwoD,
                        false, false );

                Vector< id_t > tBlockIDs( 3 );
                tBlockIDs( 0 ) = 1 ;
                tBlockIDs( 1 ) = 2 ;
                tBlockIDs( 2 ) = 3 ;
                Cell< DomainType > tBlockTypes( 3, DomainType::ThinShell );
                tBlockTypes( 0 ) = aTypeA ;
                tBlockTypes( 1 ) = aTypeB ;
                mIWG->set_blocks( tBlockIDs, tBlockTypes );

                Vector< id_t > tSideSetIDs( 1 );
                tSideSetIDs( 0 ) = 10 ;
                Cell< DomainType > tSideSetTypes( 1, DomainType::ThinShell );
                mIWG->set_sidesets( tSideSetIDs, tSideSetTypes );

                mKernel->add_equation( mIWG );

                mField = mKernel->create_field( mIWG );
                mField->create_fields( mField->iwg() );

                mField->set_solver( gDefaultSolver );
                mField->initialize();
            }

//------------------------------------------------------------------------------

            inline
            TS_TestConductorShell2D::~TS_TestConductorShell2D()
            {
                delete mKernel ;
                delete mParams ;
                delete mMesh ;
            }

//------------------------------------------------------------------------------

            inline void
            TS_TestConductorShell2D::seed_edge_field(
                    mesh::Element * aElement,
                    const real aA[ 2 ],
                    const real aB )
            {
                Vector< real > & tField = mMesh->field_data( "edge_h" );

                for ( uint e = 0; e < aElement->number_of_edges(); ++e )
                {
                    mesh::Edge * tEdge = aElement->edge( e );
                    const real p[ 2 ] = { tEdge->node( 0 )->x(), tEdge->node( 0 )->y() };
                    const real q[ 2 ] = { tEdge->node( 1 )->x(), tEdge->node( 1 )->y() };

                    tField( tEdge->index() ) =
                              aA[ 0 ] * ( q[ 0 ] - p[ 0 ] )
                            + aA[ 1 ] * ( q[ 1 ] - p[ 1 ] )
                            + aB * ( p[ 0 ] * q[ 1 ] - p[ 1 ] * q[ 0 ] );
                }
            }

//------------------------------------------------------------------------------

            inline void
            TS_TestConductorShell2D::seed_phi_uniform(
                    mesh::Element * aElement,
                    const real aA[ 2 ] )
            {
                Vector< real > & tPhi = mMesh->field_data( "phi" );
                for ( uint k = 0; k < aElement->number_of_nodes(); ++k )
                {
                    mesh::Node * tNode = aElement->node( k );
                    tPhi( tNode->index() ) = -( aA[ 0 ] * tNode->x()
                                              + aA[ 1 ] * tNode->y() );
                }
            }

//------------------------------------------------------------------------------
        } /* end namespace test */
    } /* end namespace fem */
} /* end namespace belfem */

#endif // BELFEM_CL_TS_TESTSTACK_HPP
