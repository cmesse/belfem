/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for Mesh::update_facet_nodes(): after a cut duplicate replaces a
 * node on the master element of a facet and the mesh is re-finalized, a facet
 * in a ThinShell or GeometryOnly sideset keeps the node list it was built
 * with, while a facet in a non-empty sideset of any other type whose first
 * facet has a master is re-derived from that master.
 *
 * In production the tape facet's node list is the extrusion geometry the
 * cut factory restores; re-deriving it would hand the layer elements a
 * duplicate that belongs to the volume side of the terminus. This fixture
 * runs no cut factory and reads no prism: it pins only the skip.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "en_DomainType.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
//------------------------------------------------------------------------------

    // tape terminus geometry: two TET4 sharing the z = 0 triangle, a PENTA6TS
    // layer above it, the TRI3 tape facet in sideset 10 ( master tet A face 3,
    // slave tet B face 3 ) and one plain TRI3 facet in sideset 20 on tet A
    // face 0. Connectivities are off, so the facets and their masters are the
    // ones set here and finalize() reaches update_facet_nodes() without the
    // connectivity calculator.
    //
    // Tet A carries the nodes in the order ( n0, n2, n1, n3 ), so with
    // ElementTemplate< 4, 4, 6, 4, 4 >::get_nodes_of_facet
    //   face 3 = ( mNodes[ 0 ], mNodes[ 2 ], mNodes[ 1 ] ) = ( n0, n1, n2 ),
    //            ids ( 1, 2, 3 ) -- the tape facet
    //   face 0 = ( mNodes[ 0 ], mNodes[ 1 ], mNodes[ 3 ] ) = ( n0, n2, n3 ),
    //            ids ( 1, 3, 4 ) -- the plain facet
    // Both hang on n0, the node the terminus duplicates.
    //
    // Ownership: blocks delete the volume elements, sidesets delete the
    // facets, a facet deletes its wrapper, the mesh deletes the nodes. The
    // wrappers therefore never enter a block.
    struct TapeMesh
    {
        Mesh    * mMesh ;
        Element * mTetA ;
        Element * mTetB ;
        Element * mPrism ;
        Facet   * mTapeFacet ;
        Facet   * mPlainFacet ;

        // the terminus duplicate, null until cut_terminus() has run
        Node * mDuplicate = nullptr ;

        explicit TapeMesh( const DomainType aTapeDomainType )
        {
            const real tThickness = 0.1 ;
            const real tHeight    = 0.5 ;

            mMesh = new Mesh( 3, 0, false );

            ElementFactory tFactory ;

            // triangle footprint, deliberately not equilateral
            const real tX[ 3 ] = { 0.0, 1.3, 0.9 };
            const real tY[ 3 ] = { 0.0, 0.1, 1.1 };

            Cell< Node * > & tNodes = mMesh->nodes();
            id_t tID = 1 ;

            // nodes 1-3: the tape triangle at z = 0
            for ( uint k = 0; k < 3; ++k )
            {
                tNodes.push( new Node( tID++, tX[ k ], tY[ k ], 0.0 ) );
            }

            // node 4: apex of tet A, node 5: apex of tet B
            tNodes.push( new Node( tID++, 0.7, 0.4, -tHeight ) );
            tNodes.push( new Node( tID++, 0.7, 0.4,  tHeight ) );

            // nodes 6-8: top of the layer prism
            for ( uint k = 0; k < 3; ++k )
            {
                tNodes.push( new Node( tID++, tX[ k ], tY[ k ], tThickness ) );
            }

            Cell< Element * > & tElements = mMesh->elements();

            // tet A below the tape, apex below: ( n0, n2, n1, apex )
            Block * tBlockA = new Block( 1, 1 );
            tBlockA->set_domain_type( DomainType::Air );
            mTetA = tFactory.create_element( ElementType::TET4, 1 );
            mTetA->insert_node( tNodes( 0 ), 0 );
            mTetA->insert_node( tNodes( 2 ), 1 );
            mTetA->insert_node( tNodes( 1 ), 2 );
            mTetA->insert_node( tNodes( 3 ), 3 );
            mTetA->set_block_id( 1 );
            tElements.push( mTetA );
            tBlockA->insert_element( mTetA );
            mMesh->blocks().push( tBlockA );

            // tet B above the tape, apex above: ( n0, n1, n2, apex )
            Block * tBlockB = new Block( 2, 1 );
            tBlockB->set_domain_type( DomainType::Air );
            mTetB = tFactory.create_element( ElementType::TET4, 2 );
            mTetB->insert_node( tNodes( 0 ), 0 );
            mTetB->insert_node( tNodes( 1 ), 1 );
            mTetB->insert_node( tNodes( 2 ), 2 );
            mTetB->insert_node( tNodes( 4 ), 3 );
            mTetB->set_block_id( 2 );
            tElements.push( mTetB );
            tBlockB->insert_element( mTetB );
            mMesh->blocks().push( tBlockB );

            // layer prism: bottom = tape nodes, top = nodes 6-8
            Block * tBlockL = new Block( 3, 1 );
            tBlockL->set_thickness( tThickness );
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

            // tape facet: set_master links its nodes to face 3 of tet A, the
            // slave orientation is taken from the shared nodes before the
            // terminus is cut, so the second finalize() finds it already set
            Element * tTapeSurf = tFactory.create_element( ElementType::TRI3, 100 );
            for ( uint k = 0; k < 3; ++k )
            {
                tTapeSurf->insert_node( tNodes( k ), k );
            }
            mTapeFacet = new Facet( tTapeSurf );
            mTapeFacet->set_master( mTetA, 3 );
            mTapeFacet->set_slave( mTetB, 3 );
            mTapeFacet->compute_orientation();
            mTapeFacet->set_sideset_id( 10 );

            SideSet * tTapeSet = new SideSet( 10, 1 );
            tTapeSet->set_domain_type( aTapeDomainType );
            tTapeSet->insert_facet( mTapeFacet );
            mMesh->sidesets().push( tTapeSet );

            // plain facet on face 0 of the same master, no slave
            Element * tPlainSurf = tFactory.create_element( ElementType::TRI3, 101 );
            mPlainFacet = new Facet( tPlainSurf );
            mPlainFacet->set_master( mTetA, 0 );
            mPlainFacet->set_sideset_id( 20 );

            SideSet * tPlainSet = new SideSet( 20, 1 );
            tPlainSet->set_domain_type( DomainType::Air );
            tPlainSet->insert_facet( mPlainFacet );
            mMesh->sidesets().push( tPlainSet );

            mMesh->finalize();
        }

        ~TapeMesh()
        {
            delete mMesh ;
        }

        // owning fixture: a copy or a move would double-free
        TapeMesh( const TapeMesh & ) = delete ;
        TapeMesh & operator=( const TapeMesh & ) = delete ;
        TapeMesh( TapeMesh && ) = delete ;
        TapeMesh & operator=( TapeMesh && ) = delete ;

        // a cut terminating on the tape: a duplicate of n0 takes its place on
        // the master element of both facets, then the mesh is re-finalized
        void
        cut_terminus()
        {
            Node * tOrigin = mMesh->nodes()( 0 );

            mDuplicate = new Node( 9, tOrigin->x(), tOrigin->y(), tOrigin->z() );
            mMesh->nodes().push( mDuplicate );

            mTetA->insert_node( mDuplicate, 0 );

            mMesh->unfinalize();
            mMesh->finalize();
        }
    };

//------------------------------------------------------------------------------

    // the tape facet still carries the triangle it was built with, and none of
    // its nodes is the duplicate; the relink on the master is checked too, so
    // the expectation cannot be met by a relink that did not happen
    void
    expect_tape_facet_kept( TapeMesh & aTape )
    {
        ASSERT_EQ( aTape.mTapeFacet->number_of_nodes(), 3u );

        const id_t tExpected[ 3 ] = { 1, 2, 3 };

        for ( uint k = 0; k < 3; ++k )
        {
            EXPECT_EQ( aTape.mTapeFacet->node( k )->id(), tExpected[ k ] );
            EXPECT_NE( aTape.mTapeFacet->node( k ), aTape.mDuplicate );
        }

        EXPECT_EQ( aTape.mTetA->node( 0 ), aTape.mDuplicate );
    }

//------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------

TEST( FacetNodes, ThinShellFacetKeepsItsNodesAfterRelink )
{
    TapeMesh tTape( DomainType::ThinShell );
    tTape.cut_terminus();

    expect_tape_facet_kept( tTape );
}

//------------------------------------------------------------------------------

TEST( FacetNodes, GeometryOnlyFacetKeepsItsNodesAfterRelink )
{
    TapeMesh tTape( DomainType::GeometryOnly );
    tTape.cut_terminus();

    expect_tape_facet_kept( tTape );
}

//------------------------------------------------------------------------------

TEST( FacetNodes, PlainFacetIsRederivedFromItsMaster )
{
    TapeMesh tTape( DomainType::ThinShell );
    tTape.cut_terminus();

    Cell< Node * > tNodes ;
    tTape.mTetA->get_nodes_of_facet( 0, tNodes );

    ASSERT_EQ( tNodes.size(), 3u );
    ASSERT_EQ( tTape.mPlainFacet->number_of_nodes(), 3u );

    for ( uint k = 0; k < 3; ++k )
    {
        EXPECT_EQ( tTape.mPlainFacet->node( k ), tNodes( k ) );
    }

    // face 0 of tet A is ( duplicate, n2, n3 ) after the relink
    EXPECT_EQ( tTape.mPlainFacet->node( 0 ), tTape.mDuplicate );
    EXPECT_EQ( tTape.mPlainFacet->node( 1 )->id(), static_cast< id_t >( 3 ) );
    EXPECT_EQ( tTape.mPlainFacet->node( 2 )->id(), static_cast< id_t >( 4 ) );
}

//------------------------------------------------------------------------------
