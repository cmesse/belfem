/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for mesh entities: Element factory topology catalog, Node class,
 * Element node insertion/access, facet topology, base class guards.
 * See: tests_10_mesh.md §1–§3
 *
 * Deferred: I/O, distribution, periodicity, thin shells, curves, order conversion.
 */

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "Mesh_Enums.hpp"

namespace
{
    const belfem::real tEps = 1e-12;

    // Topology descriptor for catalog test
    struct ElementCatalogEntry
    {
        belfem::ElementType mType;
        belfem::uint mN;    // number_of_nodes
        belfem::uint mC;    // number_of_corner_nodes
        belfem::uint mE;    // number_of_edges
        belfem::uint mT;    // number_of_facets
        belfem::uint mF;    // number_of_faces
        belfem::uint mDim;  // dimension
    };
}

// =============================================================================
// §1.1  Element Catalog (Data-Driven)  [semantic]
// =============================================================================

TEST( ElementCatalog, AllElementTypesTopology )
{
    // Verified 2026-03-22 against cl_Element_Factory.cpp
    // TET20, TET35 omitted: no type()/dimension() specializations
    const ElementCatalogEntry tCatalog[] = {
        { belfem::ElementType::VERTEX,   1,  1,  0, 0, 0, 0 },
        { belfem::ElementType::LINE2,    2,  2,  1, 0, 0, 1 },
        { belfem::ElementType::LINE3,    3,  2,  1, 0, 0, 1 },
        { belfem::ElementType::LINE4,    4,  2,  1, 0, 0, 1 },
        { belfem::ElementType::LINE5,    5,  2,  1, 0, 0, 1 },
        { belfem::ElementType::TRI3,     3,  3,  3, 3, 1, 2 },
        { belfem::ElementType::TRI6,     6,  3,  3, 3, 1, 2 },
        { belfem::ElementType::TRI10,   10,  3,  3, 3, 1, 2 },
        { belfem::ElementType::TRI15,   15,  3,  3, 3, 1, 2 },
        { belfem::ElementType::QUAD4,    4,  4,  4, 4, 1, 2 },
        { belfem::ElementType::QUAD8,    8,  4,  4, 4, 1, 2 },
        { belfem::ElementType::QUAD9,    9,  4,  4, 4, 1, 2 },
        { belfem::ElementType::QUAD16,  16,  4,  4, 4, 1, 2 },
        { belfem::ElementType::TET4,     4,  4,  6, 4, 4, 3 },
        { belfem::ElementType::TET10,   10,  4,  6, 4, 4, 3 },
        { belfem::ElementType::PENTA6,   6,  6,  9, 5, 5, 3 },
        { belfem::ElementType::PENTA15, 15,  6,  9, 5, 5, 3 },
        { belfem::ElementType::PENTA18, 18,  6,  9, 5, 5, 3 },
        { belfem::ElementType::PYRA5,    5,  5,  8, 5, 5, 3 },
        { belfem::ElementType::PYRA13,  13,  5,  8, 5, 5, 3 },
        { belfem::ElementType::PYRA14,  14,  5,  8, 5, 5, 3 },
        { belfem::ElementType::HEX8,     8,  8, 12, 6, 6, 3 },
        { belfem::ElementType::HEX20,   20,  8, 12, 6, 6, 3 },
        { belfem::ElementType::HEX27,   27,  8, 12, 6, 6, 3 },
        { belfem::ElementType::HEX64,   64,  8, 12, 6, 6, 3 },
        // Thin-shell variants now share facet topology with their volume
        // counterparts (Option-B canonicalization, see cl_Element_<TYPE>TS.hpp).
        // Face count (mF) remains shell-specific: 1 for the midsurface,
        // 3 for quadratic PENTA18TS (bottom/mid/top layers).
        { belfem::ElementType::QUAD4TS,   4,  4,  2, 4, 1, 2 },
        { belfem::ElementType::QUAD9TS,   9,  4,  3, 4, 1, 2 },
        { belfem::ElementType::PENTA6TS,  6,  6,  6, 5, 1, 3 },
        { belfem::ElementType::PENTA18TS,18,  6,  9, 5, 3, 3 },
    };

    belfem::mesh::ElementFactory tFactory;

    for( const auto & tEntry : tCatalog )
    {
        belfem::mesh::Element * tElem = tFactory.create_element( tEntry.mType, 1 );

        EXPECT_EQ( tElem->number_of_nodes(),        tEntry.mN )   ;
        EXPECT_EQ( tElem->number_of_corner_nodes(),  tEntry.mC )  ;
        EXPECT_EQ( tElem->number_of_edges(),         tEntry.mE )  ;
        EXPECT_EQ( tElem->number_of_facets(),        tEntry.mT )  ;
        EXPECT_EQ( tElem->number_of_faces(),         tEntry.mF )  ;
        EXPECT_EQ( tElem->dimension(),               tEntry.mDim );
        EXPECT_EQ( tElem->type(),                    tEntry.mType );

        delete tElem;
    }
}

// =============================================================================
// §1.2  Element Factory Unknown Type  [semantic]
// =============================================================================

TEST( ElementCatalog, FactoryUnknownTypeThrows )
{
    belfem::mesh::ElementFactory tFactory;

    EXPECT_THROW(
        tFactory.create_element( belfem::ElementType::UNDEFINED, 1 ),
        std::runtime_error );
}

// =============================================================================
// §2  Node Entity  [semantic]
// =============================================================================

TEST( NodeEntity, NodeConstruction )
{
    belfem::mesh::Node tNode( 1, 2.0, 3.0, 4.0 );

    EXPECT_EQ( tNode.id(), 1u );
    EXPECT_NEAR( tNode.x(), 2.0, tEps );
    EXPECT_NEAR( tNode.y(), 3.0, tEps );
    EXPECT_NEAR( tNode.z(), 4.0, tEps );
}

TEST( NodeEntity, NodeDefaultCoords )
{
    belfem::mesh::Node tNode( 42 );

    EXPECT_NEAR( tNode.x(), 0.0, tEps );
    EXPECT_NEAR( tNode.y(), 0.0, tEps );
    EXPECT_NEAR( tNode.z(), 0.0, tEps );
}

TEST( NodeEntity, NodeSetCoords3D )
{
    belfem::mesh::Node tNode( 1 );
    tNode.set_coords( 1.0, 2.0, 3.0 );

    EXPECT_NEAR( tNode.x(), 1.0, tEps );
    EXPECT_NEAR( tNode.y(), 2.0, tEps );
    EXPECT_NEAR( tNode.z(), 3.0, tEps );
}

TEST( NodeEntity, NodeEntityType )
{
    belfem::mesh::Node tNode( 1 );
    EXPECT_EQ( tNode.entity_type(), belfem::EntityType::NODE );
}

TEST( NodeEntity, NodeCoordsVector )
{
    belfem::mesh::Node tNode( 1, 5.0, 6.0, 7.0 );
    belfem::Vector< belfem::real > tC = tNode.coords();

    EXPECT_EQ( tC.length(), 3u );
    EXPECT_NEAR( tC( 0 ), 5.0, tEps );
    EXPECT_NEAR( tC( 1 ), 6.0, tEps );
    EXPECT_NEAR( tC( 2 ), 7.0, tEps );
}

TEST( NodeEntity, NodeSetCoordsVector )
{
    // Codex finding: the Vector<real> overload of set_coords was untested
    belfem::mesh::Node tNode( 1 );
    belfem::Vector< belfem::real > tC = { 10.0, 20.0, 30.0 };
    tNode.set_coords( tC );

    EXPECT_NEAR( tNode.x(), 10.0, tEps );
    EXPECT_NEAR( tNode.y(), 20.0, tEps );
    EXPECT_NEAR( tNode.z(), 30.0, tEps );
}

// =============================================================================
// §3.1  Element Node Insertion and Access  [semantic]
// =============================================================================

TEST( ElementEntity, ElementInsertAndAccessNodes )
{
    belfem::mesh::ElementFactory tFactory;
    belfem::mesh::Element * tElem = tFactory.create_element( belfem::ElementType::TRI3, 1 );

    belfem::mesh::Node tN0( 1, 0.0, 0.0 );
    belfem::mesh::Node tN1( 2, 1.0, 0.0 );
    belfem::mesh::Node tN2( 3, 0.0, 1.0 );

    tElem->insert_node( &tN0, 0 );
    tElem->insert_node( &tN1, 1 );
    tElem->insert_node( &tN2, 2 );

    EXPECT_EQ( tElem->node( 0 ), &tN0 );
    EXPECT_EQ( tElem->node( 1 ), &tN1 );
    EXPECT_EQ( tElem->node( 2 ), &tN2 );

    delete tElem;
}

TEST( ElementEntity, ElementIdStored )
{
    belfem::mesh::ElementFactory tFactory;
    belfem::mesh::Element * tElem = tFactory.create_element( belfem::ElementType::QUAD4, 42 );

    EXPECT_EQ( tElem->id(), 42u );

    delete tElem;
}

TEST( ElementEntity, ElementFlagNodes )
{
    belfem::mesh::ElementFactory tFactory;
    belfem::mesh::Element * tElem = tFactory.create_element( belfem::ElementType::TRI3, 1 );

    belfem::mesh::Node tN0( 1 );
    belfem::mesh::Node tN1( 2 );
    belfem::mesh::Node tN2( 3 );

    tElem->insert_node( &tN0, 0 );
    tElem->insert_node( &tN1, 1 );
    tElem->insert_node( &tN2, 2 );

    tElem->flag_nodes();

    EXPECT_TRUE( tN0.is_flagged() );
    EXPECT_TRUE( tN1.is_flagged() );
    EXPECT_TRUE( tN2.is_flagged() );

    delete tElem;
}

TEST( ElementEntity, ElementUnflagNodes )
{
    belfem::mesh::ElementFactory tFactory;
    belfem::mesh::Element * tElem = tFactory.create_element( belfem::ElementType::TRI3, 1 );

    belfem::mesh::Node tN0( 1 );
    belfem::mesh::Node tN1( 2 );
    belfem::mesh::Node tN2( 3 );

    tElem->insert_node( &tN0, 0 );
    tElem->insert_node( &tN1, 1 );
    tElem->insert_node( &tN2, 2 );

    tElem->flag_nodes();
    tElem->unflag_nodes();

    EXPECT_FALSE( tN0.is_flagged() );
    EXPECT_FALSE( tN1.is_flagged() );
    EXPECT_FALSE( tN2.is_flagged() );

    delete tElem;
}

// =============================================================================
// §3.2  Facet Topology (Representative Types)  [semantic]
// (idea from ChatGPT — complete per-facet node verification)
// =============================================================================

TEST( FacetTopology, TRI3FacetNodes )
{
    belfem::mesh::ElementFactory tFactory;
    belfem::mesh::Element * tElem = tFactory.create_element( belfem::ElementType::TRI3, 1 );

    belfem::mesh::Node tN0( 1 ), tN1( 2 ), tN2( 3 );
    tElem->insert_node( &tN0, 0 );
    tElem->insert_node( &tN1, 1 );
    tElem->insert_node( &tN2, 2 );

    belfem::Cell< belfem::mesh::Node * > tFacet;

    tElem->get_nodes_of_facet( 0, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN0 );
    EXPECT_EQ( tFacet( 1 ), &tN1 );

    tElem->get_nodes_of_facet( 1, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN1 );
    EXPECT_EQ( tFacet( 1 ), &tN2 );

    tElem->get_nodes_of_facet( 2, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN2 );
    EXPECT_EQ( tFacet( 1 ), &tN0 );

    delete tElem;
}

TEST( FacetTopology, QUAD4FacetNodes )
{
    belfem::mesh::ElementFactory tFactory;
    belfem::mesh::Element * tElem = tFactory.create_element( belfem::ElementType::QUAD4, 1 );

    belfem::mesh::Node tN0( 1 ), tN1( 2 ), tN2( 3 ), tN3( 4 );
    tElem->insert_node( &tN0, 0 );
    tElem->insert_node( &tN1, 1 );
    tElem->insert_node( &tN2, 2 );
    tElem->insert_node( &tN3, 3 );

    belfem::Cell< belfem::mesh::Node * > tFacet;

    tElem->get_nodes_of_facet( 0, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN0 );
    EXPECT_EQ( tFacet( 1 ), &tN1 );

    tElem->get_nodes_of_facet( 1, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN1 );
    EXPECT_EQ( tFacet( 1 ), &tN2 );

    tElem->get_nodes_of_facet( 2, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN2 );
    EXPECT_EQ( tFacet( 1 ), &tN3 );

    tElem->get_nodes_of_facet( 3, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN3 );
    EXPECT_EQ( tFacet( 1 ), &tN0 );

    delete tElem;
}

TEST( FacetTopology, TET4FacetNodes )
{
    belfem::mesh::ElementFactory tFactory;
    belfem::mesh::Element * tElem = tFactory.create_element( belfem::ElementType::TET4, 1 );

    belfem::mesh::Node tN0( 1 ), tN1( 2 ), tN2( 3 ), tN3( 4 );
    tElem->insert_node( &tN0, 0 );
    tElem->insert_node( &tN1, 1 );
    tElem->insert_node( &tN2, 2 );
    tElem->insert_node( &tN3, 3 );

    belfem::Cell< belfem::mesh::Node * > tFacet;

    // Exodus TET4 facet convention
    tElem->get_nodes_of_facet( 0, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN0 );
    EXPECT_EQ( tFacet( 1 ), &tN1 );
    EXPECT_EQ( tFacet( 2 ), &tN3 );

    tElem->get_nodes_of_facet( 1, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN1 );
    EXPECT_EQ( tFacet( 1 ), &tN2 );
    EXPECT_EQ( tFacet( 2 ), &tN3 );

    tElem->get_nodes_of_facet( 2, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN0 );
    EXPECT_EQ( tFacet( 1 ), &tN3 );
    EXPECT_EQ( tFacet( 2 ), &tN2 );

    tElem->get_nodes_of_facet( 3, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN0 );
    EXPECT_EQ( tFacet( 1 ), &tN2 );
    EXPECT_EQ( tFacet( 2 ), &tN1 );

    delete tElem;
}

TEST( FacetTopology, HEX8FacetNodes )
{
    belfem::mesh::ElementFactory tFactory;
    belfem::mesh::Element * tElem = tFactory.create_element( belfem::ElementType::HEX8, 1 );

    belfem::mesh::Node tN0( 1 ), tN1( 2 ), tN2( 3 ), tN3( 4 );
    belfem::mesh::Node tN4( 5 ), tN5( 6 ), tN6( 7 ), tN7( 8 );
    tElem->insert_node( &tN0, 0 );
    tElem->insert_node( &tN1, 1 );
    tElem->insert_node( &tN2, 2 );
    tElem->insert_node( &tN3, 3 );
    tElem->insert_node( &tN4, 4 );
    tElem->insert_node( &tN5, 5 );
    tElem->insert_node( &tN6, 6 );
    tElem->insert_node( &tN7, 7 );

    belfem::Cell< belfem::mesh::Node * > tFacet;

    // Exodus HEX8 facet convention
    tElem->get_nodes_of_facet( 0, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN0 );
    EXPECT_EQ( tFacet( 1 ), &tN1 );
    EXPECT_EQ( tFacet( 2 ), &tN5 );
    EXPECT_EQ( tFacet( 3 ), &tN4 );

    tElem->get_nodes_of_facet( 1, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN1 );
    EXPECT_EQ( tFacet( 1 ), &tN2 );
    EXPECT_EQ( tFacet( 2 ), &tN6 );
    EXPECT_EQ( tFacet( 3 ), &tN5 );

    tElem->get_nodes_of_facet( 2, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN2 );
    EXPECT_EQ( tFacet( 1 ), &tN3 );
    EXPECT_EQ( tFacet( 2 ), &tN7 );
    EXPECT_EQ( tFacet( 3 ), &tN6 );

    tElem->get_nodes_of_facet( 3, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN0 );
    EXPECT_EQ( tFacet( 1 ), &tN4 );
    EXPECT_EQ( tFacet( 2 ), &tN7 );
    EXPECT_EQ( tFacet( 3 ), &tN3 );

    tElem->get_nodes_of_facet( 4, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN0 );
    EXPECT_EQ( tFacet( 1 ), &tN3 );
    EXPECT_EQ( tFacet( 2 ), &tN2 );
    EXPECT_EQ( tFacet( 3 ), &tN1 );

    tElem->get_nodes_of_facet( 5, tFacet );
    EXPECT_EQ( tFacet( 0 ), &tN4 );
    EXPECT_EQ( tFacet( 1 ), &tN5 );
    EXPECT_EQ( tFacet( 2 ), &tN6 );
    EXPECT_EQ( tFacet( 3 ), &tN7 );

    delete tElem;
}

// =============================================================================
// §3.3  Element Base Class Guards  [semantic]
// (BELFEM_ERROR is always active — NOT debug-only)
// =============================================================================

TEST( ElementBaseGuard, BaseElementNodeThrows )
{
    // Construct base Element directly (not via factory)
    belfem::mesh::Element tElem( 1 );

    EXPECT_THROW( tElem.node( 0 ), std::runtime_error );
}

TEST( ElementBaseGuard, BaseElementTypeThrows )
{
    belfem::mesh::Element tElem( 1 );

    EXPECT_THROW( tElem.type(), std::runtime_error );
}

// =============================================================================
// §4  Basis / Vertex Counter Widths  [semantic]
// The narrow counters must refuse one past their width with a BELFEM_ERROR
// and leave the count untouched; the width below must be accepted in full.
// =============================================================================

#include "cl_Edge.hpp"
#include "cl_Face.hpp"
#include "cl_Facet.hpp"

namespace
{
    // inclusive maxima of the uint16_t and uint8_t counters under test
    const belfem::uint tMaxU16 = 65535 ;
    const belfem::uint tMaxU8  = 255 ;

    // Owns heap entities for the duration of a scope: Cell< T * > itself is
    // non-owning, so a fatal assertion or an unexpected throw would otherwise
    // skip the deletes. Empty slots stay nullptr and are safe to delete.
    template< typename T >
    struct ScopeOwner
    {
        belfem::Cell< T * > mItems ;

        explicit ScopeOwner( const belfem::uint aCount ) :
            mItems( aCount, nullptr )
        {
        }

        ~ScopeOwner()
        {
            for( belfem::uint k = 0; k < mItems.size(); ++k )
            {
                delete mItems( k );
            }
        }

        ScopeOwner( const ScopeOwner & ) = delete ;
        ScopeOwner( ScopeOwner && ) = delete ;
        ScopeOwner & operator=( const ScopeOwner & ) = delete ;
        ScopeOwner & operator=( ScopeOwner && ) = delete ;

        T * &
        operator()( const belfem::uint aIndex )
        {
            return mItems( aIndex );
        }

        belfem::uint
        size() const
        {
            return mItems.size();
        }
    };

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

    void
    create_source_nodes( ScopeOwner< belfem::mesh::Node > & aPool )
    {
        for( belfem::uint k = 0; k < aPool.size(); ++k )
        {
            aPool( k ) = new belfem::mesh::Node( k + 1 );
        }
    }

    // the first aCount pool entries, viewed as Basis pointers
    void
    view_as_bases( ScopeOwner< belfem::mesh::Node > & aPool,
                   const belfem::uint aCount,
                   belfem::Cell< belfem::mesh::Basis * > & aBases )
    {
        aBases.set_size( aCount, nullptr );
        for( belfem::uint k = 0; k < aCount; ++k )
        {
            aBases( k ) = aPool( k );
        }
    }

    void
    view_as_nodes( ScopeOwner< belfem::mesh::Node > & aPool,
                   const belfem::uint aCount,
                   belfem::Cell< belfem::mesh::Node * > & aSubset )
    {
        aSubset.set_size( aCount, nullptr );
        for( belfem::uint k = 0; k < aCount; ++k )
        {
            aSubset( k ) = aPool( k );
        }
    }
}

// --- §4.1 set_sources( Cell< Basis * >, Vector< real > ) ---

TEST( MeshEntities, SourcesBasisVectorOverloadCounterWidth )
{
    ScopeOwner< belfem::mesh::Node > tPool( tMaxU16 + 1 );
    create_source_nodes( tPool );

    belfem::Cell< belfem::mesh::Basis * > tSources ;

    {
        belfem::mesh::Node tNode( 1 );
        view_as_bases( tPool, 300, tSources );
        belfem::Vector< belfem::real > tWeights( 300, 1.0 );

        tNode.set_sources( tSources, tWeights );

        EXPECT_EQ( tNode.number_of_sources(), 300u );
        EXPECT_TRUE( tNode.is_hanging() );
        EXPECT_EQ( tNode.source( 299 ), tSources( 299 ) );
    }

    {
        belfem::mesh::Node tNode( 1 );
        view_as_bases( tPool, tMaxU16, tSources );
        belfem::Vector< belfem::real > tWeights( tMaxU16, 1.0 );

        tNode.set_sources( tSources, tWeights );

        EXPECT_EQ( tNode.number_of_sources(), tMaxU16 );
        EXPECT_EQ( tNode.source( tMaxU16 - 1 ), tSources( tMaxU16 - 1 ) );
    }

    {
        belfem::mesh::Node tNode( 1 );
        view_as_bases( tPool, tMaxU16 + 1, tSources );
        belfem::Vector< belfem::real > tWeights( tMaxU16 + 1, 1.0 );

        expect_refusal( [ & ]() { tNode.set_sources( tSources, tWeights ); },
                        "but the source counter holds at most" );

        EXPECT_EQ( tNode.number_of_sources(), 0u );
        EXPECT_FALSE( tNode.is_hanging() );
    }
}

// --- §4.2 set_sources( Cell< Basis * >, Cell< real > ) ---

TEST( MeshEntities, SourcesBasisCellOverloadCounterWidth )
{
    ScopeOwner< belfem::mesh::Node > tPool( tMaxU16 + 1 );
    create_source_nodes( tPool );

    belfem::Cell< belfem::mesh::Basis * > tSources ;

    {
        belfem::mesh::Node tNode( 1 );
        view_as_bases( tPool, 300, tSources );
        belfem::Cell< belfem::real > tWeights( 300, 1.0 );

        tNode.set_sources( tSources, tWeights );

        EXPECT_EQ( tNode.number_of_sources(), 300u );
        EXPECT_EQ( tNode.source( 299 ), tSources( 299 ) );
    }

    {
        belfem::mesh::Node tNode( 1 );
        view_as_bases( tPool, tMaxU16, tSources );
        belfem::Cell< belfem::real > tWeights( tMaxU16, 1.0 );

        tNode.set_sources( tSources, tWeights );

        EXPECT_EQ( tNode.number_of_sources(), tMaxU16 );
        EXPECT_EQ( tNode.source( tMaxU16 - 1 ), tSources( tMaxU16 - 1 ) );
    }

    {
        belfem::mesh::Node tNode( 1 );
        view_as_bases( tPool, tMaxU16 + 1, tSources );
        belfem::Cell< belfem::real > tWeights( tMaxU16 + 1, 1.0 );

        expect_refusal( [ & ]() { tNode.set_sources( tSources, tWeights ); },
                        "but the source counter holds at most" );

        EXPECT_EQ( tNode.number_of_sources(), 0u );
    }
}

// --- §4.3 set_sources( Cell< Node * >, Vector< real > ) ---

TEST( MeshEntities, SourcesNodeVectorOverloadCounterWidth )
{
    ScopeOwner< belfem::mesh::Node > tPool( tMaxU16 + 1 );
    create_source_nodes( tPool );

    belfem::Cell< belfem::mesh::Node * > tSources ;

    {
        belfem::mesh::Node tNode( 1 );
        view_as_nodes( tPool, 300, tSources );
        belfem::Vector< belfem::real > tWeights( 300, 1.0 );

        tNode.set_sources( tSources, tWeights );

        EXPECT_EQ( tNode.number_of_sources(), 300u );
        EXPECT_EQ( tNode.source_node( 299 ), tSources( 299 ) );
    }

    {
        belfem::mesh::Node tNode( 1 );
        view_as_nodes( tPool, tMaxU16, tSources );
        belfem::Vector< belfem::real > tWeights( tMaxU16, 1.0 );

        tNode.set_sources( tSources, tWeights );

        EXPECT_EQ( tNode.number_of_sources(), tMaxU16 );
        EXPECT_EQ( tNode.source_node( tMaxU16 - 1 ), tSources( tMaxU16 - 1 ) );
    }

    {
        belfem::mesh::Node tNode( 1 );
        view_as_nodes( tPool, tMaxU16 + 1, tSources );
        belfem::Vector< belfem::real > tWeights( tMaxU16 + 1, 1.0 );

        expect_refusal( [ & ]() { tNode.set_sources( tSources, tWeights ); },
                        "but the source counter holds at most" );

        EXPECT_EQ( tNode.number_of_sources(), 0u );
    }
}

// --- §4.4 allocate_source_container() + add_source() ---

TEST( MeshEntities, SourcesIncrementalPathCounterWidth )
{
    ScopeOwner< belfem::mesh::Node > tPool( tMaxU16 );
    create_source_nodes( tPool );

    {
        belfem::mesh::Node tNode( 1 );

        expect_refusal( [ & ]() { tNode.allocate_source_container( tMaxU16 + 1 ); },
                        "but the source counter holds at most" );

        EXPECT_EQ( tNode.number_of_sources(), 0u );
        EXPECT_FALSE( tNode.is_hanging() );
    }

    {
        belfem::mesh::Node tNode( 1 );

        tNode.allocate_source_container( 3 );
        tNode.add_source( tPool( 0 ), 0.25 );
        tNode.add_source( tPool( 1 ), 0.50 );
        tNode.add_source( tPool( 2 ), 0.75 );

        EXPECT_EQ( tNode.number_of_sources(), 3u );
        EXPECT_TRUE( tNode.is_hanging() );
        EXPECT_EQ( tNode.source( 2 ), tPool( 2 ) );
        EXPECT_NEAR( tNode.weight( 2 ), 0.75, tEps );
    }

    {
        // the append guard is the counter width itself
        belfem::mesh::Node tNode( 1 );

        tNode.allocate_source_container( tMaxU16 );
        for( belfem::uint k = 0; k < tMaxU16; ++k )
        {
            tNode.add_source( tPool( k ) );
        }
        EXPECT_EQ( tNode.number_of_sources(), tMaxU16 );

        expect_refusal( [ & ]() { tNode.add_source( tPool( 0 ) ); },
                        "Source counter of basis" );

        EXPECT_EQ( tNode.number_of_sources(), tMaxU16 );
        EXPECT_EQ( tNode.source( tMaxU16 - 1 ), tPool( tMaxU16 - 1 ) );
    }
}

// --- §4.5 dof counter and container ---

TEST( MeshEntities, DofContainerAllocateEmptyThenReset )
{
    belfem::mesh::Node tNode( 1 );

    tNode.allocate_dof_container();
    EXPECT_EQ( tNode.number_of_dofs(), 0u );

    tNode.reset_dof_container();
    EXPECT_EQ( tNode.number_of_dofs(), 0u );
}

TEST( MeshEntities, DofCounterRefusesAtWidth )
{
    ScopeOwner< belfem::graph::Vertex > tDofs( tMaxU8 );
    for( belfem::uint k = 0; k < tMaxU8; ++k )
    {
        tDofs( k ) = new belfem::graph::Vertex ;
        tDofs( k )->set_id( k + 1 );
    }

    {
        belfem::mesh::Node tNode( 1 );

        for( belfem::uint k = 0; k < tMaxU8; ++k )
        {
            tNode.increment_dof_counter();
        }

        // fatal on purpose: a wrapped counter would make the allocate below
        // a no-op and the first insert a write through a null pointer
        ASSERT_EQ( tNode.number_of_dofs(), tMaxU8 );

        expect_refusal( [ & ]() { tNode.increment_dof_counter(); },
                        "Dof counter of basis" );
        ASSERT_EQ( tNode.number_of_dofs(), tMaxU8 );

        tNode.allocate_dof_container();
        ASSERT_EQ( tNode.number_of_dofs(), 0u );

        for( belfem::uint k = 0; k < tMaxU8; ++k )
        {
            tNode.insert_dof( tDofs( k ) );
        }
        EXPECT_EQ( tNode.number_of_dofs(), tMaxU8 );
        EXPECT_EQ( tNode.dof( 0 ), tDofs( 0 ) );
        EXPECT_EQ( tNode.dof( tMaxU8 - 1 ), tDofs( tMaxU8 - 1 ) );

        expect_refusal( [ & ]() { tNode.insert_dof( tDofs( 0 ) ); },
                        "Dof counter of basis" );
        EXPECT_EQ( tNode.number_of_dofs(), tMaxU8 );
        EXPECT_EQ( tNode.dof( tMaxU8 - 1 ), tDofs( tMaxU8 - 1 ) );
    }
}

TEST( MeshEntities, DofContainerResetWhileUnfilled )
{
    // No-throw smoke case: the reset is keyed on the pointer, but allocate
    // has already zeroed the count and ~Basis frees the block anyway, so the
    // free itself has no observable oracle here, as in the reallocate case.
    belfem::mesh::Node tNode( 1 );

    tNode.increment_dof_counter();
    tNode.increment_dof_counter();
    tNode.increment_dof_counter();
    EXPECT_EQ( tNode.number_of_dofs(), 3u );

    tNode.allocate_dof_container();
    EXPECT_EQ( tNode.number_of_dofs(), 0u );

    tNode.reset_dof_container();
    EXPECT_EQ( tNode.number_of_dofs(), 0u );
}

TEST( MeshEntities, DofContainerReallocateAfterFreshCount )
{
    // A second allocate while the first block is still owned must free it
    // before it mallocs again. That free has no observable oracle here; this
    // case is smoke only, its value is a clean run under a leak checker.
    ScopeOwner< belfem::graph::Vertex > tDofs( 2 );
    tDofs( 0 ) = new belfem::graph::Vertex ;
    tDofs( 1 ) = new belfem::graph::Vertex ;

    {
        belfem::mesh::Node tNode( 1 );

        tNode.increment_dof_counter();
        tNode.increment_dof_counter();
        tNode.allocate_dof_container();
        ASSERT_EQ( tNode.number_of_dofs(), 0u );

        // fatal on purpose: a count that is not reset would carry the first
        // two into the second allocate, and the inserts below would write
        // past a two-slot block
        tNode.increment_dof_counter();
        tNode.increment_dof_counter();
        ASSERT_EQ( tNode.number_of_dofs(), 2u );

        tNode.allocate_dof_container();
        ASSERT_EQ( tNode.number_of_dofs(), 0u );

        tNode.insert_dof( tDofs( 0 ) );
        tNode.insert_dof( tDofs( 1 ) );
        EXPECT_EQ( tNode.number_of_dofs(), 2u );
        EXPECT_EQ( tNode.dof( 1 ), tDofs( 1 ) );
    }
}

// --- §4.6 adjacency containers: capacity width ---
// The "already allocated" precondition is an assert, so every refusal and
// every acceptance runs on its own fresh node. Three checks lock the
// acceptance side. The memory() delta across the allocate equals the stored
// capacity times the pointer size ( array_memory() charges the capacity
// fields, not the fill counts ), so a dropped capacity store goes red in
// release, where add_*() is assert-guarded and number_of_*() is only the
// counter. The fill to capacity through add_*() locks the malloc. The count
// read back locks the counter. After a refusal, memory() and the counter are
// re-read as untouched.

TEST( MeshEntities, NodeContainerCapacityWidth )
{
    {
        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        expect_refusal( [ & ]() { tNode.allocate_node_container( tMaxU8 + 1 ); },
                        "node capacity of vertex" );
        EXPECT_EQ( tNode.memory(), tMemory );
        EXPECT_EQ( tNode.number_of_nodes(), 0u );
    }
    {
        ScopeOwner< belfem::mesh::Node > tEntities( tMaxU8 );
        for( belfem::uint k = 0; k < tMaxU8; ++k )
        {
            tEntities( k ) = new belfem::mesh::Node( k + 2 );
        }

        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        ASSERT_NO_THROW( tNode.allocate_node_container( tMaxU8 ) );
        EXPECT_EQ( tNode.memory() - tMemory, tMaxU8 * sizeof( belfem::mesh::Node * ) );
        for( belfem::uint k = 0; k < tMaxU8; ++k )
        {
            tNode.add_node( tEntities( k ) );
        }
        EXPECT_EQ( tNode.number_of_nodes(), tMaxU8 );
        EXPECT_EQ( tNode.node( tMaxU8 - 1 ), tEntities( tMaxU8 - 1 ) );
    }
}

TEST( MeshEntities, EdgeContainerCapacityWidth )
{
    {
        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        expect_refusal( [ & ]() { tNode.allocate_edge_container( tMaxU8 + 1 ); },
                        "edge capacity of vertex" );
        EXPECT_EQ( tNode.memory(), tMemory );
        EXPECT_EQ( tNode.number_of_edges(), 0u );
    }
    {
        ScopeOwner< belfem::mesh::Edge > tEntities( tMaxU8 );
        for( belfem::uint k = 0; k < tMaxU8; ++k )
        {
            tEntities( k ) = new belfem::mesh::Edge ;
        }

        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        ASSERT_NO_THROW( tNode.allocate_edge_container( tMaxU8 ) );
        EXPECT_EQ( tNode.memory() - tMemory, tMaxU8 * sizeof( belfem::mesh::Edge * ) );
        for( belfem::uint k = 0; k < tMaxU8; ++k )
        {
            tNode.add_edge( tEntities( k ) );
        }
        EXPECT_EQ( tNode.number_of_edges(), tMaxU8 );
        EXPECT_EQ( tNode.edge( tMaxU8 - 1 ), tEntities( tMaxU8 - 1 ) );
    }
}

TEST( MeshEntities, FaceContainerCapacityWidth )
{
    {
        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        expect_refusal( [ & ]() { tNode.allocate_face_container( tMaxU8 + 1 ); },
                        "face capacity of vertex" );
        EXPECT_EQ( tNode.memory(), tMemory );
        EXPECT_EQ( tNode.number_of_faces(), 0u );
    }
    {
        // a Face copies its parent's nodes, so it needs a concrete element;
        // all faces share one TRI3 parent, which they do not own
        belfem::mesh::ElementFactory tFactory;
        belfem::mesh::Node tN0( 2 ), tN1( 3 ), tN2( 4 );
        ScopeOwner< belfem::mesh::Element > tParentOwner( 1 );
        tParentOwner( 0 ) = tFactory.create_element( belfem::ElementType::TRI3, 1 );
        belfem::mesh::Element * tParent = tParentOwner( 0 );
        tParent->insert_node( &tN0, 0 );
        tParent->insert_node( &tN1, 1 );
        tParent->insert_node( &tN2, 2 );

        {
            ScopeOwner< belfem::mesh::Face > tEntities( tMaxU8 );
            for( belfem::uint k = 0; k < tMaxU8; ++k )
            {
                tEntities( k ) = new belfem::mesh::Face( tParent );
            }

            belfem::mesh::Node tNode( 1 );
            const size_t tMemory = tNode.memory();
            ASSERT_NO_THROW( tNode.allocate_face_container( tMaxU8 ) );
            EXPECT_EQ( tNode.memory() - tMemory, tMaxU8 * sizeof( belfem::mesh::Face * ) );
            for( belfem::uint k = 0; k < tMaxU8; ++k )
            {
                tNode.add_face( tEntities( k ) );
            }
            EXPECT_EQ( tNode.number_of_faces(), tMaxU8 );
            EXPECT_EQ( tNode.face( tMaxU8 - 1 ), tEntities( tMaxU8 - 1 ) );
        }
    }
}

TEST( MeshEntities, FacetContainerCapacityWidth )
{
    // the facet capacity is the one 16-bit field among the five
    {
        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        expect_refusal( [ & ]() { tNode.allocate_facet_container( tMaxU16 + 1 ); },
                        "facet capacity of vertex" );
        EXPECT_EQ( tNode.memory(), tMemory );
        EXPECT_EQ( tNode.number_of_facets(), 0u );
    }
    {
        // a Facet owns and deletes its element, so each gets a bare one
        ScopeOwner< belfem::mesh::Facet > tEntities( tMaxU16 );
        for( belfem::uint k = 0; k < tMaxU16; ++k )
        {
            tEntities( k ) = new belfem::mesh::Facet( new belfem::mesh::Element( k + 1 ) );
        }

        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        ASSERT_NO_THROW( tNode.allocate_facet_container( tMaxU16 ) );
        EXPECT_EQ( tNode.memory() - tMemory, tMaxU16 * sizeof( belfem::mesh::Facet * ) );
        for( belfem::uint k = 0; k < tMaxU16; ++k )
        {
            tNode.add_facet( tEntities( k ) );
        }
        EXPECT_EQ( tNode.number_of_facets(), tMaxU16 );
        EXPECT_EQ( tNode.facet( tMaxU16 - 1 ), tEntities( tMaxU16 - 1 ) );
    }
}

TEST( MeshEntities, ElementContainerCapacityWidth )
{
    {
        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        expect_refusal( [ & ]() { tNode.allocate_element_container( tMaxU8 + 1 ); },
                        "element capacity of vertex" );
        EXPECT_EQ( tNode.memory(), tMemory );
        EXPECT_EQ( tNode.number_of_elements(), 0u );
    }
    {
        ScopeOwner< belfem::mesh::Element > tEntities( tMaxU8 );
        for( belfem::uint k = 0; k < tMaxU8; ++k )
        {
            tEntities( k ) = new belfem::mesh::Element( k + 1 );
        }

        belfem::mesh::Node tNode( 1 );
        const size_t tMemory = tNode.memory();
        ASSERT_NO_THROW( tNode.allocate_element_container( tMaxU8 ) );
        EXPECT_EQ( tNode.memory() - tMemory, tMaxU8 * sizeof( belfem::mesh::Element * ) );
        for( belfem::uint k = 0; k < tMaxU8; ++k )
        {
            tNode.add_element( tEntities( k ) );
        }
        EXPECT_EQ( tNode.number_of_elements(), tMaxU8 );
        EXPECT_EQ( tNode.element( tMaxU8 - 1 ), tEntities( tMaxU8 - 1 ) );
    }
}
