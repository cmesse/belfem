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
