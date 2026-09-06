/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Cohomology-layer tests on a programmatic annulus.
 *
 * The fixture drives the production call sequence of
 * CutFactory::compute_cohomologies with mSuggestHomologies on — flag the
 * region, build the SimplicialComplex, COREDUCE ONLY (the reduce_complex*
 * calls sit in the mSuggestHomologies==false branch, which the hardcoded
 * default never takes), construct Cohomology, replace its generators from
 * the homology, clean() — so the full clean_spfa / rectify / pocket-census
 * path runs twice, once inside the constructor and once after the update,
 * exactly as in production. Scope: this pins the Cohomology layer, not
 * CutFactory's maxwell orchestration around it.
 *
 * The one mesh-construction trap, learned the hard way: the Mesh must be
 * built with aComputeConnectivities = true (the default), or finalize_edges
 * skips connect_edges_to_elements, every edge cochain gets an empty
 * coboundary, the coreduction cannot collapse a single copair, and the
 * cohomology SNF degenerates to one generator per edge.
 */

#include <gtest/gtest.h>
#include <map>

#include "constants.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Block.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Homology.hpp"
#include "cl_Cohomology.hpp"
#include "cl_Chain.hpp"
#include "cl_Cochain.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    // one ring of elements between radii aR0 < aR1, counterclockwise.
    // node ids: 1..aN inner ring, aN+1..2*aN outer ring, both starting at
    // polar angle zero. aQuads selects one QUAD4 per sector or two TRI3.
    // The caller owns the mesh.
    Mesh *
    create_annulus( const uint aN, const bool aQuads,
        const real aR0 = 1.0, const real aR1 = 2.0 )
    {
        Mesh * tMesh = new Mesh( 2, 0, true );

        Cell< Node * > & tNodes = tMesh->nodes();

        for ( uint k = 0; k < aN; ++k )
        {
            const real tPhi = 2.0 * constant::pi * k / aN ;
            tNodes.push( new Node( k + 1,
                aR0 * std::cos( tPhi ), aR0 * std::sin( tPhi ) ) );
        }
        for ( uint k = 0; k < aN; ++k )
        {
            const real tPhi = 2.0 * constant::pi * k / aN ;
            tNodes.push( new Node( aN + k + 1,
                aR1 * std::cos( tPhi ), aR1 * std::sin( tPhi ) ) );
        }

        Block * tBlock = new Block( 1, aQuads ? aN : 2 * aN );

        ElementFactory tFactory ;
        Cell< Element * > & tElements = tMesh->elements();

        id_t tID = 1 ;
        for ( uint k = 0; k < aN; ++k )
        {
            const uint tNext = ( k + 1 ) % aN ;

            if ( aQuads )
            {
                // counterclockwise: inner k, outer k, outer k+1, inner k+1
                Element * tElement =
                    tFactory.create_element( ElementType::QUAD4, tID++ );
                tElement->insert_node( tNodes( k ),          0 );
                tElement->insert_node( tNodes( aN + k ),     1 );
                tElement->insert_node( tNodes( aN + tNext ), 2 );
                tElement->insert_node( tNodes( tNext ),      3 );
                tElement->set_block_id( 1 );
                tElements.push( tElement );
                tBlock->insert_element( tElement );
            }
            else
            {
                // two counterclockwise TRI3 per sector
                Element * tElement =
                    tFactory.create_element( ElementType::TRI3, tID++ );
                tElement->insert_node( tNodes( k ),          0 );
                tElement->insert_node( tNodes( aN + k ),     1 );
                tElement->insert_node( tNodes( aN + tNext ), 2 );
                tElement->set_block_id( 1 );
                tElements.push( tElement );
                tBlock->insert_element( tElement );

                Element * tElement2 =
                    tFactory.create_element( ElementType::TRI3, tID++ );
                tElement2->insert_node( tNodes( k ),          0 );
                tElement2->insert_node( tNodes( aN + tNext ), 1 );
                tElement2->insert_node( tNodes( tNext ),      2 );
                tElement2->set_block_id( 1 );
                tElements.push( tElement2 );
                tBlock->insert_element( tElement2 );
            }
        }

        tMesh->blocks().push( tBlock );
        tMesh->finalize();
        tMesh->create_edges( false );

        return tMesh ;
    }

    enum class Algorithm { Pellikka, CCR, PellikkaGeneralized };

    // the CutFactory::compute_cohomologies 2D recipe with
    // mSuggestHomologies on, minus the maxwell dressing
    void
    flag_region( Mesh * aMesh )
    {
        aMesh->update_node_indices();
        aMesh->update_edge_indices();
        aMesh->update_face_indices();

        aMesh->unflag_everything();

        Block * tBlock = aMesh->block( 1 );
        tBlock->flag_elements();
        for ( Element * tElement : tBlock->elements() )
        {
            tElement->flag_corner_nodes();
            tElement->flag_edges();
        }
    }

    // pairing of a 1-cochain with the inner-ring cycle, traversed
    // counterclockwise ( node ids 1..aN in ascending order ). For an annulus
    // the winding of the H^1 generator around the hole must be +-1 —
    // rectification may change the representative, never this pairing.
    int
    inner_ring_pairing( Mesh * aMesh, Cochain * aGenerator, const uint aN )
    {
        std::map< index_t, int > tCoeffs ;
        for ( const auto & [ tIndex, tCoeff ] : aGenerator->getSimplicesMap() )
        {
            tCoeffs[ tIndex ] = tCoeff ;
        }

        int tPairing = 0 ;

        for ( Edge * tEdge : aMesh->edges() )
        {
            const id_t tA = tEdge->node( 0 )->id();
            const id_t tB = tEdge->node( 1 )->id();

            // inner-ring edge: both endpoints have ids in 1..aN
            if ( tA > aN || tB > aN ) continue ;

            auto tIt = tCoeffs.find( tEdge->index() );
            if ( tIt == tCoeffs.end() ) continue ;

            // ccw direction is ring position a -> a+1 ( mod aN )
            const int tSign =
                ( tB % aN ) == ( tA % aN + 1 ) % aN ? 1 : -1 ;

            tPairing += tSign * tIt->second ;
        }

        return tPairing ;
    }

    void
    check_annulus( const Algorithm aAlgorithm, const bool aQuads )
    {
        const uint tN = 8 ;

        Mesh * tMesh = create_annulus( tN, aQuads );

        flag_region( tMesh );

        SimplicialComplex * tComplex = new SimplicialComplex( tMesh, true );

        switch ( aAlgorithm )
        {
            case Algorithm::Pellikka :
            {
                tComplex->coreduce_complexPellikka();
                break ;
            }
            case Algorithm::CCR :
            {
                tComplex->coreduce_complexCCR();
                break ;
            }
            case Algorithm::PellikkaGeneralized :
            {
                tComplex->coreduce_complexPellikkaGeneralized();
                break ;
            }
        }

        // chain-side SNF homology: betti numbers of the annulus
        Homology * tHomology = new Homology( tComplex, tMesh );
        ASSERT_EQ( tHomology->get_Generators()( 0 ).size(), 1u );
        ASSERT_EQ( tHomology->get_Generators()( 1 ).size(), 1u );
        ASSERT_EQ( tHomology->get_Generators()( 2 ).size(), 0u );

        // cochain-side SNF on the coreduced complex ( constructor runs
        // generatorsOfCohomology + clean, i.e. clean_spfa )
        Cohomology * tCohomology = new Cohomology( tComplex, tMesh );
        ASSERT_EQ( tCohomology->get_Generators()( 1 ).size(), 1u );

        // production sequence: replace the generators from the homology,
        // then rectify again
        tCohomology->updatekGeneratorsFromHomology(
            tHomology->get_Generators()( 1 ), 1 );
        tCohomology->clean();

        Cell< Cochain * > & tGenerators = tCohomology->get_Generators()( 1 );
        ASSERT_EQ( tGenerators.size(), 1u );

        // the rectification contract: every surviving coefficient is a unit
        for ( const auto & [ tIndex, tCoeff ] : tGenerators( 0 )->getSimplicesMap() )
        {
            EXPECT_TRUE( tCoeff == 1 || tCoeff == -1 )
                << "edge index " << tIndex << " carries coefficient " << tCoeff ;
        }

        // the equivalence class survives cleaning: winding around the hole
        EXPECT_EQ( std::abs( inner_ring_pairing(
            tMesh, tGenerators( 0 ), tN ) ), 1 );

        delete tCohomology ;
        delete tHomology ;
        delete tComplex ;
        delete tMesh ;
    }
}

TEST( Cohomology, AnnulusPellikkaTri )
{
    check_annulus( Algorithm::Pellikka, false );
}

// layer-tier forced repair: the six cases below always hand clean() an
// already-unit generator, so the rectification path runs but never repairs at
// this tier. Here the cleaned generator is perturbed by an EXACT node
// coboundary — the same cohomology class, so clean() must recover a unit
// representative with the winding pairing unchanged — and the perturbation
// sign is chosen so a genuine +-2 coefficient exists before the repair.
TEST( Cohomology, AnnulusForcedRepair )
{
    const uint tN = 8 ;

    Mesh * tMesh = create_annulus( tN, false );

    flag_region( tMesh );

    SimplicialComplex * tComplex = new SimplicialComplex( tMesh, true );
    tComplex->coreduce_complexPellikka();

    Homology   * tHomology   = new Homology( tComplex, tMesh );
    Cohomology * tCohomology = new Cohomology( tComplex, tMesh );

    tCohomology->updatekGeneratorsFromHomology(
        tHomology->get_Generators()( 1 ), 1 );
    tCohomology->clean();

    Cell< Cochain * > & tGenerators = tCohomology->get_Generators()( 1 );
    ASSERT_EQ( tGenerators.size(), 1u );
    Cochain * tGen = tGenerators( 0 );

    // first support edge by smallest edge index ( OrderedMap = std::map,
    // deterministic ); its coefficient is a unit after the clean above
    ASSERT_GT( tGen->getSimplicesMap().size(), 0u );
    const index_t tEdgeIndex = tGen->getSimplicesMap().begin()->first ;
    const int     tG0        = tGen->getSimplicesMap().begin()->second ;
    ASSERT_TRUE( tG0 == 1 || tG0 == -1 );

    Edge * tSupportEdge = nullptr ;
    for ( Edge * tEdge : tMesh->edges() )
    {
        if ( tEdge->index() == tEdgeIndex )
        {
            tSupportEdge = tEdge ;
            break ;
        }
    }
    ASSERT_NE( tSupportEdge, nullptr );

    // perturbation node n = node(0) of the support edge; the complex's
    // incidence convention gives delta(n)(e) = -1 at e.node(0), +1 at
    // e.node(1), so s = -g0 turns the support coefficient into 2*g0
    const index_t tNodeIndex = tSupportEdge->node( 0 )->index();
    const int     tS         = -tG0 ;

    const DynamicBitset & tDomain = tComplex->original_edges();
    for ( Edge * tEdge : tMesh->edges() )
    {
        if ( ! tDomain.test( tEdge->index() ) ) continue ;

        int tDelta = 0 ;
        if      ( tEdge->node( 0 )->index() == tNodeIndex ) tDelta = -1 ;
        else if ( tEdge->node( 1 )->index() == tNodeIndex ) tDelta =  1 ;
        if ( tDelta == 0 ) continue ;

        tGen->setCoefficient( tEdge->index(),
            tGen->getCoefficient( tEdge->index() ) + tS * tDelta );
    }

    // the forcing is real: a non-unit coefficient exists before the repair
    bool tHaveNonUnit = false ;
    for ( const auto & [ tIndex, tCoeff ] : tGen->getSimplicesMap() )
    {
        if ( tCoeff >= 2 || tCoeff <= -2 ) tHaveNonUnit = true ;
    }
    ASSERT_TRUE( tHaveNonUnit );

    // the layer-tier repair under test
    tCohomology->clean();

    // rectification contract: unit representative of the SAME class
    ASSERT_EQ( tGenerators.size(), 1u );
    for ( const auto & [ tIndex, tCoeff ] : tGenerators( 0 )->getSimplicesMap() )
    {
        EXPECT_TRUE( tCoeff == 1 || tCoeff == -1 )
            << "edge index " << tIndex << " carries coefficient " << tCoeff ;
    }
    EXPECT_EQ( std::abs( inner_ring_pairing(
        tMesh, tGenerators( 0 ), tN ) ), 1 );

    delete tCohomology ;
    delete tHomology ;
    delete tComplex ;
    delete tMesh ;
}

TEST( Cohomology, AnnulusPellikkaQuad )
{
    check_annulus( Algorithm::Pellikka, true );
}

TEST( Cohomology, AnnulusCCRTri )
{
    check_annulus( Algorithm::CCR, false );
}

TEST( Cohomology, AnnulusCCRQuad )
{
    check_annulus( Algorithm::CCR, true );
}

TEST( Cohomology, AnnulusPellikkaGeneralizedTri )
{
    check_annulus( Algorithm::PellikkaGeneralized, false );
}

TEST( Cohomology, AnnulusPellikkaGeneralizedQuad )
{
    check_annulus( Algorithm::PellikkaGeneralized, true );
}
