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
#include <memory>

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
    // aWithConductor fills the hole with block 2: node 2*aN+1 at the origin
    // and one counterclockwise TRI3 ( centre, inner k, inner k+1 ) per
    // sector, a bulk conductor for the suggested homology to bound.
    // The caller owns the mesh.
    Mesh *
    create_annulus( const uint aN, const bool aQuads,
        const real aR0 = 1.0, const real aR1 = 2.0,
        const bool aWithConductor = false )
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

        // the fan must be in place before finalize() and create_edges():
        // element indices and edge directions are assigned there
        if ( aWithConductor )
        {
            Node * tCentre = new Node( 2 * aN + 1, 0.0, 0.0 );
            tNodes.push( tCentre );

            Block * tFan = new Block( 2, aN );
            for ( uint k = 0; k < aN; ++k )
            {
                const uint tNext = ( k + 1 ) % aN ;

                Element * tElement =
                    tFactory.create_element( ElementType::TRI3, tID++ );
                tElement->insert_node( tCentre,         0 );
                tElement->insert_node( tNodes( k ),     1 );
                tElement->insert_node( tNodes( tNext ), 2 );
                tElement->set_block_id( 2 );
                tElements.push( tElement );
                tFan->insert_element( tElement );
            }
            tMesh->blocks().push( tFan );
        }

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
        for ( const auto & [ tIndex, tCoeff ] : aGenerator->get_simplices_map() )
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

    // the sense of the suggested 2-D generator on the inner ring: +1 is
    // counterclockwise ( node ids ascending ). This pins what the tree
    // produces, not the physical convention: the core's own comment names the
    // clockwise sense for +z, and which of the two holds is the homology
    // owner's call. A different value is a convention change that the owner
    // has to confirm, not a test to flip.
    constexpr int gInnerRingSense = 1 ;

    // +1 if the stored inner-ring edge aA -> aB runs counterclockwise
    // ( ring position a -> a+1 mod aN ), -1 otherwise. Read from node ids so
    // the oracle is independent of edge_direction.
    int
    inner_ring_sign( const id_t aA, const id_t aB, const uint aN )
    {
        return ( aB % aN ) == ( aA % aN + 1 ) % aN ? 1 : -1 ;
    }

    // premise of the sense tests: every fan element has positive signed
    // area, so its boundary is the counterclockwise ring
    bool
    fan_is_counterclockwise( Mesh * aMesh, const uint aN )
    {
        Block * tFan = aMesh->block( 2 );
        if ( tFan->number_of_elements() != aN ) return false ;

        for ( Element * tElement : tFan->elements() )
        {
            const real tX0 = tElement->node( 0 )->x();
            const real tY0 = tElement->node( 0 )->y();
            const real tX1 = tElement->node( 1 )->x();
            const real tY1 = tElement->node( 1 )->y();
            const real tX2 = tElement->node( 2 )->x();
            const real tY2 = tElement->node( 2 )->y();

            const real tTwiceArea =
                ( tX1 - tX0 ) * ( tY2 - tY0 ) - ( tX2 - tX0 ) * ( tY1 - tY0 );
            if ( ! ( tTwiceArea > 0.0 ) ) return false ;
        }
        return true ;
    }

    // the production prelude of CutFactory::compute_cohomologies on the
    // conductor fixture: index refresh, the suggesting Homology on the
    // { in, out } terminal pair a 2-D bulk condition sends for block 2, then
    // the reorientation. Must run BEFORE flag_region: suggest_Homology
    // unflags everything on entry and on exit. The caller owns the result.
    Homology *
    suggest_conductor_loop( Mesh * aMesh )
    {
        aMesh->update_node_indices();
        aMesh->update_edge_indices();
        aMesh->update_face_indices();

        Cell< Cell< id_t > > tTerminals = { { 2, 2 } };
        Cell< id_t > tThinShellIndices ;

        Homology * tSuggested =
            new Homology( aMesh, tTerminals, tThinShellIndices );
        tSuggested->reorient_generators();

        return tSuggested ;
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
        for ( const auto & [ tIndex, tCoeff ] : tGenerators( 0 )->get_simplices_map() )
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
    ASSERT_GT( tGen->get_simplices_map().size(), 0u );
    const index_t tEdgeIndex = tGen->get_simplices_map().begin()->first ;
    const int     tG0        = tGen->get_simplices_map().begin()->second ;
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

        tGen->set_coefficient( tEdge->index(),
            tGen->get_coefficient( tEdge->index() ) + tS * tDelta );
    }

    // the forcing is real: a non-unit coefficient exists before the repair
    bool tHaveNonUnit = false ;
    for ( const auto & [ tIndex, tCoeff ] : tGen->get_simplices_map() )
    {
        if ( tCoeff >= 2 || tCoeff <= -2 ) tHaveNonUnit = true ;
    }
    ASSERT_TRUE( tHaveNonUnit );

    // the layer-tier repair under test
    tCohomology->clean();

    // rectification contract: unit representative of the SAME class
    ASSERT_EQ( tGenerators.size(), 1u );
    for ( const auto & [ tIndex, tCoeff ] : tGenerators( 0 )->get_simplices_map() )
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

// the SIGNED sense of the suggested 2-D generator: gamma is the inner ring
// with coefficient gInnerRingSense * ( +1 on a stored edge that runs
// counterclockwise, -1 otherwise ), and nothing else survives — the radial
// fan edges cancel between neighbours.
TEST( Cohomology, SuggestedGeneratorIsTheCcwInnerRing )
{
    const uint tN = 8 ;

    std::unique_ptr< Mesh > tMesh( create_annulus( tN, false, 1.0, 2.0, true ) );
    ASSERT_TRUE( fan_is_counterclockwise( tMesh.get(), tN ) );

    std::unique_ptr< Homology > tSuggested( suggest_conductor_loop( tMesh.get() ) );
    ASSERT_EQ( tSuggested->get_Generators()( 1 ).size(), 1u );
    Chain * tGamma = tSuggested->get_Generators()( 1 )( 0 );

    // support is exactly the inner ring: one entry per ring edge
    ASSERT_EQ( tGamma->get_simplices_map().size(), tN );

    uint tRingEdges = 0 ;
    for ( Edge * tEdge : tMesh->edges() )
    {
        const id_t tA = tEdge->node( 0 )->id();
        const id_t tB = tEdge->node( 1 )->id();
        if ( tA > tN || tB > tN ) continue ;
        ++tRingEdges ;

        EXPECT_EQ( tGamma->get_coefficient( tEdge->index() ),
            gInnerRingSense * inner_ring_sign( tA, tB, tN ) )
            << "inner-ring edge " << tA << " -> " << tB ;
    }
    EXPECT_EQ( tRingEdges, tN );
}

// the production handoff driven with the SUGGESTED chain: after
// updatekGeneratorsFromHomology + clean the cut cochain pairs with the
// counterclockwise ring to gInnerRingSense, and, independently of that
// convention, with gamma itself to +1 ( Pellikka et al. 2013, Eq. 5.3 ).
TEST( Cohomology, SuggestedGeneratorPairsWithTheCleanedCochain )
{
    const uint tN = 8 ;

    std::unique_ptr< Mesh > tMesh( create_annulus( tN, false, 1.0, 2.0, true ) );
    ASSERT_TRUE( fan_is_counterclockwise( tMesh.get(), tN ) );

    std::unique_ptr< Homology > tSuggested( suggest_conductor_loop( tMesh.get() ) );
    ASSERT_EQ( tSuggested->get_Generators()( 1 ).size(), 1u );
    Chain * tGamma = tSuggested->get_Generators()( 1 )( 0 );

    // the phi region is the annulus alone; the edge list is untouched, so
    // the second index refresh inside leaves gamma's edge indices valid
    flag_region( tMesh.get() );

    std::unique_ptr< SimplicialComplex > tComplex(
        new SimplicialComplex( tMesh.get(), true ) );
    tComplex->coreduce_complexPellikkaGeneralized();

    std::unique_ptr< Cohomology > tCohomology(
        new Cohomology( tComplex.get(), tMesh.get() ) );
    ASSERT_EQ( tCohomology->get_Generators()( 1 ).size(), 1u );

    // unit coupling: the update multiplies by F^T, which inverts F only
    // when its single invariant factor is 1
    const int tCoupling =
        ( *tCohomology->get_Generators()( 1 )( 0 ) )( tGamma );
    ASSERT_TRUE( tCoupling == 1 || tCoupling == -1 )
        << "< c, gamma > = " << tCoupling ;

    tCohomology->updatekGeneratorsFromHomology(
        tSuggested->get_Generators()( 1 ), 1 );
    tCohomology->clean();

    // the update deletes the old generators: reacquire the pointer
    Cell< Cochain * > & tGenerators = tCohomology->get_Generators()( 1 );
    ASSERT_EQ( tGenerators.size(), 1u );
    Cochain * tCut = tGenerators( 0 );

    EXPECT_EQ( inner_ring_pairing( tMesh.get(), tCut, tN ), gInnerRingSense );
    EXPECT_EQ( ( *tCut )( tGamma ), 1 );
}
