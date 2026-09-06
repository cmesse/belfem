/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Periodic-quotient cohomology tests.
 *
 * The reference topology is the infinite tape (the tapestack3d / sidecoating
 * setup): a periodic domain whose end planes are meshed as copies, where the
 * cohomology carries one generator per encircled conductor PLUS one free
 * axial generator from the periodicity itself. Two miniatures of that:
 *
 *   PeriodicBar     — 1x1xN HEX8 bar, z-fold. Quotient is a solid torus:
 *                     betti = ( 1, 1, 0 ). The single generator is the free
 *                     axial one, paired +-1 with the axial cycle through
 *                     the fold.
 *   PeriodicBarHole — 3x3xN HEX8 bar with the center column removed,
 *                     z-fold. Quotient is ( annulus x S1 ), the skeleton of
 *                     "periodic cylinder minus one wire": betti_1 = 2
 *                     ( encircling + axial ), betti_2 = 1. Generators are
 *                     basis-dependent, so the invariant assertion is the
 *                     2x2 pairing matrix against the two reference cycles
 *                     being unimodular.
 *
 * The periodicity is hand-wired — distinct slave-plane nodes at z = L,
 * symmetric set_periodic() links, and the Periodicity master/slave
 * node/edge/face pair lists populated by id-matching — bypassing the
 * geometric PeriodicityFactory on purpose: the subject under test is the
 * periodic fold inside SimplicialComplex, plus the clean() pass the
 * Cohomology constructor runs over it ( a no-op repair on these meshes —
 * forced repairs live in the solver-tier SPFA tests ), not the plane
 * matcher. As in the annulus suite, the mesh must be built with
 * aComputeConnectivities = true.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <map>
#include <vector>

#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Block.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Mesh_Periodicity.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Homology.hpp"
#include "cl_Cohomology.hpp"
#include "cl_Chain.hpp"
#include "cl_Cochain.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    // grid of ( aM+1 ) x ( aM+1 ) nodes per layer, aN element rings along z
    // plus a distinct slave layer at z = aN. aHole removes the center
    // element column ( requires aM == 3 ). node id: layer * tNPL + j * ( aM
    // + 1 ) + i + 1. The slave-layer nodes pair to layer 0.
    struct PeriodicBar
    {
        Mesh *        mMesh ;
        Periodicity * mPeriodicity ;   // owned by the mesh
        uint          mM ;
        uint          mN ;
        uint          mNPL ;           // nodes per layer

        id_t
        node_id( const uint i, const uint j, const uint aLayer ) const
        {
            return aLayer * mNPL + j * ( mM + 1 ) + i + 1 ;
        }

        PeriodicBar( const uint aM, const uint aN, const bool aHole )
            : mM( aM ), mN( aN ), mNPL( ( aM + 1 ) * ( aM + 1 ) )
        {
            mMesh = new Mesh( 3, 0, true );

            Cell< Node * > & tNodes = mMesh->nodes();

            // layers 0 .. aN, layer aN = slave copies of layer 0 at z = aN
            for ( uint l = 0; l <= aN; ++l )
            {
                for ( uint j = 0; j <= aM; ++j )
                {
                    for ( uint i = 0; i <= aM; ++i )
                    {
                        tNodes.push( new Node( this->node_id( i, j, l ),
                            ( real ) i, ( real ) j, ( real ) l ) );
                    }
                }
            }

            // elements: all columns, minus the center for aHole
            uint tCount = 0 ;
            for ( uint j = 0; j < aM; ++j )
            {
                for ( uint i = 0; i < aM; ++i )
                {
                    if ( aHole && i == aM / 2 && j == aM / 2 ) continue ;
                    ++tCount ;
                }
            }

            Block * tBlock = new Block( 1, tCount * aN );

            ElementFactory tFactory ;
            Cell< Element * > & tElements = mMesh->elements();

            id_t tID = 1 ;
            for ( uint l = 0; l < aN; ++l )
            {
                for ( uint j = 0; j < aM; ++j )
                {
                    for ( uint i = 0; i < aM; ++i )
                    {
                        if ( aHole && i == aM / 2 && j == aM / 2 ) continue ;

                        Element * tElement =
                            tFactory.create_element( ElementType::HEX8, tID++ );

                        // node id - 1 is the push position; the node map the
                        // by-id accessor needs does not exist before finalize
                        // gmsh convention: bottom quad ccw, then top
                        tElement->insert_node( tNodes( this->node_id( i,     j,     l     ) - 1 ), 0 );
                        tElement->insert_node( tNodes( this->node_id( i + 1, j,     l     ) - 1 ), 1 );
                        tElement->insert_node( tNodes( this->node_id( i + 1, j + 1, l     ) - 1 ), 2 );
                        tElement->insert_node( tNodes( this->node_id( i,     j + 1, l     ) - 1 ), 3 );
                        tElement->insert_node( tNodes( this->node_id( i,     j,     l + 1 ) - 1 ), 4 );
                        tElement->insert_node( tNodes( this->node_id( i + 1, j,     l + 1 ) - 1 ), 5 );
                        tElement->insert_node( tNodes( this->node_id( i + 1, j + 1, l + 1 ) - 1 ), 6 );
                        tElement->insert_node( tNodes( this->node_id( i,     j + 1, l + 1 ) - 1 ), 7 );

                        tElement->set_block_id( 1 );
                        tElements.push( tElement );
                        tBlock->insert_element( tElement );
                    }
                }
            }

            mMesh->blocks().push( tBlock );
            mMesh->finalize();
            mMesh->create_edges( false );
            mMesh->create_faces( false );

            this->wire_periodicity();
        }

        ~PeriodicBar()
        {
            delete mMesh ;   // owns nodes, elements, blocks, periodicity
        }

        // symmetric node links plus id-matched edge and face pairs, then
        // register with the mesh so the complex sees has_periodicity()
        void
        wire_periodicity()
        {
            mPeriodicity = new Periodicity( mMesh );

            // node pairs: slave layer mN <-> master layer 0
            for ( uint j = 0; j <= mM; ++j )
            {
                for ( uint i = 0; i <= mM; ++i )
                {
                    Node * tMaster = mMesh->node( this->node_id( i, j, 0 ) );
                    Node * tSlave  = mMesh->node( this->node_id( i, j, mN ) );

                    tMaster->set_periodic( tSlave );
                    tSlave->set_periodic( tMaster );

                    mPeriodicity->master_nodes().push( tMaster );
                    mPeriodicity->slave_nodes().push( tSlave );
                }
            }

            const id_t tSlaveFirst = this->node_id( 0, 0, mN );

            // an entity is in the slave plane if ALL its nodes are slave
            // nodes; its master is found by mapping node ids down one layer
            auto tMapDown = [ & ]( const id_t aID ) -> id_t
            {
                return aID - mN * mNPL ;
            } ;

            // edge pairs by sorted mapped node ids
            std::map< std::pair< id_t, id_t >, Edge * > tMasterEdges ;
            for ( Edge * tEdge : mMesh->edges() )
            {
                const id_t a = tEdge->node( 0 )->id();
                const id_t b = tEdge->node( 1 )->id();
                if ( a < tSlaveFirst && b < tSlaveFirst
                    && a <= mNPL && b <= mNPL )
                {
                    tMasterEdges[ { std::min( a, b ), std::max( a, b ) } ] = tEdge ;
                }
            }
            for ( Edge * tEdge : mMesh->edges() )
            {
                const id_t a = tEdge->node( 0 )->id();
                const id_t b = tEdge->node( 1 )->id();
                if ( a < tSlaveFirst || b < tSlaveFirst ) continue ;

                const id_t am = tMapDown( a );
                const id_t bm = tMapDown( b );
                auto tIt = tMasterEdges.find(
                    { std::min( am, bm ), std::max( am, bm ) } );
                ASSERT_TRUE( tIt != tMasterEdges.end() )
                    << "no master partner for slave edge " << a << "-" << b ;

                tEdge->set_periodic( tIt->second );
                tIt->second->set_periodic( tEdge );
                mPeriodicity->master_edges().push( tIt->second );
                mPeriodicity->slave_edges().push( tEdge );
            }

            // face pairs by sorted mapped node ids
            std::map< std::vector< id_t >, Face * > tMasterFaces ;
            for ( Face * tFace : mMesh->faces() )
            {
                std::vector< id_t > tIDs ;
                bool tOnMasterPlane = true ;
                for ( uint k = 0; k < tFace->number_of_nodes(); ++k )
                {
                    const id_t tNodeID = tFace->node( k )->id();
                    if ( tNodeID > mNPL ) { tOnMasterPlane = false ; break ; }
                    tIDs.push_back( tNodeID );
                }
                if ( ! tOnMasterPlane ) continue ;
                std::sort( tIDs.begin(), tIDs.end() );
                tMasterFaces[ tIDs ] = tFace ;
            }
            for ( Face * tFace : mMesh->faces() )
            {
                std::vector< id_t > tIDs ;
                bool tOnSlavePlane = true ;
                for ( uint k = 0; k < tFace->number_of_nodes(); ++k )
                {
                    const id_t tNodeID = tFace->node( k )->id();
                    if ( tNodeID < tSlaveFirst ) { tOnSlavePlane = false ; break ; }
                    tIDs.push_back( tMapDown( tNodeID ) );
                }
                if ( ! tOnSlavePlane ) continue ;
                std::sort( tIDs.begin(), tIDs.end() );
                auto tIt = tMasterFaces.find( tIDs );
                ASSERT_TRUE( tIt != tMasterFaces.end() )
                    << "no master partner for a slave face" ;

                tFace->set_periodic( tIt->second );
                tIt->second->set_periodic( tFace );
                mPeriodicity->master_faces().push( tIt->second );
                mPeriodicity->slave_faces().push( tFace );
            }

            // production convention ( PeriodicityFactory::fix_face_slaves,
            // private, replicated here ): a slave-plane boundary face is
            // re-based from master to SLAVE role — create_complex reads
            // face->slave() for the fold orientation and would null-deref
            // on the FaceFactory default ( master = element, slave = null )
            Cell< Node * > tCorners ;
            for ( Face * tFace : mPeriodicity->slave_faces() )
            {
                Element * tElement = tFace->master() ;
                tFace->flag_nodes() ;

                uint tOrientation = BELFEM_UINT_MAX ;
                for ( uint f = 0 ; f < tElement->number_of_faces() ; ++f )
                {
                    uint tHits = 0 ;
                    tElement->get_corner_nodes_of_facet( f, tCorners ) ;
                    for ( Node * tNode : tCorners )
                    {
                        if ( ! tNode->is_flagged() ) break ;
                        ++tHits ;
                    }
                    if ( tHits != tFace->number_of_corner_nodes() ) continue ;

                    const id_t tFirst = tCorners.first()->id();
                    for ( uint k = 0 ; k < tHits ; ++k )
                    {
                        if ( tFirst == tFace->node( k )->id() )
                        {
                            tOrientation = k + 1 ;
                            tFace->set_master( nullptr, BELFEM_UINT_MAX ) ;
                            tFace->set_slave( tElement, f, tOrientation ) ;
                            break ;
                        }
                    }
                    if ( tOrientation != BELFEM_UINT_MAX ) break ;
                }
                ASSERT_TRUE( tOrientation != BELFEM_UINT_MAX )
                    << "could not orient slave face" ;
                tFace->unflag_nodes() ;
            }

            mMesh->set_periodicity( mPeriodicity );
        }
    };

    // the CutFactory::compute_cohomologies 3D recipe with
    // mSuggestHomologies on, minus the maxwell dressing
    void
    flag_region_3d( Mesh * aMesh )
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
            tElement->flag_faces();
        }
    }

    // pairing of a 1-cochain with a directed node path given as node ids;
    // consecutive ids must be joined by a mesh edge. NO implicit wrap: a
    // cycle through the fold closes only in the quotient ( its last node is
    // the slave copy of its first ), so callers list every real edge step
    // explicitly — a plain closed loop repeats its first node at the end
    int
    cycle_pairing( Mesh * aMesh, Cochain * aGenerator,
        const std::vector< id_t > & aCycle )
    {
        std::map< index_t, int > tCoeffs ;
        for ( const auto & [ tIndex, tCoeff ] : aGenerator->getSimplicesMap() )
        {
            tCoeffs[ tIndex ] = tCoeff ;
        }

        // edge lookup by unordered endpoint ids
        std::map< std::pair< id_t, id_t >, Edge * > tEdgeMap ;
        for ( Edge * tEdge : aMesh->edges() )
        {
            const id_t a = tEdge->node( 0 )->id();
            const id_t b = tEdge->node( 1 )->id();
            tEdgeMap[ { std::min( a, b ), std::max( a, b ) } ] = tEdge ;
        }

        int tPairing = 0 ;
        const uint n = aCycle.size();
        for ( uint k = 0; k + 1 < n; ++k )
        {
            const id_t tA = aCycle[ k ];
            const id_t tB = aCycle[ k + 1 ];

            auto tIt = tEdgeMap.find( { std::min( tA, tB ), std::max( tA, tB ) } );
            EXPECT_TRUE( tIt != tEdgeMap.end() )
                << "cycle step " << tA << "->" << tB << " is not a mesh edge" ;
            if ( tIt == tEdgeMap.end() ) continue ;

            Edge * tEdge = tIt->second ;

            auto tCoeff = tCoeffs.find( tEdge->index() );
            if ( tCoeff == tCoeffs.end() ) continue ;

            const int tSign = tEdge->node( 0 )->id() == tA ? 1 : -1 ;
            tPairing += tSign * tCoeff->second ;
        }
        return tPairing ;
    }

    void
    check_unit_coefficients( Cochain * aGenerator )
    {
        for ( const auto & [ tIndex, tCoeff ] : aGenerator->getSimplicesMap() )
        {
            EXPECT_TRUE( tCoeff == 1 || tCoeff == -1 )
                << "edge index " << tIndex << " carries coefficient " << tCoeff ;
        }
    }
}

// solid torus: the free axial generator alone
TEST( CohomologyPeriodic, PeriodicBar )
{
    const uint tN = 4 ;
    PeriodicBar tBar( 1, tN, false );
    Mesh * tMesh = tBar.mMesh ;

    flag_region_3d( tMesh );

    SimplicialComplex * tComplex = new SimplicialComplex( tMesh, true );
    tComplex->coreduce_complexPellikkaGeneralized();

    Homology * tHomology = new Homology( tComplex, tMesh );
    ASSERT_EQ( tHomology->get_Generators()( 0 ).size(), 1u );
    ASSERT_EQ( tHomology->get_Generators()( 1 ).size(), 1u );
    EXPECT_EQ( tHomology->get_Generators()( 2 ).size(), 0u );

    // constructor = SNF generators + clean(), i.e. clean_spfa over the
    // periodic fold — the subject under test. updatekGeneratorsFromHomology is
    // deliberately NOT called here: its 3D contract expects the suggested
    // homology as terminal in/out PAIRS ( [in0, out0, ...], one condition
    // per pair ), which raw SNF generators are not; that path is covered by
    // the 2D annulus suite, whose 2D contract is one-per-condition
    Cohomology * tCohomology = new Cohomology( tComplex, tMesh );

    Cell< Cochain * > & tGenerators = tCohomology->get_Generators()( 1 );
    ASSERT_EQ( tGenerators.size(), 1u );
    check_unit_coefficients( tGenerators( 0 ) );

    // axial cycle through the fold: corner column, slave id last
    std::vector< id_t > tAxial ;
    for ( uint l = 0; l < tN; ++l )
    {
        tAxial.push_back( tBar.node_id( 0, 0, l ) );
    }
    tAxial.push_back( tBar.node_id( 0, 0, tN ) );   // folds onto layer 0

    EXPECT_EQ( std::abs( cycle_pairing(
        tMesh, tGenerators( 0 ), tAxial ) ), 1 );

    delete tCohomology ;
    delete tHomology ;
    delete tComplex ;
}

// ( annulus x S1 ): encircling generator + free axial generator — the
// infinite-tape skeleton ( "periodic cylinder minus one wire" )
TEST( CohomologyPeriodic, PeriodicBarHole )
{
    const uint tN = 4 ;
    PeriodicBar tBar( 3, tN, true );
    Mesh * tMesh = tBar.mMesh ;

    flag_region_3d( tMesh );

    SimplicialComplex * tComplex = new SimplicialComplex( tMesh, true );
    tComplex->coreduce_complexPellikkaGeneralized();

    Homology * tHomology = new Homology( tComplex, tMesh );
    ASSERT_EQ( tHomology->get_Generators()( 0 ).size(), 1u );
    ASSERT_EQ( tHomology->get_Generators()( 1 ).size(), 2u );
    EXPECT_EQ( tHomology->get_Generators()( 2 ).size(), 1u );

    // ctor generators only — see the note in PeriodicBar: the 3D update
    // contract wants terminal in/out pairs, not raw SNF generators
    Cohomology * tCohomology = new Cohomology( tComplex, tMesh );

    Cell< Cochain * > & tGenerators = tCohomology->get_Generators()( 1 );
    ASSERT_EQ( tGenerators.size(), 2u );
    check_unit_coefficients( tGenerators( 0 ) );
    check_unit_coefficients( tGenerators( 1 ) );

    // reference cycle 1: axial through the fold, at the outer corner
    std::vector< id_t > tAxial ;
    for ( uint l = 0; l < tN; ++l )
    {
        tAxial.push_back( tBar.node_id( 0, 0, l ) );
    }
    tAxial.push_back( tBar.node_id( 0, 0, tN ) );

    // reference cycle 2: around the removed center column, in the master
    // plane: the four nodes bounding the hole, first repeated to close
    std::vector< id_t > tLoop = {
        tBar.node_id( 1, 1, 0 ), tBar.node_id( 2, 1, 0 ),
        tBar.node_id( 2, 2, 0 ), tBar.node_id( 1, 2, 0 ),
        tBar.node_id( 1, 1, 0 ) } ;

    // generators are basis-dependent; the pairing matrix must be unimodular
    const int a = cycle_pairing( tMesh, tGenerators( 0 ), tAxial );
    const int b = cycle_pairing( tMesh, tGenerators( 0 ), tLoop );
    const int c = cycle_pairing( tMesh, tGenerators( 1 ), tAxial );
    const int d = cycle_pairing( tMesh, tGenerators( 1 ), tLoop );

    EXPECT_EQ( std::abs( a * d - b * c ), 1 )
        << "pairing matrix [ " << a << " " << b << " ; "
        << c << " " << d << " ] is not unimodular" ;

    delete tCohomology ;
    delete tHomology ;
    delete tComplex ;
}
