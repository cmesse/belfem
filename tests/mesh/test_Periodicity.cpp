/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit tests for mesh::PeriodicityFactory::create_periodicity() on hand-built
 * HEX8 bars: the node bijection and its positional pairing contract, the
 * absolute tolerance at the facet and node stages, the candidate named by the
 * k-d tree on a rejection, and conductor edges tied through the master-element
 * fallback when the cap wrappers carry no edge container.
 *
 * The mesh main forces the throwing reaction, so every BELFEM_ERROR on the
 * refusal paths is observable as a std::runtime_error in both builds. A mesh
 * is never reused after a refusal: the factory leaves node links, flags and
 * indices dirty when it throws.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Set.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "cl_Mesh_Periodicity.hpp"
#include "cl_Mesh_PeriodicityFactory.hpp"
#include "meshtools.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    const real tTiny = 1e-12 ;

    // runs aCall, requires a std::runtime_error whose what() carries the
    // literal part of the BELFEM_ERROR text
    template< typename F >
    void
    expect_refusal( F && aCall, const std::string & aSubstring )
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

    // wraps a fresh QUAD4 element into a facet whose nodes are taken from
    // facet aIndex of aMaster; the facet owns the wrapper, the sideset owns
    // the facet, so the wrapper is never pushed into a block
    Facet *
    create_cap_facet( ElementFactory & aFactory,
                      const id_t aID,
                      Element * aMaster,
                      const uint aIndex,
                      const id_t aSideSetID )
    {
        Element * tSurf = aFactory.create_element( ElementType::QUAD4, aID );
        Facet * tFacet = new Facet( tSurf );
        tFacet->set_master( aMaster, aIndex );
        tFacet->set_sideset_id( aSideSetID );
        return tFacet ;
    }

    // M x M x N HEX8 bar in gmsh corner order, one block. Grid point
    // ( i, j, l ) has id l * ( M + 1 )^2 + j * ( M + 1 ) + i + 1 and sits at
    // ( i - c, j - c, l ), c = M / 2 when centred. Sideset 1 is the z = 0 cap
    // ( master side ), sideset 2 the z = N cap ( slave side ); facet ids are
    // disjoint from the brick ids. Coordinates may be edited before
    // finalize(). The plane triples are members so the factory's non-const
    // Cell< Node * > & overload binds to construction pointers.
    struct PeriodicBrick
    {
        Mesh * mMesh ;
        uint   mM ;
        uint   mN ;
        uint   mNPL ;
        real   mShift ;

        Cell< Node * > mMasterTriple ;
        Cell< Node * > mSlaveTriple ;

        id_t
        node_id( const uint i, const uint j, const uint l ) const
        {
            return l * mNPL + j * ( mM + 1 ) + i + 1 ;
        }

        // node id - 1 is the push position; the by-id accessor needs the map
        // that finalize() builds
        Node *
        node( const uint i, const uint j, const uint l )
        {
            return mMesh->nodes()( this->node_id( i, j, l ) - 1 );
        }

        void
        perturb( const uint i, const uint j, const uint l,
                 const real aDx, const real aDy, const real aDz )
        {
            Node * tNode = this->node( i, j, l );
            tNode->set_coords( tNode->x() + aDx, tNode->y() + aDy, tNode->z() + aDz );
        }

        PeriodicBrick( const uint aM, const uint aN, const bool aCentred ) :
            mM( aM ), mN( aN ), mNPL( ( aM + 1 ) * ( aM + 1 ) ),
            mShift( aCentred ? 0.5 * ( real ) aM : 0.0 )
        {
            mMesh = new Mesh( 3, 0, true );

            Cell< Node * > & tNodes = mMesh->nodes();

            for ( uint l = 0; l <= aN; ++l )
            {
                for ( uint j = 0; j <= aM; ++j )
                {
                    for ( uint i = 0; i <= aM; ++i )
                    {
                        tNodes.push( new Node( this->node_id( i, j, l ),
                            ( real ) i - mShift, ( real ) j - mShift, ( real ) l ) );
                    }
                }
            }

            Block * tBlock = new Block( 1, aM * aM * aN );

            ElementFactory tFactory ;
            Cell< Element * > & tElements = mMesh->elements();

            id_t tID = 1 ;
            for ( uint l = 0; l < aN; ++l )
            {
                for ( uint j = 0; j < aM; ++j )
                {
                    for ( uint i = 0; i < aM; ++i )
                    {
                        Element * tElement =
                            tFactory.create_element( ElementType::HEX8, tID++ );

                        // gmsh convention: bottom quad ccw, then top
                        tElement->insert_node( this->node( i,     j,     l     ), 0 );
                        tElement->insert_node( this->node( i + 1, j,     l     ), 1 );
                        tElement->insert_node( this->node( i + 1, j + 1, l     ), 2 );
                        tElement->insert_node( this->node( i,     j + 1, l     ), 3 );
                        tElement->insert_node( this->node( i,     j,     l + 1 ), 4 );
                        tElement->insert_node( this->node( i + 1, j,     l + 1 ), 5 );
                        tElement->insert_node( this->node( i + 1, j + 1, l + 1 ), 6 );
                        tElement->insert_node( this->node( i,     j + 1, l + 1 ), 7 );

                        tElement->set_block_id( 1 );
                        tElements.push( tElement );
                        tBlock->insert_element( tElement );
                    }
                }
            }
            mMesh->blocks().push( tBlock );

            // brick ( i, j ) of layer l is element l * M^2 + j * M + i
            SideSet * tBottom = new SideSet( 1, aM * aM );
            SideSet * tTop    = new SideSet( 2, aM * aM );
            for ( uint j = 0; j < aM; ++j )
            {
                for ( uint i = 0; i < aM; ++i )
                {
                    const uint tColumn = j * aM + i ;

                    tBottom->insert_facet( create_cap_facet( tFactory,
                        1000 + tColumn + 1,
                        tElements( tColumn ),
                        bottom_facet_index( ElementType::HEX8 ), 1 ) );

                    tTop->insert_facet( create_cap_facet( tFactory,
                        2000 + tColumn + 1,
                        tElements( ( aN - 1 ) * aM * aM + tColumn ),
                        top_facet_index( ElementType::HEX8 ), 2 ) );
                }
            }
            mMesh->sidesets().push( tBottom );
            mMesh->sidesets().push( tTop );
        }

        ~PeriodicBrick()
        {
            delete mMesh ;   // owns nodes, elements, blocks, sidesets, edges, periodicity
        }

        // owning fixture: a copy or a move would double-free
        PeriodicBrick( const PeriodicBrick & ) = delete ;
        PeriodicBrick & operator=( const PeriodicBrick & ) = delete ;
        PeriodicBrick( PeriodicBrick && ) = delete ;
        PeriodicBrick & operator=( PeriodicBrick && ) = delete ;

        void
        finalize()
        {
            mMesh->finalize();
        }

        // edges on the bricks only: the cap wrappers stay without an edge
        // container, so the factory must reach the edges through the master
        void
        create_conductor_edges()
        {
            Vector< id_t > tBlocks = { 1 };
            Vector< id_t > tSideSets ;
            mMesh->create_edges( false, tBlocks, tSideSets, false );
        }

        // P = e_x, Q = e_y, N = e_z on both caps: pure translation along z
        void
        set_translation_planes( PeriodicityFactory & aFactory )
        {
            mMasterTriple = { this->node( 0, 0, 0 ), this->node( 1, 0, 0 ), this->node( 0, 1, 0 ) };
            mSlaveTriple  = { this->node( 0, 0, mN ), this->node( 1, 0, mN ), this->node( 0, 1, mN ) };
            aFactory.set_master_plane( mMasterTriple );
            aFactory.set_slave_plane( mSlaveTriple );
        }

        // slave frame rotated 90 degrees about z through the origin
        // ( P = e_y, Q = -e_x ): the partner of ( x, y, 0 ) is ( -y, x, N ).
        // Needs the centred cap, because the in-plane coordinates are taken
        // from absolute positions
        void
        set_twist_planes( PeriodicityFactory & aFactory )
        {
            const uint tCentre = mM / 2 ;
            mMasterTriple = { this->node( tCentre,     tCentre,     0 ),
                              this->node( tCentre + 1, tCentre,     0 ),
                              this->node( tCentre,     tCentre + 1, 0 ) };
            mSlaveTriple  = { this->node( tCentre,     tCentre,     mN ),
                              this->node( tCentre,     tCentre + 1, mN ),
                              this->node( tCentre - 1, tCentre,     mN ) };
            aFactory.set_master_plane( mMasterTriple );
            aFactory.set_slave_plane( mSlaveTriple );
        }

        // the mesh takes ownership of the periodicity
        Periodicity *
        create( PeriodicityFactory & aFactory )
        {
            Periodicity * tPeriodicity = aFactory.create_periodicity();
            mMesh->set_periodicity( tPeriodicity );
            return tPeriodicity ;
        }

        // partner id under the twist: ( i, j ) -> ( M - j, i ) on layer N
        id_t
        twisted_partner_id( const id_t aMasterID ) const
        {
            const uint i = ( aMasterID - 1 ) % ( mM + 1 );
            const uint j = ( aMasterID - 1 ) / ( mM + 1 );
            return this->node_id( mM - j, i, mN );
        }
    };

    // the pair lists are a bijection between the two caps, ids are untouched,
    // node indices are restored, and slave( k ) is the partner of master( k )
    void
    expect_node_bijection( Mesh * aMesh, Periodicity * aPeriodicity, const index_t aNumCapNodes )
    {
        Cell< Node * > & tMasters = aPeriodicity->master_nodes();
        Cell< Node * > & tSlaves  = aPeriodicity->slave_nodes();

        EXPECT_TRUE( aPeriodicity->is_flagged( EntityType::NODE ) );
        ASSERT_EQ( tMasters.size(), aNumCapNodes );
        ASSERT_EQ( tSlaves.size(), aNumCapNodes );

        Set< Node * > tPartners ;
        for ( index_t k = 0; k < aNumCapNodes; ++k )
        {
            Node * P = tMasters( k );
            Node * Q = tSlaves( k );
            ASSERT_NE( P, nullptr );
            ASSERT_NE( Q, nullptr );
            EXPECT_EQ( P->periodic(), Q );
            EXPECT_EQ( Q->periodic(), P );
            EXPECT_NE( P, Q );
            tPartners.insert( Q );
        }
        EXPECT_EQ( tPartners.size(), ( std::size_t ) aNumCapNodes );

        Cell< Node * > & tNodes = aMesh->nodes();
        for ( index_t k = 0; k < tNodes.size(); ++k )
        {
            EXPECT_EQ( tNodes( k )->index(), k );
            EXPECT_EQ( tNodes( k )->id(), ( id_t ) k + 1 );
        }
    }

    // every master edge is tied to a slave edge whose endpoints are the
    // partners of its own, in the same order; aNumReversed counts the slave
    // edges the factory had to flip out of their creation order
    void
    expect_edges_tied( Periodicity * aPeriodicity,
                       const index_t aNumEdges,
                       const index_t aNumReversed )
    {
        Cell< Edge * > & tMasters = aPeriodicity->master_edges();
        Cell< Edge * > & tSlaves  = aPeriodicity->slave_edges();

        EXPECT_TRUE( aPeriodicity->is_flagged( EntityType::EDGE ) );
        ASSERT_EQ( tMasters.size(), aNumEdges );
        ASSERT_EQ( tSlaves.size(), aNumEdges );

        index_t tReversed = 0 ;
        for ( index_t k = 0; k < aNumEdges; ++k )
        {
            Edge * E = tMasters( k );
            Edge * F = tSlaves( k );
            ASSERT_NE( E, nullptr );
            ASSERT_NE( F, nullptr );
            EXPECT_EQ( E->periodic(), F );
            EXPECT_EQ( F->periodic(), E );

            for ( uint i = 0; i < 2; ++i )
            {
                EXPECT_EQ( F->node( i )->original(), E->node( i )->periodic() );
            }
            if ( F->node( 0 )->id() > F->node( 1 )->id() )
            {
                ++tReversed ;
            }
        }
        EXPECT_EQ( tReversed, aNumReversed );
    }
}

// =============================================================================
// node bijection
// =============================================================================

TEST( Periodicity, TranslationBijectionPreservesIds )
{
    PeriodicBrick tBar( 2, 2, false );
    tBar.finalize();

    PeriodicityFactory tFactory( tBar.mMesh );
    tBar.set_translation_planes( tFactory );
    Periodicity * tPeriodicity = tBar.create( tFactory );

    expect_node_bijection( tBar.mMesh, tPeriodicity, 9 );
    EXPECT_EQ( tPeriodicity->master_facets().size(), ( index_t ) 4 );

    // the partner of a cap node is its copy N layers up
    const id_t tOffset = tBar.mN * tBar.mNPL ;
    for ( Node * P : tPeriodicity->master_nodes() )
    {
        EXPECT_NEAR( P->z(), 0.0, tTiny );
        EXPECT_EQ( P->periodic()->id(), P->id() + tOffset );
    }
}

TEST( Periodicity, TwistedCapMatchesRotatedPartner )
{
    PeriodicBrick tBar( 2, 2, true );
    tBar.finalize();

    PeriodicityFactory tFactory( tBar.mMesh );
    tBar.set_twist_planes( tFactory );
    Periodicity * tPeriodicity = tBar.create( tFactory );

    expect_node_bijection( tBar.mMesh, tPeriodicity, 9 );

    for ( Node * P : tPeriodicity->master_nodes() )
    {
        Node * Q = P->periodic();
        EXPECT_EQ( Q->id(), tBar.twisted_partner_id( P->id() ) );
        EXPECT_NEAR( Q->x(), -P->y(), tTiny );
        EXPECT_NEAR( Q->y(),  P->x(), tTiny );
        EXPECT_NEAR( Q->z(), ( real ) tBar.mN, tTiny );
    }
}

// =============================================================================
// tolerance boundary
// =============================================================================

TEST( Periodicity, FacetPerturbationIsRefused )
{
    PeriodicBrick tBar( 2, 2, false );

    // the shared centre node moves every slave centroid by 2.5e-7 m
    tBar.perturb( 1, 1, 2, 1e-6, 0.0, 0.0 );
    tBar.finalize();

    PeriodicityFactory tFactory( tBar.mMesh );
    tBar.set_translation_planes( tFactory );

    expect_refusal( [&]() { tFactory.create_periodicity(); },
        "Could not find a corresponding periodic facet" );
}

TEST( Periodicity, SubToleranceMismatchIsAccepted )
{
    PeriodicBrick tBar( 2, 2, false );

    tBar.perturb( 1, 1, 2, 1e-10, 0.0, 0.0 );
    Node * tMoved = tBar.node( 1, 1, 2 );
    tBar.finalize();

    PeriodicityFactory tFactory( tBar.mMesh );
    tBar.set_translation_planes( tFactory );
    Periodicity * tPeriodicity = tBar.create( tFactory );

    expect_node_bijection( tBar.mMesh, tPeriodicity, 9 );
    EXPECT_EQ( tBar.node( 1, 1, 0 )->periodic(), tMoved );
    EXPECT_NEAR( tMoved->x(), 1.0 + 1e-10, tTiny );
}

TEST( Periodicity, NodeMismatchIsRefusedInEveryBuild )
{
    PeriodicBrick tBar( 1, 2, false );

    // opposite in-plane shifts keep the cap centroid, so the facet stage
    // accepts the pair and the node stage is left with two orphans
    tBar.perturb( 0, 0, 2,  1e-6, 0.0, 0.0 );
    tBar.perturb( 1, 1, 2, -1e-6, 0.0, 0.0 );
    tBar.finalize();

    PeriodicityFactory tFactory( tBar.mMesh );
    tBar.set_translation_planes( tFactory );

    expect_refusal( [&]() { tFactory.create_periodicity(); },
        "matched 2 of 4 target nodes of the periodic pair within" );
}

// =============================================================================
// k-d tree rejection
// =============================================================================

namespace
{
    // unit HEX8 bricks that share no nodes: three masters with a bottom cap
    // ( sideset 1, z = 0 ) and three slaves with a top cap ( sideset 2,
    // z = 1 ). The first master cap sits in-plane at ( 3, 0 ) and no slave
    // cap is within the tolerance of any master cap
    struct LooseBricks
    {
        Mesh * mMesh ;
        id_t   mNodeID = 1 ;
        id_t   mElementID = 1 ;

        Cell< Node * > mMasterTriple ;
        Cell< Node * > mSlaveTriple ;

        Cell< Facet * > mSlaveFacets ;
        Cell< real >    mSlaveX ;
        Cell< real >    mSlaveY ;

        ElementFactory     mFactory ;
        Cell< Element * >  mBricks ;

        Element *
        add_brick( const real aX, const real aY, Block * aBlock )
        {
            Cell< Node * > & tNodes = mMesh->nodes();
            Element * tElement = mFactory.create_element( ElementType::HEX8, mElementID++ );

            const real tX[ 4 ] = { aX - 0.5, aX + 0.5, aX + 0.5, aX - 0.5 };
            const real tY[ 4 ] = { aY - 0.5, aY - 0.5, aY + 0.5, aY + 0.5 };

            for ( uint l = 0; l < 2; ++l )
            {
                for ( uint k = 0; k < 4; ++k )
                {
                    Node * tNode = new Node( mNodeID++, tX[ k ], tY[ k ], ( real ) l );
                    tNodes.push( tNode );
                    tElement->insert_node( tNode, 4 * l + k );
                }
            }
            tElement->set_block_id( 1 );
            mMesh->elements().push( tElement );
            aBlock->insert_element( tElement );
            mBricks.push( tElement );
            return tElement ;
        }

        LooseBricks()
        {
            mMesh = new Mesh( 3, 0, true );
            Block * tBlock = new Block( 1, 6 );

            SideSet * tBottom = new SideSet( 1, 3 );
            const real tMasterX[ 3 ] = {  3.0,  20.0, -20.0 };
            const real tMasterY[ 3 ] = {  0.0,  20.0, -20.0 };
            for ( uint k = 0; k < 3; ++k )
            {
                Element * tBrick = this->add_brick( tMasterX[ k ], tMasterY[ k ], tBlock );
                tBottom->insert_facet( create_cap_facet( mFactory, 1001 + k, tBrick,
                    bottom_facet_index( ElementType::HEX8 ), 1 ) );
            }
            mMasterTriple = { mBricks( 0 )->node( 0 ), mBricks( 0 )->node( 1 ), mBricks( 0 )->node( 3 ) };

            SideSet * tTop = new SideSet( 2, 3 );
            const real tSlaveX[ 3 ] = { -1.0,   0.0, 3.0 };
            const real tSlaveY[ 3 ] = {  0.0, 100.0, 5.0 };
            mSlaveFacets.set_size( 3, nullptr );
            mSlaveX.set_size( 3, 0.0 );
            mSlaveY.set_size( 3, 0.0 );
            for ( uint k = 0; k < 3; ++k )
            {
                Element * tBrick = this->add_brick( tSlaveX[ k ], tSlaveY[ k ], tBlock );
                mSlaveFacets( k ) = create_cap_facet( mFactory, 2001 + k, tBrick,
                    top_facet_index( ElementType::HEX8 ), 2 );
                mSlaveX( k ) = tSlaveX[ k ];
                mSlaveY( k ) = tSlaveY[ k ];
                tTop->insert_facet( mSlaveFacets( k ) );
            }
            mSlaveTriple = { mBricks( 3 )->node( 4 ), mBricks( 3 )->node( 5 ), mBricks( 3 )->node( 7 ) };

            mMesh->blocks().push( tBlock );
            mMesh->sidesets().push( tBottom );
            mMesh->sidesets().push( tTop );
            mMesh->finalize();
        }

        ~LooseBricks()
        {
            delete mMesh ;
        }

        LooseBricks( const LooseBricks & ) = delete ;
        LooseBricks & operator=( const LooseBricks & ) = delete ;
        LooseBricks( LooseBricks && ) = delete ;
        LooseBricks & operator=( LooseBricks && ) = delete ;

        // brute-force nearest slave cap to an in-plane point
        Facet *
        nearest_slave_facet( const real aX, const real aY ) const
        {
            Facet * tBest = nullptr ;
            real tBestDist = BELFEM_REAL_MAX ;
            for ( index_t k = 0; k < mSlaveFacets.size(); ++k )
            {
                const real tDist = std::sqrt(
                    ( aX - mSlaveX( k ) ) * ( aX - mSlaveX( k ) ) +
                    ( aY - mSlaveY( k ) ) * ( aY - mSlaveY( k ) ) );
                if ( tDist < tBestDist )
                {
                    tBestDist = tDist ;
                    tBest = mSlaveFacets( k );
                }
            }
            return tBest ;
        }
    };
}

TEST( Periodicity, ClosestCandidateIsNamedByLinearDistance )
{
    LooseBricks tBricks ;

    PeriodicityFactory tFactory( tBricks.mMesh );
    tFactory.set_master_plane( tBricks.mMasterTriple );
    tFactory.set_slave_plane( tBricks.mSlaveTriple );

    // the tree roots at ( 0, 100 ) and descends to ( 3, 5 ) first; the far
    // side ( -1, 0 ) is 4 m away and must not be pruned by the 3 m split
    Facet * tNearest = tBricks.nearest_slave_facet( 3.0, 0.0 );
    ASSERT_NE( tNearest, nullptr );
    EXPECT_EQ( tNearest, tBricks.mSlaveFacets( 0 ) );

    const std::string tExpected =
        "for facet 1001 (master: 1), closest candidate: "
        + std::to_string( tNearest->id() )
        + " (master: " + std::to_string( tNearest->master()->id() ) + ")";

    expect_refusal( [&]() { tFactory.create_periodicity(); }, tExpected );
}

// =============================================================================
// conductor edges
// =============================================================================

TEST( Periodicity, ConductorEdgesAreTiedThroughTheMasterFallback )
{
    // translation: creation order ( lower id first ) is preserved by the
    // layer shift, so no slave edge needs flipping
    {
        PeriodicBrick tBar( 2, 2, false );
        tBar.finalize();
        tBar.create_conductor_edges();

        PeriodicityFactory tFactory( tBar.mMesh );
        tBar.set_translation_planes( tFactory );
        Periodicity * tPeriodicity = tBar.create( tFactory );

        expect_node_bijection( tBar.mMesh, tPeriodicity, 9 );
        expect_edges_tied( tPeriodicity, 12, 0 );
    }

    // twist: the partner ids of a y-edge ( i, j ) -> ( i, j + 1 ) descend,
    // so each of the six y-edges is flipped and the six x-edges are not
    {
        PeriodicBrick tBar( 2, 2, true );
        tBar.finalize();
        tBar.create_conductor_edges();

        PeriodicityFactory tFactory( tBar.mMesh );
        tBar.set_twist_planes( tFactory );
        Periodicity * tPeriodicity = tBar.create( tFactory );

        expect_node_bijection( tBar.mMesh, tPeriodicity, 9 );
        expect_edges_tied( tPeriodicity, 12, 6 );
    }
}
