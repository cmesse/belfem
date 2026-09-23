/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * What these cases lock:
 *
 *  - CutData refuses a coefficient that is not a unit and names the edge it
 *    sits on, instead of recording a silent zero weight in its bitsets,
 *  - the CutData bitsets reproduce the cochain edge by edge, sign included,
 *  - clean() repairs a generator perturbed by an exact node coboundary to
 *    unit coefficients in the same class, on no more edges than the
 *    unperturbed generator acted on,
 *  - a disk with k holes yields k generators whose pairing matrix against the
 *    k hole rings is unimodular, whichever coreduction produced them
 *    ( Pellikka et al. 2013, Eq. 5.3 ),
 *  - updatekGeneratorsFromHomology refuses more conditions than generators,
 *    and conditions that couple to the generators with too low a rank, each
 *    with its own message,
 *  - clean() refuses a generator that has no unit representative on the mesh
 *    as a thick cut, and still rectifies one that has.
 *
 * These fixtures are not the CutFactory::mSuggestHomologies sequence. They
 * build Homology( SimplicialComplex *, Mesh * ), the chain-side Smith form,
 * where the factory builds Homology( Mesh *, ... ) and reorients its
 * generators first; reorient_generators() is never called here. The Mesh must
 * be built with aComputeConnectivities = true, or finalize_edges skips
 * connect_edges_to_elements and the coreduction cannot collapse a copair.
 *
 * Scope: only the two refusal cases and the bitset case construct CutData.
 * The cleaned-support case and the k-hole cases drive Cohomology::clean and
 * the homology of the complex, and construct none. CutProcessor, and with it
 * the thick to thin conversion and the 3D face pruning, is not constructed at
 * all.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <unistd.h>

#include "constants.hpp"
#include "cl_Matrix.hpp"
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
#include "cl_CutData.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    enum class Algorithm { Pellikka, CCR, PellikkaGeneralized };

    // one ring of TRI3 elements between radii aR0 < aR1, counterclockwise.
    // node ids: 1..aN inner ring, aN+1..2*aN outer ring, both starting at
    // polar angle zero. The caller owns the mesh.
    Mesh *
    create_annulus( const uint aN,
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

        Block * tBlock = new Block( 1, 2 * aN );

        ElementFactory tFactory ;
        Cell< Element * > & tElements = tMesh->elements();

        id_t tID = 1 ;
        for ( uint k = 0; k < aN; ++k )
        {
            const uint tNext = ( k + 1 ) % aN ;

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

        tMesh->blocks().push( tBlock );

        tMesh->finalize();
        tMesh->create_edges( false );

        return tMesh ;
    }

    // a ( 3k+1 ) x 4 grid of unit quads at integer lattice positions, each
    // split into two counterclockwise TRI3, with k rectangular holes: the
    // 2 x 2 quad block at columns 3h+1, 3h+2 and rows 1, 2 is left out, which
    // leaves a one-cell frame around the strip and a one-cell bridge between
    // neighbouring holes. The result is a disk with k holes, so b1 = k. The
    // node at the centre of a hole would carry no element and is not created.
    // aRings receives, per hole, the counterclockwise walk of the eight
    // boundary nodes of the removed block as node ids, the first repeated at
    // the end to close the loop. The caller owns the mesh.
    Mesh *
    create_punched_strip( const uint aK,
        std::vector< std::vector< id_t > > & aRings )
    {
        const uint tNumCols = 3 * aK + 2 ;
        const uint tNumRows = 5 ;

        auto tIsHoleCentre = [ aK ]( const uint aI, const uint aJ ) -> bool
        {
            if ( aJ != 2 ) return false ;
            for ( uint h = 0; h < aK; ++h )
            {
                if ( aI == 3 * h + 2 ) return true ;
            }
            return false ;
        };

        auto tIsHoleQuad = [ aK ]( const uint aI, const uint aJ ) -> bool
        {
            if ( aJ < 1 || aJ > 2 ) return false ;
            for ( uint h = 0; h < aK; ++h )
            {
                if ( aI == 3 * h + 1 || aI == 3 * h + 2 ) return true ;
            }
            return false ;
        };

        Mesh * tMesh = new Mesh( 2, 0, true );

        Cell< Node * > & tNodes = tMesh->nodes();

        // lattice position to node, nullptr at the hole centres
        std::vector< Node * > tLattice( tNumCols * tNumRows, nullptr );

        id_t tNodeID = 1 ;
        for ( uint j = 0; j < tNumRows; ++j )
        {
            for ( uint i = 0; i < tNumCols; ++i )
            {
                if ( tIsHoleCentre( i, j ) ) continue ;

                Node * tNode =
                    new Node( tNodeID++, ( real ) i, ( real ) j );
                tNodes.push( tNode );
                tLattice[ i + j * tNumCols ] = tNode ;
            }
        }

        Block * tBlock = new Block( 1, 16 * aK + 8 );

        ElementFactory tFactory ;
        Cell< Element * > & tElements = tMesh->elements();

        id_t tID = 1 ;
        for ( uint j = 0; j + 1 < tNumRows; ++j )
        {
            for ( uint i = 0; i + 1 < tNumCols; ++i )
            {
                if ( tIsHoleQuad( i, j ) ) continue ;

                Node * tA = tLattice[ i     + j       * tNumCols ];
                Node * tB = tLattice[ i + 1 + j       * tNumCols ];
                Node * tC = tLattice[ i + 1 + ( j+1 ) * tNumCols ];
                Node * tD = tLattice[ i     + ( j+1 ) * tNumCols ];

                Element * tElement =
                    tFactory.create_element( ElementType::TRI3, tID++ );
                tElement->insert_node( tA, 0 );
                tElement->insert_node( tB, 1 );
                tElement->insert_node( tC, 2 );
                tElement->set_block_id( 1 );
                tElements.push( tElement );
                tBlock->insert_element( tElement );

                Element * tElement2 =
                    tFactory.create_element( ElementType::TRI3, tID++ );
                tElement2->insert_node( tA, 0 );
                tElement2->insert_node( tC, 1 );
                tElement2->insert_node( tD, 2 );
                tElement2->set_block_id( 1 );
                tElements.push( tElement2 );
                tBlock->insert_element( tElement2 );
            }
        }

        tMesh->blocks().push( tBlock );

        tMesh->finalize();
        tMesh->create_edges( false );

        aRings.clear();
        for ( uint h = 0; h < aK; ++h )
        {
            const uint tLeft  = 3 * h + 1 ;
            const uint tRight = 3 * h + 3 ;

            const std::vector< std::pair< uint, uint > > tWalk = {
                { tLeft,      1 }, { tLeft + 1,  1 }, { tRight,     1 },
                { tRight,     2 }, { tRight,     3 }, { tLeft + 1,  3 },
                { tLeft,      3 }, { tLeft,      2 }, { tLeft,      1 } };

            std::vector< id_t > tRing ;
            for ( const auto & [ tI, tJ ] : tWalk )
            {
                tRing.push_back( tLattice[ tI + tJ * tNumCols ]->id() );
            }
            aRings.push_back( tRing );
        }

        return tMesh ;
    }

    // the region flagging CutFactory::compute_cohomologies does in 2D before
    // the complex is built, minus the maxwell dressing
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

    void
    coreduce( SimplicialComplex * aComplex, const Algorithm aAlgorithm )
    {
        switch ( aAlgorithm )
        {
            case Algorithm::Pellikka :
            {
                aComplex->coreduce_complexPellikka();
                break ;
            }
            case Algorithm::CCR :
            {
                aComplex->coreduce_complexCCR();
                break ;
            }
            case Algorithm::PellikkaGeneralized :
            {
                aComplex->coreduce_complexPellikkaGeneralized();
                break ;
            }
        }
    }

    // +1 if the stored inner-ring edge aA -> aB runs counterclockwise
    // ( ring position a -> a+1 mod aN ), -1 otherwise. Read from node ids so
    // the oracle is independent of edge_direction.
    int
    inner_ring_sign( const id_t aA, const id_t aB, const uint aN )
    {
        return ( aB % aN ) == ( aA % aN + 1 ) % aN ? 1 : -1 ;
    }

    // pairing of a 1-cochain with the inner-ring cycle of the annulus,
    // traversed counterclockwise ( node ids 1..aN in ascending order )
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

            tPairing += inner_ring_sign( tA, tB, aN ) * tIt->second ;
        }

        return tPairing ;
    }

    // pairing of a 1-cochain with a directed node path given as node ids;
    // consecutive ids must be joined by a mesh edge, and a closed loop
    // repeats its first node at the end
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
            const id_t tA = tEdge->node( 0 )->id();
            const id_t tB = tEdge->node( 1 )->id();
            tEdgeMap[ { std::min( tA, tB ), std::max( tA, tB ) } ] = tEdge ;
        }

        int tPairing = 0 ;
        const uint tN = aCycle.size();
        for ( uint k = 0; k + 1 < tN; ++k )
        {
            const id_t tFrom = aCycle[ k ];
            const id_t tTo   = aCycle[ k + 1 ];

            auto tIt = tEdgeMap.find(
                { std::min( tFrom, tTo ), std::max( tFrom, tTo ) } );
            EXPECT_TRUE( tIt != tEdgeMap.end() )
                << "cycle step " << tFrom << "->" << tTo
                << " is not a mesh edge" ;
            if ( tIt == tEdgeMap.end() ) continue ;

            Edge * tEdge = tIt->second ;

            auto tCoeff = tCoeffs.find( tEdge->index() );
            if ( tCoeff == tCoeffs.end() ) continue ;

            const int tSign = tEdge->node( 0 )->id() == tFrom ? 1 : -1 ;
            tPairing += tSign * tCoeff->second ;
        }
        return tPairing ;
    }

    // cofactor expansion along the first row, over the rows and columns still
    // in play. Exact in int: the unimodularity of the pairing matrix is an
    // integer statement, and a floating inverse would not decide it
    int
    minor_determinant( const Matrix< int > & aF,
        const std::vector< uint > & aRows, const std::vector< uint > & aCols )
    {
        const uint n = aRows.size();
        if ( n == 1 )
        {
            return aF( aRows[ 0 ], aCols[ 0 ] );
        }

        int tDeterminant = 0 ;
        int tSign = 1 ;

        for ( uint j = 0; j < n; ++j )
        {
            std::vector< uint > tRows( aRows.begin() + 1, aRows.end() );

            std::vector< uint > tCols ;
            for ( uint l = 0; l < n; ++l )
            {
                if ( l != j ) tCols.push_back( aCols[ l ] );
            }

            tDeterminant += tSign * aF( aRows[ 0 ], aCols[ j ] )
                * minor_determinant( aF, tRows, tCols );

            tSign = -tSign ;
        }

        return tDeterminant ;
    }

    int
    integer_determinant( const Matrix< int > & aF, const uint aK )
    {
        std::vector< uint > tRows ;
        std::vector< uint > tCols ;
        for ( uint i = 0; i < aK; ++i )
        {
            tRows.push_back( i );
            tCols.push_back( i );
        }
        return minor_determinant( aF, tRows, tCols );
    }

    std::string
    print_matrix( const Matrix< int > & aF, const uint aK )
    {
        std::string tText = "[" ;
        for ( uint i = 0; i < aK; ++i )
        {
            for ( uint h = 0; h < aK; ++h )
            {
                tText += " " + std::to_string( aF( i, h ) );
            }
            tText += i + 1 < aK ? " ;" : " ]" ;
        }
        return tText ;
    }

    // number of edges the generator actually acts on
    index_t
    support_size( Cochain * aGenerator )
    {
        index_t tCount = 0 ;
        for ( const auto & [ tIndex, tCoeff ] : aGenerator->getSimplicesMap() )
        {
            if ( tCoeff != 0 ) ++tCount ;
        }
        return tCount ;
    }

    Edge *
    find_edge_by_index( Mesh * aMesh, const index_t aIndex )
    {
        for ( Edge * tEdge : aMesh->edges() )
        {
            if ( tEdge->index() == aIndex ) return tEdge ;
        }
        return nullptr ;
    }

    // runs aCall and requires a std::runtime_error whose what() carries both
    // literal fragments of the BELFEM_ERROR text
    template< typename F >
    void
    expect_refusal( F && aCall,
        const std::string & aFirst, const std::string & aSecond )
    {
        bool tThrown = false ;
        try
        {
            aCall();
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true ;

            const std::string tMessage( aError.what() );
            EXPECT_NE( tMessage.find( aFirst ), std::string::npos )
                << "message was: " << tMessage ;
            EXPECT_NE( tMessage.find( aSecond ), std::string::npos )
                << "message was: " << tMessage ;
        }
        EXPECT_TRUE( tThrown )
            << "expected a refusal carrying \"" << aFirst
            << "\" and \"" << aSecond << "\"" ;
    }

    // the annulus coreduced, its cohomology generators replaced from the
    // homology and cleaned, owned as a unit: members are destroyed in reverse
    // declaration order, so the cohomology dies before the complex and the
    // mesh it points into
    struct CleanedAnnulus
    {
        std::unique_ptr< Mesh >              mMesh ;
        std::unique_ptr< SimplicialComplex > mComplex ;
        std::unique_ptr< Homology >          mHomology ;
        std::unique_ptr< Cohomology >        mCohomology ;

        CleanedAnnulus( const uint aN, const Algorithm aAlgorithm )
        {
            mMesh.reset( create_annulus( aN ) );

            flag_region( mMesh.get() );

            mComplex.reset( new SimplicialComplex( mMesh.get(), true ) );
            coreduce( mComplex.get(), aAlgorithm );

            mHomology.reset( new Homology( mComplex.get(), mMesh.get() ) );
            mCohomology.reset( new Cohomology( mComplex.get(), mMesh.get() ) );

            mCohomology->updatekGeneratorsFromHomology(
                mHomology->get_Generators()( 1 ), 1 );
            mCohomology->clean();
        }

        Cochain *
        generator()
        {
            return mCohomology->get_Generators()( 1 )( 0 );
        }
    };

    // hand the cut layer a coefficient that is not a unit, past clean(), and
    // require the refusal to name the edge it sits on
    void
    check_non_unit_is_refused( const int aCoefficient )
    {
        const uint tN = 8 ;

        CleanedAnnulus tFixture( tN, Algorithm::PellikkaGeneralized );

        ASSERT_EQ( tFixture.mCohomology->get_Generators()( 1 ).size(), 1u );

        Cochain * tGenerator = tFixture.generator();
        ASSERT_GT( tGenerator->getSimplicesMap().size(), 0u );

        // first support edge by smallest edge index ( OrderedMap = std::map,
        // deterministic )
        const index_t tEdgeIndex = tGenerator->getSimplicesMap().begin()->first ;
        const int     tG0        = tGenerator->getSimplicesMap().begin()->second ;
        ASSERT_TRUE( tG0 == 1 || tG0 == -1 );

        Edge * tEdge = find_edge_by_index( tFixture.mMesh.get(), tEdgeIndex );
        ASSERT_NE( tEdge, nullptr );

        tGenerator->setCoefficient( tEdgeIndex, aCoefficient );

        CutData tData( tFixture.mMesh.get(), tFixture.mCohomology.get(), 0 );
        tData.collect_edges();

        const index_t tNumEdges = tFixture.mMesh->edges().size();

        expect_refusal( [ & ]() { tData.collect_coefficients( tNumEdges ); },
            "edge " + std::to_string( tEdge->id() ) + " carries",
            "not a unit" );
    }

    // k generators against k hole rings: both are bases of the same free
    // group, so their pairing matrix is unimodular whichever representatives
    // the coreduction and the rectification picked. The test reads |det|, so
    // it is blind to a global sign for every k
    void
    check_punched_strip( const uint aK, const Algorithm aAlgorithm )
    {
        std::vector< std::vector< id_t > > tRings ;

        std::unique_ptr< Mesh > tMesh( create_punched_strip( aK, tRings ) );

        ASSERT_EQ( tRings.size(), aK );
        ASSERT_EQ( tMesh->nodes().size(), 5 * ( 3 * aK + 2 ) - aK );
        ASSERT_EQ( tMesh->elements().size(), 16 * aK + 8 );

        flag_region( tMesh.get() );

        std::unique_ptr< SimplicialComplex > tComplex(
            new SimplicialComplex( tMesh.get(), true ) );
        coreduce( tComplex.get(), aAlgorithm );

        // premise: the strip is a connected disk with k holes
        std::unique_ptr< Homology > tHomology(
            new Homology( tComplex.get(), tMesh.get() ) );
        ASSERT_EQ( tHomology->get_Generators()( 0 ).size(), 1u );
        ASSERT_EQ( tHomology->get_Generators()( 1 ).size(), aK );
        ASSERT_EQ( tHomology->get_Generators()( 2 ).size(), 0u );

        std::unique_ptr< Cohomology > tCohomology(
            new Cohomology( tComplex.get(), tMesh.get() ) );
        ASSERT_EQ( tCohomology->get_Generators()( 1 ).size(), aK );

        tCohomology->updatekGeneratorsFromHomology(
            tHomology->get_Generators()( 1 ), 1 );
        tCohomology->clean();

        Cell< Cochain * > & tGenerators = tCohomology->get_Generators()( 1 );
        ASSERT_EQ( tGenerators.size(), aK );

        for ( uint i = 0; i < aK; ++i )
        {
            for ( const auto & [ tIndex, tCoeff ] :
                tGenerators( i )->getSimplicesMap() )
            {
                EXPECT_TRUE( tCoeff == 1 || tCoeff == -1 )
                    << "generator " << i << ", edge index " << tIndex
                    << " carries coefficient " << tCoeff ;
            }
        }

        Matrix< int > tF( aK, aK, 0 );
        for ( uint i = 0; i < aK; ++i )
        {
            for ( uint h = 0; h < aK; ++h )
            {
                tF( i, h ) =
                    cycle_pairing( tMesh.get(), tGenerators( i ), tRings[ h ] );
            }
        }

        EXPECT_EQ( std::abs( integer_determinant( tF, aK ) ), 1 )
            << "pairing matrix " << print_matrix( tF, aK )
            << " is not unimodular" ;
    }
}

TEST( CutData, PlusTwoCoefficientIsRefusedNamingTheEdge )
{
    check_non_unit_is_refused( 2 );
}

TEST( CutData, MinusTwoCoefficientIsRefusedNamingTheEdge )
{
    check_non_unit_is_refused( -2 );
}

// the signed pairing the cut cases consume: the plus and minus bitsets carry
// the cochain's signs, not only its support
TEST( CutData, CutDataBitsetsReproduceTheSignedPairing )
{
    const uint tN = 8 ;

    CleanedAnnulus tFixture( tN, Algorithm::PellikkaGeneralized );

    ASSERT_EQ( tFixture.mCohomology->get_Generators()( 1 ).size(), 1u );
    Cochain * tGenerator = tFixture.generator();

    CutData tData( tFixture.mMesh.get(), tFixture.mCohomology.get(), 0 );
    tData.collect_edges();
    tData.collect_coefficients( tFixture.mMesh->edges().size() );

    // the bitsets against the cochain on every edge, not only where the cut
    // acts: the cochain is keyed by the edge's slot in the mesh container,
    // which is how CutData reads it, and an edge outside the support must
    // weigh zero
    Cell< Edge * > & tEdges = tFixture.mMesh->edges();
    for ( index_t tSlot = 0; tSlot < tEdges.size(); ++tSlot )
    {
        EXPECT_EQ( tData.weight( tEdges( tSlot ) ),
            tGenerator->getCoefficient( tSlot ) )
            << "edge " << tEdges( tSlot )->id() ;
    }

    int tPairing = 0 ;
    for ( Edge * tEdge : tFixture.mMesh->edges() )
    {
        const id_t tA = tEdge->node( 0 )->id();
        const id_t tB = tEdge->node( 1 )->id();

        // inner-ring edge: both endpoints have ids in 1..tN
        if ( tA > tN || tB > tN ) continue ;

        tPairing += inner_ring_sign( tA, tB, tN ) * tData.weight( tEdge );
    }

    EXPECT_EQ( tPairing,
        inner_ring_pairing( tFixture.mMesh.get(), tGenerator, tN ) );
    EXPECT_EQ( std::abs( tPairing ), 1 );
}

// rectification repairs a perturbed generator without paying for it in
// support: the cleaned representative acts on no more edges than the
// generator did before it was perturbed, and it is a unit representative of
// the same class
TEST( CutData, CleanedSupportIsNotLarger )
{
    const uint tN = 8 ;

    CleanedAnnulus tFixture( tN, Algorithm::Pellikka );

    Cell< Cochain * > & tGenerators =
        tFixture.mCohomology->get_Generators()( 1 );
    ASSERT_EQ( tGenerators.size(), 1u );
    Cochain * tGenerator = tGenerators( 0 );

    ASSERT_GT( tGenerator->getSimplicesMap().size(), 0u );
    const index_t tEdgeIndex = tGenerator->getSimplicesMap().begin()->first ;
    const int     tG0        = tGenerator->getSimplicesMap().begin()->second ;
    ASSERT_TRUE( tG0 == 1 || tG0 == -1 );

    Edge * tSupportEdge = find_edge_by_index( tFixture.mMesh.get(), tEdgeIndex );
    ASSERT_NE( tSupportEdge, nullptr );

    // the baseline is the generator before the perturbation, not after it
    const index_t tSupportBefore = support_size( tGenerator );

    // perturbation by an exact node coboundary: the same cohomology class.
    // The complex's incidence convention gives delta(n)(e) = -1 at e.node(0)
    // and +1 at e.node(1), so s = -g0 turns the support coefficient into 2*g0
    const index_t tNodeIndex = tSupportEdge->node( 0 )->index();
    const int     tS         = -tG0 ;

    const DynamicBitset & tDomain = tFixture.mComplex->original_edges();
    for ( Edge * tEdge : tFixture.mMesh->edges() )
    {
        if ( ! tDomain.test( tEdge->index() ) ) continue ;

        int tDelta = 0 ;
        if      ( tEdge->node( 0 )->index() == tNodeIndex ) tDelta = -1 ;
        else if ( tEdge->node( 1 )->index() == tNodeIndex ) tDelta =  1 ;
        if ( tDelta == 0 ) continue ;

        tGenerator->setCoefficient( tEdge->index(),
            tGenerator->getCoefficient( tEdge->index() ) + tS * tDelta );
    }

    // the forcing is real: a non-unit coefficient exists before the repair
    bool tHaveNonUnit = false ;
    for ( const auto & [ tIndex, tCoeff ] : tGenerator->getSimplicesMap() )
    {
        if ( tCoeff >= 2 || tCoeff <= -2 ) tHaveNonUnit = true ;
    }
    ASSERT_TRUE( tHaveNonUnit );

    tFixture.mCohomology->clean();

    // clean() mutates the generators in place, so the Cell still holds them
    ASSERT_EQ( tGenerators.size(), 1u );
    EXPECT_LE( support_size( tGenerators( 0 ) ), tSupportBefore );

    // a support that shrank by discarding the class does not count as a
    // repair: the result is a unit cochain still pairing to a unit with the
    // inner ring
    for ( const auto & [ tIndex, tCoeff ] :
        tGenerators( 0 )->getSimplicesMap() )
    {
        EXPECT_TRUE( tCoeff == 1 || tCoeff == -1 )
            << "edge index " << tIndex << " carries coefficient " << tCoeff ;
    }
    EXPECT_EQ( std::abs( inner_ring_pairing(
        tFixture.mMesh.get(), tGenerators( 0 ), tN ) ), 1 );
}

TEST( CutData, KHoleDiskPellikkaTwoHoles )
{
    check_punched_strip( 2, Algorithm::Pellikka );
}

TEST( CutData, KHoleDiskPellikkaThreeHoles )
{
    check_punched_strip( 3, Algorithm::Pellikka );
}

TEST( CutData, KHoleDiskCCRTwoHoles )
{
    check_punched_strip( 2, Algorithm::CCR );
}

TEST( CutData, KHoleDiskCCRThreeHoles )
{
    check_punched_strip( 3, Algorithm::CCR );
}

TEST( CutData, KHoleDiskPellikkaGeneralizedTwoHoles )
{
    check_punched_strip( 2, Algorithm::PellikkaGeneralized );
}

TEST( CutData, KHoleDiskPellikkaGeneralizedThreeHoles )
{
    check_punched_strip( 3, Algorithm::PellikkaGeneralized );
}

//------------------------------------------------------------------------------
// the guards of updatekGeneratorsFromHomology
//------------------------------------------------------------------------------

// the annulus has one generator; two conditions on it are one too many. The
// two suggested generators are the same chain, so without this guard the
// duplicate columns would reach the rank guard instead, with its own text
TEST( Cohomology, MoreConditionsThanGeneratorsIsRefused )
{
    CleanedAnnulus tFixture( 8, Algorithm::PellikkaGeneralized );

    ASSERT_EQ( tFixture.mCohomology->get_Generators()( 1 ).size(), 1u );
    ASSERT_GE( tFixture.mHomology->get_Generators()( 1 ).size(), 1u );

    Chain * tSuggested = tFixture.mHomology->get_Generators()( 1 )( 0 );
    Cell< Chain * > tTwo = { tSuggested, tSuggested };

    expect_refusal( [ & ]()
        {
            tFixture.mCohomology->updatekGeneratorsFromHomology( tTwo, 1 );
        },
        "defines 2 current or voltage conditions",
        "only provides 1 cohomology generator," );
}

// one condition on a chain that pairs to zero with the generator: few enough
// conditions, but the coupling matrix has rank 0
TEST( Cohomology, ConditionOnANullChainIsRefusedAsDependent )
{
    CleanedAnnulus tFixture( 8, Algorithm::PellikkaGeneralized );

    ASSERT_EQ( tFixture.mCohomology->get_Generators()( 1 ).size(), 1u );

    std::unique_ptr< Chain > tEmpty( new Chain( 1, tFixture.mMesh.get(), false, false ) );
    ASSERT_EQ( ( *tFixture.generator() )( tEmpty.get() ), 0 );

    Cell< Chain * > tOne = { tEmpty.get() };

    expect_refusal( [ & ]()
        {
            tFixture.mCohomology->updatekGeneratorsFromHomology( tOne, 1 );
        },
        "The 1 current or voltage conditions are not independent",
        "with rank 0" );
}

//------------------------------------------------------------------------------
// clean() on a generator without a unit representative
//------------------------------------------------------------------------------

namespace
{
    // a fresh directory under the system temp path that is the working
    // directory for the lifetime of this object. clean() writes error.exo
    // relative to the working directory, overwriting any file of that name,
    // so the refusal case runs in here. The previous working directory is
    // restored and the directory removed on scope exit
    class ScratchWorkDir
    {
        std::filesystem::path mPrevious ;
        std::filesystem::path mPath ;

    public:

        ScratchWorkDir() :
            mPrevious( std::filesystem::current_path() ),
            mPath( std::filesystem::temp_directory_path()
                / ( "belfem_test_cutdata_" + std::to_string( ::getpid() ) ) )
        {
            std::filesystem::remove_all( mPath );
            std::filesystem::create_directory( mPath );
            std::filesystem::current_path( mPath );
        }

        ~ScratchWorkDir()
        {
            // a failed restore would leave every later case of this binary in
            // a deleted directory, so it is reported, not swallowed
            std::error_code tRestore ;
            std::filesystem::current_path( mPrevious, tRestore );
            if ( tRestore )
            {
                ADD_FAILURE() << "could not restore the working directory "
                              << mPrevious << ": " << tRestore.message() ;
            }
            std::error_code tRemove ;
            std::filesystem::remove_all( mPath, tRemove );
        }

        const std::filesystem::path &
        path() const
        {
            return mPath ;
        }
    };

    // multiplies every coefficient of the annulus generator by aScale after
    // checking the premise the feasibility argument rests on: all aN edges
    // of the inner ring are constraint edges of clean(), which builds its
    // arcs from the complex's edges BEFORE coreduction. Returns the pairing
    // of the scaled generator with the inner ring
    int
    scale_generator( CleanedAnnulus & aFixture, const uint aN, const int aScale )
    {
        const DynamicBitset & tDomain = aFixture.mComplex->original_edges();

        uint tRingEdges = 0 ;
        for ( Edge * tEdge : aFixture.mMesh->edges() )
        {
            if ( tEdge->node( 0 )->id() > aN || tEdge->node( 1 )->id() > aN ) continue ;
            EXPECT_TRUE( tDomain.test( tEdge->index() ) )
                << "inner-ring edge " << tEdge->id() << " is not a constraint edge" ;
            ++tRingEdges ;
        }
        EXPECT_EQ( tRingEdges, aN );

        Cochain * tGenerator = aFixture.generator();

        // copied first: setCoefficient writes into the map being walked
        Cell< std::pair< index_t, int > > tEntries ;
        for ( const auto & tEntry : tGenerator->getSimplicesMap() )
        {
            tEntries.push( tEntry );
        }

        for ( const auto & [ tIndex, tCoeff ] : tEntries )
        {
            tGenerator->setCoefficient( tIndex, aScale * tCoeff );
        }

        return inner_ring_pairing( aFixture.mMesh.get(), tGenerator, aN );
    }
}

// A unit representative exists only if |< c, z >| <= length( z ) on every
// closed loop z of constraint edges; clean() tests this as a negative cycle
// of its difference constraints ( src/homology/doc/thin_cut_nonunit_rectification.md,
// sec. 2-3 ). The inner ring has 8 edges, and the cleaned generator pairs to
// +-1 with it, so nine times the generator pairs to +-9 and the ring is a
// negative cycle of weight 8 - 9 = -1. That makes the system infeasible
// whatever else the complex holds; which negative cycle clean() reports is
// its own choice, and only its form is checked below
TEST( Cohomology, InfeasibleGeneratorIsRefusedAsAThickCut )
{
#ifndef BELFEM_EXODUS
    GTEST_SKIP() << "the refusal writes its diagnostic mesh through Exodus" ;
#endif
    const uint tN = 8 ;

    CleanedAnnulus tFixture( tN, Algorithm::PellikkaGeneralized );
    ASSERT_EQ( tFixture.mCohomology->get_Generators()( 1 ).size(), 1u );

    ASSERT_EQ( std::abs( scale_generator( tFixture, tN, 9 ) ), 9 );

    ScratchWorkDir tWorkDir ;

    std::string tMessage ;
    try
    {
        tFixture.mCohomology->clean();
    }
    catch( const std::runtime_error & aError )
    {
        tMessage = aError.what();
    }
    ASSERT_FALSE( tMessage.empty() ) << "clean() accepted a generator with no unit representative" ;

    // the first generator is numbered 0
    EXPECT_NE( tMessage.find( "found cohomology generator 0," ), std::string::npos )
        << "message was: " << tMessage ;
    EXPECT_NE( tMessage.find( "only -1, 0, and 1" ), std::string::npos )
        << "message was: " << tMessage ;

    // the certificate: a loop length and the node pairs of its edges
    const std::size_t tLoop = tMessage.find( "The loop of " );
    ASSERT_NE( tLoop, std::string::npos ) << "message was: " << tMessage ;
    EXPECT_TRUE( std::isdigit( static_cast< unsigned char >(
        tMessage[ tLoop + std::string( "The loop of " ).size() ] ) ) )
        << "message was: " << tMessage ;
    EXPECT_NE( tMessage.find( "between nodes (" ), std::string::npos )
        << "message was: " << tMessage ;
    EXPECT_NE( tMessage.find( ")->(" ), std::string::npos )
        << "message was: " << tMessage ;

    // the diagnostic mesh landed in the scratch directory
    EXPECT_TRUE( std::filesystem::exists( tWorkDir.path() / "error.exo" ) );
}

// the control: twice the generator pairs to +-2 <= 8 with every loop around
// the hole, so a unit representative exists ( same note, same criterion ). It is two thin cuts, a unit
// cochain in the class of twice the generator, and the pairing with the ring
// is a class invariant, so it stays 2
TEST( Cohomology, DoubledGeneratorIsStillRectified )
{
    const uint tN = 8 ;

    CleanedAnnulus tFixture( tN, Algorithm::PellikkaGeneralized );
    ASSERT_EQ( tFixture.mCohomology->get_Generators()( 1 ).size(), 1u );

    ASSERT_EQ( std::abs( scale_generator( tFixture, tN, 2 ) ), 2 );

    ASSERT_NO_THROW( tFixture.mCohomology->clean() );

    Cochain * tGenerator = tFixture.generator();
    for ( const auto & [ tIndex, tCoeff ] : tGenerator->getSimplicesMap() )
    {
        EXPECT_TRUE( tCoeff == 1 || tCoeff == -1 )
            << "edge index " << tIndex << " carries coefficient " << tCoeff ;
    }
    EXPECT_EQ( std::abs( inner_ring_pairing( tFixture.mMesh.get(), tGenerator, tN ) ), 2 );
}
