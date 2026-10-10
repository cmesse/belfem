/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * See LICENSE for the full text of the license.
 */

/**
 * Unit test for the order in which CutProcessor appends the nodes it creates
 * to the mesh: the abstract nodes first, then every cut copy in ascending
 * id. The test pins that postcondition on one run; the bucket order of the
 * cut sets' hash maps is not under its control.
 *
 * The fixture is the punched strip of test_CutData.cpp with two holes, so
 * the cohomology has two generators and the copies fall into at least two
 * cut sets ( two duplicate maps ). The strip is coreduced, the cohomology
 * generators replaced from the homology and cleaned, exactly as the
 * k-hole tests do; nothing here reads a chain or a cochain.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <map>
#include <memory>
#include <vector>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Block.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Homology.hpp"
#include "cl_Cohomology.hpp"
#include "cl_CutProcessor.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    // a ( 3k+2 ) x 5 node lattice, so ( 3k+1 ) x 4 unit quads, each split
    // into two counterclockwise TRI3, with k rectangular holes: the 2 x 2 quad block at
    // columns 3h+1, 3h+2 and rows 1, 2 is left out. A disk with k holes,
    // b1 = k. The node at the centre of a hole is not created. The caller owns
    // the mesh ( the builder of test_CutData.cpp, without the hole rings )
    Mesh *
    create_punched_strip( const uint aK )
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

        Cell< Node * > tLattice( tNumCols * tNumRows, nullptr );

        id_t tNodeID = 1 ;
        for ( uint j = 0; j < tNumRows; ++j )
        {
            for ( uint i = 0; i < tNumCols; ++i )
            {
                if ( tIsHoleCentre( i, j ) ) continue ;

                Node * tNode = new Node( tNodeID++, ( real ) i, ( real ) j );
                tNodes.push( tNode );
                tLattice( i + j * tNumCols ) = tNode ;
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

                Node * tA = tLattice( i     + j       * tNumCols );
                Node * tB = tLattice( i + 1 + j       * tNumCols );
                Node * tC = tLattice( i + 1 + ( j+1 ) * tNumCols );
                Node * tD = tLattice( i     + ( j+1 ) * tNumCols );

                Element * tElement = tFactory.create_element( ElementType::TRI3, tID++ );
                tElement->insert_node( tA, 0 );
                tElement->insert_node( tB, 1 );
                tElement->insert_node( tC, 2 );
                tElement->set_block_id( 1 );
                tElements.push( tElement );
                tBlock->insert_element( tElement );

                Element * tElement2 = tFactory.create_element( ElementType::TRI3, tID++ );
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

        return tMesh ;
    }

    // the region flagging CutFactory::compute_cohomologies does in 2D
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

    // the strip with its cleaned cohomology; members die in reverse order,
    // so the cohomology goes before the complex and the mesh it points into
    struct PunchedStrip
    {
        std::unique_ptr< Mesh >              mMesh ;
        std::unique_ptr< SimplicialComplex > mComplex ;
        std::unique_ptr< Homology >          mHomology ;
        std::unique_ptr< Cohomology >        mCohomology ;

        PunchedStrip( const uint aK )
        {
            mMesh.reset( create_punched_strip( aK ) );
            flag_region( mMesh.get() );
            mComplex.reset( new SimplicialComplex( mMesh.get(), true ) );
            mComplex->coreduce_complexPellikkaGeneralized();
            mHomology.reset( new Homology( mComplex.get(), mMesh.get() ) );
            mCohomology.reset( new Cohomology( mComplex.get(), mMesh.get() ) );
            mCohomology->updatekGeneratorsFromHomology(
                mHomology->get_Generators()( 1 ), 1 );
            mCohomology->clean();
        }
    };
}

//------------------------------------------------------------------------------

// the processor appends the abstract nodes, then every cut copy in ascending
// id, so a copy's position in the mesh follows its allocation and not the
// bucket order of a cut set's hash map. The two-hole strip has two
// generators, and its copies span at least two cut sets; the periodic branch
// of the allocation is outside this fixture
TEST( CutProcessor, CopiesAreAppendedInAscendingIdAfterTheAbstractNodes )
{
    const uint tK = 2 ;

    PunchedStrip tFixture( tK );
    ASSERT_EQ( tFixture.mCohomology->get_Generators()( 1 ).size(), tK );

    Mesh * tMesh = tFixture.mMesh.get();

    const index_t tNumOriginal = tMesh->nodes().size();
    const id_t    tMaxID       = tMesh->max_node_id();

    // the originals as they sit before the processor runs
    Cell< Node * > tOriginals( tNumOriginal, nullptr );
    for ( index_t k = 0; k < tNumOriginal; ++k )
    {
        tOriginals( k ) = tMesh->nodes()( k );
    }

    Vector< id_t > tPhiBlocks( 1 );
    tPhiBlocks( 0 ) = 1 ;
    Vector< id_t > tNone ;

    CutProcessor tProcessor( tMesh, tFixture.mCohomology.get(),
        tPhiBlocks, tNone, tNone, tNone, tNone );

    Cell< Node * > & tNodes = tMesh->nodes();

    ASSERT_EQ( tProcessor.abstract_nodes().size(), tK );
    ASSERT_GT( tNodes.size(), tNumOriginal + tK )
        << "no cut copy was appended" ;

    // the originals keep their slots
    for ( index_t k = 0; k < tNumOriginal; ++k )
    {
        EXPECT_EQ( tNodes( k ), tOriginals( k ) ) << "slot " << k ;
    }

    // the abstract nodes come first, in cut order
    for ( uint c = 0; c < tK; ++c )
    {
        EXPECT_EQ( tNodes( tNumOriginal + c ), tProcessor.abstract_nodes()( c ) )
            << "abstract node " << c ;
        EXPECT_EQ( tNodes( tNumOriginal + c )->id(), tMaxID + 1 + c );
    }

    // every copy after them, in strictly ascending id, each hanging on one
    // original plus the abstract nodes of its pattern; the copies are
    // counted per pattern
    std::map< std::vector< id_t >, uint > tPatterns ;

    for ( index_t k = tNumOriginal + tK; k < tNodes.size(); ++k )
    {
        Node * tCopy = tNodes( k );

        EXPECT_GT( tCopy->id(), tMaxID + tK ) << "slot " << k ;

        if ( k > tNumOriginal + tK )
        {
            EXPECT_GT( tCopy->id(), tNodes( k - 1 )->id() )
                << "copy " << tCopy->id() << " at slot " << k
                << " follows copy " << tNodes( k - 1 )->id() ;
        }

        ASSERT_TRUE( tCopy->is_hanging() ) << "copy " << tCopy->id() ;

        std::vector< id_t > tPattern ;
        uint tNumOriginals = 0 ;

        for ( uint s = 0; s < tCopy->number_of_sources(); ++s )
        {
            const id_t tID = tCopy->source( s )->id();

            if ( tID > tMaxID && tID <= tMaxID + tK )
            {
                tPattern.push_back( tID );
            }
            else
            {
                EXPECT_LE( tID, tMaxID ) << "copy " << tCopy->id()
                    << " hangs on " << tID << ", which is no original" ;
                ++tNumOriginals ;
            }
        }

        EXPECT_EQ( tNumOriginals, 1u ) << "copy " << tCopy->id() ;
        EXPECT_FALSE( tPattern.empty() ) << "copy " << tCopy->id() ;

        std::sort( tPattern.begin(), tPattern.end() );
        ++tPatterns[ tPattern ];
    }

    EXPECT_GE( tPatterns.size(), 2u )
        << "the copies span fewer than two cut sets, so the test does not "
           "exercise the order across duplicate maps" ;

    // a duplicate map with one entry iterates in order whatever its bucket
    // order; the fixture must hold a set with at least two copies
    uint tMaxCopies = 0 ;
    for ( const auto & tPair : tPatterns )
    {
        tMaxCopies = std::max( tMaxCopies, tPair.second );
    }
    EXPECT_GE( tMaxCopies, 2u )
        << "no cut set holds two copies, so no duplicate map could iterate "
           "out of order" ;
}
