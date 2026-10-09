/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * ThinShellFactory::update_node_facet_tables on two hand-built TRI3 facets
 * whose nodes carry duplicates. Every original and every duplicate must end
 * up with exactly the facets that touch its original, in facet order, in a
 * container of exactly that size, with flag slot 2 clear. The node list is
 * given once as originals and once with duplicates in place of their
 * originals; both must give the same tables.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_Facet.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_ThinShellFactory.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    // originals 1-4; node 2 has the duplicates 12 and 22, node 3 has the
    // duplicate 13. Facet 101 holds 1-2-3, facet 102 holds 12-13-4, so the
    // facet lookup must go through the originals
    struct Patch
    {
        Cell< Node * > mO ;
        Node * mD2a = nullptr ;
        Node * mD2b = nullptr ;
        Node * mD3  = nullptr ;
        Cell< Facet * > mFacets ;

        Patch()
        {
            for ( index_t k = 0; k < 4; ++k )
            {
                Node * tNode = new Node( k + 1, ( real ) k, 0.0, 0.0 );
                tNode->set_index( k );
                mO.push( tNode );
            }

            mD2a = this->duplicate( mO( 1 ), 12, 4 );
            mD2b = this->duplicate( mO( 1 ), 22, 5 );
            mD3  = this->duplicate( mO( 2 ), 13, 6 );

            mO( 1 )->allocate_duplicate_container( 2 );
            mO( 1 )->add_duplicate( mD2a );
            mO( 1 )->add_duplicate( mD2b );
            mO( 2 )->allocate_duplicate_container( 1 );
            mO( 2 )->add_duplicate( mD3 );

            mFacets.push( this->facet( 101, { mO( 0 ), mO( 1 ), mO( 2 ) } ) );
            mFacets.push( this->facet( 102, { mD2a, mD3, mO( 3 ) } ) );
        }

        ~Patch()
        {
            for ( Facet * tFacet : mFacets )
            {
                delete tFacet ;
            }
            delete mD2a ;
            delete mD2b ;
            delete mD3 ;
            for ( Node * tNode : mO )
            {
                delete tNode ;
            }
        }

        Patch( const Patch & ) = delete ;
        Patch & operator=( const Patch & ) = delete ;

        Node *
        duplicate( Node * aOriginal, const id_t aID, const index_t aIndex )
        {
            Node * tNode = new Node( aID, aOriginal->x(), aOriginal->y(), aOriginal->z() );
            tNode->set_original( aOriginal );
            tNode->set_index( aIndex );
            return tNode ;
        }

        Facet *
        facet( const id_t aID, const Cell< Node * > & aNodes )
        {
            ElementFactory tFactory ;
            Element * tTri = tFactory.create_element( ElementType::TRI3, aID );
            for ( uint k = 0; k < 3; ++k )
            {
                tTri->insert_node( aNodes( k ), k );
            }
            return new Facet( tTri );
        }

        void
        expect_facets( Node * aNode, const Cell< id_t > & aIDs )
        {
            ASSERT_EQ( aNode->number_of_facets(), ( uint ) aIDs.size() ) << "node " << aNode->id() ;
            ASSERT_EQ( aNode->facet_capacity(), ( uint16_t ) aIDs.size() ) << "node " << aNode->id() ;
            for ( index_t k = 0; k < aIDs.size(); ++k )
            {
                EXPECT_EQ( aNode->facet( k )->id(), aIDs( k ) ) << "node " << aNode->id() << " slot " << k ;
            }
            EXPECT_FALSE( aNode->is_flagged( 2 ) ) << "node " << aNode->id() ;
        }

        void
        expect_tables()
        {
            this->expect_facets( mO( 0 ), { 101 } );
            this->expect_facets( mO( 1 ), { 101, 102 } );
            this->expect_facets( mD2a, { 101, 102 } );
            this->expect_facets( mD2b, { 101, 102 } );
            this->expect_facets( mO( 2 ), { 101, 102 } );
            this->expect_facets( mD3, { 101, 102 } );
            this->expect_facets( mO( 3 ), { 102 } );
        }
    };
}

//------------------------------------------------------------------------------

TEST( ThinShellNodeFacetTables, OriginalNodes )
{
    Patch tPatch ;
    Cell< Node * > tNodes = { tPatch.mO( 0 ), tPatch.mO( 1 ), tPatch.mO( 2 ), tPatch.mO( 3 ) };

    EXPECT_NO_THROW( ThinShellFactory::update_node_facet_tables( tNodes, tPatch.mFacets ) );
    tPatch.expect_tables();

    for ( index_t k = 0; k < 4; ++k )
    {
        EXPECT_EQ( tPatch.mO( k )->index(), k );
    }
}

//------------------------------------------------------------------------------

// the duplicates 12 and 13 stand in for their originals
TEST( ThinShellNodeFacetTables, DuplicateNodes )
{
    Patch tPatch ;
    Cell< Node * > tNodes = { tPatch.mO( 0 ), tPatch.mD2a, tPatch.mD3, tPatch.mO( 3 ) };

    EXPECT_NO_THROW( ThinShellFactory::update_node_facet_tables( tNodes, tPatch.mFacets ) );
    tPatch.expect_tables();

    EXPECT_EQ( tPatch.mD2a->index(), ( index_t ) 4 );
    EXPECT_EQ( tPatch.mD3->index(), ( index_t ) 6 );
}
