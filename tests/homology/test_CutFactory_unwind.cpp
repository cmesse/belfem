/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * A CutFactory that dies after CutFactory::duplicate_and_relink_facets()
 * and before the thin-shell sidesets are restored, as on an exception
 * between the two, must free every facet exactly once: the original facets
 * through its temporary sideset, the duplicates through the mesh sideset.
 * The test asserts the state it builds; a double free is what Valgrind or
 * the allocator reports when the factory and then the mesh are destroyed.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <memory>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "cl_Topology.hpp"
#include "cl_Protoshell.hpp"
#include "cl_CutFactory.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    // two TET4 conductors, 1 = ( 1, 2, 3, 4 ) in block 1 and
    // 2 = ( 1, 2, 3, 5 ) in block 2, share the face 1-2-3, which is the one
    // facet of the thin-shell sideset 11 with element 1 as master and
    // element 2 as slave
    struct Shell
    {
        std::unique_ptr< Mesh > mMesh ;
        std::unique_ptr< Topology > mTopology ;
        Cell< Node * > mP ;

        Shell() :
            mMesh( new Mesh( 3, 0, false ) )
        {
            const real tX[ 5 ] = { 0.0, 1.0, 0.0, 0.0,  0.0 };
            const real tY[ 5 ] = { 0.0, 0.0, 1.0, 0.0,  0.0 };
            const real tZ[ 5 ] = { 0.0, 0.0, 0.0, 1.0, -1.0 };

            for ( index_t k = 0; k < 5; ++k )
            {
                Node * tNode = new Node( k + 1, tX[ k ], tY[ k ], tZ[ k ] );
                mMesh->nodes().push( tNode );
                mP.push( tNode );
            }

            Block * tB1 = new Block( 1, 1 );
            Block * tB2 = new Block( 2, 1 );
            tB1->set_domain_type( DomainType::Conductor );
            tB2->set_domain_type( DomainType::Conductor );

            Element * tMaster = this->tet( 1, tB1, { 0, 1, 2, 3 } );
            Element * tSlave  = this->tet( 2, tB2, { 0, 1, 2, 4 } );

            mMesh->blocks().push( tB1 );
            mMesh->blocks().push( tB2 );

            ElementFactory tFactory ;
            Facet * tFacet = new Facet( tFactory.create_element( ElementType::TRI3, 101 ) );
            tFacet->set_master( tMaster, this->facet_index( tMaster ) );
            tFacet->set_slave( tSlave, this->facet_index( tSlave ) );
            tFacet->set_sideset_id( 11 );

            SideSet * tSideSet = new SideSet( 11, 1 );
            tSideSet->insert_facet( tFacet );
            tSideSet->set_domain_type( DomainType::ThinShell );
            mMesh->sidesets().push( tSideSet );

            mMesh->finalize();

            mTopology.reset( new Topology( mMesh.get() ) );
            mTopology->run();
        }

        Element *
        tet( const id_t aID, Block * aBlock, const Cell< index_t > & aSlots )
        {
            ElementFactory tFactory ;
            Element * tTet = tFactory.create_element( ElementType::TET4, aID );
            for ( uint k = 0; k < 4; ++k )
            {
                tTet->insert_node( mP( aSlots( k ) ), k );
            }
            tTet->set_block_id( aBlock->id() );
            aBlock->insert_element( tTet );
            mMesh->elements().push( tTet );
            return tTet ;
        }

        // index of the facet of aElement whose corners are the nodes 1-2-3
        uint
        facet_index( Element * aElement )
        {
            Cell< Node * > tCorners ;
            for ( uint f = 0; f < aElement->number_of_facets(); ++f )
            {
                aElement->get_corner_nodes_of_facet( f, tCorners );
                Cell< id_t > tIDs ;
                for ( Node * tNode : tCorners )
                {
                    tIDs.push( tNode->id() );
                }
                std::sort( tIDs.begin(), tIDs.end() );
                if ( tIDs.size() == 3 && tIDs( 0 ) == 1 && tIDs( 1 ) == 2 && tIDs( 2 ) == 3 )
                {
                    return f ;
                }
            }
            ADD_FAILURE() << "no facet of element " << aElement->id() << " is 1-2-3" ;
            return 0 ;
        }
    };
}

//------------------------------------------------------------------------------

TEST( CutFactoryUnwind, DestroyedBeforeThinShellRestore )
{
    Shell tShell ;
    ASSERT_EQ( tShell.mTopology->groups( DomainType::ThinShell ).length(), ( index_t ) 1 );

    Facet * tOriginal = tShell.mMesh->sideset( 11 )->facets()( 0 );

    {
        Cell< Protoshell * > tProtoshells ;
        CutFactory tFactory( tShell.mMesh.get(), tShell.mTopology.get(), tProtoshells,
                             CutAlgorithm::PellikkaGeneralized, false );

        tFactory.duplicate_and_relink_facets();

        // the mesh sideset holds the master- and slave-side duplicates, the
        // factory's temporary sideset holds the original facet
        SideSet * tSideSet = tShell.mMesh->sideset( 11 );
        ASSERT_EQ( tSideSet->facets().size(), ( index_t ) 2 );
        EXPECT_NE( tSideSet->facets()( 0 ), tOriginal );
        EXPECT_NE( tSideSet->facets()( 1 ), tOriginal );
        EXPECT_EQ( tSideSet->facets()( 0 )->master()->id(), ( id_t ) 1 );
        EXPECT_EQ( tSideSet->facets()( 1 )->master()->id(), ( id_t ) 2 );
    }

    tShell.mMesh.reset();
}
