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
 * Unit tests for the interface duplicates of InterfaceProcessor on a
 * hand-built TET4 mesh. Two conductor tets C1 and C2 share the face
 * ( n, m, q ); the edge n-m lies in the plane z = 1, where an air tet sits on
 * the top face of each conductor tet. The nodes n and m are on the line where
 * both conductor blocks meet the air: they must get exactly one conductor-side
 * copy, referenced by every conductor block that touches them, while the air
 * keeps the original and the copy hangs on it with weight 1. The face shared
 * by C1 and C2 must then reference identical nodes from both sides.
 *
 * Layouts: one air sideset per conductor block; one air sideset over both
 * blocks; a third conductor block that touches n only through the edge n-q of
 * the shared face and owns no air facet; a ferro block in place of C2, whose
 * copy must differ from the conductor copy; a periodic variant whose node
 * pairs n-m and a1-a2 lie in one set and across two sets; two ferro blocks,
 * whose shared face keeps its split; and a buffer block in place of C1, which
 * keeps private copies.
 *
 * The facets carry the conductor ( or, against the air, the ferro ) element
 * as master and the other side as slave, as fix_facet_masters would leave
 * them. The homology main forces the throwing reaction, so the refusals are
 * observable as std::runtime_error in both builds.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "cl_Mesh_Periodicity.hpp"
#include "cl_Topology.hpp"
#include "cl_InterfaceProcessor.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    enum class Layout
    {
        SideSetPerBlock,
        OneSideSet,
        EdgeOnlyBlock,
        Ferro,
        FerroFerro,
        Buffer
    };

    // node slots: n = 0, m = 1, a1 = 2, a2 = 3, q = 4, t = 5, b1 = 6, b2 = 7.
    // C1 = ( n, m, a1, q ), C2 = ( n, m, a2, q ), A1 = ( n, m, a1, t ),
    // A2 = ( n, m, a2, t ), C3 = ( n, q, b1, b2 ). Blocks: C1 -> 1
    // ( conductor, or ferro, or buffer ), C2 -> 2 ( conductor, or ferro ),
    // A1 and A2 -> 3 ( air ), C3 -> 4 ( conductor ).
    // The periodic variant moves a1 into y = 0 and a2 into y = 1, the planes
    // of n and m
    struct Junction
    {
        std::unique_ptr< Mesh > mMesh ;
        std::unique_ptr< Topology > mTopology ;
        Cell< Node * > mP ;
        Element * mC1 = nullptr ;
        Element * mC2 = nullptr ;
        Element * mA1 = nullptr ;
        Element * mA2 = nullptr ;
        Element * mC3 = nullptr ;
        Cell< Node * > mAbstractNodes ;
        index_t mNumOriginalNodes = 0 ;

        Junction( const Layout aLayout, const bool aPeriodic = false ) :
            mMesh( new Mesh( 3, 0, false ) )
        {
            const real tA1y = aPeriodic ? 0.0 : 0.5 ;
            const real tA2y = aPeriodic ? 1.0 : 0.5 ;

            const real tX[ 8 ] = { 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, -0.5, 0.5 };
            const real tY[ 8 ] = { 0.0, 1.0, tA1y, tA2y, 0.5, 0.5, -1.0, -1.0 };
            const real tZ[ 8 ] = { 1.0, 1.0, 1.0, 1.0, 0.0, 2.0, 0.5, 0.5 };

            for ( index_t k = 0; k < 8; ++k )
            {
                Node * tNode = new Node( k + 1, tX[ k ], tY[ k ], tZ[ k ] );
                mMesh->nodes().push( tNode );
                mP.push( tNode );
            }
            mNumOriginalNodes = mP.size() ;

            const bool tEdgeOnly = aLayout == Layout::EdgeOnlyBlock ;

            Block * tB1 = new Block( 1, 1 );
            Block * tB2 = new Block( 2, 1 );
            Block * tB3 = new Block( 3, 2 );
            tB1->set_domain_type( aLayout == Layout::FerroFerro ? DomainType::Ferro
                                : aLayout == Layout::Buffer ? DomainType::Buffer : DomainType::Conductor );
            tB2->set_domain_type( aLayout == Layout::Ferro || aLayout == Layout::FerroFerro
                                ? DomainType::Ferro : DomainType::Conductor );
            tB3->set_domain_type( DomainType::Air );

            mC1 = this->tet( 1, tB1, { 0, 1, 2, 4 } );
            mC2 = this->tet( 2, tB2, { 0, 1, 3, 4 } );
            mA1 = this->tet( 3, tB3, { 0, 1, 2, 5 } );
            mA2 = this->tet( 4, tB3, { 0, 1, 3, 5 } );

            mMesh->blocks().push( tB1 );
            mMesh->blocks().push( tB2 );
            mMesh->blocks().push( tB3 );

            if ( tEdgeOnly )
            {
                Block * tB4 = new Block( 4, 1 );
                tB4->set_domain_type( DomainType::Conductor );
                mC3 = this->tet( 5, tB4, { 0, 4, 6, 7 } );
                mMesh->blocks().push( tB4 );
            }

            Facet * tTop1 = this->facet( 101, mC1, mA1, { 0, 1, 2 } );
            Facet * tTop2 = this->facet( 102, mC2, mA2, { 0, 1, 3 } );

            if ( aLayout == Layout::OneSideSet )
            {
                this->sideset( 11, { tTop1, tTop2 } );
            }
            else
            {
                this->sideset( 11, { tTop1 } );
                this->sideset( 12, { tTop2 } );
            }

            if ( aLayout == Layout::Ferro || aLayout == Layout::FerroFerro )
            {
                this->sideset( 13, { this->facet( 103, mC1, mC2, { 0, 1, 4 } ) } );
            }

            mMesh->finalize();
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

        // index of the facet of aElement whose corners are the slots aSlots
        uint
        facet_index( Element * aElement, const Cell< index_t > & aSlots )
        {
            Cell< id_t > tWanted ;
            for ( index_t s : aSlots )
            {
                tWanted.push( mP( s )->id() );
            }
            std::sort( tWanted.begin(), tWanted.end() );

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
                if ( tIDs.vector_data() == tWanted.vector_data() )
                {
                    return f ;
                }
            }
            ADD_FAILURE() << "no facet of element " << aElement->id() << " matches" ;
            return 0 ;
        }

        Facet *
        facet( const id_t aID, Element * aMaster, Element * aSlave, const Cell< index_t > & aSlots )
        {
            ElementFactory tFactory ;
            Facet * tFacet = new Facet( tFactory.create_element( ElementType::TRI3, aID ) );
            tFacet->set_master( aMaster, this->facet_index( aMaster, aSlots ) );
            tFacet->set_slave( aSlave, this->facet_index( aSlave, aSlots ) );
            return tFacet ;
        }

        void
        sideset( const id_t aID, const Cell< Facet * > & aFacets )
        {
            SideSet * tSideSet = new SideSet( aID, aFacets.size() );
            for ( Facet * tFacet : aFacets )
            {
                tFacet->set_sideset_id( aID );
                tSideSet->insert_facet( tFacet );
            }
            mMesh->sidesets().push( tSideSet );
        }

        void
        run()
        {
            mTopology.reset( new Topology( mMesh.get() ) );
            mTopology->run();

            id_t tMaxNodeID = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                tMaxNodeID = std::max( tMaxNodeID, tNode->id() );
            }

            InterfaceProcessor tProcessor( mMesh.get(), mTopology.get(), mAbstractNodes,
                                           mMesh->number_of_sidesets(), tMaxNodeID );
        }

        // the node aElement holds where it was built with the original aSlot
        Node *
        held( Element * aElement, const index_t aSlot )
        {
            for ( uint k = 0; k < aElement->number_of_nodes(); ++k )
            {
                Node * tNode = aElement->node( k );
                if ( tNode == mP( aSlot ) || ( tNode->is_duplicate() && tNode->original() == mP( aSlot ) ) )
                {
                    return tNode ;
                }
            }
            ADD_FAILURE() << "element " << aElement->id() << " holds neither node " << mP( aSlot )->id() << " nor a copy" ;
            return nullptr ;
        }

        // the copies appended by the processor, each exactly once
        index_t
        number_of_new_nodes()
        {
            Cell< Node * > & tNodes = mMesh->nodes() ;
            Cell< Node * > tNew ;
            for ( index_t k = mNumOriginalNodes; k < tNodes.size(); ++k )
            {
                tNew.push( tNodes( k ) );
            }
            std::sort( tNew.begin(), tNew.end() );
            EXPECT_TRUE( std::adjacent_find( tNew.begin(), tNew.end() ) == tNew.end() )
                << "a copy is appended to the mesh more than once" ;
            return tNew.size();
        }

        // aCopy is the one duplicate of aSlot and hangs on it with weight 1
        void
        expect_hanging_copy( Node * aCopy, const index_t aSlot )
        {
            Node * tOriginal = mP( aSlot );
            ASSERT_NE( aCopy, tOriginal );
            EXPECT_EQ( aCopy->original(), tOriginal );
            ASSERT_EQ( aCopy->number_of_sources(), ( uint ) 1 );
            EXPECT_EQ( aCopy->source( 0 ), tOriginal );
            EXPECT_EQ( aCopy->weight( 0 ), 1.0 );
        }

        // the triple-line node aSlot: one conductor copy held by all
        // conductor elements, the air on the original
        void
        expect_one_conductor_copy( const index_t aSlot, const Cell< Element * > & aConductors )
        {
            Node * tOriginal = mP( aSlot );
            Node * tCopy = this->held( aConductors( 0 ), aSlot );
            ASSERT_NE( tCopy, nullptr );
            EXPECT_NE( tCopy, tOriginal ) << "conductor element " << aConductors( 0 )->id() << " kept the original" ;

            for ( Element * tElement : aConductors )
            {
                EXPECT_EQ( this->held( tElement, aSlot ), tCopy )
                    << "conductor element " << tElement->id() << " holds another node for " << tOriginal->id() ;
            }

            EXPECT_EQ( this->held( mA1, aSlot ), tOriginal );
            EXPECT_EQ( this->held( mA2, aSlot ), tOriginal );

            EXPECT_EQ( tOriginal->number_of_duplicates(), ( uint ) 1 )
                << "node " << tOriginal->id() << " has more than one copy" ;
            this->expect_hanging_copy( tCopy, aSlot );
        }

        // C1 and C2 reference identical nodes on their shared face
        void
        expect_conforming_shared_face()
        {
            for ( index_t tSlot : { ( index_t ) 0, ( index_t ) 1, ( index_t ) 4 } )
            {
                EXPECT_EQ( this->held( mC1, tSlot )->id(), this->held( mC2, tSlot )->id() )
                    << "shared face splits at node " << mP( tSlot )->id() ;
            }
        }
    };

    // the Hesse forms of a hand-wired Periodicity, which only the
    // PeriodicityFactory may set otherwise. Access to the protected accessor
    // is checked where the pointer to member is formed, inside this derived
    // class, so calling it through a Periodicity & is well-formed
    struct PeriodicityPlanes : public Periodicity
    {
        static void
        set( Periodicity & aPeriodicity, const Vector< real > & aMaster, const Vector< real > & aSlave )
        {
            Vector< real > & ( Periodicity::*tMaster )() = &PeriodicityPlanes::master_hesse ;
            Vector< real > & ( Periodicity::*tSlave )() = &PeriodicityPlanes::slave_hesse ;
            ( aPeriodicity.*tMaster )() = aMaster ;
            ( aPeriodicity.*tSlave )() = aSlave ;
        }
    };

    // periodic planes y = 0 ( master ) and y = 1 ( slave ); aLinks are
    // ( from, to ) slot pairs, applied as Node::set_periodic( from -> to )
    Periodicity *
    wire_periodicity( Junction & aJunction, const Cell< std::pair< index_t, index_t > > & aLinks )
    {
        Periodicity * tPeriodicity = new Periodicity( aJunction.mMesh.get() );

        Vector< real > tMaster( 4, 0.0 );
        tMaster( 1 ) = 1.0 ;
        Vector< real > tSlave( 4, 0.0 );
        tSlave( 1 ) = 1.0 ;
        tSlave( 3 ) = 1.0 ;
        PeriodicityPlanes::set( *tPeriodicity, tMaster, tSlave );

        for ( auto tLink : aLinks )
        {
            aJunction.mP( tLink.first )->set_periodic( aJunction.mP( tLink.second ) );
        }

        aJunction.mMesh->set_periodicity( tPeriodicity );
        return tPeriodicity ;
    }

    bool
    contains( const std::string & aText, const std::string & aPart )
    {
        return aText.find( aPart ) != std::string::npos ;
    }
}

//------------------------------------------------------------------------------

// one air sideset per conductor block: two interface sets share n and m
TEST( InterfaceProcessor, SideSetPerBlockSharesTheConductorCopy )
{
    Junction tJunction( Layout::SideSetPerBlock );
    tJunction.run();

    const Cell< Element * > tConductors = { tJunction.mC1, tJunction.mC2 };
    tJunction.expect_one_conductor_copy( 0, tConductors );
    tJunction.expect_one_conductor_copy( 1, tConductors );
    tJunction.expect_conforming_shared_face();

    // n, m, a1 and a2 get one copy each
    EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 4 );
}

//------------------------------------------------------------------------------

// one air sideset over both conductor blocks, C1 owning its first facet
TEST( InterfaceProcessor, OneSideSetRelinksEveryMasterBlock )
{
    Junction tJunction( Layout::OneSideSet );
    tJunction.run();

    const Cell< Element * > tConductors = { tJunction.mC1, tJunction.mC2 };
    tJunction.expect_one_conductor_copy( 0, tConductors );
    tJunction.expect_one_conductor_copy( 1, tConductors );
    tJunction.expect_conforming_shared_face();

    // the visualization split exists on the second block too
    EXPECT_NE( tJunction.held( tJunction.mC2, 3 ), tJunction.mP( 3 ) );
    EXPECT_EQ( tJunction.held( tJunction.mA2, 3 ), tJunction.mP( 3 ) );

    EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 4 );
}

//------------------------------------------------------------------------------

// C3 touches n only through the edge n-q of the shared face and owns no air
// facet; it must hold the conductor copy of n all the same
TEST( InterfaceProcessor, EdgeOnlyConductorBlockGetsTheCopy )
{
    Junction tJunction( Layout::EdgeOnlyBlock );
    tJunction.run();

    tJunction.expect_one_conductor_copy( 0, { tJunction.mC1, tJunction.mC2, tJunction.mC3 } );
    tJunction.expect_one_conductor_copy( 1, { tJunction.mC1, tJunction.mC2 } );
    tJunction.expect_conforming_shared_face();

    // q is no interface node
    EXPECT_EQ( tJunction.held( tJunction.mC3, 4 ), tJunction.mP( 4 ) );
    EXPECT_EQ( tJunction.held( tJunction.mC1, 4 ), tJunction.mP( 4 ) );

    EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 4 );
}

//------------------------------------------------------------------------------

// a ferro block in place of C2: the conductor and the ferro side keep
// distinct copies of n, both distinct from the air original
TEST( InterfaceProcessor, ConductorAndFerroCopiesStayDistinct )
{
    Junction tJunction( Layout::Ferro );
    tJunction.run();

    for ( index_t tSlot : { ( index_t ) 0, ( index_t ) 1 } )
    {
        Node * tOriginal = tJunction.mP( tSlot );
        Node * tConductor = tJunction.held( tJunction.mC1, tSlot );
        Node * tFerro = tJunction.held( tJunction.mC2, tSlot );

        EXPECT_NE( tConductor, tOriginal );
        EXPECT_NE( tFerro, tOriginal );
        EXPECT_NE( tConductor, tFerro );
        EXPECT_EQ( tJunction.held( tJunction.mA1, tSlot ), tOriginal );
        EXPECT_EQ( tJunction.held( tJunction.mA2, tSlot ), tOriginal );

        EXPECT_EQ( tOriginal->number_of_duplicates(), ( uint ) 2 )
            << "node " << tOriginal->id() << " has a copy per interface set" ;
        tJunction.expect_hanging_copy( tConductor, tSlot );
        tJunction.expect_hanging_copy( tFerro, tSlot );
    }

    // the conductor|ferro sideset copies q on the conductor side only
    EXPECT_NE( tJunction.held( tJunction.mC1, 4 ), tJunction.mP( 4 ) );
    EXPECT_EQ( tJunction.held( tJunction.mC2, 4 ), tJunction.mP( 4 ) );

    // conductor: n, m, a1, q ; ferro: n, m, a2
    EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 7 );
}

//------------------------------------------------------------------------------

// n-m is a pair inside both sets, a1-a2 a pair across the two sets: each pair
// of copies is tied once and backed up once
TEST( InterfaceProcessor, PeriodicPairsAreTiedOnce )
{
    Junction tJunction( Layout::SideSetPerBlock, true );
    Periodicity * tPeriodicity = wire_periodicity( tJunction, { { 0, 1 }, { 1, 0 }, { 2, 3 }, { 3, 2 } } );

    EXPECT_NO_THROW( tJunction.run() );

    const Cell< Element * > tConductors = { tJunction.mC1, tJunction.mC2 };
    tJunction.expect_one_conductor_copy( 0, tConductors );
    tJunction.expect_one_conductor_copy( 1, tConductors );

    Node * tN = tJunction.held( tJunction.mC1, 0 );
    Node * tM = tJunction.held( tJunction.mC2, 1 );
    Node * tA1 = tJunction.held( tJunction.mC1, 2 );
    Node * tA2 = tJunction.held( tJunction.mC2, 3 );

    EXPECT_EQ( tN->periodic(), tM );
    EXPECT_EQ( tM->periodic(), tN );
    EXPECT_EQ( tA1->periodic(), tA2 );
    EXPECT_EQ( tA2->periodic(), tA1 );

    // the originals keep their own pairing
    EXPECT_EQ( tJunction.mP( 0 )->periodic(), tJunction.mP( 1 ) );
    EXPECT_EQ( tJunction.mP( 2 )->periodic(), tJunction.mP( 3 ) );

    ASSERT_TRUE( tPeriodicity->has_node_pair_backup() );
    tPeriodicity->restore_node_pairs();
    EXPECT_EQ( tPeriodicity->master_nodes().size(), ( index_t ) 2 )
        << "a pair of copies is backed up more than once" ;
}

//------------------------------------------------------------------------------

// n and a2 are partners; a1 names a2 as its partner too. The copy of a2 is
// tied to the copy of n inside the second set, so pairing the copy of a1 with
// it is a conflict
TEST( InterfaceProcessor, PeriodicPartnerTiedElsewhereIsRefused )
{
    Junction tJunction( Layout::SideSetPerBlock, true );
    wire_periodicity( tJunction, { { 0, 3 }, { 3, 0 }, { 2, 3 } } );

    std::string tMessage ;
    try
    {
        tJunction.run();
    }
    catch ( const std::runtime_error & aError )
    {
        tMessage = aError.what();
    }

    EXPECT_TRUE( contains( tMessage, "already tied to another node" ) ) << tMessage ;
}

//------------------------------------------------------------------------------

// two ferro blocks under the air, with a ferro|ferro sideset between them: the
// ferro copies stay private to their sets, so the face keeps its split
TEST( InterfaceProcessor, FerroFerroFaceKeepsItsSplit )
{
    Junction tJunction( Layout::FerroFerro );
    EXPECT_NO_THROW( tJunction.run() );

    for ( index_t tSlot : { ( index_t ) 0, ( index_t ) 1, ( index_t ) 4 } )
    {
        EXPECT_NE( tJunction.held( tJunction.mC1, tSlot ), tJunction.mP( tSlot ) )
            << "the master side of the ferro|ferro face kept original " << tJunction.mP( tSlot )->id() ;
        EXPECT_NE( tJunction.held( tJunction.mC1, tSlot )->id(), tJunction.held( tJunction.mC2, tSlot )->id() )
            << "the ferro|ferro face is no longer split at node " << tJunction.mP( tSlot )->id() ;
    }

    for ( index_t tSlot : { ( index_t ) 0, ( index_t ) 1 } )
    {
        EXPECT_EQ( tJunction.held( tJunction.mA1, tSlot ), tJunction.mP( tSlot ) );
        EXPECT_EQ( tJunction.held( tJunction.mA2, tSlot ), tJunction.mP( tSlot ) );
    }
}

//------------------------------------------------------------------------------

// a buffer block beside a conductor block under the air: the buffer|air set
// is processed without a refusal and keeps private copies
TEST( InterfaceProcessor, BufferMasterKeepsPrivateCopies )
{
    Junction tJunction( Layout::Buffer );
    EXPECT_NO_THROW( tJunction.run() );

    for ( index_t tSlot : { ( index_t ) 0, ( index_t ) 1 } )
    {
        Node * tOriginal = tJunction.mP( tSlot );
        Node * tBuffer = tJunction.held( tJunction.mC1, tSlot );
        Node * tConductor = tJunction.held( tJunction.mC2, tSlot );

        EXPECT_NE( tBuffer, tOriginal );
        EXPECT_NE( tConductor, tOriginal );
        EXPECT_NE( tBuffer, tConductor );
        EXPECT_EQ( tJunction.held( tJunction.mA1, tSlot ), tOriginal );
        EXPECT_EQ( tJunction.held( tJunction.mA2, tSlot ), tOriginal );
        EXPECT_EQ( tOriginal->number_of_duplicates(), ( uint ) 2 );
        tJunction.expect_hanging_copy( tBuffer, tSlot );
        tJunction.expect_hanging_copy( tConductor, tSlot );
    }

    // buffer: n, m, a1 ; conductor: n, m, a2
    EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 6 );
}
