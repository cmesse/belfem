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
 * whose shared face keeps its split; a buffer block in place of C1, which
 * keeps private copies; a coil block in place of C2, whose air set is
 * decoupled; a ferro block in place of C2 under ONE sideset, which mixes
 * master domain types and is refused; and that sideset with an air|air
 * facet first or last, refused in either order; a sideset of an air|air
 * and a conductor|conductor facet, which this processor never duplicates at
 * and therefore leaves alone; two ferro blocks under ONE air sideset,
 * refused because private copies serve one block; a coil that is only the
 * slave of a later facet, refused; a conductor|conductor or a boundary facet
 * beside the conductor|air facet in one sideset, refused in either order
 * because an interface sideset holds interface facets only; and two ferro
 * sidesets with swapped masters on their common face, kept apart. The
 * periodic variant on the
 * buffer and the coil layout puts the pair a1-a2 across sets of different
 * master domain type or treatment, which the processor refuses by name. The
 * nodes can be pushed into the mesh in reverse order, which must not move a
 * copy's id or position.
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
        Buffer,
        Coil,
        MixedSideSet,
        MixedAirFirst,
        MixedAirLast,
        InactiveSideSet,
        MultiBlockFerro,
        CoilSlaveLater,
        ConductorFacet,
        BoundaryFacet,
        OverlappingFerro
    };

    // node slots: n = 0, m = 1, a1 = 2, a2 = 3, q = 4, t = 5, b1 = 6, b2 = 7.
    // C1 = ( n, m, a1, q ), C2 = ( n, m, a2, q ), A1 = ( n, m, a1, t ),
    // A2 = ( n, m, a2, t ), C3 = ( n, q, b1, b2 ). Blocks: C1 -> 1
    // ( conductor, or ferro, or buffer ), C2 -> 2 ( conductor, ferro or coil ),
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

        Junction( const Layout aLayout, const bool aPeriodic = false, const bool aReverseNodes = false,
                  const bool aSwapFacets = false ) :
            mMesh( new Mesh( 3, 0, false ) )
        {
            const real tA1y = aPeriodic ? 0.0 : 0.5 ;
            const real tA2y = aPeriodic ? 1.0 : 0.5 ;

            const real tX[ 8 ] = { 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, -0.5, 0.5 };
            const real tY[ 8 ] = { 0.0, 1.0, tA1y, tA2y, 0.5, 0.5, -1.0, -1.0 };
            const real tZ[ 8 ] = { 1.0, 1.0, 1.0, 1.0, 0.0, 2.0, 0.5, 0.5 };

            for ( index_t k = 0; k < 8; ++k )
            {
                mP.push( new Node( k + 1, tX[ k ], tY[ k ], tZ[ k ] ) );
            }
            // the mesh may hold the nodes in reverse order; ids stay
            for ( index_t k = 0; k < 8; ++k )
            {
                mMesh->nodes().push( mP( aReverseNodes ? 7 - k : k ) );
            }
            mNumOriginalNodes = mP.size() ;

            const bool tMixed = aLayout == Layout::MixedSideSet
                             || aLayout == Layout::MixedAirFirst
                             || aLayout == Layout::MixedAirLast ;

            const bool tEdgeOnly = aLayout == Layout::EdgeOnlyBlock ;

            Block * tB1 = new Block( 1, 1 );
            Block * tB2 = new Block( 2, 1 );
            Block * tB3 = new Block( 3, 2 );
            tB1->set_domain_type( aLayout == Layout::FerroFerro || aLayout == Layout::MultiBlockFerro || aLayout == Layout::OverlappingFerro ? DomainType::Ferro
                                : aLayout == Layout::Buffer ? DomainType::Buffer : DomainType::Conductor );
            tB2->set_domain_type( aLayout == Layout::Ferro || aLayout == Layout::FerroFerro || aLayout == Layout::MultiBlockFerro || aLayout == Layout::OverlappingFerro || tMixed ? DomainType::Ferro
                                : aLayout == Layout::Coil || aLayout == Layout::CoilSlaveLater ? DomainType::Coil : DomainType::Conductor );
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

            if ( aLayout == Layout::OneSideSet || aLayout == Layout::MixedSideSet || aLayout == Layout::MultiBlockFerro )
            {
                this->sideset( 11, { tTop1, tTop2 } );
            }
            else if ( aLayout == Layout::CoilSlaveLater || aLayout == Layout::ConductorFacet )
            {
                // the face C1|C2 beside the conductor|air facet, in either order
                Facet * tShared = this->facet( 104, mC1, mC2, { 0, 1, 4 } );
                this->sideset( 11, aSwapFacets ? Cell< Facet * >( { tShared, tTop1 } ) : Cell< Facet * >( { tTop1, tShared } ) );
                this->sideset( 12, { tTop2 } );
            }
            else if ( aLayout == Layout::OverlappingFerro )
            {
                // each ferro sideset holds its air facet and the common face
                // with itself as master: same key blocks, different masters
                this->sideset( 11, { tTop1, this->facet( 103, mC1, mC2, { 0, 1, 4 } ) } );
                this->sideset( 12, { tTop2, this->facet( 106, mC2, mC1, { 0, 1, 4 } ) } );
            }
            else if ( aLayout == Layout::BoundaryFacet )
            {
                // the face ( m, a1, q ) of C1 has no neighbor: a boundary facet
                Facet * tBoundary = this->boundary_facet( 105, mC1, { 1, 2, 4 } );
                this->sideset( 11, aSwapFacets ? Cell< Facet * >( { tBoundary, tTop1 } ) : Cell< Facet * >( { tTop1, tBoundary } ) );
                this->sideset( 12, { tTop2 } );
            }
            else if ( aLayout == Layout::MixedAirFirst || aLayout == Layout::MixedAirLast )
            {
                // A1 and A2 share the face ( n, m, t ): an air|air facet
                Facet * tAirAir = this->facet( 103, mA1, mA2, { 0, 1, 5 } );
                if ( aLayout == Layout::MixedAirFirst )
                {
                    this->sideset( 11, { tAirAir, tTop1 } );
                }
                else
                {
                    this->sideset( 11, { tTop1, tAirAir } );
                }
                this->sideset( 12, { tTop2 } );
            }
            else if ( aLayout == Layout::InactiveSideSet )
            {
                // neither facet is duplicated at: air|air and conductor|conductor
                this->sideset( 11, { this->facet( 103, mA1, mA2, { 0, 1, 5 } ),
                                     this->facet( 104, mC1, mC2, { 0, 1, 4 } ) } );
                this->sideset( 12, { tTop1 } );
                this->sideset( 13, { tTop2 } );
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

        // a facet with a master only
        Facet *
        boundary_facet( const id_t aID, Element * aMaster, const Cell< index_t > & aSlots )
        {
            ElementFactory tFactory ;
            Facet * tFacet = new Facet( tFactory.create_element( ElementType::TRI3, aID ) );
            tFacet->set_master( aMaster, this->facet_index( aMaster, aSlots ) );
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

        // position of a node in the mesh's node container
        index_t
        index_in_mesh( Node * aNode )
        {
            Cell< Node * > & tNodes = mMesh->nodes() ;
            for ( index_t k = 0; k < tNodes.size(); ++k )
            {
                if ( tNodes( k ) == aNode ) return k ;
            }
            ADD_FAILURE() << "node " << aNode->id() << " is not in the mesh" ;
            return 0 ;
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

    // the message of the refusal aJunction.run() raises, empty if none
    std::string
    refusal( Junction & aJunction )
    {
        try
        {
            aJunction.run();
        }
        catch ( const std::runtime_error & aError )
        {
            return aError.what();
        }
        return "" ;
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

    const std::string tMessage = refusal( tJunction );
    EXPECT_TRUE( contains( tMessage, "already tied to another node" ) ) << tMessage ;
}

//------------------------------------------------------------------------------

// a1 lies on the buffer|air set only and a2 on the conductor|air set only:
// the pair crosses a change of master domain type and is refused by name,
// with both node ids and both types
TEST( InterfaceProcessor, PeriodicPairAcrossDomainTypesIsRefused )
{
    Junction tJunction( Layout::Buffer, true );
    wire_periodicity( tJunction, { { 0, 1 }, { 1, 0 }, { 2, 3 }, { 3, 2 } } );

    const std::string tMessage = refusal( tJunction );
    EXPECT_TRUE( contains( tMessage, "change of domain type" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "node 3: Buffer" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "node 4: Conductor" ) ) << tMessage ;
    EXPECT_FALSE( contains( tMessage, "not duplicated in any interface set" ) ) << tMessage ;
}

//------------------------------------------------------------------------------

// a1 lies on the tie-weight conductor|air set only and a2 on the decoupled
// coil|air set only: the pair crosses a change of treatment and is refused
// by name, with the label of each side
TEST( InterfaceProcessor, PeriodicPairAcrossTreatmentsIsRefused )
{
    Junction tJunction( Layout::Coil, true );
    wire_periodicity( tJunction, { { 0, 1 }, { 1, 0 }, { 2, 3 }, { 3, 2 } } );

    const std::string tMessage = refusal( tJunction );
    EXPECT_TRUE( contains( tMessage, "interfaces of different treatment" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "move the periodic plane away from the coil" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "node 3: tie-weight" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "node 4: decoupled" ) ) << tMessage ;
    EXPECT_FALSE( contains( tMessage, "not duplicated in any interface set" ) ) << tMessage ;
}

//------------------------------------------------------------------------------

// a1 names q as its partner, and q lies on no interface set at all: this is
// the one case that keeps the plain "not duplicated" refusal
TEST( InterfaceProcessor, PeriodicPartnerOutsideEveryInterfaceSetIsRefused )
{
    Junction tJunction( Layout::SideSetPerBlock, true );
    wire_periodicity( tJunction, { { 2, 4 } } );

    const std::string tMessage = refusal( tJunction );
    EXPECT_TRUE( contains( tMessage, "periodic partner 5 of interface node 3 was not duplicated in any interface set" ) ) << tMessage ;
}

//------------------------------------------------------------------------------

// two ferro blocks under the air: both ferro|air sets have the same master
// type and private copies, so a1-a2 is tied across the sets; n-m is tied
// inside each of the three sets, the ferro|ferro one included
TEST( InterfaceProcessor, PeriodicPairsAcrossFerroSetsAreTied )
{
    Junction tJunction( Layout::FerroFerro, true );
    Periodicity * tPeriodicity = wire_periodicity( tJunction, { { 0, 1 }, { 1, 0 }, { 2, 3 }, { 3, 2 } } );

    EXPECT_NO_THROW( tJunction.run() );

    Node * tA1 = tJunction.held( tJunction.mC1, 2 );
    Node * tA2 = tJunction.held( tJunction.mC2, 3 );
    ASSERT_NE( tA1, nullptr );
    ASSERT_NE( tA2, nullptr );
    EXPECT_NE( tA1, tJunction.mP( 2 ) );
    EXPECT_NE( tA2, tJunction.mP( 3 ) );
    EXPECT_EQ( tA1->periodic(), tA2 );
    EXPECT_EQ( tA2->periodic(), tA1 );

    ASSERT_TRUE( tPeriodicity->has_node_pair_backup() );
    tPeriodicity->restore_node_pairs();
    EXPECT_EQ( tPeriodicity->master_nodes().size(), ( index_t ) 4 )
        << "expected one backup entry per pair of copies: n-m in three sets and a1-a2 across two" ;
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

//------------------------------------------------------------------------------

// one sideset over a conductor|air and a ferro|air facet: refused by name,
// before any copy is appended. The accepted multi-block case is
// OneSideSetRelinksEveryMasterBlock above, one sideset over two conductor
// blocks
TEST( InterfaceProcessor, MixedMasterSideSetIsRefused )
{
    Junction tJunction( Layout::MixedSideSet );

    const std::string tMessage = refusal( tJunction );
    EXPECT_TRUE( contains( tMessage, "interface sideset 11 mixes master domain types" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "facet 101 ( element 1 ) has a Conductor ( 5 ) master" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "facet 102 ( element 2 ) a Ferro ( 3 ) master" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "one master domain type" ) ) << tMessage ;
    EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 0 );
}

//------------------------------------------------------------------------------

// an air|air facet beside the conductor|air facet in one sideset: refused
// whether the air facet comes first ( the sideset would otherwise be dropped
// by its first facet ) or last
TEST( InterfaceProcessor, MixedSideSetWithAirFacetIsRefusedInEitherOrder )
{
    for ( Layout tLayout : { Layout::MixedAirFirst, Layout::MixedAirLast } )
    {
        Junction tJunction( tLayout );

        const std::string tMessage = refusal( tJunction );
        EXPECT_TRUE( contains( tMessage, "interface sideset 11 mixes master domain types" ) ) << tMessage ;
        EXPECT_TRUE( contains( tMessage, "Air ( 1 )" ) ) << tMessage ;
        EXPECT_TRUE( contains( tMessage, "Conductor ( 5 )" ) ) << tMessage ;
        EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 0 );
    }
}

//------------------------------------------------------------------------------

// copies are numbered set by set and in ascending id of the original, and
// appended in that order: the ids and the mesh positions must not move when
// the mesh holds the nodes in another order, and on the forward build they
// are the ones the rule predicts. Reversing the node order changes only the
// insertion order of the hash maps; it cannot force a particular hash order,
// so for the copies the absolute ids are the check, not the agreement of the
// two builds. SideSetPerBlock: sets "5" ( sideset 11:
// n, m, a1 ) then "6" ( sideset 12: a2 ). Ferro: sets "3" ( sideset 13:
// n, m, q, shared conductor copies ), "5" ( sideset 11: a1 ), "6" ( sideset
// 12: private ferro copies of n, m, a2 ). The periodic layout also backs up
// its pairs in that order
TEST( InterfaceProcessor, CopyNumberingIsIndependentOfNodeOrder )
{
    // SideSetPerBlock
    for ( bool tReverse : { false, true } )
    {
        Junction tJunction( Layout::SideSetPerBlock, false, tReverse );
        tJunction.run();

        EXPECT_EQ( tJunction.held( tJunction.mC1, 0 )->id(), ( id_t ) 9 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC1, 1 )->id(), ( id_t ) 10 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC1, 2 )->id(), ( id_t ) 11 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC2, 3 )->id(), ( id_t ) 12 ) << "reverse " << tReverse ;
        for ( id_t tID : { ( id_t ) 9, ( id_t ) 10, ( id_t ) 11, ( id_t ) 12 } )
        {
            EXPECT_EQ( tJunction.mMesh->nodes()( tID - 1 )->id(), tID ) << "reverse " << tReverse ;
        }
    }

    // Ferro: conductor copies first, in set order, then the private ferro copies
    for ( bool tReverse : { false, true } )
    {
        Junction tJunction( Layout::Ferro, false, tReverse );
        tJunction.run();

        EXPECT_EQ( tJunction.held( tJunction.mC1, 0 )->id(), ( id_t ) 9 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC1, 1 )->id(), ( id_t ) 10 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC1, 4 )->id(), ( id_t ) 11 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC1, 2 )->id(), ( id_t ) 12 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC2, 0 )->id(), ( id_t ) 13 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC2, 1 )->id(), ( id_t ) 14 ) << "reverse " << tReverse ;
        EXPECT_EQ( tJunction.held( tJunction.mC2, 3 )->id(), ( id_t ) 15 ) << "reverse " << tReverse ;
        for ( id_t tID = 9; tID <= 15; ++tID )
        {
            EXPECT_EQ( tJunction.mMesh->nodes()( tID - 1 )->id(), tID ) << "reverse " << tReverse ;
        }
        // every held copy sits at the position its id says
        for ( index_t tSlot : { ( index_t ) 0, ( index_t ) 1, ( index_t ) 2, ( index_t ) 4 } )
        {
            Node * tCopy = tJunction.held( tJunction.mC1, tSlot );
            EXPECT_EQ( tJunction.index_in_mesh( tCopy ), ( index_t ) tCopy->id() - 1 ) << "reverse " << tReverse ;
        }
        for ( index_t tSlot : { ( index_t ) 0, ( index_t ) 1, ( index_t ) 3 } )
        {
            Node * tCopy = tJunction.held( tJunction.mC2, tSlot );
            EXPECT_EQ( tJunction.index_in_mesh( tCopy ), ( index_t ) tCopy->id() - 1 ) << "reverse " << tReverse ;
        }
    }

    // periodic: the backup lists the pairs in the same order on both builds,
    // and on the forward build it is the rule's order — the n-m copies tied
    // inside set "5" ( n on the master plane y = 0 ), then the a1-a2 copies
    // across sets "5" and "6"
    Cell< id_t > tForward ;
    for ( bool tReverse : { false, true } )
    {
        Junction tJunction( Layout::SideSetPerBlock, true, tReverse );
        Periodicity * tPeriodicity = wire_periodicity( tJunction, { { 0, 1 }, { 1, 0 }, { 2, 3 }, { 3, 2 } } );
        tJunction.run();

        ASSERT_TRUE( tPeriodicity->has_node_pair_backup() );
        tPeriodicity->restore_node_pairs();

        Cell< id_t > tPairs ;
        for ( index_t k = 0; k < tPeriodicity->master_nodes().size(); ++k )
        {
            tPairs.push( tPeriodicity->master_nodes()( k )->id() );
            tPairs.push( tPeriodicity->slave_nodes()( k )->id() );
        }
        if ( ! tReverse )
        {
            tForward = tPairs ;
            const Cell< id_t > tExpected = { 9, 10, 11, 12 };
            EXPECT_EQ( tPairs.vector_data(), tExpected.vector_data() ) << "the periodic backup is not in the rule's order" ;
        }
        else
        {
            EXPECT_EQ( tPairs.vector_data(), tForward.vector_data() ) << "the periodic backup order depends on the node order" ;
        }
    }
}

//------------------------------------------------------------------------------

// a sideset of an air|air and a conductor|conductor facet mixes master types
// but is no interface sideset: this processor duplicates at neither facet,
// so the sideset is left alone and the junction is processed as usual
TEST( InterfaceProcessor, SideSetWithoutInterfaceFacetIsLeftAlone )
{
    Junction tJunction( Layout::InactiveSideSet );
    EXPECT_NO_THROW( tJunction.run() );

    const Cell< Element * > tConductors = { tJunction.mC1, tJunction.mC2 };
    tJunction.expect_one_conductor_copy( 0, tConductors );
    tJunction.expect_one_conductor_copy( 1, tConductors );
    tJunction.expect_conforming_shared_face();
    EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 4 );
}

//------------------------------------------------------------------------------

// two ferro blocks under one air sideset: the copies of a ferro set are
// private to one block, so the sideset is refused by name, nothing appended
TEST( InterfaceProcessor, MultiBlockFerroSideSetIsRefused )
{
    Junction tJunction( Layout::MultiBlockFerro );

    const std::string tMessage = refusal( tJunction );
    EXPECT_TRUE( contains( tMessage, "interface sideset 11 has Ferro ( 3 ) masters on more than one block" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "facet 101 ( element 1 ) on block 1, facet 102 ( element 2 ) on block 2" ) ) << tMessage ;
    EXPECT_TRUE( contains( tMessage, "a sideset of its own" ) ) << tMessage ;
    EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 0 );
}

//------------------------------------------------------------------------------

// C2 is a coil and the face C1|C2 sits beside the conductor|air facet in one
// sideset: with the air facet first the processor refuses the mixed
// treatment; with the coil facet first the topology refuses the contact
// before the processor runs. Either way the run stops
TEST( InterfaceProcessor, CoilSlaveOnLaterFacetIsRefused )
{
    {
        Junction tJunction( Layout::CoilSlaveLater, false, false, false );
        const std::string tMessage = refusal( tJunction );
        EXPECT_TRUE( contains( tMessage, "interface sideset 11 mixes facets that touch a coil block with facets that do not" ) ) << tMessage ;
        EXPECT_TRUE( contains( tMessage, "facet 101 ( element 1 ) touches no coil, facet 104 ( element 1 ) touches a coil" ) ) << tMessage ;
        EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 0 );
    }
    {
        Junction tJunction( Layout::CoilSlaveLater, false, false, true );
        const std::string tMessage = refusal( tJunction );
        EXPECT_TRUE( contains( tMessage, "a conductor and a coil must not touch" ) ) << tMessage ;
        EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 0 );
    }
}

//------------------------------------------------------------------------------

// the face C1|C2 ( both conductor ) beside the conductor|air facet in one
// sideset, in either order: refused, an interface sideset holds interface
// facets only
TEST( InterfaceProcessor, ConductorFacetInInterfaceSideSetIsRefused )
{
    for ( bool tSwap : { false, true } )
    {
        Junction tJunction( Layout::ConductorFacet, false, false, tSwap );
        const std::string tMessage = refusal( tJunction );
        EXPECT_TRUE( contains( tMessage, "interface sideset 11 mixes interface facets with facets that are no interface" ) ) << "swap " << tSwap << ": " << tMessage ;
        EXPECT_TRUE( contains( tMessage, "facet 101 ( element 1 ) is duplicated at, facet 104 ( element 1 ) is not ( neither side is air or ferro )" ) ) << tMessage ;
        EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 0 ) << "swap " << tSwap ;
    }
}

//------------------------------------------------------------------------------

// a boundary facet of C1 beside the conductor|air facet in one sideset, in
// either order: refused the same way
TEST( InterfaceProcessor, BoundaryFacetInInterfaceSideSetIsRefused )
{
    for ( bool tSwap : { false, true } )
    {
        Junction tJunction( Layout::BoundaryFacet, false, false, tSwap );
        const std::string tMessage = refusal( tJunction );
        EXPECT_TRUE( contains( tMessage, "interface sideset 11 mixes interface facets with facets that are no interface" ) ) << "swap " << tSwap << ": " << tMessage ;
        EXPECT_TRUE( contains( tMessage, "facet 105 ( element 1 ) is not ( a boundary facet without a neighbor )" ) ) << tMessage ;
        EXPECT_EQ( tJunction.number_of_new_nodes(), ( index_t ) 0 ) << "swap " << tSwap ;
    }
}

//------------------------------------------------------------------------------

// two ferro sidesets over the same blocks with swapped masters on their
// common face: keyed apart by their master blocks, each keeps its private
// copies and the face stays split
TEST( InterfaceProcessor, OverlappingFerroSideSetsStayApart )
{
    Junction tJunction( Layout::OverlappingFerro );
    EXPECT_NO_THROW( tJunction.run() );

    for ( index_t tSlot : { ( index_t ) 0, ( index_t ) 1, ( index_t ) 4 } )
    {
        EXPECT_NE( tJunction.held( tJunction.mC1, tSlot ), tJunction.mP( tSlot ) );
        EXPECT_NE( tJunction.held( tJunction.mC2, tSlot ), tJunction.mP( tSlot ) );
        EXPECT_NE( tJunction.held( tJunction.mC1, tSlot )->id(), tJunction.held( tJunction.mC2, tSlot )->id() )
            << "the ferro|ferro face is no longer split at node " << tJunction.mP( tSlot )->id() ;
    }
}
