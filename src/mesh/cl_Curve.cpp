/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "commtools.hpp"
#include "cl_Curve.hpp"

#include "cl_Element_Factory.hpp"
#include "cl_Mesh.hpp"
#include "fn_reverse.hpp"

namespace belfem
{
    namespace mesh
    {
        Curve::Curve( const id_t aID, SideSet * aSideSetA, SideSet * aSideSetB  ) :
            mID( aID ),
            mSideSetA( aSideSetA ),
            mSideSetB( aSideSetB ),
            mElementType( interpolation_order_numeric( mSideSetA->element_type() ) == 2 &&
                   interpolation_order_numeric( mSideSetB->element_type() ) == 2 ?
                   ElementType::LINE3 : ElementType::LINE2 )
        {

        }

        Curve::Curve( const id_t aID, const ElementType aElementType ) :
            mID( aID ),
            mSideSetA( nullptr ),
            mSideSetB( nullptr ),
            mElementType( aElementType )
        {

        }

        Curve::~Curve()
        {
            for ( Segment * tSegment : mSegments )
            {
                delete tSegment ;
            }
        }

        void
        Curve::save( const std::string & aFileName )
        {
            if ( comm_rank() == 0 )
            {
                // create a temporary mesh
                Mesh * tMesh = new Mesh( 3 );

                // collect nodes
                index_t tCount = 0 ;
                for ( Segment * tSegment : mSegments )
                {
                    tCount += tSegment->number_of_nodes() ;
                }

                Cell< Node * > tOriginalNodes( tCount, nullptr );
                tCount = 0 ;
                for ( Segment * tSegment : mSegments )
                {
                    for ( uint k=0; k<tSegment->number_of_nodes(); ++k )
                    {
                        tOriginalNodes( tCount++ ) = tSegment->node( k );
                    }
                }
                unique( tOriginalNodes );


                // copy the nodes
                Cell< Node * > & tNodes = tMesh->nodes() ;
                tNodes.set_size( tOriginalNodes.size(), nullptr );

                tCount = 0 ;
                for ( Node * tNode : tOriginalNodes )
                {
                    tNodes( tCount++ ) = new Node( tNode->id(), tNode->x(), tNode->y(), tNode->z() );
                }

                // create a node map
                Map< id_t, Node * > tMap ;
                for ( Node * tNode : tNodes )
                {
                    tMap[ tNode->id() ] = tNode ;
                }

                // create a new block
                Block * tBlock = new Block( 1, mSegments.size() );
                Cell< Element * > & tElements = tBlock->elements() ;

                // create the elements and link them
                ElementFactory tFactory ;

                tCount = 0 ;
                for ( Segment * tSegment : mSegments )
                {
                    Element * tElement = tFactory.create_element( this->element_type(),  tSegment->id() );

                    for ( uint k=0; k<tSegment->number_of_nodes(); ++k )
                    {
                        tElement->insert_node( tMap( tSegment->node( k )->id()) , k );
                    }
                    tElements( tCount++ ) = tElement ;
                }
                tMesh->add_block( tBlock );

                tMesh->finalize();

                // save the mesh
                tMesh->save( aFileName );

                delete tMesh ;
            }
        }

        void
        Curve::reverse()
        {
            belfem::reverse( mSegments );
            belfem::reverse( mNodes );
            mS = belfem::reverse( mS );

            BELFEM_ASSERT( mElementType == ElementType::LINE2 || mElementType == ElementType::LINE3,
                "Higher order lines are not supported" );

            for ( Segment * tSegment : mSegments )
            {
                Node * tSwap = tSegment->node( 0 );
                tSegment->element()->insert_node( tSegment->element()->node( 1 ), 0 );
                tSegment->element()->insert_node( tSwap, 1 );
            }
        }

        void Curve::sideset_a( SideSet *aSideset )
        {
           mSideSetA = aSideset ;
        }

        void Curve::sideset_b( SideSet *aSideset )
        {
            mSideSetB = aSideset ;
        }

        void
        Curve::assign_edges()
        {
            BELFEM_ASSERT( mEdges.size() == 0, "Edges are already assigned" );
            mEdges.set_size( mSegments.size(), nullptr );

            index_t tCount = 0 ;

            for ( Segment * tSegment : mSegments )
            {
                tSegment->clear_edge() ;

                Node * tA = tSegment->node( 0 );
                Node * tB = tSegment->node( 1 );

                Edge * tEdge = nullptr ;
                for ( uint e=0; e<tA->number_of_edges(); ++e )
                {
                    tEdge = tA->edge( e );

                    if ( tEdge->node( 0 )->id() == tB->id() || tEdge->node( 1 )->id() == tB->id() )
                    {
                        tSegment->insert_edge(  tEdge );
                        mEdges( tCount++ ) = tEdge ;
                        break ;
                    }
                    tEdge = nullptr ;
                }
                BELFEM_ASSERT(  tEdge != nullptr, "Edge not found for segment %lu", ( long unsigned int ) tSegment->id() );

            }

        }

    }
}
