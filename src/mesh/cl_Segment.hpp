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

#ifndef CL_SEGMENT_HPP
#define CL_SEGMENT_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    namespace mesh
    {
        class Segment : public Vertex
        {
            // wrapped 1D element
            Element * mElement;

            Edge * mEdge = nullptr ;

        public:

            // bring parent edge(uint) overloads into scope
            using Vertex::edge;

            Segment( Element * aElement );

            ~Segment() override;

            id_t
            id() const;

            ElementType
            type() const;

            uint
            number_of_nodes() const override;

            Node *
            node( const uint aIndex ) override;

            const Node *
            node( const uint aIndex ) const override ;

            Edge *
            edge();

            const Edge *
            edge() const ;

            void
            insert_edge( Edge * aEdge );

            void
            clear_edge();

            int
            edge_direction() const ;

//------------------------------------------------------------------------------

            /**
             * expose the element pointer
             */
            Element *
            element();

//------------------------------------------------------------------------------

            /**
             * inherited, but forbidden function
             */
            Element *
            element( const uint aIndex ) override;

//------------------------------------------------------------------------------

            /**
             * inherited, but forbidden function
             */
            const Element *
            element( const uint aIndex ) const override ;

//------------------------------------------------------------------------------

            size_t
            memory() const;

//------------------------------------------------------------------------------
        };

        inline id_t
        Segment::id() const
        {
            return mElement->id();
        }

//------------------------------------------------------------------------------

        inline ElementType
        Segment::type() const
        {
            return mElement->type();
        }

//------------------------------------------------------------------------------

        inline uint
        Segment::number_of_nodes() const
        {
            return mElement->number_of_nodes();
        }

//------------------------------------------------------------------------------

        inline Node *
        Segment::node( const uint aIndex )
        {
            return mElement->node( aIndex );
        }

        inline const Node *
        Segment::node( const uint aIndex ) const
        {
            return mElement->node( aIndex );
        }

//------------------------------------------------------------------------------

        inline Element *
        Segment::element()
        {
            return mElement;
        }

        inline Edge *
        Segment::edge()
        {
            BELFEM_ASSERT( mEdge != nullptr, "Edge has not been assigned to segment %lu", ( long unsigned int ) mElement->id() );
            return mEdge;
        }

        inline const Edge *
        Segment::edge() const
        {
            BELFEM_ASSERT( mEdge != nullptr, "Edge hast not been assigned to segment %lu", ( long unsigned int ) mElement->id() );
            return mEdge;
        }

        inline void
        Segment::insert_edge( Edge * aEdge )
        {
            mEdge = aEdge;
        }

        inline void Segment::clear_edge()
        {
            mEdge = nullptr;
        }

        inline int Segment::edge_direction() const
        {
            BELFEM_ASSERT( mEdge != nullptr, "Edge hast not been assigned to segment %lu", ( long unsigned int ) mElement->id() );
            return mEdge->node( 0 )->id() == mElement->node( 0 )->id() ? 1 : -1;
        }

        inline size_t Segment::memory() const
        {
            return sizeof( Segment )
                + this->number_of_vertices() * sizeof( Vertex * );
        }
//------------------------------------------------------------------------------
    }
}

#endif //CL_SEGMENT_HPP
