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

#ifndef BELFEM_CL_EDGE_HPP
#define BELFEM_CL_EDGE_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"
#include "cl_Vertex.hpp"

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------

        /**
         * \brief Special Edge class for NEDELEC-Type elements
         */
        class Edge : public Vertex
        {
            Edge * mPeriodic = nullptr;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            Edge();

//-----------------------------------------------------------------------------

            ~Edge() override;

//-----------------------------------------------------------------------------

            EntityType
            entity_type() const override ;

//-----------------------------------------------------------------------------
            size_t
            memory() const ;

//------------------------------------------------------------------------------

            void
            set_periodic( Edge * aEdge ) ;

            Edge *
            periodic() ;

            const Edge *
            periodic() const ;

            bool
            is_periodic() const ;

//-----------------------------------------------------------------------------
        };
//-----------------------------------------------------------------------------

        inline EntityType
        Edge::entity_type() const
        {
            return EntityType::EDGE ;
        }

 //-----------------------------------------------------------------------------

        inline size_t
        Edge::memory() const
        {
            return sizeof( Edge ) + this->array_memory() ;
        }

//-----------------------------------------------------------------------------

        inline void
        Edge::set_periodic( Edge * aEdge )
        {
            mPeriodic = aEdge ;
        }

        inline Edge *
        Edge::periodic()
        {
            return mPeriodic ;
        }

        inline const Edge *
        Edge::periodic() const
        {
            return mPeriodic ;
        }

        inline bool
        Edge::is_periodic() const
        {
            return mPeriodic != nullptr ;
        }

 //-----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EDGE_HPP
