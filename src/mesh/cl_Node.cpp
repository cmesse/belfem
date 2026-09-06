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

#include "cl_Node.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Node::Node(
                const id_t & aID,
                const real aX,
                const real aY,
                const real aZ ) :
                Vertex()
        {
            // set the id
            this->set_id( aID );

            // copy node coordinates into coordinate vector
            mCoords[ 0 ] = aX;
            mCoords[ 1 ] = aY;
            mCoords[ 2 ] = aZ;
        }

//------------------------------------------------------------------------------

        Node::~Node()
        {
            this->reset_duplicate_container();
            this->delete_containers();
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const Vector< real > & aCoords )
        {
            std::copy( aCoords.data(),
                       aCoords.data()
                        + aCoords.length(), mCoords );
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const real aX, const real aY )
        {
            mCoords[ 0 ] = aX;
            mCoords[ 1 ] = aY;
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const real aX, const real aY, const real & aZ )
        {
            mCoords[ 0 ] = aX;
            mCoords[ 1 ] = aY;
            mCoords[ 2 ] = aZ;
        }

//------------------------------------------------------------------------------
    }
}