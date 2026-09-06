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

#include "cl_Segment.hpp"
#include "assert.hpp"
namespace belfem
{
    namespace mesh
    {
        Segment::Segment( Element *aElement ) :
            mElement( aElement )
        {

        }

//------------------------------------------------------------------------------

        Segment::~Segment()
        {
            delete mElement;
        }

//------------------------------------------------------------------------------

        inline Element *
        Segment::element( const uint aIndex )
        {
            BELFEM_ERROR( false, "forbidden call to element( const uint aIndex ) for mesh::Facet");
            return nullptr ;
        }

//------------------------------------------------------------------------------

        inline const Element *
        Segment::element( const uint aIndex ) const
        {
            BELFEM_ERROR( false, "invalid call to element( const uint aIndex ) const for mesh::Facet");
            return nullptr ;
        }

//------------------------------------------------------------------------------
    }
}