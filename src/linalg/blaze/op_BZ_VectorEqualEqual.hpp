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

#ifndef BELFEM_OP_BZ_VECTOREQUALEQUAL_HPP
#define BELFEM_OP_BZ_VECTOREQUALEQUAL_HPP

#include "cl_BZ_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    template < typename T >
    bool
    operator==( const Vector< T > & aA,
                const Vector< T > & aB )
    {
        return aA.vector_data() == aB.vector_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename T >
    bool
    operator==( const T           & aA,
                const Vector< T > & aB )
    {
        return aA == aB.vector_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename T >
    bool
    operator==( const Vector< T > & aA,
                const T           & aB )
    {
        return aA.vector_data() == aB;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_OP_BZ_VECTOREQUALEQUAL_HPP
