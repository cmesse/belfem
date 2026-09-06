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

#ifndef BELFEM_OP_BZ_MATRIXEQUALEQUAL_HPP
#define BELFEM_OP_BZ_MATRIXEQUALEQUAL_HPP

#include "cl_BZ_Matrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template < typename T >
    bool
    operator==( const Matrix< T > & aA,
                const Matrix< T > & aB )
    {
        return aA.matrix_data() == aB.matrix_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename T >
    bool
    operator==( const T           & aA,
                const Matrix< T > & aB )
    {
        return aA == aB.matrix_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template < typename T >
    bool
    operator==( const Matrix< T > & aA,
                const T           & aB )
    {
        return aA.matrix_data() == aB;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_OP_BZ_MATRIXEQUALEQUAL_HPP
