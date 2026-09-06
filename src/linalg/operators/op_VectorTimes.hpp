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

#ifndef BELFEM_OP_VECTORTIMES_HPP
#define BELFEM_OP_VECTORTIMES_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename A, typename B >
    inline auto
    operator*( const Vector< A > & aA,
               const           B  & aB )
        -> decltype( aA.vector_data() * aB )
    {

        return aA.vector_data() * aB ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename A, typename B >
    inline auto
    operator*( const              A & aA,
               const Vector< B > & aB )
        -> decltype( aA * aB.vector_data() )
    {

        return aA * aB.vector_data() ;
    }

//------------------------------------------------------------------------------

}
#endif //BELFEM_OP_VECTORTIMES_HPP
