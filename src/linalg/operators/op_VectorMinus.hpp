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

#ifndef BELFEM_OP_VECTORMINUS_HPP
#define BELFEM_OP_VECTORMINUS_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    inline auto
    operator-( const Vector<T> & aA,
               const Vector<T> & aB )
        -> decltype( aA.vector_data() - aB.vector_data() )
    {
        BELFEM_ASSERT( aA.length() == aB.length(),
                      "Length of vectors does not match ( %lu and %lu ).",
                      ( long unsigned int ) aA.length(),
                      ( long unsigned int ) aB.length() );

        return aA.vector_data() - aB.vector_data();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename A, typename B >
    inline auto
    operator-( const Vector< A > & aA,
               const            B  & aB )
    -> decltype( aA.vector_data() - aB )
    {

        return aA.vector_data() - aB ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< typename A, typename B >
    inline auto
    operator-( const              A & aA,
               const Vector< B > & aB )
        -> decltype( aA - aB.vector_data() )
    {

        return aA - aB.vector_data() ;
    }

//------------------------------------------------------------------------------

}
#endif //BELFEM_OP_VECTORMINUS_HPP
