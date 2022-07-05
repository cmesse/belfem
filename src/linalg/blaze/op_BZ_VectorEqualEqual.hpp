//
// Created by Christian Messe on 03.09.19.
//

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
