//
// Created by Christian Messe on 04.09.19.
//

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
