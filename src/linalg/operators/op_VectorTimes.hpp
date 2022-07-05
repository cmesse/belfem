//
// Created by Christian Messe on 03.09.19.
//

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
