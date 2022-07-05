//
// Created by Christian Messe on 03.09.19.
//

#ifndef BELFEM_OP_VECTORDIV_HPP
#define BELFEM_OP_VECTORDIV_HPP

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    inline auto
    operator/( const Vector< T > & aA,
               const             T  & aB )
        -> decltype( aA.vector_data() / aB )
    {

        return aA.vector_data() / aB ;
    }

//------------------------------------------------------------------------------

}
#endif //BELFEM_OP_VECTORDIV_HPP
