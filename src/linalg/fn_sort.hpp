//
// Created by Christian Messe on 2019-01-21.
//

#ifndef BELFEM_FN_SORT_HPP
#define BELFEM_FN_SORT_HPP

#include "cl_Vector.hpp"

#ifdef BELFEM_BLAZE
#include <memory>
#include <vector>
#include <algorithm>
#endif

namespace belfem
{
//------------------------------------------------------------------------------

    template < typename T >
    void
    sort( Vector< T > & aVector )
    {
#ifdef BELFEM_ARMADILLO
        // call armadillo interface
        aVector.vector_data() = sort( aVector.vector_data() );
#elif  BELFEM_BLAZE
        // sort data
        std::sort( aVector.data(), aVector.data() + aVector.length() );
#endif
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_SORT_HPP
