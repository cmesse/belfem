//
// Created by Christian Messe on 2018-12-25.
//

#ifndef BELFEM_FN_UNIQUE_HPP
#define BELFEM_FN_UNIQUE_HPP

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
    unique( Vector< T > & aVector )
    {
#ifdef BELFEM_ARMADILLO
        // call armadillo interface
        aVector.vector_data() = unique( aVector.vector_data() );
#elif  BELFEM_BLAZE

        // get length of vector
        std::size_t tLength = aVector.length();

        // get pointer to raw data
        T * tData = aVector.data();

        // sort data
        std::sort( tData, tData + tLength );

        // make vector unique and resize
        aVector.vector_data().resize(
                std::distance( tData,
                    std::unique( tData, tData + tLength ) ),
                true );
#endif
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_UNIQUE_HPP
