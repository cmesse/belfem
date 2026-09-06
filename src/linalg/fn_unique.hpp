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

/**
 * @file
 * @brief Sorts a vector in place and removes duplicate entries.
 * @ingroup grp_linalg
 *
 */

#ifndef BELFEM_FN_UNIQUE_HPP
#define BELFEM_FN_UNIQUE_HPP

#ifdef BELFEM_BLAZE
#include <memory>
#include <vector>
#include <algorithm>
#endif

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Sorts a vector in place and drops duplicate entries.
     * @ingroup grp_linalg
     *
     * The vector is **sorted as a side effect**, on both backends -- this is not a
     * duplicate filter that preserves the original order. The vector shrinks to the
     * number of distinct values.
     *
     * @param aVector sorted, deduplicated and resized in place
     */
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
