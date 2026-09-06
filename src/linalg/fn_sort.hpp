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
 * @brief Sorts a vector in place, ascending.
 * @ingroup grp_linalg
 *
 */

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

    /**
     * @brief Sorts a vector in place, in ascending order.
     * @ingroup grp_linalg
     * @param aVector sorted in place; its length does not change
     */
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
