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
 * @brief Converts a Cell into a Vector.
 * @ingroup grp_linalg
 *
 * Copies the entries into a new Vector. Order and length are preserved and no storage is
 * shared with the source.
 *
 */

#ifndef BELFEM_FN_TO_VECTOR_HPP
#define BELFEM_FN_TO_VECTOR_HPP

#include "cl_Vector.hpp"
#include "cl_Cell.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * convert a Cell< T > into a Vector< T >
     */
    template < typename T >
    Vector< T >
    to_vector( const Cell< T > & aCell )
    {
        // get the number of entries
        size_t tN = aCell.size();

        // allocate the output vector
        Vector< T > aResult( tN );

        // raw data pointers (valid for both Armadillo and Blaze backends)
        T * tTarget = aResult.data();
        const T * tSource = aCell.data();

        for( size_t k = 0; k < tN; ++k )
        {
            tTarget[ k ] = tSource[ k ];
        }

        return aResult;
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_TO_VECTOR_HPP
