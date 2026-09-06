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
 * @brief Appends one vector onto another, in place.
 * @ingroup grp_linalg
 *
 */

#ifndef FN_APPEND_HPP
#define FN_APPEND_HPP

#ifdef BELFEM_BLAZE
#include <algorithm>   // for std::copy
#endif
#include "cl_Vector.hpp"

namespace belfem
{
    /**
     * @brief Appends one vector onto the end of another, in place.
     * @ingroup grp_linalg
     *
     * @warning @p aA and @p aB must be different objects. In assertion-enabled builds the
     * Blaze path checks this with BELFEM_ASSERT; the Armadillo path has no such check,
     * and in a release build neither does. Self-appending would silently produce a
     * doubled vector. Do not rely on either backend's aliasing behaviour.
     *
     * The overload taking @p aB by const reference carries exactly the same hazard.
     *
     * @param aA grown by the length of @p aB, with @p aB's entries copied onto its end
     * @param aB the vector to append; not modified, despite the non-const reference
     */
    template< typename T >
    void
    append( Vector< T > & aA,  Vector< T > & aB )
{
#ifdef BELFEM_ARMADILLO
    aA.vector_data() = arma::join_cols( aA.vector_data(), aB.vector_data() );

#elif BELFEM_BLAZE
    // grow in place ( resize with preserve keeps the existing entries ),
    // then copy aB into the new tail; DynamicVector has no insert()
    BELFEM_ASSERT( &aA != &aB, "append: aA and aB must not be the same vector" );

    const size_t tOldLength = aA.length();

    aA.vector_data().resize( tOldLength + aB.length(), true );

    std::copy( aB.data(), aB.data() + aB.length(), aA.data() + tOldLength );

#else
    // Store original data in temporary container
    Vector< T > tTempA = aA;
    aA.set_size( tTempA.length() + aB.length() );

    index_t tCount = 0;
    for ( T tVal : tTempA )
    {
        aA( tCount++ ) = tVal;
    }
    for ( T tVal : aB )
    {
        aA( tCount++ ) = tVal;
    }
#endif
}

template< typename T >
void
append( Vector< T > & aA,  const Vector< T > & aB )
{
#ifdef BELFEM_ARMADILLO
    aA.vector_data() = arma::join_cols( aA.vector_data(), aB.vector_data() );

#elif BELFEM_BLAZE
    // grow in place ( resize with preserve keeps the existing entries ),
    // then copy aB into the new tail; DynamicVector has no insert()
    BELFEM_ASSERT( &aA != &aB, "append: aA and aB must not be the same vector" );

    const size_t tOldLength = aA.length();

    aA.vector_data().resize( tOldLength + aB.length(), true );

    std::copy( aB.data(), aB.data() + aB.length(), aA.data() + tOldLength );

#else
    // Store original data in temporary container
    Vector< T > tTempA = aA;
    aA.set_size( tTempA.length() + aB.length() );

    index_t tCount = 0;
    for ( T tVal : tTempA )
    {
        aA( tCount++ ) = tVal;
    }
    for ( T tVal : aB )
    {
        aA( tCount++ ) = tVal;
    }
#endif
}

} // namespace belfem

#endif //FN_APPEND_HPP
