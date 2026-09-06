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
 * @brief Cross product of two vectors.
 * @ingroup grp_linalg
 *
 * The header selects the Armadillo or Blaze implementation at compile time. Overloads
 * templated on `ET` accept unevaluated backend vector expressions directly, without
 * materialising a temporary `Vector`.
 */

/**
 * @fn template<typename T> auto belfem::cross( const Vector<T> & aA, const Vector<T> & aB )
 * @brief Cross product of two length-3 vectors.
 * @ingroup grp_linalg
 * @param aA left input vector; must have length 3
 * @param aB right input vector; must have length 3
 * @return the right-hand-rule vector @p aA x @p aB, with magnitude equal to the area of
 *         the parallelogram they span
 */

#ifndef BELFEM_FN_CROSS_HPP
#define BELFEM_FN_CROSS_HPP

#include "assert.hpp"
#ifdef BELFEM_ARMADILLO
#include "fn_AR_cross.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_cross.hpp"
#endif


namespace belfem
{
    template < typename T >
    auto
    cross( const Vector< T > & aA, const Vector< T > & aB )
        -> decltype( cross( aA.vector_data(), aB.vector_data() ) )
    {
        BELFEM_ASSERT( aA.length() == 3 && aB.length() == 3,
            "Both vectors must have a length of 3");

        return cross( aA.vector_data(), aB.vector_data() );
    }
}
#endif //BELFEM_FN_CROSS_HPP
