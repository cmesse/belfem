/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_HUST_HPP
#define BELFEM_FN_HUST_HPP

#include <cmath>
#include "cl_Vector.hpp"

namespace belfem
{
    namespace material
    {
        /**
         *
         * This function is a semi-empirical model used to describe the intrinsic
         * thermal resistivity or electrical resistivity of metals,
         * particularly at cryogenic temperatures.
         *
         * The function is designed to bridge the gap between low-temperature power-law behavior
         * (where resistivity often scales as $T^n$) and high-temperature limits.
         *
         * Hust, J. G., & Lankford, A. B. (1984)
         *
         * Update of thermal conductivity and electrical resistivity of electrolytic iron,
         * tungsten, and stainless steel (NBS Special Publication 260-90).
         * National Bureau of Standards.
         *
         * @param P
         * @param T
         * @return
         */
        template< typename value > value
        hust( const Vector< value > & P, const value T )
        {
            return P(1) * std::pow( T, P(2)) /
                 ( 1.0 +
                     P( 1 ) * P(3 )
                     * std::pow( T, P(2 ) + P(4) )
                     * std::exp( - std::pow( P( 5 ) / T, P( 6 ) ) ) );
        }

    }
}
#endif //BELFEM_FN_HUST_HPP