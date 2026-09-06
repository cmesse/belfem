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


#ifndef CL_HASH_HPP
#define CL_HASH_HPP

#include <functional>
#include <cstddef>
#include "typedefs.hpp"

namespace belfem
{
    /**
     * @brief Incremental hash computation.
     *
     * @ingroup grp_core
     * @see @ref core_core_usage_guide
     */
    class Hash
    {
        std::size_t mValue = 0 ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Hash() = default ;

//------------------------------------------------------------------------------

        ~Hash() = default ;

//------------------------------------------------------------------------------

        inline void
        reset()
        {
            mValue = 0 ;
        }

//------------------------------------------------------------------------------

        inline std::size_t
        value() const
        {
            return mValue;
        }

//------------------------------------------------------------------------------

        inline void
        set_value( const std::size_t aValue )
        {
            mValue = aValue ;
        }

//------------------------------------------------------------------------------

        template< typename T >
        Hash &
        operator += ( const T aValue )
        {
            mValue +=
                // Compute and add the hash of the new value.
                std::hash<T>{}(aValue)
                // Add an entropy constant (derived from the golden ratio) for better distribution.
                + 0x9e3779b97f4a7c15ULL
                // Mix in the current hash by shifting left (emphasizes high-order bits).
                + (mValue << 6)
                // Mix in the current hash by shifting right (brings in low-order bits).
                + (mValue >> 2);

            return *this;
        }

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------
}
#endif //CL_HASH_HPP
