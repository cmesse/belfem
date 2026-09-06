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

#ifndef BELFEM_CL_BITSET_HPP
#define BELFEM_CL_BITSET_HPP


#include <bitset>

#include "typedefs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Compile-time fixed-size bitset.
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    template < index_t N >
    class Bitset
    {
//------------------------------------------------------------------------------

        // the wrapped standard bitset
        std::bitset< N > mBitset;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------


        /**
         * trivial constructor
         */
        Bitset() = default;

//------------------------------------------------------------------------------

        /**
         * copy constructor
         */
        Bitset( const Bitset< N > & aBitset ) :
                mBitset( aBitset.mBitset )
        {
        }

//------------------------------------------------------------------------------

        /**
         * move constructor
         */
        Bitset( Bitset< N > && aBitset ) :
                mBitset( std::move( aBitset.mBitset ) )
        {
        }

//------------------------------------------------------------------------------

        /**
         * tivial destructor
         */
        ~Bitset() = default;

//------------------------------------------------------------------------------

        /**
         * the size of the bitset
         */
        inline index_t
        size() const
        {
            return N;
        }

//------------------------------------------------------------------------------

        /**
         * set a bit to true
         */
        inline void
        set( const index_t aIndex )
        {
            mBitset.set( aIndex );
        }

//------------------------------------------------------------------------------

        /**
         * set a bit to false
         */
        inline void
        reset( const index_t aIndex )
        {
            mBitset.reset( aIndex );
        }


//------------------------------------------------------------------------------

        /**
         * set all bits to false
         */
        inline void
        reset()
        {
            mBitset.reset();
        }


//------------------------------------------------------------------------------

        /**
         * invert a bit
         */
        inline void
        flip( const index_t aIndex )
        {
            mBitset.flip( aIndex );
        }

//------------------------------------------------------------------------------

        /**
         * test if a bit is set
         */
        inline bool
        test( const index_t aIndex ) const
        {
            return mBitset.test( aIndex );
        }

//------------------------------------------------------------------------------

        /**
         * count the number of true bits
         */
        inline index_t
        count() const
        {
            return mBitset.count();
        }

//------------------------------------------------------------------------------

        /**
         * expose data container
         */
        inline std::bitset< N > &
        data()
        {
            return mBitset ;
        }

//------------------------------------------------------------------------------

        /**
         * expose data container (const)
         */
        inline const std::bitset< N > &
        data() const
        {
            return mBitset ;
        }

//------------------------------------------------------------------------------

        /**
         * copy assignment operator
         */
        Bitset< N > &
        operator=( const Bitset< N > & aBitset )
        {
            if( this != & aBitset )
            {
                mBitset = aBitset.mBitset ;
            }
            return *this ;
        }

//------------------------------------------------------------------------------

        /**
         * move assignment operator
         */
        Bitset< N > &
        operator=( Bitset< N > && aBitset )
        {
            if( this != & aBitset )
            {
                mBitset = std::move( aBitset.mBitset ) ;
            }
            return *this ;
        }

//------------------------------------------------------------------------------

        /**
         * equality operator
         */
        inline bool
        operator==( const Bitset< N > & aBitset ) const
        {
            return mBitset == aBitset.mBitset ;
        }

//------------------------------------------------------------------------------

        /**
         * inequality operator
         */
        inline bool
        operator!=( const Bitset< N > & aBitset ) const
        {
            return mBitset != aBitset.mBitset ;
        }

//------------------------------------------------------------------------------
    };
}

#endif //BELFEM_CL_BITSET_HPP
