//
// Created by Christian Messe on 09.09.19.
//

#ifndef BELFEM_CL_BITSET_HPP
#define BELFEM_CL_BITSET_HPP


#include <bitset>

#include "typedefs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

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
    };
}

#endif //BELFEM_CL_BITSET_HPP
