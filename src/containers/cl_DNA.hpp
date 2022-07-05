//
// Created by Christian Messe on 02.05.21.
//

#ifndef BELFEM_CL_DNA_HPP
#define BELFEM_CL_DNA_HPP

#include <bitset>
#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{

    template < index_t N >
    class DNA
    {
//------------------------------------------------------------------------------

        // the wrapped standard bitset which represents the DNA
        std::bitset< N * sizeof( float ) * CHAR_BIT > mBitset;

        // flag telling if this gene is alive
        bool mAlive = true;

        // fitness value
        real mFitness = BELFEM_REAL_MAX;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        DNA() = default;

//------------------------------------------------------------------------------

        /**
         * tivial destructor
         */
        ~DNA() = default;

//------------------------------------------------------------------------------

        /**
         * sets the values of the bitset
         */
        void
        set_values( const Vector <real> & aValues )
        {
            // make sure that memory sizes match
            BELFEM_ASSERT( sizeof( float ) == sizeof( uint ),
                          "Wrong sizes of float and unsigned int ( %u vs %u )",
                          ( unsigned int ) sizeof( float ),
                          ( unsigned int ) sizeof( uint ));

            // make sure that vector has the right length
            BELFEM_ASSERT( aValues.length() == N, "Vector has length %u, but expect %u",
                          ( unsigned int ) aValues.length(),
                          ( unsigned int ) N );

            // set alive flag
            mAlive = true;

            // convert values to float
            float tValues[N];

            for ( index_t k = 0; k < N; ++k )
            {
                tValues[ k ] = float( aValues( k ));
            }


            // determine size of temporary bitset
            const index_t tSize = sizeof( float ) * CHAR_BIT;

            // offset
            index_t tCount = 0;

            // temporary value
            uint tSwap;

            for ( index_t k = 0; k < 3; ++k )
            {
                // convert value to int
                memcpy( &tSwap, &tValues[ k ], sizeof( float ));

                // create a temporary bitset
                std::bitset< tSize > tBitset( tSwap );

                for ( index_t i = 0; i < tSize; ++i )
                {
                    if ( tBitset.test( i ))
                    {
                        mBitset.set( tCount++ );
                    }
                    else
                    {
                        mBitset.reset( tCount++ );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * return one specific value of the bitset
         * @param aIndex
         * @return
         */
        real
        get_value( const index_t & aIndex ) const
        {
            // make sure that vector has the right length
            BELFEM_ASSERT( aIndex < N, " requested index out of bounds ( is %u, but must be < %u )",
                          ( unsigned int ) aIndex,
                          ( unsigned int ) N );

            // create a temporary bitset
            std::bitset< 32 > tBitset;

            // counter for dna
            uint tCount = aIndex * sizeof( float ) * CHAR_BIT;

            for ( index_t k = 0; k < 32; ++k )
            {
                if ( mBitset.test( tCount++ ))
                {
                    tBitset.set( k );
                }
                else
                {
                    tBitset.reset( k );
                }
            }

            // convert bitset to unsigned int
            uint tSwap = uint( tBitset.to_ulong());
            float tValue;
            memcpy( &tValue, &tSwap, sizeof( float ));

            // return value
            return real( tValue );
        }
//------------------------------------------------------------------------------

        /**
         * the length of the genome
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
         * deactivate the alive flag for this bit and kill fitness
         */
        inline void
        kill()
        {
            mAlive = false ;
            mFitness = BELFEM_REAL_MAX ;
        }
//------------------------------------------------------------------------------

        /**
         * activate alive flag
         */
        inline void
        resurrect()
        {
            mAlive = true;
            mFitness = 0.0 ;
        }

//------------------------------------------------------------------------------

        /**
         * tell if this dna is alive
         */
        inline bool
        alive() const
        {
            return mAlive;
        }

//------------------------------------------------------------------------------

        /**
         * set the fitness value
         */
        inline void
        set_fitness( const real & aFitness )
        {
            mFitness = aFitness;
        }

//------------------------------------------------------------------------------

        /**
         * get the fitness value
         */
        const real &
        fitness() const
        {
            return mFitness;
        }

//------------------------------------------------------------------------------

        /**
         * expose data container
         */
        inline std::bitset< N * sizeof( float ) * CHAR_BIT > &
        data()
        {
            return mBitset;
        }

    };

//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_DNA_HPP
