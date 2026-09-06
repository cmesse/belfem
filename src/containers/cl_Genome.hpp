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

#ifndef BELFEM_CL_DNA_HPP
#define BELFEM_CL_DNA_HPP

#include <algorithm>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Bitset.hpp"
#include "random.hpp"

namespace belfem
{
    // B: number of bits per chromosome
    // N: number of chromosomes
    template  < size_t B, size_t N >
    class Genome
    {
        const Vector< real > & mMinVals ; //<-- no need to copy every time,
        const Vector< real > & mMaxVals ; // they are the same for all genomes
        const Bitset< N >    & mTypes ; // false: linear, true : log scale
        Bitset< N * B > mDNA ; // the dna string
        real mFitness = BELFEM_REAL_MAX;

    public:

        Genome( const Vector< real > & aMinVals,
                const Vector< real > & aMaxVals,
                const Bitset< N >    & aTypes ) :
                mMinVals( aMinVals ),
                mMaxVals( aMaxVals ),
                mTypes( aTypes )
        {

        }

        void
        set_values( const Vector< real > & aValues )
        {
            mDNA.reset();
            for( size_t i=0; i<N; ++i )
            {
                real tValue = aValues( i );
                real tMin = mMinVals( i );
                real tMax = mMaxVals( i );

                if ( std::abs( tMin - tMax ) < BELFEM_EPSILON )
                {
                    continue;
                }

                // check if we store in logarithmic scaling
                if ( mTypes.test( i ) )
                {
                    tValue = std::log( tValue );
                    tMin = std::log( tMin );
                    tMax = std::log( tMax );
                }

                // clamp value
                tValue = std::max( std::min( tValue, tMax ), tMin );

                // normalize to [0, 1]
                tValue = ( tValue - tMin ) / ( tMax - tMin );

                // scale to [0, 2^B - 1] and encode
                encode( i, std::round( tValue * ( ( 1u << B ) - 1 ) ) );
            }
        }

        void
        get_values( Vector< real > & aValues ) const
        {
            for( size_t i=0; i<N; ++i )
            {
                // scale to [tMin, tMax]
                real tMin = mMinVals( i );
                real tMax = mMaxVals( i );

                if ( std::abs( tMin - tMax ) < BELFEM_EPSILON )
                {
                    aValues( i ) = tMin ;
                    continue;
                }

                // decode uint from bitset
                uint tBits = decode( i );

                // normalize to [0, 1]
                real tValue = real( tBits ) / real( ( 1u << B ) - 1 );



                // check if we stored in logarithmic scaling
                if ( mTypes.test( i ) )
                {
                    tMin = std::log( tMin );
                    tMax = std::log( tMax );

                    tValue = std::exp( tMin + tValue * ( tMax - tMin ) );

                }
                else
                {
                    tValue = tMin + tValue * ( tMax - tMin );
                }

                // log and exp do not round-trip exactly -- exp( log( 100 ) )
                // overshoots by a few ulp -- so the top of the range decodes to
                // marginally more than the caller's maximum. set_values() and
                // randomize() clamp on the way in; clamping here as well makes
                // the bounds an invariant the caller can rely on.
                aValues( i ) = std::max( std::min( tValue, mMaxVals( i ) ),
                                                   mMinVals( i ) );
            }
        }

        void
        randomize()
        {
            mDNA.reset();
            for( size_t i=0; i<N; ++i )
            {
                real tMin = mMinVals( i );
                real tMax = mMaxVals( i );

                if ( std::abs( tMin - tMax ) < BELFEM_EPSILON )
                {
                    // do nothing
                    continue;
                }
                // for log-scale params, work in log space
                if ( mTypes.test( i ) )
                {
                    tMin = std::log( tMin );
                    tMax = std::log( tMax );
                }

                // Gaussian distribution centered at midpoint
                real tMean = 0.5 * ( tMin + tMax );
                real tStdDev = ( tMax - tMin ) / 6.0;  // ±3σ covers ~99.7% of  range

                // Box-Muller transform to generate Gaussian from uniform
                real u1 = rand();
                real u2 = rand();

                // avoid log(0)
                u1 = std::max( u1, BELFEM_EPSILON );

                real tZ = std::sqrt( -2.0 * std::log( u1 ) ) * std::cos( 2.0 * constant::pi * u2 );

                real tValue = tMean + tStdDev * tZ ;

                // clamp to [tMin, tMax]
                tValue = std::max( std::min( tValue, tMax ), tMin );

                // normalize to [0, 1]
                tValue = ( tValue - tMin ) / ( tMax - tMin );

                // encode
                encode( i, std::round( tValue * ( ( 1u << B ) - 1 ) ) );
            }
        }

        void
        inherit( const Genome * aMom, const Genome * aDad )
        {
            // total number of bits in the genome
            size_t tNumBits = N * B;

            // random split point for crossover
            size_t tSplit = size_t( rand() * real( tNumBits ) );

            // random bit to mutate
            size_t tMutate = size_t( rand() * real( tNumBits ) );

            // start with a clean slate
            mDNA.reset();

            // inherit bits [0, tSplit) from mom
            for( size_t k=0; k<tSplit; ++k )
            {
                if( aMom->mDNA.test( k ) )
                {
                    mDNA.set( k );
                }
            }

            // inherit bits [tSplit, tNumBits) from dad
            for( size_t k=tSplit; k<tNumBits; ++k )
            {
                if( aDad->mDNA.test( k ) )
                {
                    mDNA.set( k );
                }
            }

            // mutate one random bit
            mDNA.flip( tMutate );
        }

        real
        fitness() const
        {
            return mFitness;
        }

        void
        set_fitness( const real aFitness )
        {
            mFitness = aFitness;
        }

        void
        kill()
        {
            mDNA.reset();
            mFitness = BELFEM_REAL_MAX;
        }

        bool
        is_alive() const
        {
            return mFitness != BELFEM_REAL_MAX;
        }

    private:

        uint
        decode( const size_t aIndex ) const
        {


            uint aValue = 0 ;
            uint tExp = 1 ;
            size_t tCount = aIndex * B ;
            for( size_t i=0; i<B; ++i )
            {
                if ( mDNA.test( tCount++ ) )
                {
                    aValue += tExp ;
                }
                tExp *= 2 ;
            }
            return aValue ;
        }

        void
        encode( const size_t aIndex, uint aValue )
        {
            size_t tCount = aIndex * B;
            for(size_t i=0; i<B; ++i)
            {
                if (aValue & 1)  // bitwise test
                {
                    mDNA.set(tCount++);
                } else
                {
                    mDNA.reset(tCount++);
                }
                aValue >>= 1; // shift right
            }
        }
    };

    template < size_t B, size_t N >
    struct opGenomeSort
    {
        bool
        operator()( const Genome< B, N > * aA, const Genome< B, N > * aB )
        {
            return aA->fitness() < aB->fitness();
        }
    };

}
#endif //BELFEM_CL_DNA_HPP