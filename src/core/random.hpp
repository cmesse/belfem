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

#ifndef BELFEM_RANDOM_HPP
#define BELFEM_RANDOM_HPP

#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <fstream>
#include <random>
#include "typedefs.hpp"

#include "cl_Communicator.hpp"
#ifdef BELFEM_MPI
extern belfem::Communicator gComm;
#else
namespace belfem
{
    template < typename T >
    void
    random_seed( T & aSeed )
    {
        static_assert( std::is_integral< T >::value && std::is_unsigned< T >::value,
            "random_seed: the seed must be an unsigned integer" );

        std::ifstream tStream ( "/dev/urandom", std::ios::binary );

        // test if stream exists
        if( tStream )
        {
            char tMemblock[ sizeof( T ) ];
            tStream.read( tMemblock, sizeof( T ) );
            const bool tOk = tStream.gcount() == std::streamsize( sizeof( T ) );
            tStream.close();

            if ( tOk )
            {
                std::memcpy( &aSeed, tMemblock, sizeof( T ) );
                return ;
            }
        }

        // no random device, or a short read: use the clock
        aSeed = ( T ) time( NULL ) ;
    }
}
#endif
namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * seed the random generator. MPI builds seed the communicator's engine,
     * the serial branch seeds std::rand(), which rand() uses in that build.
     */
    inline
    void
    random_seed()
    {
#ifdef BELFEM_MPI
        std::random_device rd;
        gComm.random().seed(rd());  // Seed once using random device
#else
        unsigned int tSeed;
        random_seed( tSeed );
        std::srand( tSeed );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * a random number between 0 and 1
     * must call seed first
     */
    inline real
    rand()
    {
#ifdef BELFEM_MPI
        std::uniform_real_distribution<real> distribution(0.0, 1.0);
        return distribution( gComm.random() ); // Generates a real number between 0 and 1
#else
        return ( ( real ) std::rand() )/ RAND_MAX ;
#endif
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_RANDOM_HPP
