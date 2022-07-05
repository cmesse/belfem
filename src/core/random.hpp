//
// Created by Christian Messe on 27.10.19.
//

#ifndef BELFEM_RANDOM_HPP
#define BELFEM_RANDOM_HPP

#include <fstream>
#include "typedefs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template < typename T >
    void
    random_seed( T & aSeed )
    {
        std::ifstream tStream ( "/dev/urandom", std::ios::binary );

        // test if stream exists
        if( tStream )
        {
            // create memory block
            char tMemblock[ sizeof( T ) ];

            // read random data from stream
            tStream.read( tMemblock, sizeof( T ) );

            // close stream
            tStream.close();

            // cast stream to data
            aSeed = *reinterpret_cast< T * >( tMemblock );
        }
        else
        {
            // use the clock
            aSeed = ( T ) time( NULL ) ;
        }
    }

//------------------------------------------------------------------------------

    /**
     * seed the C++ random generator
     */
    inline
    void
    random_seed()
    {
        unsigned int tSeed;
        random_seed( tSeed );
        std::srand( tSeed );
    }

//------------------------------------------------------------------------------

    /**
     * a random number between 0 and 1
     * must call seed first
     */
    inline real
    rand()
    {
        return ( ( real ) std::rand() )/ RAND_MAX ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_RANDOM_HPP
