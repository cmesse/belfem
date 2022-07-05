//
// Created by Christian Messe on 04.09.19.
//

#ifndef BELFEM_FN_BZ_LINSPACE_HPP
#define BELFEM_FN_BZ_LINSPACE_HPP

#include "cl_BZ_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    template< typename T >
    void
    linspace( const T & aStart,
              const T & aEnd,
              const belfem::size_t & aN,
              Vector <T> & aValues )
    {
        // yes, this can be done a lot better using std::for_each

        aValues.set_size( aN );

        belfem::size_t tN = aN - 1;

        // stepsize
        T tDX = ( aEnd - aStart ) / ( aN - 1 );

        // first step
        aValues( 0 ) = aStart;

        // intermediate steps
        for ( uint k = 1; k < tN; ++k )
        {
            aValues( k ) = aValues( k - 1 ) + tDX;
        }

        // last step
        aValues( tN ) = aEnd;
    }

//------------------------------------------------------------------------------

    template< typename T >
    Vector <T>
    linspace( const T & aStart,
              const T & aEnd,
              const belfem::size_t & aN )
    {
        Vector <real> aValues;
        linspace( aStart, aEnd, aN, aValues );
        return aValues;
    }

//------------------------------------------------------------------------------

}
#endif //BELFEM_FN_BZ_LINSPACE_HPP
