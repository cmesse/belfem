//
// Created by Christian Messe on 15.09.19.
//

#ifndef BELFEM_FN_POLYVAL_HPP
#define BELFEM_FN_POLYVAL_HPP

#ifdef BELFEM_ARMADILLO
#include "armadillo.hpp"
#endif

#include "typedefs.hpp"
#include "assert.h"
#include "cl_Vector.hpp"

namespace belfem
{

//------------------------------------------------------------------------------

    template < typename T >
    T
    polyval( const Vector< T > & aCoeffs, const T & aX )
    {
        const index_t tN = aCoeffs.length();
        T aResult = aCoeffs( 0 );

        for( index_t k=1; k<tN; ++k )
        {
            aResult *= aX;
            aResult += aCoeffs( k );
        }

        return aResult;
    }

//------------------------------------------------------------------------------


    template < typename T >
    void
    polyval( const Vector< T > & aCoeffs, const Vector< T > & aX, Vector< T > & aY )
    {
#ifdef BELFEM_ARMADILLO
        aY.vector_data() = arma::polyval( aCoeffs.vector_data(), aX.vector_data() );
#else
        index_t tN = aX.vector_data();
        aY.set_size( tN );
        for( uint k=0; k<tN; ++k )
        {
            aY( k ) = polyval( aCoeffs, aX( k ) );
        }
#endif
    }



//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_POLYVAL_HPP
