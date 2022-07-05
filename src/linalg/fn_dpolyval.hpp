//
// Created by Christian Messe on 02.07.20.
//

#ifndef BELFEM_FN_DPOLYVAL_HPP
#define BELFEM_FN_DPOLYVAL_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    /**
     * computes the derivative of a polynomial.
     */
    template < typename T >
    T
    dpolyval( const Vector< T > & aCoeffs, const T & aX )
    {
        const index_t tN = aCoeffs.length() - 1 ;
        T tPow = ( T ) tN ;
        T aResult = tPow * aCoeffs( 0 );

        for( index_t k=1; k<tN; ++k )
        {
            aResult *= aX;
            tPow -= 1.0 ;
            aResult += tPow * aCoeffs( k );
        }

        return aResult;
    }
}
#endif //BELFEM_FN_DPOLYVAL_HPP
