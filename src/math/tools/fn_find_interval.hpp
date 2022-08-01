//
// Created by Christian Messe on 02.12.19.
//

#ifndef BELFEM_FN_FIND_INTERVAL_HPP
#define BELFEM_FN_FIND_INTERVAL_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    inline void
    find_interval( const Vector< real > & aData, const real aValue, index_t & aIndex, real & aXi )
    {
        index_t i = 0;
        index_t k = aData.length() - 1;
        index_t j = ( i + k ) / 2 ;

        if( aValue <= aData( 1 ) )
        {
            aIndex = 0;

            aXi = ( aValue - aData( 0 ) ) / ( aData( 1 ) - aData( 0 ) );
        }
        else if ( aValue >= aData( k - 1 ) )
        {
            aIndex = k - 1;
            aXi = ( aValue - aData( k-1 ) ) / ( aData( k ) - aData( k - 1 ) ) ;
        }
        else
        {
            while ( true )
            {
                if (( aData( i ) <= aValue ) && ( aData( j ) >= aValue ))
                {
                    // this is the correct interval
                    k = j;
                }
                else
                {
                    // the other one is the correct interval
                    i = j;
                }
                if (( k - i ) == 1 )
                {
                    break;
                }
                else
                {
                    j = ( i + k ) / 2 + ( i + k ) % 2;
                }
                aIndex = j;

                aXi = ( aValue - aData( j )) / ( aData( j + 1 ) - aData( j ));

            }
        }
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_FIND_INTERVAL_HPP
