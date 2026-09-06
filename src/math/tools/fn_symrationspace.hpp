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

#ifndef BELFEM_FN_SYMRATIONSPACE_HPP
#define BELFEM_FN_SYMRATIONSPACE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "assert.hpp"

namespace belfem
{
    template < typename  T >
    void
    symratiospace( const T & aXmin, const T & aXmax, const T & aRatio, const index_t aN, Vector< T > & aX )
    {
        // check input
        BELFEM_ASSERT( aN % 2 == 1 , "aN must be odd" );

        // set size of vector
        aX.set_size( aN );

        // number of steps per height
        index_t tN = 0.5 * ( aN - 1 ) ;

        // stepper
        T tDeltaX = 1.0 ;
        T tX  = 0.0 ;
        for( index_t k=0; k<tN; ++k )
        {
            tX = tX + tDeltaX ;
            tDeltaX *= aRatio ;
        }

        // populate beginning, end and middle point
        aX( 0 )      = aXmin ;
        aX( aN - 1 ) = aXmax ;
        aX( tN ) = aXmin + 0.5 * ( aXmax - aXmin );

        // adapt scale
        tDeltaX = aX( tN ) / tX ;

        // populate the rest
        for( index_t k=1; k<tN; ++k )
        {
            aX( k ) = aX( k - 1 ) + tDeltaX ;
            aX( aN - k - 1 ) = aX( aN - k ) - tDeltaX ;
            tDeltaX *= aRatio ;
        }
    }

}
#endif //BELFEM_FN_SYMRATIONSPACE_HPP
