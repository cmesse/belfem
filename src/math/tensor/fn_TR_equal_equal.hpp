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

#ifndef BELFEM_FN_TR_EQUAL_EQUAL_HPP
#define BELFEM_FN_TR_EQUAL_EQUAL_HPP
#include "typedefs.hpp"
namespace belfem
{
    namespace tensor
    {
//----------------------------------------------------------------------------

        template < typename T >
        inline bool
        equal_equal( const T * A, const T * B, const index_t aCapacity )
        {
            for( index_t k=0; k<aCapacity; ++k )
            {
                if( A[ k ] != B[ k ] )
                {
                    return false ;
                }
            }

            return true ;
        }

//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */
#endif //BELFEM_FN_TR_EQUAL_EQUAL_HPP
