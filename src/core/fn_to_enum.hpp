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


#ifndef BELFEM_FN_TO_ENUM_HPP
#define BELFEM_FN_TO_ENUM_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "stringtools.hpp"

namespace belfem
{

//------------------------------------------------------------------------------

    template < typename T >
    void
    to_enum( const string & aString, T & aEnum )
    {
        int tNumEntries = ( int ) T::UNDEFINED ;

        aEnum = T::UNDEFINED ;

        string tString = string_to_lower( aString );

        for( int k=0; k<tNumEntries; ++k )
        {
            aEnum = ( T ) k ;
            string tEnum = string_to_lower( to_string( aEnum ) );

            if( tString == tEnum )
            {
                return ;
            }
        }

        // overwrite enum with false value if no value was detected
        aEnum = T::UNDEFINED ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_TO_ENUM_HPP
