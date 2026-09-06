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

#include <cstring>

#include "cl_StringList.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    StringList::StringList( const uint aNumberOfStrings ) :
        mMemory( aNumberOfStrings ),
        mCount( 0 )
    {
        // allocate data container
        mData = ( char** ) malloc( sizeof( *mData ) * ( aNumberOfStrings + 1 ) );

        BELFEM_ERROR( mData != nullptr,
            "Failed to allocate string list for %u entries",
            ( unsigned int ) aNumberOfStrings );
    }

//------------------------------------------------------------------------------

    StringList::~StringList()
    {
        for( uint k=0; k<mCount; ++k )
        {
            free( mData[ k ] );
        }

        free( mData );
    }

//------------------------------------------------------------------------------

    void
    StringList::push( const string & aString )
    {
        // always-active: a release-mode overflow would write past the container
        BELFEM_ERROR( mCount < mMemory,
                "Error adding '%s' to stringlist: memory full.",
                aString.c_str() );

        mData[ mCount ] = ( char * ) malloc( sizeof( char ) * ( aString.length() + 1 ) );

        BELFEM_ERROR( mData[ mCount ] != nullptr,
                "Failed to allocate string list entry for '%s'",
                aString.c_str() );

        std::strcpy( mData[ mCount++ ], aString.c_str() );
    }

//------------------------------------------------------------------------------

}