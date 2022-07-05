//
// Created by Christian Messe on 2019-07-30.
//
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
        BELFEM_ASSERT( mCount < mMemory,
                "Error adding '%s' to stringlist: memory full.",
                aString.c_str() );

        mData[ mCount ] = ( char * ) malloc( sizeof( char ) * ( aString.length() + 1 ) );
        std::strcpy( mData[ mCount++ ], aString.c_str() );
    }

//------------------------------------------------------------------------------

}