//
// Created by Christian Messe on 01.09.19.
//

#include "cl_Arguments.hpp"
#include "commtools.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Arguments::Arguments( int & argc, char * argv[] )
    {
        mArguments.clear();

        // if ( comm_rank() == 0 )
        {
            for ( int k = 0; k < argc; ++k )
            {
                mArguments.push( string( argv[ k ] ));
            }

            // string for communicator
            string tArgString = "";
            for ( int k = 1; k < argc; ++k )
            {
                tArgString += string(  argv[ k ] );
                if( k < argc-1 )
                {
                    tArgString += " ";
                }
            }

            // write string into communicator
            gComm.set_arguments( tArgString );
        }
    }

//------------------------------------------------------------------------------

    const Cell< string > &
    Arguments::data() const
    {
        return mArguments ;
    }

//------------------------------------------------------------------------------

    const string &
    Arguments::data( const index_t aIndex ) const
    {
        return mArguments( aIndex );
    }

//------------------------------------------------------------------------------
}