//
// Created by Christian Messe on 01.09.19.
//

#include <cctype>

#include "cl_Arguments.hpp"
#include "cl_Logger.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace arguments
    {
        // returns true if the string is a plain unsigned integer
        static bool
        is_uint( const string & aString )
        {
            if ( aString.empty() )
            {
                return false;
            }
            for ( const char tChar : aString )
            {
                if ( ! std::isdigit( static_cast< unsigned char >( tChar ) ) )
                {
                    return false;
                }
            }
            return true;
        }
    }

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

            /*
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

            // write string into communicator (obsolete)
            gComm.set_arguments( tArgString ); */
        }

        this->set_verbosity_from_arguments();
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

    void
    Arguments::set_verbosity_from_arguments()
    {
        const index_t tNumArgs = mArguments.size();

        for ( index_t k = 1; k < tNumArgs; ++k )
        {
            const string & tArg = mArguments( k );

            if ( tArg == "-v" || tArg == "--verbose" )
            {
                // an integer in the next argument selects the level,
                // a bare flag means maximum verbosity
                if ( k + 1 < tNumArgs && arguments::is_uint( mArguments( k + 1 ) ) )
                {
                    gLog.set_info_level( std::stoi( mArguments( k + 1 ) ) );
                    ++k;
                }
                else
                {
                    gLog.set_info_level( InfoLevel::Everything );
                }
            }
            else if ( tArg.rfind( "--verbose=", 0 ) == 0 )
            {
                const string tValue = tArg.substr( 10 );
                BELFEM_ERROR( arguments::is_uint( tValue ),
                        "--verbose expects an unsigned integer, got '%s'",
                        tValue.c_str() );
                gLog.set_info_level( std::stoi( tValue ) );
            }
            else if ( tArg.rfind( "-v", 0 ) == 0
                      && arguments::is_uint( tArg.substr( 2 ) ) )
            {
                // attached short form, e.g. -v3
                gLog.set_info_level( std::stoi( tArg.substr( 2 ) ) );
            }
        }
    }

//------------------------------------------------------------------------------
}
