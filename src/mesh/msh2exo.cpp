//
// Created by christian on 9/1/26.
//
//
// Created by christian on 8/28/26.
//

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

#include <cctype>
#include <cmath>
#include <iostream>

#include "assert.hpp"
#include "banner.hpp"
#include "cl_Communicator.hpp"
#include "cl_HDF5.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Arguments.hpp"
#include "cl_Mesh.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );

namespace
{
//------------------------------------------------------------------------------

    void
    print_help()
    {
        std::cout << "Usage: msh2exo <gmshfile>" << std::endl;
        std::cout << "  converts a gmsh file into an exodus file of the same name" << std::endl;
        std::cout << "  -h, --help                                print this help message" << std::endl;
        std::cout << "  -v, --verbose [level]                     set the info level" << std::endl;
    }

//------------------------------------------------------------------------------

    bool
    is_uint( const string & aString )
    {
        if ( aString.empty() )
        {
            return false ;
        }
        for ( const char tChar : aString )
        {
            if ( ! std::isdigit( static_cast< unsigned char >( tChar ) ) )
            {
                return false ;
            }
        }
        return true ;
    }

//------------------------------------------------------------------------------

    /**
     * the shared -v / --verbose flags are read by Arguments but not removed
     * from the list, so the file is the first argument that is neither a flag
     * nor the operand of one. Returns an empty string if there is none
     */
    string
    input_file_from_arguments( const Arguments & aArguments )
    {
        const Cell< string > & tArgs = aArguments.data() ;

        index_t tNumArgs = tArgs.size() ;

        for ( index_t k=1; k<tNumArgs; ++k )
        {
            const string & tArg = tArgs( k ) ;

            if ( tArg == "-v" || tArg == "--verbose" )
            {
                // an unsigned integer behind the bare flag is its operand
                if ( k+1 < tNumArgs && is_uint( tArgs( k+1 ) ) )
                {
                    ++k ;
                }
                continue ;
            }

            if ( ! tArg.empty() && tArg.front() == '-' )
            {
                continue ;
            }

            return tArg ;
        }

        return "" ;
    }

//------------------------------------------------------------------------------

    bool
    help_requested( const Arguments & aArguments )
    {
        for ( const string & tArg : aArguments.data() )
        {
            if ( tArg == "-h" || tArg == "--help" )
            {
                return true ;
            }
        }
        return false ;
    }

//------------------------------------------------------------------------------

    /**
     * replaces the extension of the file name, not of the path : a directory
     * in the path may carry a dot while the file name does not, and the
     * leading dot of a hidden file belongs to its name rather than marking an
     * extension
     */
    string
    exodus_path( const string & aPath )
    {
        std::size_t tSlash = aPath.find_last_of( "/" ) ;

        // first character of the file name
        std::size_t tBegin = ( tSlash == string::npos ) ? 0 : tSlash + 1 ;

        std::size_t tDot = aPath.find_last_of( "." ) ;

        // no extension : either no dot in the file name at all, or the only
        // dot is the one that makes the file hidden
        if ( tDot == string::npos || tDot <= tBegin )
        {
            return aPath + ".exo" ;
        }

        return aPath.substr( 0, tDot ) + ".exo" ;
    }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
}

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // create Arguments ( this also reads the shared -v / --verbose flags )
    Arguments tArguments( argc, argv );

    if ( gComm.rank() == 0 )
    {
        string tPath = input_file_from_arguments( tArguments );

        if ( tPath.empty() || help_requested( tArguments ) )
        {
            print_help();
            return gComm.finalize();
        }

        print_banner( "msh2exo" );

        string tName = exodus_path( tPath );

        Mesh tMesh( tPath );
        tMesh.save( tName ) ;
    }

    return gComm.finalize();
}
