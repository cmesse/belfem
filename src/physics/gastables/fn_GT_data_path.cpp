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

#include "filetools.hpp"
#include "globals.hpp"
#include "fn_GT_data_path.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        //! subdirectory of the data path that holds the gas tables
        constexpr char gGastablesSubdir[] = "/fluid";

        //! a directory counts as the gas table directory if this file is in it
        constexpr char gGastablesMarker[] = "/gasdata.inp";

//------------------------------------------------------------------------------

        /**
         * return aPath if it holds the gas tables, otherwise an empty string
         */
        static string
        check_path( const string & aPath )
        {
            if( aPath.size() == 0 )
            {
                return "";
            }

            return file_exists( aPath + gGastablesMarker ) ? aPath : "";
        }

//------------------------------------------------------------------------------

        string
        data_path()
        {
            string tPath;

            // gBelfemDataPath points at the share directory. It is set from
            // $BELFEM_DATA by Communicator::set_globals(), and a code may also set
            // it itself, for instance from a config file. It is returned without
            // checking, so that a wrong path produces a clear error when the file
            // is opened rather than silently falling through to a different one.
            if( gBelfemDataPath.size() > 0 )
            {
                return gBelfemDataPath + gGastablesSubdir;
            }

            // relative to the working directory, so that running from a build
            // directory finds the tables without any configuration
            for( const char * tRelative : { "share/fluid",
                                            "../share/fluid",
                                            "../../share/fluid",
                                            "../../../share/fluid" } )
            {
                tPath = check_path( tRelative );
                if( tPath.size() > 0 ) return tPath;
            }

            return "";
        }

//------------------------------------------------------------------------------

        bool
        data_available()
        {
            return check_path( data_path() ).size() > 0;
        }

//------------------------------------------------------------------------------
    }
}
