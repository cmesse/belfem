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

#include <fstream>
#include <sstream>

#include "fn_embed_python_guide.hpp"
#include "assert.hpp"
#include "cl_HDF5.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "filetools.hpp"
#include "globals.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        //! location of the reference reader below the data directory
        constexpr char gPythonGuideSubpath[] = "/python/database/howtoread.py";

//------------------------------------------------------------------------------

        bool
        embed_python_guide( HDF5 & aFile )
        {
            BELFEM_ASSERT( comm_rank() == 0,
                "embed_python_guide() must only run on rank 0 -- the save "
                "routines that call it guard the rank before opening the file" );

            // without a data path there is nothing to look for
            if ( gBelfemDataPath.size() == 0 )
            {
                return false ;
            }

            const string tPath = gBelfemDataPath + gPythonGuideSubpath ;

            if ( ! file_exists( tPath ) )
            {
                return false ;
            }

            // one-off setup I/O: read the whole guide into a string
            std::ifstream tStream( tPath );
            if ( ! tStream.is_open() )
            {
                return false ;
            }
            std::stringstream tText ;
            tText << tStream.rdbuf() ;

            // a guide that EXISTS but cannot be read completely must not be
            // embedded half-way: a truncated reference reader is worse than
            // none ( audit finding, 2026-08-28 ). Absent stays silent by
            // design; broken gets a visible warning.
            if ( tStream.bad() )
            {
                tStream.close() ;
                message( InfoLevel::Default,
                    "    Warning: could not read %s completely -- the python "
                    "guide is not embedded.", tPath.c_str() );
                return false ;
            }
            tStream.close() ;

            // the HDF5 string writer rejects empty strings; a directory that
            // passes file_exists() also lands here with an empty read
            if ( tText.str().empty() )
            {
                return false ;
            }

            aFile.create_group( "python" );
            aFile.save_data( "howtoread.py", tText.str() );
            aFile.close_active_group() ;

            return true ;
        }

//------------------------------------------------------------------------------
    }
}
