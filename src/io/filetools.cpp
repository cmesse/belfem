//
// Created by Christian Messe on 2018-12-25.
//

#include "cl_Communicator.hpp"
#include "filetools.hpp"

// Externally Defined Global Communicator
extern belfem::Communicator gComm;

namespace belfem
{
//------------------------------------------------------------------------------

    bool
    file_exists( const std::string & aPath )
    {

        // test if file exists
        std::ifstream tFile( aPath );

        // save result into output variable
        bool aFileExists;

        if( tFile )
        {
            // close file
            tFile.close();
            aFileExists = true;
        }
        else
        {
            aFileExists = false;
        }

        return aFileExists;
    }

//------------------------------------------------------------------------------

    /**
     * this function takes a path and makes it parallel
     */
    std::string
    make_path_parallel( const std::string & aPath )
    {

        // test if running in parallel mode
        if ( gComm.size() > 1 )
        {
            // get file extesion
            std::string tFileExt = filetype( aPath );

            // get base path
            std::string tBasePath = aPath.substr( 0, aPath.find_last_of(".") );

            // add proc number to path
            std::string aParallelPath = tBasePath + "_"
                                        +  std::to_string( gComm.size() ) + "."
                                        +  std::to_string( gComm.rank() ) + "."
                                        + tFileExt;
            return aParallelPath;
        }
        else
        {
            // do not modify path
            return aPath;
        }
    }

//------------------------------------------------------------------------------
}