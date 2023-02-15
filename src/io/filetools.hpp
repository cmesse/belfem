//
// Created by Christian Messe on 2018-12-25.
//

#ifndef BELFEM_FILE_TOOLS_HPP
#define BELFEM_FILE_TOOLS_HPP

#include <cstdio>
#include <fstream>
#include <string>

#include "stringtools.hpp"

//------------------------------------------------------------------------------
namespace belfem
{
//------------------------------------------------------------------------------

    enum class FileMode
    {
        NEW,
        OPEN_RDONLY,
        OPEN_RDONLY_PARALLEL,
        OPEN_RDWR
    };

//------------------------------------------------------------------------------

    /**
     * this function tests if a file exists
     * @param aPath
     * @return
     */
    bool
    file_exists( const std::string & aPath );

//------------------------------------------------------------------------------

    /**
     * this function takes a path and makes it parrallel
     */
    std::string
    make_path_parallel( const std::string & aPath );

//------------------------------------------------------------------------------
} /* namespace belfem */
//------------------------------------------------------------------------------
#endif //BELFEM_FILE_TOOLS_HPP
