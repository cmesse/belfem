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

    /**
     * Resolve a data file named in an input file, looking below aSubDirectory
     * of the shared data directory.
     *
     * The run directory always wins, so a local copy overrides the shared
     * database. Failing that the file is looked up below the data directory,
     * first under the same relative path and then by name alone, so that
     * "MatData/bhdata.hdf5" also resolves for a run that has no MatData
     * directory of its own.
     *
     * If nothing is found, aFile is returned UNCHANGED. That is load-bearing
     * in two ways: whoever opens the file reports the name the user wrote,
     * and the dlopen search path stays intact for plugin libraries, which
     * resolve through the loader and need not exist as a file at all.
     * Note that dlopen only consults the platform loader path
     * ( $LD_LIBRARY_PATH on Linux, $DYLD_LIBRARY_PATH on macOS ) for a name
     * carrying no slash -- for a slashed name the two fallbacks below are the
     * only lookups beyond the run directory.
     *
     * The root is gBelfemDataPath, NOT $BELFEM_DATA directly. It is filled by
     * Communicator::set_globals() from the environment, and on an installed
     * tree from BELFEM_INSTALL_DATADIR when the environment is silent, so an
     * unset $BELFEM_DATA does not imply an empty root. A code may also set it
     * itself. An empty root skips both fallbacks; an empty aSubDirectory
     * searches the data directory itself.
     *
     * gastables::data_path() deliberately does NOT use this. It is a
     * directory resolver with marker-file validation and a relative-path
     * ladder, and it intentionally skips validation when the root is set
     * (fn_GT_data_path.cpp:51-55). That is not duplication to be collapsed.
     */
    std::string
    search_data_file( const std::string & aFile,
                      const std::string & aSubDirectory );

//------------------------------------------------------------------------------
} /* namespace belfem */
//------------------------------------------------------------------------------
#endif //BELFEM_FILE_TOOLS_HPP
