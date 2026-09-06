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

#ifndef BELFEM_FN_MATERIAL_DATA_PATH_HPP
#define BELFEM_FN_MATERIAL_DATA_PATH_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace material
    {
        /**
         * Return the directory holding the shared material databases, that is
         * the material subdirectory of gBelfemDataPath. Returns an empty string
         * when that global is empty.
         *
         * Note that an unset $BELFEM_DATA does not imply an empty global:
         * Communicator::set_globals() fills it from the environment, and on an
         * installed tree falls back to the compiled-in install data directory
         * when the environment is silent.
         */
        string
        data_path();

        /**
         * Resolve a material data file named in an input file.
         *
         * The run directory always wins, so a local copy overrides the shared
         * database. If the file is not there, it is looked up below data_path(),
         * first under the same relative path and then by name alone, so that
         * "MatData/bhdata.hdf5" also resolves for a run that has no MatData
         * directory of its own.
         *
         * If nothing is found, aFile is returned unchanged, so that whoever opens
         * it reports the name the user wrote. This also keeps the loader search
         * intact for the plugin libraries, which need not exist as a file at all.
         *
         * The search order itself lives in belfem::search_data_file()
         * ( io/filetools.hpp ); this is a forwarder supplying the material
         * subdirectory. It is shared with the source-function plugin loader in
         * numerics/sources, which sits below physics and cannot include this.
         */
        string
        data_file( const string & aFile );

    }
}
#endif //BELFEM_FN_MATERIAL_DATA_PATH_HPP
