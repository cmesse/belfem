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

#ifndef BELFEM_FN_GET_PATH_HPP
#define BELFEM_FN_GET_PATH_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace gastables
    {
        /**
         * Return the directory holding the gas tables.
         *
         * If gBelfemDataPath is set, which Communicator::set_globals() does from
         * $BELFEM_DATA and a config file may also do, its fluid subdirectory is
         * returned without being checked, so a wrong path fails where the file is
         * opened rather than falling through to a different directory. Otherwise
         * share/fluid is searched relative to the working directory. Returns an
         * empty string if neither holds the tables.
         */
        string
        data_path();

        /**
         * Test whether the gas tables can be found. Never aborts, so a caller can
         * degrade gracefully rather than fail to construct.
         */
        bool
        data_available();

    }
}
#endif //BELFEM_FN_GET_PATH_HPP
