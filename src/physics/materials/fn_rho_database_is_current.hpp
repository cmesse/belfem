/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_FN_RHO_DATABASE_IS_CURRENT_HPP
#define BELFEM_FN_RHO_DATABASE_IS_CURRENT_HPP

#include "typedefs.hpp"
#include "hdf5_tools.hpp"

namespace belfem
{
    /**
     * Test whether a cached rho database was written by the current format.
     *
     * The file name encodes only the material label and the RRR value, so a file
     * left over from an older BELFEM is indistinguishable by name and would be
     * loaded rather than rebuilt. Older files store the ratio as a lowercase
     * "rrr" ( alongside a "lambda" dataset ); the current writer stores "RRR".
     * The presence of that marker is therefore the format test.
     *
     * Returns true when the file may be loaded, false when it must be rebuilt.
     * An unreadable file counts as "must rebuild" rather than an error — the
     * caller regenerates, which is the recoverable outcome.
     */
    inline bool
    rho_database_is_current( const string & aPath )
    {
#ifdef BELFEM_HDF5
        hid_t tFile = H5Fopen( aPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

        if ( tFile < 0 )
        {
            return false;
        }

        // H5Lexists returns a tri-state: > 0 exists, 0 absent, < 0 error.
        // Only a definite hit counts as current -- a probe ERROR must land
        // on the rebuild side, not be swallowed by a bool conversion
        htri_t tExists = H5Lexists( tFile, "RRR", H5P_DEFAULT );

        H5Fclose( tFile );

        return tExists > 0;
#else
        // without HDF5 the database path cannot work at all; keep the
        // pre-existing behavior rather than forcing a rebuild here
        return true;
#endif
    }
}
#endif //BELFEM_FN_RHO_DATABASE_IS_CURRENT_HPP
