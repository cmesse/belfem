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
     * Content revision of the cached rho databases ( <label>_RRR<n>.hdf5 ).
     *
     * The file name encodes only the material label and the integer RRR, so a
     * table sampled from an older law is indistinguishable by name. Bump this
     * value whenever anything sampled into the table changes: a Kohler law or
     * its switch points, rho_0( T ), the alloy homogenization, or the grid.
     * Every cache written under another value is rebuilt once.
     */
    constexpr uint gRhoDatabaseRevision = 1 ;

    /**
     * Test whether a cached rho database may be loaded: it must carry the
     * "RRR" dataset and a one-element "revision" dataset equal to
     * gRhoDatabaseRevision. Files without the revision were written before it
     * existed, by a law that may have changed since.
     *
     * Returns true when the file may be loaded, false when it must be rebuilt.
     * Never throws: the caller probes on the master only and broadcasts the
     * verdict, so an unreadable or malformed file counts as "must rebuild",
     * the recoverable outcome.
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
        bool tIsCurrent = H5Lexists( tFile, "RRR", H5P_DEFAULT ) > 0
                       && H5Lexists( tFile, "revision", H5P_DEFAULT ) > 0 ;

        if ( tIsCurrent )
        {
            tIsCurrent = false ;

            hid_t tDataSet = H5Dopen2( tFile, "revision", H5P_DEFAULT );

            if ( tDataSet >= 0 )
            {
                hid_t tDataSpace = H5Dget_space( tDataSet );

                hid_t tDataType = H5Dget_type( tDataSet );

                if ( tDataSpace >= 0 && tDataType >= 0 )
                {
                    // the read below fills exactly one value, and only from an
                    // integer: HDF5 would convert a stored 1.5 to 1
                    if (    H5Sget_simple_extent_npoints( tDataSpace ) == 1
                         && H5Tget_class( tDataType ) == H5T_INTEGER )
                    {
                        unsigned int tRevision = 0 ;

                        if ( H5Dread( tDataSet, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
                                      H5P_DEFAULT, &tRevision ) >= 0 )
                        {
                            tIsCurrent = tRevision == gRhoDatabaseRevision ;
                        }
                    }
                }
                if ( tDataType >= 0 )
                {
                    H5Tclose( tDataType );
                }
                if ( tDataSpace >= 0 )
                {
                    H5Sclose( tDataSpace );
                }
                H5Dclose( tDataSet );
            }
        }

        H5Fclose( tFile );

        return tIsCurrent;
#else
        // without HDF5 the database path cannot work at all; keep the
        // pre-existing behavior rather than forcing a rebuild here
        return true;
#endif
    }
}
#endif //BELFEM_FN_RHO_DATABASE_IS_CURRENT_HPP
