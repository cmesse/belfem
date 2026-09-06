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

#ifndef BELFEM_FN_CIRCUIT_STATE_IO_HPP
#define BELFEM_FN_CIRCUIT_STATE_IO_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_ShiftRegister.hpp"
#include "hdf5_types.hpp"
#ifdef BELFEM_HDF5
#include "hdf5_tools.hpp"
#endif

namespace belfem
{
    namespace electronics
    {
//------------------------------------------------------------------------------

        /**
         * write one time-stepping history register into the restart group,
         * as `aLabel` ( samples, newest first ) and `aLabel_cap` ( capacity )
         */
        inline void
        save_shift_register(
                hid_t                       & aGroup,
                const string                & aLabel,
                const ShiftRegister< real > & aRegister )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;

            uint tCapacity = ( uint ) aRegister.capacity() ;
            hdf5::save_scalar_to_file( aGroup, aLabel + "_cap", tCapacity, tStatus );

            Vector< real > tSamples( aRegister.size() );
            for ( index_t k = 0; k < aRegister.size(); ++k )
            {
                tSamples( k ) = aRegister( k );
            }
            hdf5::save_vector_to_file( aGroup, aLabel, tSamples, tStatus );
#endif
        }

//------------------------------------------------------------------------------

        /**
         * restore one history register saved by save_shift_register().
         *
         * The rebuild is clear() + oldest-first push(): the visible window is
         * identical to the dumped one. The revert flag after a FULL rebuild is
         * CanRevert rather than the live CanRevertFull ( the backup slot is not
         * refilled ) — calling shift_back() before the first post-restore
         * shift() is therefore FORBIDDEN; every production reject path shifts
         * first.
         */
        inline void
        load_shift_register(
                hid_t                 & aGroup,
                const string          & aLabel,
                ShiftRegister< real > & aRegister )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;

            uint tCapacity = 0 ;
            hdf5::load_scalar_from_file( aGroup, aLabel + "_cap", tCapacity, tStatus );
            BELFEM_ERROR( tCapacity == ( uint ) aRegister.capacity(),
                "Register %s in the restart file has capacity %u, the circuit expects %u ( BDF order changed? ). Delete the memdump to restart with a cold circuit.",
                aLabel.c_str(),
                tCapacity,
                ( uint ) aRegister.capacity() );

            Vector< real > tSamples ;
            hdf5::load_vector_from_file( aGroup, aLabel, tSamples, tStatus );
            BELFEM_ERROR( tSamples.length() <= aRegister.capacity(),
                "Register %s in the restart file holds %lu samples, more than its stated capacity %u. Delete the memdump to restart with a cold circuit.",
                aLabel.c_str(),
                ( long unsigned int ) tSamples.length(),
                tCapacity );

            aRegister.clear() ;
            for ( index_t k = tSamples.length(); k > 0; --k )
            {
                aRegister.push( tSamples( k - 1 ) );
            }
#endif
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_CIRCUIT_STATE_IO_HPP
