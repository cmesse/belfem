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

#include "cl_Switch.hpp"
#ifdef BELFEM_HDF5
#include "hdf5_tools.hpp"
#endif

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        Switch::Switch( const bool aIsClosed, const real aSwitchTime, Cell < ElectricNode* > & aTerminals, string aLabel ) :
                TwoTerminals(aTerminals,aLabel), mIsClosed(aIsClosed), mSwitchTime(aSwitchTime),
                mPrevIsClosed(aIsClosed)
        {}

//----------------------------------------------------------------------------

        bool
        Switch::is_closed() const
        {
            return mIsClosed ;
        }

//------------------------------------------------------------------------------

        void
        Switch::switch_state()
        {
            mIsClosed = !mIsClosed ;
        }

//------------------------------------------------------------------------------

        void
        Switch::shift( const real aTime, const real aDeltaTime )
        {
            // snapshot for a one-deep shift_back(), matching the L/C shift registers
            mPrevIsClosed = mIsClosed ;
            mPrevIsSwitched = mIsSwitched ;

            if (aTime >= mSwitchTime && !mIsSwitched)
            {
                mIsSwitched = true;
                this->switch_state() ;
            }
        }

//------------------------------------------------------------------------------

        void
        Switch::shift_back()
        {
            // undo a latch that fired during a rejected timestep attempt
            mIsClosed = mPrevIsClosed ;
            mIsSwitched = mPrevIsSwitched ;
        }

//------------------------------------------------------------------------------

        void
        Switch::save_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            hdf5::save_bool_to_file( aGroup, aPrefix + "is_closed", mIsClosed, tStatus );
            hdf5::save_bool_to_file( aGroup, aPrefix + "is_switched", mIsSwitched, tStatus );
#endif
        }

//------------------------------------------------------------------------------

        void
        Switch::load_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            hdf5::load_bool_from_file( aGroup, aPrefix + "is_closed", mIsClosed, tStatus );
            hdf5::load_bool_from_file( aGroup, aPrefix + "is_switched", mIsSwitched, tStatus );

            // the rollback snapshots are attempt-local; a restored switch has
            // no attempt to revert, so they mirror the restored state
            mPrevIsClosed = mIsClosed ;
            mPrevIsSwitched = mIsSwitched ;
#endif
        }
    }
}
