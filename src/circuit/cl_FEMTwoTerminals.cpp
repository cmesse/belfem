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

#include "cl_FEMTwoTerminals.hpp"
#include "fn_circuit_state_io.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        FEMTwoTerminals::FEMTwoTerminals(Cell < ElectricNode* > & aTerminals, string aLabel) :
                TwoTerminals(aTerminals,aLabel), mH(ShiftRegister<real>(1)), mI(ShiftRegister<real>(1)), mV(ShiftRegister<real>(1))
        {}

//----------------------------------------------------------------------------

        void
        FEMTwoTerminals::compute_current()
        {
            this->set_current(mIn) ;
        }

//------------------------------------------------------------------------------

        void
        FEMTwoTerminals::shift( const real aTime, const real aDeltaTime )
        {
            mH.push(aDeltaTime) ;
            mI.push(mIn) ;
            mV.push(mVn) ;
        }

//------------------------------------------------------------------------------

        void
        FEMTwoTerminals::shift_back()
        {
            mH.revert() ;
            mI.revert() ;
            mV.revert() ;
        }

//------------------------------------------------------------------------------

        void
        FEMTwoTerminals::set_IV( const real I, const real V )
        {
            mIn = I ;
            mVn = V ;
            this->set_current(mIn) ;
        }

//------------------------------------------------------------------------------

        void
        FEMTwoTerminals::save_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            save_shift_register( aGroup, aPrefix + "h", mH );
            save_shift_register( aGroup, aPrefix + "i", mI );
            save_shift_register( aGroup, aPrefix + "v", mV );
            hdf5::save_scalar_to_file( aGroup, aPrefix + "in", mIn, tStatus );
            hdf5::save_scalar_to_file( aGroup, aPrefix + "vn", mVn, tStatus );
#endif
        }

//------------------------------------------------------------------------------

        void
        FEMTwoTerminals::load_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            load_shift_register( aGroup, aPrefix + "h", mH );
            load_shift_register( aGroup, aPrefix + "i", mI );
            load_shift_register( aGroup, aPrefix + "v", mV );
            hdf5::load_scalar_from_file( aGroup, aPrefix + "in", mIn, tStatus );
            hdf5::load_scalar_from_file( aGroup, aPrefix + "vn", mVn, tStatus );

            // a terminal pair is not an unknown-current dof: the Jacobian
            // reads get_current(), which update_components() cannot restore
            this->set_current( mIn );
#endif
        }

//------------------------------------------------------------------------------

        /*real
        FEMTwoTerminals::dIdV() const
        {
            return 1.0/mRFEM ;
        }*/

//------------------------------------------------------------------------------
    }
}