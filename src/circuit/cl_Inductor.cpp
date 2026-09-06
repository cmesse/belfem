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

#include "cl_Inductor.hpp"
#include "fn_circuit_state_io.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        Inductor::Inductor(const real aValue, const uint aOrder, Cell < ElectricNode* > & aTerminals, string aLabel) :
                TwoTerminals(aTerminals,aLabel), mValue(aValue), mH(ShiftRegister<real>(aOrder)), mI(ShiftRegister<real>(aOrder)), mOrder(aOrder)
        {
            mBDF = new ode::BDF(mH) ;
        }

//----------------------------------------------------------------------------

        Inductor::~Inductor()
        {
            delete mBDF ;
        }

//----------------------------------------------------------------------------

        void
        Inductor::compute_current()
        {
            //Computing the current of the inductor for time discretization
            this->set_current((mTerminals[0]->get_voltage()-mTerminals[1]->get_voltage())/mRL + miLh);
        }

//------------------------------------------------------------------------------

        void
        Inductor::shift( const real aTime, const real aDeltaTime )
        {
            //Update the time step
            mDeltaTime = aDeltaTime;
            mH.push(aDeltaTime);
            mI.push(this->get_current());

            this->update_companions() ;

            // Reset the current
            this->set_current(0.0) ;
        }

//------------------------------------------------------------------------------

        void
        Inductor::update_companions()
        {
            BELFEM_ASSERT( ! mH.empty(), "update_companions() called on an empty history register" );

            //Compute the coefficients
            mBDF->compute_coefficients();

            //Update the discretized resistance
            mRL =  (mValue * mBDF->coefficients()(0))/mH(0);

            //Update the discretized current source
            miLh = 0.0;
            for (uint i = 1; i < mBDF->coefficients().length(); ++i)
            {
                miLh -= mBDF->coefficients()(i) * mI(i-1);
            }
            miLh /= mBDF->coefficients()(0);
        }

//------------------------------------------------------------------------------

        void
        Inductor::shift_back(  )
        {
            mH.revert() ;
            mI.revert() ;

            // rebuild the ACCEPTED step's companions from the reverted
            // registers, so the circuit walk's compute_current() reconstructs
            // the accepted current exactly before the retry pushes it back
            // into the history. A first-step rejection reverts to empty
            // registers: keep the set_timestep() seed and the zero source
            if ( ! mH.empty() )
            {
                this->update_companions() ;
            }
        }

//------------------------------------------------------------------------------

        void
        Inductor::set_timestep( const real aDeltaTime )
        {
            mDeltaTime = aDeltaTime ;
            mRL = mValue/mDeltaTime ;
        }

//------------------------------------------------------------------------------

        void
        Inductor::save_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            save_shift_register( aGroup, aPrefix + "h", mH );
            save_shift_register( aGroup, aPrefix + "i", mI );
            real tCurrent = this->get_current() ;
            hdf5::save_scalar_to_file( aGroup, aPrefix + "current", tCurrent, tStatus );
#endif
        }

//------------------------------------------------------------------------------

        void
        Inductor::load_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            load_shift_register( aGroup, aPrefix + "h", mH );
            load_shift_register( aGroup, aPrefix + "i", mI );
            BELFEM_ERROR( mH.size() == mI.size(),
                "Registers %sh and %si in the restart file have different lengths ( corrupt dump ). Delete the memdump to restart with a cold circuit.",
                aPrefix.c_str(), aPrefix.c_str() );
            real tCurrent = 0.0 ;
            hdf5::load_scalar_from_file( aGroup, aPrefix + "current", tCurrent, tStatus );
            this->set_current( tCurrent );

            // rebuild the companions so a stamp without a prior shift is coherent
            if ( ! mH.empty() )
            {
                this->update_companions() ;
            }
#endif
        }

//------------------------------------------------------------------------------

        real
        Inductor::get_value() const
        {
            return mValue ;
        }

//------------------------------------------------------------------------------

        real
        Inductor::get_discretized_resistance() const
        {
            return mRL ;
        }

//------------------------------------------------------------------------------

        real
        Inductor::get_discretized_source() const
        {
            return miLh ;
        }

//------------------------------------------------------------------------------
    }
}