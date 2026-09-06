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

#include "cl_Capacitor.hpp"
#include "fn_circuit_state_io.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        Capacitor::Capacitor(const real aValue, const uint aOrder, Cell < ElectricNode* > & aTerminals, string aLabel) :
                TwoTerminals(aTerminals,aLabel), mValue(aValue), mH(ShiftRegister<real>(aOrder)), mV(ShiftRegister<real>(aOrder)), mOrder(aOrder)
        {
            mBDF = new ode::BDF(mH) ;
        }

//----------------------------------------------------------------------------

        Capacitor::~Capacitor()
        {
            delete mBDF ;
        }

//----------------------------------------------------------------------------

        void
        Capacitor::compute_current()
        {
            //Computing the current of the capacitor for time discretization
            this->set_current((mTerminals[0]->get_voltage()-mTerminals[1]->get_voltage())/mRC + miCh);
        }

//------------------------------------------------------------------------------

        void
        Capacitor::shift( const real aTime, const real aDeltaTime )
        {

            //Update the time step
            mDeltaTime = aDeltaTime ;
            mH.push(aDeltaTime) ;
            mV.push(this->get_node_plus()->get_voltage() - this->get_node_minus()->get_voltage()) ;

            this->update_companions() ;

            // Reset the current
            this->set_current(0.0) ;
        }

//------------------------------------------------------------------------------

        void
        Capacitor::update_companions()
        {
            BELFEM_ASSERT( ! mH.empty(), "update_companions() called on an empty history register" );

            //Compute the coefficients
            mBDF->compute_coefficients() ;

            //Update the discretized resistance
            mRC = mH(0)/(mValue*mBDF->coefficients()(0)) ;

            //Update the discretized current source
            miCh = 0.0 ;
            for (uint i = 1 ; i < mBDF->coefficients().length(); ++i )
            {
                miCh += mBDF->coefficients()(i) * mV(i-1) ;
            }
            miCh *= mValue/mH(0) ;
        }

//------------------------------------------------------------------------------

        void
        Capacitor::shift_back(  )
        {
            mH.revert() ;
            mV.revert() ;

            // rebuild the ACCEPTED step's companions from the reverted
            // registers, so the circuit walk's compute_current() restores
            // the accepted current ( latent for the capacitor — its history
            // stores voltages — but keeps the rollback state coherent ).
            // A first-step rejection reverts to empty registers: keep the
            // rejected shift's BDF1 seed companions ( a zero source only on
            // a zero-voltage cold start — the pushed quantity is a voltage )
            if ( ! mH.empty() )
            {
                this->update_companions() ;
            }
        }

//------------------------------------------------------------------------------

        void
        Capacitor::set_timestep( const real aDeltaTime )
        {
            mDeltaTime = aDeltaTime ;
            mRC = mDeltaTime/mValue ;
        }

//------------------------------------------------------------------------------

        void
        Capacitor::save_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            save_shift_register( aGroup, aPrefix + "h", mH );
            save_shift_register( aGroup, aPrefix + "v", mV );
            real tCurrent = this->get_current() ;
            hdf5::save_scalar_to_file( aGroup, aPrefix + "current", tCurrent, tStatus );
#endif
        }

//------------------------------------------------------------------------------

        void
        Capacitor::load_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            load_shift_register( aGroup, aPrefix + "h", mH );
            load_shift_register( aGroup, aPrefix + "v", mV );
            BELFEM_ERROR( mH.size() == mV.size(),
                "Registers %sh and %sv in the restart file have different lengths ( corrupt dump ). Delete the memdump to restart with a cold circuit.",
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
        Capacitor::get_value() const
        {
            return mValue ;
        }

//------------------------------------------------------------------------------

        real
        Capacitor::get_discretized_resistance() const
        {
            return mRC ;
        }

//------------------------------------------------------------------------------

        real
        Capacitor::get_discretized_source() const
        {
            return miCh ;
        }

//------------------------------------------------------------------------------
    }

}
