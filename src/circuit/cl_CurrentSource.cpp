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

#include "cl_CurrentSource.hpp"
#include "constants.hpp"
#ifdef BELFEM_HDF5
#include "hdf5_tools.hpp"
#endif


namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        CurrentSource::CurrentSource( SourceFunction * aFunction, Cell < ElectricNode* > & aTerminals, string aLabel) :
                TwoTerminals(aTerminals,aLabel), mFunction(aFunction)
        {
            //Initialize the current of the current source
            this->set_current(mFunction->compute(0.0)) ;
        }

//----------------------------------------------------------------------------

        CurrentSource::~CurrentSource()
        {
            delete mFunction ;
        }

//----------------------------------------------------------------------------

        real
        CurrentSource::get_value() const
        {
            return this->get_current() ;
        }

//-----------------------------------------------------------------------------

        void
        CurrentSource::shift( const real aTime, const real aDeltaTime )
        {
            this->set_current(mFunction->compute(aTime)) ;
        }

//------------------------------------------------------------------------------

        void
        CurrentSource::save_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            real tCurrent = this->get_current() ;
            hdf5::save_scalar_to_file( aGroup, aPrefix + "current", tCurrent, tStatus );
#endif
        }

//------------------------------------------------------------------------------

        void
        CurrentSource::load_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            // a stamp without a prior shift() reads the current; the next
            // shift() overwrites it from the source function
            herr_t tStatus = 0 ;
            real tCurrent = 0.0 ;
            hdf5::load_scalar_from_file( aGroup, aPrefix + "current", tCurrent, tStatus );
            this->set_current( tCurrent );
#endif
        }

//------------------------------------------------------------------------------
    }
}