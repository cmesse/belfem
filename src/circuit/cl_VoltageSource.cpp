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

#include "cl_VoltageSource.hpp"
#include "constants.hpp"
#ifdef BELFEM_HDF5
#include "hdf5_tools.hpp"
#endif

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        VoltageSource::VoltageSource( SourceFunction * aFunction, Cell < ElectricNode * > & aTerminals, string aLabel) :
                TwoTerminals(aTerminals, aLabel), mFunction(aFunction)
        {}

//----------------------------------------------------------------------------

        VoltageSource::~VoltageSource()
        {
            delete mFunction ;
        }

//----------------------------------------------------------------------------

        real
        VoltageSource::get_value() const
        {
            return mValueTime ;
        }

//-----------------------------------------------------------------------------

        void
        VoltageSource::set_value( const real aValue )
        {
            mValueTime = aValue ;
        }

//-----------------------------------------------------------------------------

        void
        VoltageSource::shift( const real aTime, const real aDeltaTime )
        {
            mValueTime = mFunction->compute( aTime ) ;
        }

//------------------------------------------------------------------------------

        void
        VoltageSource::save_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            hdf5::save_scalar_to_file( aGroup, aPrefix + "value_time", mValueTime, tStatus );
#endif
        }

//------------------------------------------------------------------------------

        void
        VoltageSource::load_state( hid_t aGroup, const string & aPrefix )
        {
#ifdef BELFEM_HDF5
            // a stamp without a prior shift() reads mValueTime; the next
            // shift() overwrites it from the source function
            herr_t tStatus = 0 ;
            hdf5::load_scalar_from_file( aGroup, aPrefix + "value_time", mValueTime, tStatus );
#endif
        }

//------------------------------------------------------------------------------
    }
}
