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

#include "assert.hpp"
#include "constants.hpp"
#include "cl_GM_HelmholtzTransport.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        // the model check must run BEFORE the downcast is dereferenced
        // in the initializer list
        static Helmholtz *
        checked_helmholtz_eos( Gas & aParent )
        {
            BELFEM_ERROR( aParent.gas_model() == GasModel::HELMHOLTZ,
                         "HelmholtzTransport needs a gas with a Helmholtz EoS" );

            return static_cast< Helmholtz * >( aParent.eos() );
        }

//----------------------------------------------------------------------------

        HelmholtzTransport::HelmholtzTransport( Gas & aParent ) :
            mEoS( * checked_helmholtz_eos( aParent ) )
        {
            BELFEM_ERROR( aParent.number_of_components() == 1,
                         "Parent must have only one gas" );
        }

//----------------------------------------------------------------------------


        real
        HelmholtzTransport::mu( const real T, const real p ) const
        {
            BELFEM_ERROR( false, "mu function is not implemented for this gas" );
            return BELFEM_QUIET_NAN;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport::lambda( const real T, const real p ) const
        {
            BELFEM_ERROR( false, "lambda function is not implemented for this gas" );
            return BELFEM_QUIET_NAN;
        }

//----------------------------------------------------------------------------
    }
}