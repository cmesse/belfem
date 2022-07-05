//
// Created by Christian Messe on 23.11.20.
//
#include "assert.hpp"
#include "constants.hpp"
#include "cl_GM_HelmholtzTransport.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        HelmholtzTransport::HelmholtzTransport( Gas & aParent ) :
            mEoS( *reinterpret_cast< Helmholtz * >( aParent.eos() ) ),
            mRefgas( * aParent.component( 0 ) )
        {

            BELFEM_ERROR( aParent.number_of_components() == 1,
                         "Parent must have only one gas" );
        }

//----------------------------------------------------------------------------


        real
        HelmholtzTransport::mu( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "mu function is not implemented for this gas" );
            return BELFEM_QUIET_NAN;
        }

//----------------------------------------------------------------------------

        real
        HelmholtzTransport::lambda( const real & aT, const real & aP )
        {
            BELFEM_ERROR( false, "lambda function is not implemented for this gas" );
            return BELFEM_QUIET_NAN;
        }

//----------------------------------------------------------------------------
    }
}