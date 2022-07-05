//
// Created by Christian Messe on 23.11.20.
//

#ifndef BELFEM_CL_GM_HelmholtzTransport_HPP
#define BELFEM_CL_GM_HelmholtzTransport_HPP

#include "typedefs.hpp"
#include "cl_Gas.hpp"
#include "cl_GT_RefGas.hpp"

#include "cl_GM_EoS.hpp"
#include "cl_GM_Helmholtz.hpp"
namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        /**
         * a special class for trans properties of liqid gases
         */
        class HelmholtzTransport
        {
//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            Helmholtz         & mEoS ;

            gastables::RefGas & mRefgas ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            HelmholtzTransport( Gas & aParent );

            virtual ~HelmholtzTransport() = default ;

//----------------------------------------------------------------------------

            virtual real
            mu( const real & aT, const real & aP );

//----------------------------------------------------------------------------

            virtual real
            lambda( const real & aT, const real & aP );

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_HelmholtzTransport_HPP
