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
         * Base class for the transport properties of a Helmholtz fluid.
         *
         * The default implementation does NOT return a value: mu() and
         * lambda() raise BELFEM_ERROR, so asking a fluid without a transport
         * correlation for one aborts the run rather than handing back a
         * silently wrong number. A species with a real correlation derives
         * from this class: HelmholtzTransport_Methane,
         * HelmholtzTransport_Hydrogen and HelmholtzTransport_LemmonJacobsen
         * ( nitrogen, oxygen ).
         *
         * Every Helmholtz fluid Gas can build currently carries one of them;
         * the base is only reached by a fluid added without a correlation.
         */
        class HelmholtzTransport
        {
//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            Helmholtz         & mEoS ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            HelmholtzTransport( Gas & aParent );

            virtual ~HelmholtzTransport() = default ;

//----------------------------------------------------------------------------

            virtual real
            mu( const real T, const real p ) const;

//----------------------------------------------------------------------------

            virtual real
            lambda( const real T, const real p ) const;

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_HelmholtzTransport_HPP
