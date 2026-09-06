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
#include "cl_FEM_DofMgr_Parameters.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"

namespace belfem
{
    namespace fem
    {
        namespace dofmgr
        {
//-----------------------------------------------------------------------------

            Parameters::Parameters( Kernel * aKernel  ):
                mBlockIntegrationOrder( 0 ),
                mSideSetIntegrationOrder( 0 ),
                mIntegrationScheme( aKernel->params()->integration_scheme() )
            {

            }

//-----------------------------------------------------------------------------

            void
            Parameters::print()
            {
                std::cout << "Block Integration Order " << mBlockIntegrationOrder << std::endl;
                std::cout << "Sideset Integration Order " << mSideSetIntegrationOrder << std::endl;
                std::cout << "Integration Scheme " << to_string( mIntegrationScheme ) << std::endl;
            }

//-----------------------------------------------------------------------------
        }
    }
}
