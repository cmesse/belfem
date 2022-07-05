//
// Created by christian on 7/7/21.
//

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

            Parameters::Parameters( Kernel * aKernel,
                                    const index_t aFieldIndex  ):
                mBlockIntegrationOrder(
                        aKernel->params()->block_integration_order( aFieldIndex  ) ),
                mSideSetIntegrationOrder(
                        aKernel->params()->sideset_integration_order( aFieldIndex ) ),
                mIntegrationScheme( aKernel->params()->integration_scheme() ),
                mEnforceLinear( aKernel->params()->linear_enforcement_flag( aFieldIndex ) )
            {

            }

//-----------------------------------------------------------------------------

            void
            Parameters::print()
            {
                std::cout << "Block Integration Order " << mBlockIntegrationOrder << std::endl;
                std::cout << "Sideset Integration Order " << mSideSetIntegrationOrder << std::endl;
                std::cout << "Integration Scheme " << to_string( mIntegrationScheme ) << std::endl;
                std::cout << "Enforce Linear  " << this->enforce_linear() << std::endl;
            }

//-----------------------------------------------------------------------------
        }
    }
}
