//
// Created by christian on 7/7/21.
//

#ifndef BELFEM_CL_FEM_DOFMGR_PARAMETERS_HPP
#define BELFEM_CL_FEM_DOFMGR_PARAMETERS_HPP

#include "typedefs.hpp"
#include "en_IntegrationScheme.hpp"
namespace belfem
{
    class Mesh ;

    namespace fem
    {
        class Kernel ;

        namespace dofmgr
        {
//-----------------------------------------------------------------------------

            /**
             * this is a parameter object that contains the relevant information
             * for managing the DOFs
             */
            class Parameters
            {
                const uint & mBlockIntegrationOrder ;
                const uint & mSideSetIntegrationOrder ;
                const IntegrationScheme & mIntegrationScheme ;
                const uint & mEnforceLinear ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                Parameters( Kernel * aKernel,
                            const index_t aFieldIndex );

//------------------------------------------------------------------------------

                ~Parameters() = default ;

//------------------------------------------------------------------------------

                uint
                block_integration_order() const ;

//------------------------------------------------------------------------------

                uint
                sideset_integration_order() const ;

//------------------------------------------------------------------------------

                IntegrationScheme
                integration_scheme() const ;

//------------------------------------------------------------------------------

                bool
                enforce_linear() const ;

//------------------------------------------------------------------------------

                void
                print();

//------------------------------------------------------------------------------
            };

//------------------------------------------------------------------------------

            inline uint
            Parameters::block_integration_order() const
            {
                return mBlockIntegrationOrder ;
            }

//------------------------------------------------------------------------------

            inline uint
            Parameters::sideset_integration_order() const
            {
                return mSideSetIntegrationOrder ;
            }

//------------------------------------------------------------------------------

            inline IntegrationScheme
            Parameters::integration_scheme() const
            {
                return mIntegrationScheme ;
            }

//------------------------------------------------------------------------------

            inline bool
            Parameters::enforce_linear() const
            {
                return mEnforceLinear > 0 ;
            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */
#endif //BELFEM_CL_FEM_DOFMGR_PARAMETERS_HPP
