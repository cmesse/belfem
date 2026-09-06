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
                const uint mBlockIntegrationOrder ;
                const uint mSideSetIntegrationOrder ;
                const IntegrationScheme & mIntegrationScheme ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                Parameters( Kernel * aKernel );

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
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */
#endif //BELFEM_CL_FEM_DOFMGR_PARAMETERS_HPP
