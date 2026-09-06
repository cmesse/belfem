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

#ifndef BELFEM_CL_IWG_POISSON_HPP
#define BELFEM_CL_IWG_POISSON_HPP

#include "cl_IWG_Timestep.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Poisson : public IWG_Timestep
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             *
             * @param aModelDimensionality The IWG has no direct access to the
             * mesh, so the dimension information must be passed here
             *
             * @param aType usually as default, but children can override the type
             *              from their constructor
             *
             * @param aMode direct or iterative assembly mode
             */
            IWG_Poisson (
                    const ModelDimensionality  aModelDimensionality,
                    const IwgType aType=IwgType::Poisson,
                    const IwgMode aMode=IwgMode::Direct );

//------------------------------------------------------------------------------

            ~IWG_Poisson() override = default ;

//------------------------------------------------------------------------------

            void
            compute_jacobian(  Element        * aElement,
                               Matrix< real > & aJacobian ) override;

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_IWG_POISSON_HPP
