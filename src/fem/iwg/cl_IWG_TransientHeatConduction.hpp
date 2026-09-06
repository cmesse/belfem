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

#ifndef BELFEM_CL_IWG_TRANSIENTHEATCONDUCTION_HPP
#define BELFEM_CL_IWG_TRANSIENTHEATCONDUCTION_HPP

#include "cl_IWG_Timestep.hpp"

/**
 * look at :
 *     o   src/fem/maxwell/matrices/mt_maxwell_phi.cpp
 *     o   cl_IWG_StaticHeatConduction.cpp
 */
namespace belfem
{
    namespace fem
    {
        class IWG_TransientHeatConduction : public IWG_Timestep
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_TransientHeatConduction( const ModelDimensionality  aModelDimensionality,
                    const IwgType aType=IwgType::TransientHeatConduction,
                    const IwgMode aMode=IwgMode::Iterative );

            ~IWG_TransientHeatConduction() override = default;

//------------------------------------------------------------------------------

        protected:

            void
            compute_mkf( Element * aElement );

        private:


        };
    }
}
#endif //BELFEM_CL_IWG_TRANSIENTHEATCONDUCTION_HPP