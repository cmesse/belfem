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

#ifndef BELFEM_CL_IWG_MAXWELLTHERMAL_HPP
#define BELFEM_CL_IWG_MAXWELLTHERMAL_HPP

#include "cl_IWG_TransientHeatConduction.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief Transient heat conduction coupled to the Maxwell solution.
         *
         * @ingroup grp_fem_thermal
         * @see @ref fem_thermal_index
         */
        class IWG_MaxwellThermal : public IWG_TransientHeatConduction
        {

            // link to function; set per group by link_to_group(). Never
            // called while nullptr: link_to_group() guards the assembled
            // case with an always-active error
            void
            ( *mFunMKF )
                      ( Calculator * aCalc, TimestepMatrices * aMatrices) = nullptr;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

            IWG_MaxwellThermal( ModelDimensionality  aModelDimensionality,
                    const IwgType aType=IwgType::MaxwellThermal,
                    const IwgMode aMode=IwgMode::Iterative );

            ~IWG_MaxwellThermal() override = default ;

            void
            link_to_group( Group  * aGroup ) override;

            void
            create_custom_vectors_and_matrices( Calculator * aCalc ) override ;

//------------------------------------------------------------------------------
            protected:
//------------------------------------------------------------------------------

                void
                compute_mkf( Element * aElement ) override;

//------------------------------------------------------------------------------

        };
    }
}

#endif //BELFEM_CL_IWG_MAXWELLTHERMAL_HPP