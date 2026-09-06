//
// Created by christian on 10/27/25.
//

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