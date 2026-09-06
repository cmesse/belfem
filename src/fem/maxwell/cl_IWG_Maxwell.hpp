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

#ifndef CL_IWG_MAXWELL_HPP
#define CL_IWG_MAXWELL_HPP

#include "cl_IWG_Timestep.hpp"
#include "cl_Maxwell_FieldList.hpp"
#include "cl_FEM_Dof.hpp"
#include "en_Maxwell_Formulations.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief Main electromagnetic integral weak form (h-phi formulation).
         *
         * @ingroup grp_fem_maxwell
         * @see @ref fem_maxwell_maxwell_usage_guide
         */
        class IWG_Maxwell : public IWG_Timestep
        {
            // defines which formulation we have
            const maxwell::Formulation mFormulation ;

            //! contains the number of spatial dimensions
            const uint mNumberOfDimensions;

            const bool mHigherOrder ;

            const bool mUseEnrichment ;

            bool mUseEdges = false ;

            //! list with dofs per entity
            maxwell::FieldList mFields ;

            // link to function
            void
            ( *mFunMKF )
                      ( Calculator * aCalc,
                        TimestepMatrices * aMatrices) = nullptr ;

            bool mHaveThermal = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell(
            const maxwell::Formulation aFormulation,
            const ModelDimensionality aDimensionality,
            const bool aHigherOrder = false,
            const bool aUseEnrichment = false );

            ~IWG_Maxwell() override = default;

            void
            set_currents( Vector< real > & aI );

            void
            initialize() override ;

            void
            link_to_group( Group  * aGroup ) override ;

//------------------------------------------------------------------------------

            bool
            has_edge_dofs() const override ;

            void
            create_custom_vectors_and_matrices( Calculator * aCalc ) override ;

            void
            custom_postprocess() override;

            void
            collect_abstract_node_dofs() override ;

            bool
            have_thermal() const ;

//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            init_activation_maps() override;

            void
            compute_mkf( Element * aElement ) override ;

        };

        inline bool
        IWG_Maxwell::has_edge_dofs() const
        {
            return mUseEdges ;
        }

        inline bool
        IWG_Maxwell::have_thermal() const
        {
            return mHaveThermal ;
        }


    }
}
#endif // CL_IWG_MAXWELL_HPP
