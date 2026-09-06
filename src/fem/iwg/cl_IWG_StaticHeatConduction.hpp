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

#ifndef BELFEM_CL_IWG_STATICHEATCONDUCTION_HPP
#define BELFEM_CL_IWG_STATICHEATCONDUCTION_HPP

#include "cl_IWG_Timestep.hpp"
namespace belfem
{
    namespace fem
    {
        class IWG_StaticHeatConduction : public IWG_Timestep
        {

            // link to function
            void
            ( IWG_StaticHeatConduction::*mFunComputeJacobian )
                    (       Element        * aElement,
                            Matrix< real > & aJacobian,
                            Vector< real > & aRHS );

        protected:

            // matrix for thermal conductivity
            Matrix< real > mLambda ;

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
             */
            IWG_StaticHeatConduction (
                    const ModelDimensionality  aModelDimensionality,
                    const IwgType aType=IwgType::StaticHeatConduction );

//------------------------------------------------------------------------------

            ~IWG_StaticHeatConduction () override = default;

//------------------------------------------------------------------------------
// Functions called by Field during assembly
//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) override;

//------------------------------------------------------------------------------

            void
            link_to_group( Group * aGroup ) override;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_conduction(
                    Element        * aElement,
                    Matrix< real > & aK,
                    Vector< real > & aQ );

            void
            compute_convection(
                    Element        * aElement,
                    Vector< real > & aConvection ) override ;

//------------------------------------------------------------------------------

            void
            compute_alpha_boundary_condition(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) override;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline void
        IWG_StaticHeatConduction::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            ( this->*mFunComputeJacobian )( aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_IWG_STATICHEATCONDUCTION_HPP
