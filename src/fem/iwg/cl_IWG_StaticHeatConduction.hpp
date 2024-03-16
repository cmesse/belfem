//
// Created by Christian Messe on 24.07.20.
//

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
             * @param aNumberOfDimensions The IWG has no direct access to the
             * mesh, so the dimension information must be passed here
             *
             * @param aType usally as default, but children can override the type
             *              from their constructor
             */
            IWG_StaticHeatConduction (
                    const ModelDimensionality  aModelDimensionality,
                    const IwgType aType=IwgType::StaticHeatConduction );

//------------------------------------------------------------------------------

            virtual ~IWG_StaticHeatConduction () = default;

//------------------------------------------------------------------------------
// Functions called by Field during assembly
//------------------------------------------------------------------------------

            virtual void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            virtual void
            link_to_group( Group * aGroup );

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
                    Vector< real > & aConvection ) ;

//------------------------------------------------------------------------------

            void
            compute_alpha_boundary_condition(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

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
