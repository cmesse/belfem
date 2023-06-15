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
                    const uint    aNumberOfDimensions,
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

            void
            compute_convection(
                    Element        * aElement,
                    Vector< real > & aConvection ) ;

//------------------------------------------------------------------------------

            virtual void
            compute_alpha_boundary_condition(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            virtual void
            allocate_work_matrices( Group    * aGroup );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_IWG_STATICHEATCONDUCTION_HPP
