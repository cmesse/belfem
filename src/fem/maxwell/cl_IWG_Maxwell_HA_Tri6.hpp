//
// Created by christian on 3/23/22.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HA_TRI6_HPP
#define BELFEM_CL_IWG_MAXWELL_HA_TRI6_HPP

#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_HA_Tri6 : public IWG_Maxwell
        {
            // link to function
            void
            ( IWG_Maxwell_HA_Tri6::*mFunJacobian )
                    (       Element        * aElement,
                            Matrix< real > & aJacobian,
                            Vector< real > & aRHS );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_HA_Tri6();

//------------------------------------------------------------------------------

            ~IWG_Maxwell_HA_Tri6();

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            link_jacobian_function( Group * aGroup );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_superconductor(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );
//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_ferro(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_air(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_coil(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_interface(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            ( this->*mFunJacobian ) ( aElement, aJacobian, aRHS ) ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_superconductor(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_superconductor(
                    aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_ferro(
                    aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_interface(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->compute_interface_ha_tri6( aElement, aJacobian, mGroup->work_K() );

            // compute RHS
            aRHS = aJacobian * this->collect_q0_ha_2d( aElement );

            // compone matrices
            aJacobian += mGroup->work_K() * this->timestep() ;
        }

//------------------------------------------------------------------------------
    }
}


#endif //BELFEM_CL_IWG_MAXWELL_HA_TRI6_HPP
