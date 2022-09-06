//
// Created by christian on 12/17/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HPHI_TET10_HPP
#define BELFEM_CL_IWG_MAXWELL_HPHI_TET10_HPP

#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_HPhi_Tet10 : public IWG_Maxwell
        {
            // link to function
            void
            ( IWG_Maxwell_HPhi_Tet10::*mFunJacobian )
            (       Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_HPhi_Tet10();

//------------------------------------------------------------------------------

            ~IWG_Maxwell_HPhi_Tet10() ;

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
            compute_jacobian_and_rhs_air(
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
            compute_jacobian_and_rhs_scair(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_scfm(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_fmair(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_cut(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            ( this->*mFunJacobian ) ( aElement, aJacobian, aRHS ) ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_superconductor(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // call the function from the parent class
            IWG_Maxwell::compute_jacobian_and_rhs_superconductor(
                    aElement,
                    aJacobian,
                    aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_air_phi_higher_order(
                    aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_ferro(
                    aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_scfm(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->compute_interface_ha_tet10( aElement, aJacobian, mGroup->work_K() );

            // compute RHS
            aRHS = aJacobian * this->collect_q0_ha_3d( aElement );

            // combine matrices
            aJacobian += mGroup->work_K() * this->timestep() ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_cut(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_cut( aElement, aJacobian, aRHS );
        }
//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IWG_MAXWELL_HPHI_TET10_HPP
