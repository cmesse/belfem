//
// Created by christian on 12/22/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HPHI_TRI6_HPP
#define BELFEM_CL_IWG_MAXWELL_HPHI_TRI6_HPP

#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_HPhi_Tri6 : public IWG_Maxwell
        {
            // work vector for thin shell stiffnesses
            real mWork[ 16 ];

            // tangential components of H per layer
            Vector< real > mHt ;

            // curl operator for thin shell
            Matrix< real > mC ;

            // link to function
            void
            ( IWG_Maxwell_HPhi_Tri6::*mFunJacobian )
            (       Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );


            // matrix for element values that need to be integrated
            // over thin shells
            Matrix< real > mLayerData ;

            // work vector for normal component of H
            Vector< real > mH ;

            // edge interpolation vector
            Vector< real > mE ;

            // edge index vector
            Vector< uint > mU ;

            Vector< int > mPivot ;
            Vector< real > mAlpha;
            Vector< real > mBeta;
            Vector< real > mGamma ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_HPhi_Tri6();

//------------------------------------------------------------------------------

            ~IWG_Maxwell_HPhi_Tri6();

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
                    Vector< real > & aRHS ) ;
//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_ferro(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_air(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_scair(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_scfm(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_fmair(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_wave(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_farfield(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_cut(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_thinshell(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_layer_mass(
                    const uint aLayer,
                    const real aE0,
                    const real aE1,
                    Matrix< real > & aM );
//------------------------------------------------------------------------------

            void
            compute_layer_stiffness(
                    const uint aLayer,
                    const uint aIntPoint,
                    const real aE0,
                    const real aE1,
                    const Vector< real > & aHt,
                    const real aHn,
                    const real adHndx,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_thinshell( Element * aElement );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            ( this->*mFunJacobian ) ( aElement, aJacobian, aRHS ) ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_superconductor(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_superconductor(
                aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_air_phi_higher_order(
                    aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_ferro(
            Element        * aElement,
            Matrix< real > & aJacobian,
            Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_ferro(
                    aElement, aJacobian, aRHS );
        }

 //------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_scfm(
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

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_cut(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_cut( aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IWG_MAXWELL_HPHI_TRI6_HPP
