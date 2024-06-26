//
// Created by Christian Messe on 12.01.22.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HPHI_TRI3_HPP
#define BELFEM_CL_IWG_MAXWELL_HPHI_TRI3_HPP

#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_HPhi_Tri3 : public IWG_Maxwell
        {

            //! link to function
            void
            ( IWG_Maxwell_HPhi_Tri3::*mFunJacobian )
            (       Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------
                public:
//------------------------------------------------------------------------------

                    IWG_Maxwell_HPhi_Tri3();

//------------------------------------------------------------------------------

                    ~IWG_Maxwell_HPhi_Tri3();

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
            compute_jacobian_and_rhs_fmfm(
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
            compute_jacobian_and_rhs_wave(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_symmetry_air(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_antisymmetry_air(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_symmetry_ferro(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_antisymmetry_ferro(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_cut(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_thinshell(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_layer_stiffness( Element * aElement,
                                     const uint aLayer,
                                     const Vector < real > & aHt,
                                     const real              aHn,
                                     Matrix< real > & aK );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_thinshell( Element * aElement );

//------------------------------------------------------------------------------

        };
//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            ( this->*mFunJacobian ) ( aElement, aJacobian, aRHS ) ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_superconductor(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_superconductor(
                    aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_ferro(
                    aElement, aJacobian, aRHS );
        }


//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->compute_jacobian_and_rhs_air_phi_linear(
                    aElement, aJacobian, aRHS );
        }


//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_scfm(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->compute_interface_ha_tri3( aElement, aJacobian, mGroup->work_K() );

            // compute RHS
            aRHS = aJacobian * this->collect_q0_ha_2d( aElement );

            // compone matrices
            aJacobian += mGroup->work_K() * this->timestep() ;
         }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_cut(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_cut( aElement, aJacobian, aRHS );
        }

//-----------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_fmfm(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->compute_interface_aa_tri3( aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IWG_MAXWELL_HPHI_TRI3_HPP
