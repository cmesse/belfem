//
// Created by christian on 10/18/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HPHI2D_HPP
#define BELFEM_CL_IWG_MAXWELL_HPHI2D_HPP

#include "cl_IWG_Maxwell_Old.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_HPhi2D : public IWG_Maxwell_Old
        {

            void
            ( IWG_Maxwell_HPhi2D::*mFunJacobian )
                    (       Element        * aElement,
                            Matrix< real > & aJacobian,
                            Vector< real > & aRHS );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_HPhi2D();

//------------------------------------------------------------------------------

            ~IWG_Maxwell_HPhi2D() ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            shift_fields();

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            allocate_work_matrices( Group    * aGroup );

//------------------------------------------------------------------------------

            void
            link_jacobian_function( Group * aGroup ) ;

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
            compute_jacobian_and_rhs_interface_hphi(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_interface_aphi(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_interface_ha(
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

            void
            compute_jacobian_and_rhs_interface(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_farfield(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_symmetry(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_boundary(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );


//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            return ( this->*mFunJacobian )( aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_IWG_MAXWELL_HPHI2D_HPP
