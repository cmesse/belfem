//
// Created by christian on 10/27/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HPHIL2B2D_HPP
#define BELFEM_CL_IWG_MAXWELL_HPHIL2B2D_HPP

#include "cl_IWG_Maxwell_Old.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_HPhiL2B2D : public IWG_Maxwell_Old
        {
            void
            ( IWG_Maxwell_HPhiL2B2D::*mFunJacobian )
                    ( Element * aElement,
                      Matrix< real > & aJacobian );

            void
            ( IWG_Maxwell_HPhiL2B2D::*mFunRHS )
                    ( Element * aElement,
                      Vector< real > & aRHS );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_HPhiL2B2D();

//------------------------------------------------------------------------------

            ~IWG_Maxwell_HPhiL2B2D() = default ;

//------------------------------------------------------------------------------

            void
            compute_jacobian(
                    Element        * aElement,
                    Matrix< real > & aJacobian ) ;

//------------------------------------------------------------------------------

            void
            compute_rhs(
                    Element        * aElement,
                    Vector< real > & aRHS ) ;

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
            compute_rhs_superconductor(
                    Element        * aElement,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_rhs_air(
                    Element        * aElement,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_rhs_ferro(
                    Element        * aElement,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_rhs_cut(
                    Element        * aElement,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_block (
                    Element        * aElement,
                    Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

            void
            compute_jacobian_cut (
                    Element        * aElement,
                    Matrix< real > & aJacobian );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhiL2B2D ::compute_jacobian(
                Element        * aElement,
                Matrix< real > & aJacobian  )
        {
            return ( this->*mFunJacobian )( aElement, aJacobian );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhiL2B2D ::compute_rhs(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            return ( this->*mFunRHS )( aElement, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhiL2B2D::compute_jacobian_block(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            this->projection_jacobian_b2d( aElement, aJacobian );
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_IWG_MAXWELL_HPHIL2B2D_HPP
