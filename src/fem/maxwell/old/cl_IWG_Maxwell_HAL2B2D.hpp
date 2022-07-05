//
// Created by christian on 10/14/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HAL2B2D_HPP
#define BELFEM_CL_IWG_MAXWELL_HAL2B2D_HPP

#include "cl_IWG_Maxwell_Old.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_HAL2B2D : public IWG_Maxwell_Old
        {
            void
            ( IWG_Maxwell_HAL2B2D::*mFunRHS )
                    ( Element * aElement,
                      Vector< real > & aRHS );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_HAL2B2D();

//------------------------------------------------------------------------------

            ~IWG_Maxwell_HAL2B2D() = default ;

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
        };

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HAL2B2D ::compute_rhs(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            return ( this->*mFunRHS )( aElement, aRHS );
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_IWG_MAXWELL_HAL2B2D_HPP
