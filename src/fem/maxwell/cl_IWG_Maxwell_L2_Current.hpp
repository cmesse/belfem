//
// Created by christian on 3/18/22.
//

#ifndef BELFEM_CL_IWG_MAXWELL_L2_CURRENT_HPP
#define BELFEM_CL_IWG_MAXWELL_L2_CURRENT_HPP

#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_L2_Current : public IWG_Maxwell
        {
            // alpha parameter for J = M + alpha * K
            const real mAlpha ;

            // link to function
            void
            ( IWG_Maxwell_L2_Current::*mFunLink ) ( Group * aGroup );


            // link to function
            void
            ( IWG_Maxwell_L2_Current::*mFunJacobian )
                    (       Element        * aElement,
                            Matrix< real > & aJacobian );

            // link to function
            void
            ( IWG_Maxwell_L2_Current::*mFunRHS )
                    (       Element        * aElement,
                            Vector< real > & aRHS );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_L2_Current(
                    const ElementType aElementType,
                    const real            aAlpha = 0.0 );

//------------------------------------------------------------------------------

            ~IWG_Maxwell_L2_Current();

//------------------------------------------------------------------------------

            void
            compute_jacobian(
                    Element        * aElement,
                    Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

            void
            compute_rhs(
                    Element        * aElement,
                    Vector< real > & aRHS );


//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            link_jacobian_function( Group * aGroup );

//------------------------------------------------------------------------------

            void
            allocate_work_matrices( Group * aGroup );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            link_jacobian_2d( Group * aGroup );

//------------------------------------------------------------------------------

            void
            jacobian_2d_alpha0( Element * aElement, Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            jacobian_2d_alpha1( Element * aElement, Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            rhs_2d( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------
        } ;
//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_L2_Current::link_jacobian_function( Group * aGroup )
        {
            ( this->*mFunLink )( aGroup );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_L2_Current::compute_jacobian(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            ( this->*mFunJacobian )( aElement, aJacobian );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_L2_Current::compute_rhs(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            ( this->*mFunRHS )( aElement, aRHS );
        }

//------------------------------------------------------------------------------
    }
}


#endif //BELFEM_CL_IWG_MAXWELL_L2_CURRENT_HPP
