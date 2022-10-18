//
// Created by christian on 3/22/22.
//

#ifndef BELFEM_CL_IWG_MAXWELL_L2_MAGFIELD_HPP
#define BELFEM_CL_IWG_MAXWELL_L2_MAGFIELD_HPP


#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        enum class Magfield_L2_Mode
        {
            HA   = 0,
            HPHI = 1,
            UNDEFINED
        };

        class IWG_Maxwell_L2_Magfield : public IWG_Maxwell
        {
            // alpha parameter for J = M + alpha * K
            const real mAlpha ;

            // link to function
            void
            ( IWG_Maxwell_L2_Magfield::*mFunLink ) ( Group * aGroup );


            // link to function
            void
            ( IWG_Maxwell_L2_Magfield::*mFunJacobian )
                    (       Element        * aElement,
                            Matrix< real > & aJacobian );

            // link to function
            void
            ( IWG_Maxwell_L2_Magfield::*mFunRHS )
                    (       Element        * aElement,
                            Vector< real > & aRHS );

            Vector< real > mH ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_L2_Magfield(
                    const ElementType aElementType,
                    const Magfield_L2_Mode aMode,
                    const real            aAlpha = 0.0 );

//------------------------------------------------------------------------------

            ~IWG_Maxwell_L2_Magfield();

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
            link_jacobian_2d_hphi( Group * aGroup );

//------------------------------------------------------------------------------

            void
            link_jacobian_2d_ha( Group * aGroup );

//------------------------------------------------------------------------------

            void
            compute_jacobian_2d_alpha0( Element * aElement, Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

            void
            compute_jacobian_2d_alpha1( Element * aElement, Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

            void
            compute_jacobian_cut( Element * aElement,
                                  Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

            void
            compute_rhs_2d_h( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_rhs_2d_phi( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_rhs_2d_phi_ferro( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_rhs_2d_a( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_rhs_cut( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------
        } ;
//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_L2_Magfield::link_jacobian_function( Group * aGroup )
        {
            ( this->*mFunLink )( aGroup );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_L2_Magfield::compute_jacobian(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            ( this->*mFunJacobian )( aElement, aJacobian );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_L2_Magfield::compute_rhs(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            ( this->*mFunRHS )( aElement, aRHS );
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IWG_MAXWELL_L2_MAGFIELD_HPP
