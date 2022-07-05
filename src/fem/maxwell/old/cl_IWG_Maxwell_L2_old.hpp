//
// Created by Christian Messe on 03.02.22.
//

#ifndef BELFEM_CL_IWG_MAXWELL_L2_OLD_HPP
#define BELFEM_CL_IWG_MAXWELL_L2_OLD_HPP

#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        enum class Maxwell_L2_Mode
        {
            B2H,
            H2B,
            H2J
        };

        class IWG_Maxwell_L2 : public IWG_Maxwell
        {

            // link to function
            void
            ( IWG_Maxwell_L2::*mFunJacobian )
            (       Element        * aElement,
                    Matrix< real > & aJacobian );

            // link to function
            void
            ( IWG_Maxwell_L2::*mFunRHS )
            (   Element        * aElement,
                Vector< real > & aRHS );

            // alpha parameter for J = M ] alpha * K
            real mAlpha ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_L2( const ElementType aElementType,
                            const Maxwell_L2_Mode aMode ) ;

//------------------------------------------------------------------------------

            ~IWG_Maxwell_L2() ;

//------------------------------------------------------------------------------

            void
            link_jacobian_function( Group * aGroup );

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

            void
            set_block(
                    const id_t aBlockID,
                    const DomainType aDomainType );

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            allocate_work_matrices( Group    * aGroup );

//------------------------------------------------------------------------------

            void
            l2_node_tri3_scalar( Element * aElement,
                                Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_node_tri6_scalar( Element * aElement,
                                 Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_node_tet4_scalar( Element * aElement,
                                 Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_node_tet10_scalar( Element * aElement,
                                 Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_node_tri3_vector( Element * aElement,
                                 Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_node_tri6_vector( Element * aElement,
                                 Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_node_tet4_vector( Element * aElement,
                                 Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_node_tet10_vector( Element * aElement,
                                  Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_edge( Element * aElement,
                          Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            l2_b2h( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            l2_h2b( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            l2_h2jz( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            l2_h2j( Element * aElement, Vector< real > & aRHS );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            set_mode_tri3( const Maxwell_L2_Mode aMode );

            void
            set_mode_tri6( const Maxwell_L2_Mode aMode );

            void
            set_mode_tet4( const Maxwell_L2_Mode aMode );

            void
            set_mode_tet10( const Maxwell_L2_Mode aMode );

//------------------------------------------------------------------------------
            const Matrix< real > &
            J_tri_straight( Element * aElement );

//------------------------------------------------------------------------------

            const Matrix< real > &
            J_tet_straight( Element * aElement );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_L2::compute_jacobian(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            (this->*mFunJacobian)( aElement, aJacobian );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_L2::compute_rhs(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            (this->*mFunRHS)(aElement, aRHS );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IWG_MAXWELL_L2_OLD_HPP
