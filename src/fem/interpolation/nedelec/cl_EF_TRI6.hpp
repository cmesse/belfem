//
// Created by christian on 12/2/21.
//

#ifndef BELFEM_CL_EF_TRI6_HPP
#define BELFEM_CL_EF_TRI6_HPP

#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        class EF_TRI6 : public EdgeFunction
        {
            // node coordinates in matrix
            Matrix< real > mNodeCoords ;

            // node coordinates in vectors
            Vector< real > mX;
            Vector< real > mY;

            real mS[ 3 ];

            real mNablaXi[ 2 ];
            real mNablaEta[ 2 ];

            //! derivatives of node functions
            Matrix< real > mNxi ;
            Matrix< real > mNeta ;

            // help coefficients for shape function
            Matrix< real > mG ;
            Matrix< real > mGxi ;
            Matrix< real > mGeta ;

            Matrix< real > mH ;
            Matrix< real > mHxi ;
            Matrix< real > mHeta ;

            real mExi[ 16 ] ;
            real mEeta[ 16 ] ;

            // additional work vector
            real mW[ 16 ];

            // contracted second derivatives of the quadratic geometry map,
            // S^(pq)_m ( IF-Voigt rows xx, yy, xy ), constant per element;
            // filled on the CURVED branch of link(), NaN-filled on the
            // straight branch ( the hess channel of G must never read it
            // there — the affine map has S = 0 and the channel is skipped )
            real mCurv[ 3 ][ 2 ];

            uint mLastJ = BELFEM_UINT_MAX ;
            uint mLastNabla = BELFEM_UINT_MAX ;
            void
            ( EF_TRI6::*mFunInterpolation )( const uint aIndex );

            void
            ( EF_TRI6::*mFunCurl )( const uint aIndex );

            void
            ( EF_TRI6::*mFunGrad )( const uint aIndex );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_TRI6();

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~EF_TRI6() override = default;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi ) override;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element * aElement ) override;

//------------------------------------------------------------------------------

            /**
             * interpolation operator for h-field
             */
            const Matrix< real > &
            E( const uint aIndex ) override;

//------------------------------------------------------------------------------

            /**
              * curl operator for h-field
              */
            const Matrix< real > &
            C( const uint aIndex = 0 ) override;

//------------------------------------------------------------------------------

            /**
              * gradient operator for h-field
              * ( straight path: signed term1 with frozen nablas;
              *   curved path: + the Q-folded inverse-map Hessian channel;
              *   layout contract in EdgeFunction::mGrad )
              */
            const Matrix< real > &
            G( const uint aIndex = 0 ) override ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            E_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            C_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            C_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_J( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_nabla( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_E( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_E_xi( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_E_eta( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_C();

//------------------------------------------------------------------------------

            /**
             * signed term1 of the gradient: reference partials chained with
             * the current nablas ( the partial tables of this class are
             * UNSIGNED — the signs are applied here, the compute_C way )
             */
            void
            G_term1( const uint aIndex );

//------------------------------------------------------------------------------

            void
            G_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            G_curved( const uint aIndex );

//------------------------------------------------------------------------------
        };


//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6::E( const uint aIndex )
        {
            ( this->*mFunInterpolation )( aIndex );
            return mE ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6::C( const uint aIndex )
        {
            ( this->*mFunCurl )( aIndex );
            return mC ;
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_EF_TRI6_HPP
