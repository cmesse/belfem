//
// Created by christian on 12/2/21.
//

#ifndef BELFEM_CL_EF_TRI6_HPP
#define BELFEM_CL_EF_TRI6_HPP

#include "cl_EF_EdgeFunction.hpp"
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

            real mdNablaXidXi[ 2 ];
            real mdNablaXidEta[ 2 ];
            real mdNablaEtadXi[ 2 ];
            real mdNablaEtadEta[ 2 ];

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

            uint mLastJ = BELFEM_UINT_MAX ;
            uint mLastNabla = BELFEM_UINT_MAX ;
            uint mLastNablaDeriv = BELFEM_UINT_MAX ;

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
            ~EF_TRI6() = default;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi );

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element * aElement,
                  const bool aOperatorCurlH,
                  const bool aOperatorGradPhi,
                  const bool aOperatorCurlA );

//------------------------------------------------------------------------------

            /**
             * interpolation operator for h-field
             */
            const Matrix< real > &
            E( const uint aIndex );

//------------------------------------------------------------------------------

            /**
              * curl operator for h-field
              */
            const Matrix< real > &
            C( const uint aIndex = 0 );

//------------------------------------------------------------------------------

            /**
             * gradient operator for phi-field
             */
            const Matrix< real > &
            B( const uint aIndex = 0 ) ;

//------------------------------------------------------------------------------

            /**
              * curl operator for a-field
              */
            const Matrix< real > &
            CA( const uint aIndex = 0 ) ;

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
            B_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            B_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_J( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_nabla( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_nabla_derivatives( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_E( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_E_xi_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_E_xi_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_E_eta_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_E_eta_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_C();

//------------------------------------------------------------------------------


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

        inline void
        EF_TRI6::B_curved( const uint aIndex )
        {
            this->compute_J( aIndex );
            mB = mInvJ * mGroup->integration()->dNdXi( aIndex );
        }

//------------------------------------------------------------------------------

        inline void
        EF_TRI6::B_straight( const uint aIndex )
        {
            mB = mInvJ * mGroup->integration()->dNdXi( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6::B( const uint aIndex )
        {
           ( this->*mFunGrad )( aIndex );
            return mB ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6::CA( const uint aIndex )
        {
            ( this->*mFunGrad )( aIndex );

            mCurlA( 0, 0 ) =  mB( 1, 0 );
            mCurlA( 1, 0 ) = -mB( 0, 0 );
            mCurlA( 0, 1 ) =  mB( 1, 1 );
            mCurlA( 1, 1 ) = -mB( 0, 1 );
            mCurlA( 0, 2 ) =  mB( 1, 2 );
            mCurlA( 1, 2 ) = -mB( 0, 2 );
            mCurlA( 0, 3 ) =  mB( 1, 3 );
            mCurlA( 1, 3 ) = -mB( 0, 3 );
            mCurlA( 0, 4 ) =  mB( 1, 4 );
            mCurlA( 1, 4 ) = -mB( 0, 4 );
            mCurlA( 0, 5 ) =  mB( 1, 5 );
            mCurlA( 1, 5 ) = -mB( 0, 5 );

            return mCurlA ;
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_EF_TRI6_HPP
