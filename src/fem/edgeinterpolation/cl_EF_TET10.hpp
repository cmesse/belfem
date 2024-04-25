//
// Created by christian on 12/1/21.
//

#ifndef BELFEM_CL_EF_TET10_HPP
#define BELFEM_CL_EF_TET10_HPP

#include "cl_EF_EdgeFunction.hpp"

namespace belfem
{
    namespace mesh
    {
        class Element ;
    }

    namespace fem
    {
        class Element ;

        class EF_TET10 : public EdgeFunction
        {
            //! node coordinates as Matrix
            Matrix< real > mNodeCoords ;

            //! node coordinates as vector
            Vector< real > mX ;
            Vector< real > mY ;
            Vector< real > mZ ;

            //! derivatives of node functions
            Matrix< real > mNxi ;
            Matrix< real > mNeta ;
            Matrix< real > mNzeta ;

            //! edge directions
            real mS[ 6 ];

            //! face orientations
            uint mT[ 4 ];

            //! values for the face function
            Matrix< real > mF ;

            //! precomputed values for edges
            Matrix< real > mG ;
            Matrix< real > mH ;

            //! precomputed values for faces
            Matrix< real > mU ;
            Matrix< real > mV ;
            Matrix< real > mW ;


            //! nabla values
            real mNablaXi[ 3 ];
            real mNablaEta[ 3 ];
            real mNablaZeta[ 3 ];
            real mNablaTau[ 3 ];

            real mdNablaXidXi[ 3 ];
            real mdNablaXidEta[ 3 ];
            real mdNablaXidZeta[ 3 ];

            real mdNablaEtadXi[ 3 ];
            real mdNablaEtadEta[ 3 ];
            real mdNablaEtadZeta[ 3 ];

            real mdNablaZetadXi[ 3 ];
            real mdNablaZetadEta[ 3 ];
            real mdNablaZetadZeta[ 3 ];

            real mdNablaTaudXi[ 3 ];
            real mdNablaTaudEta[ 3 ];
            real mdNablaTaudZeta[ 3 ];

            // edge shape functions
            Matrix< real > mExi ;
            Matrix< real > mEeta ;
            Matrix< real > mEzeta ;

            // face shape functions
            Matrix< real > mFxi ;
            Matrix< real > mFeta ;
            Matrix< real > mFzeta ;

            // containers with precomputed factors for edges
            Matrix< real > mGxi ;
            Matrix< real > mGeta ;
            Matrix< real > mGzeta ;

            // containers with precomputed factors for edges
            Matrix< real > mHxi ;
            Matrix< real > mHeta ;
            Matrix< real > mHzeta ;

            // containers with precomputed factors for faces
            Matrix< real > mUxi ;
            Matrix< real > mUeta ;
            Matrix< real > mUzeta ;

            Matrix< real > mVxi ;
            Matrix< real > mVeta ;
            Matrix< real > mVzeta ;

            Matrix< real > mWxi ;
            Matrix< real > mWeta ;
            Matrix< real > mWzeta ;

            // additional work vector
            real mM [ 76 ];

            uint mLastJ = BELFEM_UINT_MAX ;
            uint mLastNabla = BELFEM_UINT_MAX ;
            uint mLastNablaDeriv = BELFEM_UINT_MAX ;

            void
            ( EF_TET10::*mFunInterpolation )( const uint aIndex );

            void
            ( EF_TET10::*mFunDerivatives )( const uint aIndex );

            void
            ( EF_TET10::*mFunGrad )( const uint aIndex );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_TET10();

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~EF_TET10() = default;

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

            void
            precompute( const Matrix< real > & aXi );

//------------------------------------------------------------------------------

            // compute the edge function
            const Matrix< real > &
            E( const uint aIndex );

//------------------------------------------------------------------------------

            // compute the curl function
            const Matrix< real > &
            C( const uint aIndex = 0 );

//------------------------------------------------------------------------------

            // compute the B-Operator
            const Matrix< real > &
            B( const uint aIndex = 0 ) ;

//------------------------------------------------------------------------------

            // compute the Curl operator for the a-field
            const Matrix< real > &
            CA( const uint aIndex = 0 ) ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            E_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            E_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            E_xi_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            E_xi_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            C_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            B_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            B_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_jacobian( const uint aIndex ) ;

//------------------------------------------------------------------------------

            void
            compute_nabla( const uint aIndex ) ;

//------------------------------------------------------------------------------

            void
            compute_nabla_derivatives( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_edge_functions( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_face_functions( const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * Initially E only contains edge based functions.
             * The face functions depend on the face orientation and are
             * properly arranged in this additional step
             *
             * @param aE edge functions + space for face functions
             * @param aF computed face functions
             */
            void
            combine_functions( Matrix< real > & aE, const Matrix< real > & aF );

//------------------------------------------------------------------------------

            void
            compute_edge_derivatives_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_edge_derivatives_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_face_derivatives_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_face_derivatives_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            print_memory();

//------------------------------------------------------------------------------
        } ;

//------------------------------------------------------------------------------

        // compute the edge function
        inline const Matrix< real > &
        EF_TET10::E( const uint aIndex )
        {
            // compute the matrix
            ( this->*mFunInterpolation )( aIndex );

            // return the result
            return mE ;
        }

//------------------------------------------------------------------------------

        // compute the B-function
        inline const Matrix< real > &
        EF_TET10::B( const uint aIndex )
        {
            // compute the B-matrix
            ( this->*mFunGrad )( aIndex );

            // return the result
            return mB ;
        }

//------------------------------------------------------------------------------

        inline void
        EF_TET10::B_curved( const uint aIndex )
        {
            this->compute_jacobian( aIndex );
            mB = mInvJ * mGroup->integration()->dNdXi( aIndex );
        }

//------------------------------------------------------------------------------

        inline void
        EF_TET10::B_straight( const uint aIndex )
        {
            mB = mInvJ * mGroup->integration()->dNdXi( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TET10::CA( const uint aIndex )
        {
            // compute the B-Matrix
            ( this->*mFunGrad )( aIndex );

            // expand the B-Matrix
            mCurlA( 1,  0 ) =  mB( 2, 0 ) ;
            mCurlA( 2,  0 ) = -mB( 1, 0 ) ;
            mCurlA( 0,  1 ) = -mB( 2, 0 ) ;
            mCurlA( 2,  1 ) =  mB( 0, 0 ) ;
            mCurlA( 0,  2 ) =  mB( 1, 0 ) ;
            mCurlA( 1,  2 ) = -mB( 0, 0 ) ;
            mCurlA( 1,  3 ) =  mB( 2, 1 ) ;
            mCurlA( 2,  3 ) = -mB( 1, 1 ) ;
            mCurlA( 0,  4 ) = -mB( 2, 1 ) ;
            mCurlA( 2,  4 ) =  mB( 0, 1 ) ;
            mCurlA( 0,  5 ) =  mB( 1, 1 ) ;
            mCurlA( 1,  5 ) = -mB( 0, 1 ) ;
            mCurlA( 1,  6 ) =  mB( 2, 2 ) ;
            mCurlA( 2,  6 ) = -mB( 1, 2 ) ;
            mCurlA( 0,  7 ) = -mB( 2, 2 ) ;
            mCurlA( 2,  7 ) =  mB( 0, 2 ) ;
            mCurlA( 0,  8 ) =  mB( 1, 2 ) ;
            mCurlA( 1,  8 ) = -mB( 0, 2 ) ;
            mCurlA( 1,  9 ) =  mB( 2, 3 ) ;
            mCurlA( 2,  9 ) = -mB( 1, 3 ) ;
            mCurlA( 0, 10 ) = -mB( 2, 3 ) ;
            mCurlA( 2, 10 ) =  mB( 0, 3 ) ;
            mCurlA( 0, 11 ) =  mB( 1, 3 ) ;
            mCurlA( 1, 11 ) = -mB( 0, 3 ) ;
            mCurlA( 1, 12 ) =  mB( 2, 4 ) ;
            mCurlA( 2, 12 ) = -mB( 1, 4 ) ;
            mCurlA( 0, 13 ) = -mB( 2, 4 ) ;
            mCurlA( 2, 13 ) =  mB( 0, 4 ) ;
            mCurlA( 0, 14 ) =  mB( 1, 4 ) ;
            mCurlA( 1, 14 ) = -mB( 0, 4 ) ;
            mCurlA( 1, 15 ) =  mB( 2, 5 ) ;
            mCurlA( 2, 15 ) = -mB( 1, 5 ) ;
            mCurlA( 0, 16 ) = -mB( 2, 5 ) ;
            mCurlA( 2, 16 ) =  mB( 0, 5 ) ;
            mCurlA( 0, 17 ) =  mB( 1, 5 ) ;
            mCurlA( 1, 17 ) = -mB( 0, 5 ) ;
            mCurlA( 1, 18 ) =  mB( 2, 6 ) ;
            mCurlA( 2, 18 ) = -mB( 1, 6 ) ;
            mCurlA( 0, 19 ) = -mB( 2, 6 ) ;
            mCurlA( 2, 19 ) =  mB( 0, 6 ) ;
            mCurlA( 0, 20 ) =  mB( 1, 6 ) ;
            mCurlA( 1, 20 ) = -mB( 0, 6 ) ;
            mCurlA( 1, 21 ) =  mB( 2, 7 ) ;
            mCurlA( 2, 21 ) = -mB( 1, 7 ) ;
            mCurlA( 0, 22 ) = -mB( 2, 7 ) ;
            mCurlA( 2, 22 ) =  mB( 0, 7 ) ;
            mCurlA( 0, 23 ) =  mB( 1, 7 ) ;
            mCurlA( 1, 23 ) = -mB( 0, 7 ) ;
            mCurlA( 1, 24 ) =  mB( 2, 8 ) ;
            mCurlA( 2, 24 ) = -mB( 1, 8 ) ;
            mCurlA( 0, 25 ) = -mB( 2, 8 ) ;
            mCurlA( 2, 25 ) =  mB( 0, 8 ) ;
            mCurlA( 0, 26 ) =  mB( 1, 8 ) ;
            mCurlA( 1, 26 ) = -mB( 0, 8 ) ;
            mCurlA( 1, 27 ) =  mB( 2, 9 ) ;
            mCurlA( 2, 27 ) = -mB( 1, 9 ) ;
            mCurlA( 0, 28 ) = -mB( 2, 9 ) ;
            mCurlA( 2, 28 ) =  mB( 0, 9 ) ;
            mCurlA( 0, 29 ) =  mB( 1, 9 ) ;
            mCurlA( 1, 29 ) = -mB( 0, 9 ) ;

            return mCurlA ;
        }

//------------------------------------------------------------------------------

    } /* end namespace fem */
}  /* end namespace belfem */
#endif //BELFEM_CL_EF_TET10_HPP