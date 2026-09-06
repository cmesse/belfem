//
// Created by christian on 12/1/21.
//

#ifndef BELFEM_CL_EF_TET10_HPP
#define BELFEM_CL_EF_TET10_HPP

#include "nedelec/cl_EF_EdgeFunction.hpp"

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

            // matrix for second derivatives
            // N_xi2, N_eta2, N_zeta2, N_etazeta, N_xizeta, N_xieta
            Matrix< real > mD ;

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

            // contracted second derivatives of the quadratic geometry map,
            // S^(pq)_m ( IF-Voigt rows xx, yy, zz, yz, xz, xy — the mD
            // member IS the d2N table, reused ), constant per element;
            // filled on the CURVED branch of link(), NaN-filled on the
            // straight branch ( the hess channel of G must never read it
            // there — the affine map has S = 0 and the channel is skipped )
            real mCurv[ 6 ][ 3 ];

            uint mLastJ = BELFEM_UINT_MAX ;
            uint mLastNabla = BELFEM_UINT_MAX ;

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
            ~EF_TET10() override = default;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element * aElement ) override ;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi ) override ;

//------------------------------------------------------------------------------

            // compute the edge function
            const Matrix< real > &
            E( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            // compute the curl function
            const Matrix< real > &
            C( const uint aIndex = 0 ) override ;

//------------------------------------------------------------------------------

            /**
              * gradient operator for h-field
              * ( straight path: term1 from the pre-signed combined
              *   partials with frozen nablas; curved path: + the Q-folded
              *   inverse-map Hessian channel;
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
            E_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            E_xi_curved( const uint aIndex );

//------------------------------------------------------------------------------

            void
            E_xi_straight( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_nabla( const uint aIndex ) ;

//------------------------------------------------------------------------------

            void
            compute_edge_functions( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_face_functions( const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * term1 of the gradient: reference partials chained with the
             * current nablas, in C()'s exact derivative + combine sequence
             * ( the tables of this class are pre-signed and
             * orientation-combined )
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
            compute_edge_derivatives( const uint aIndex );

//------------------------------------------------------------------------------

            void
            compute_face_derivatives( const uint aIndex );

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

    } /* end namespace fem */
}  /* end namespace belfem */
#endif //BELFEM_CL_EF_TET10_HPP