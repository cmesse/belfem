//
// Created by christian on 8/11/22.
//

#ifndef BELFEM_CL_EF_TRI6_TS_HPP
#define BELFEM_CL_EF_TRI6_TS_HPP

#include "cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        /**
         * edge function for thin shell
         */
        class EF_TRI6_TS : public EdgeFunction
        {
            // polynomial for first dof on edge
            Vector< real > mEdgePoly0 ;

            // polynomial for second dof on edge
            Vector< real > mEdgePoly1 ;

            // needed for computation of Nablas for each edge
            Cell< Matrix< real > > mNxi0 ;
            Cell< Matrix< real > > mNxi1 ;
            Cell< Matrix< real > > mNxi2 ;

            // values for coordinates of connected element master
            Vector< real > mTau ;

            // coordinates of master
            Matrix< real > mX ;

            // Work Vector
            real mWork[ 10 ];

            // function for edge function
            const Matrix< real >  &
            ( EF_TRI6_TS :: * mComputeE )( const uint aIndex );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_TRI6_TS();

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~EF_TRI6_TS() = default;

//------------------------------------------------------------------------------

            /**
             * @param aXi
             */
            void
            precompute( const Matrix< real > & aXi ) ;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             * ( here, we don't care about the switches )
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

            /**
             * help function for precomputing N_xi
             */
             void
             compute_N_xi( const real aXi, const real aEta, Matrix< real > & aNxi );

//--------------------------------------------------------------------------

            const Matrix< real > &
            compute_E_edge0_curved( const uint aIndex );

//--------------------------------------------------------------------------

            const Matrix< real > &
            compute_E_edge1_curved( const uint aIndex );

//--------------------------------------------------------------------------

            const Matrix< real > &
            compute_E_edge2_curved( const uint aIndex );

//--------------------------------------------------------------------------

            const Matrix< real > &
            compute_E_edge_straight( const uint aIndex );

//--------------------------------------------------------------------------

            void
            compute_J_edge0( const uint aIndex );

//--------------------------------------------------------------------------

            void
            compute_J_edge1( const uint aIndex );

//--------------------------------------------------------------------------

            void
            compute_J_edge2( const uint aIndex );

//--------------------------------------------------------------------------

            void
            compute_r( const uint aIndex );

//--------------------------------------------------------------------------

            void
            compute_nabla();

//--------------------------------------------------------------------------

            void
            compute_edgepoly_edge0() ;

//--------------------------------------------------------------------------

            void
            compute_edgepoly_edge1() ;

//--------------------------------------------------------------------------

            void
            compute_edgepoly_edge2() ;

//--------------------------------------------------------------------------

        };

//--------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6_TS::E( const uint aIndex )
        {
            return ( this->*mComputeE )( aIndex );
        }

//--------------------------------------------------------------------------

        inline void
        EF_TRI6_TS::compute_J_edge0( const uint aIndex )
        {
            mJ = mNxi0( aIndex ) * mX ;
        }

//--------------------------------------------------------------------------

        inline void
        EF_TRI6_TS::compute_J_edge1( const uint aIndex )
        {
            mJ = mNxi1( aIndex ) * mX ;
        }

//--------------------------------------------------------------------------

        inline void
        EF_TRI6_TS::compute_J_edge2( const uint aIndex )
        {
            mJ = mNxi2( aIndex ) * mX ;
        }

//--------------------------------------------------------------------------

        inline void
        EF_TRI6_TS::compute_r( const uint aIndex )
        {
            mWork[ 4 ] = mWork[ 0 ] * mTau( aIndex ) + mWork[ 1 ];
            mWork[ 5 ] = mWork[ 2 ] * mTau( aIndex ) + mWork[ 3 ];
            mAbsDetJ = std::sqrt( mWork[ 4 ]*mWork[ 4 ] + mWork[ 5 ]*mWork[ 5 ] );
        }

//--------------------------------------------------------------------------


        inline void
        EF_TRI6_TS::compute_nabla()
        {
            mDetJ =   mJ( 0, 0 ) * mJ( 1, 1 )
                    - mJ( 0, 1 ) * mJ( 1, 0 );

            // k0 * mAbsDetJ * mDetJ
            mWork[ 6 ] =   mJ(1,1)*mWork[4]
                         - mJ(1,0)*mWork[5] ;

            // k1 * mAbsDetJ * mDetJ
            mWork[ 7 ] =mJ(0,0)*mWork[5]
                      - mJ(0,1)*mWork[4] ;
        }

//--------------------------------------------------------------------------

        inline void
        EF_TRI6_TS::compute_edgepoly_edge0()
        {
            mEdgePoly0( 0 ) =  mWork[ 6 ] + mWork[ 7 ] ;
            mEdgePoly0( 1 ) =  0.5 * mWork[ 6 ] - mWork[ 7 ];
            mEdgePoly0( 2 ) = -0.5 * mWork[ 6 ];

            mEdgePoly1( 0 ) = -mWork[ 6 ]-mWork[ 7 ];
            mEdgePoly1( 1 ) = 0.5*mWork[ 7 ]-mWork[ 6 ];
            mEdgePoly1( 2 ) = 0.5*mWork[ 7 ] ;
        }

//--------------------------------------------------------------------------

        inline void
        EF_TRI6_TS::compute_edgepoly_edge1()
        {
            mEdgePoly0( 0 ) = -mWork[ 6 ] ;
            mEdgePoly0( 1 ) =  mWork[ 6 ] + 1.5* mWork[ 7 ] ;
            mEdgePoly0( 2 ) = -0.5*mWork[ 7 ] ;

            mEdgePoly1( 0 ) =  mWork[ 6 ] ;
            mEdgePoly1( 1 ) =  0.5*mWork[ 6 ] + 1.5*mWork[ 7 ];
            mEdgePoly1( 2 ) = -0.5*(mWork[ 6 ]+mWork[ 7 ]);
        }

//--------------------------------------------------------------------------

        inline void
        EF_TRI6_TS::compute_edgepoly_edge2()
        {
            mEdgePoly0( 0 ) = -mWork[ 7 ] ;
            mEdgePoly0( 1 ) = -1.5*mWork[ 6 ] - 0.5*mWork[ 7 ] ;
            mEdgePoly0( 2 ) =  0.5*(mWork[ 6 ]+mWork[ 7 ]) ;
            mEdgePoly1( 0 ) = mWork[ 7 ] ;
            mEdgePoly1( 1 ) = 1.5*mWork[ 6 ] + mWork[ 7 ] ;
            mEdgePoly1( 2 ) = 0.5*mWork[ 6 ] ;
        }

//--------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6_TS::compute_E_edge0_curved( const uint aIndex )
        {
            this->compute_J_edge0( aIndex );
            this->compute_r( aIndex );
            this->compute_nabla();
            this->compute_edgepoly_edge0() ;

            mE( 0, 0 ) = polyval( mEdgePoly0, mTau( aIndex ) );
            mE( 0, 1 ) = polyval( mEdgePoly1, mTau( aIndex ) );
            mE /= mDetJ * mAbsDetJ ;
            return mE ;
        }

//--------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6_TS::compute_E_edge1_curved( const uint aIndex )
        {
            this->compute_J_edge1( aIndex );
            this->compute_r( aIndex );
            this->compute_nabla();
            this->compute_edgepoly_edge1() ;

            mE( 0, 0 ) = polyval( mEdgePoly0, mTau( aIndex ) );
            mE( 0, 1 ) = polyval( mEdgePoly1, mTau( aIndex ) );
            mE /= mDetJ * mAbsDetJ ;
            return mE ;
        }

//--------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6_TS::compute_E_edge2_curved( const uint aIndex )
        {
            this->compute_J_edge2( aIndex );
            this->compute_r( aIndex );
            this->compute_nabla();
            this->compute_edgepoly_edge2() ;

            mE( 0, 0 ) = polyval( mEdgePoly0, mTau( aIndex ) );
            mE( 0, 1 ) = polyval( mEdgePoly1, mTau( aIndex ) );
            mE /= mDetJ * mAbsDetJ ;
            return mE ;
        }

//--------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TRI6_TS::compute_E_edge_straight( const uint aIndex )
        {
            mE( 0, 0 ) = polyval( mEdgePoly0, mTau( aIndex ) );
            mE( 0, 1 ) = polyval( mEdgePoly1, mTau( aIndex ) );
            return mE ;
        }

//--------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EF_TRI6_TS_HPP
