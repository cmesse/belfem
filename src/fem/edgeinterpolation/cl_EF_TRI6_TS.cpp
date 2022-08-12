//
// Created by christian on 8/11/22.
//

#include "cl_EF_TRI6_TS.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        EF_TRI6_TS::EF_TRI6_TS()
        {
            mE.set_size( 1, 2 );
            mX.set_size(  6, 2 );
            mJ.set_size( 2, 2 );
            mEdgePoly0.set_size( 3 );
            mEdgePoly1.set_size( 3 );
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6_TS::precompute( const Matrix< real > & aXi )
        {
            // number of integration points
            uint tNumIntPoints = aXi.n_cols() ;

            mTau.set_size( tNumIntPoints );
            for( uint k = 0; k<tNumIntPoints; ++k )
            {
                mTau( k ) = aXi( 0, k );
            }

            // allocate memory
            mNxi0.set_size( tNumIntPoints, Matrix< real >( 2, 6 ) );
            mNxi1.set_size( tNumIntPoints, Matrix< real >( 2, 6 ) );
            mNxi2.set_size( tNumIntPoints, Matrix< real >( 2, 6 ) );

            // compute values for edge 0
            for( uint k = 0; k<tNumIntPoints; ++k )
            {
                real tEta = mTau( k ) ;
                real tXi  = - tEta ;

                this->compute_N_xi( tXi, tEta, mNxi0( k ) );
            }

            // compute values for edge 1
            for( uint k = 0; k<tNumIntPoints; ++k )
            {
                real tZeta = mTau( k ) ;
                real tEta  = - tZeta ;
                real tXi   = -1. - tEta - tZeta ;
                this->compute_N_xi( tXi, tEta, mNxi1( k ) );
            }

            // compute values for edge 2
            for( uint k = 0; k<tNumIntPoints; ++k )
            {
                real tXi  = mTau( k ) ;
                real tZeta = - tXi ;
                real tEta = -1. - tXi - tZeta ;
                this->compute_N_xi( tXi, tEta, mNxi2( k ) );
            }

        }

//------------------------------------------------------------------------------

        void
        EF_TRI6_TS::link(
                Element * aElement,
                const bool aOperatorCurlH,
                const bool aOperatorGradPhi,
                const bool aOperatorCurlA )
        {
            // grab facet master
            mesh::Element * tElement = aElement->facet()->master() ;

            BELFEM_ASSERT( tElement->type() == ElementType::TRI6 ,
                "Can only link to TRI6 element" );

            // collect node coords
            mX( 0, 0 ) = tElement->node( 0 )->x() ;
            mX( 1, 0 ) = tElement->node( 1 )->x() ;
            mX( 2, 0 ) = tElement->node( 2 )->x() ;
            mX( 3, 0 ) = tElement->node( 3 )->x() ;
            mX( 4, 0 ) = tElement->node( 4 )->x() ;
            mX( 5, 0 ) = tElement->node( 5 )->x() ;
            mX( 0, 1 ) = tElement->node( 0 )->y() ;
            mX( 1, 1 ) = tElement->node( 1 )->y() ;
            mX( 2, 1 ) = tElement->node( 2 )->y() ;
            mX( 3, 1 ) = tElement->node( 3 )->y() ;
            mX( 4, 1 ) = tElement->node( 4 )->y() ;
            mX( 5, 1 ) = tElement->node( 5 )->y() ;

            //if( aElement->master()->element()->is_curved() )
            if( true )
            {
                switch ( aElement->facet()->master_index())
                {
                    case ( 0 ) :
                    {
                        // link jacobian function
                        mComputeE = & EF_TRI6_TS::compute_E_edge0_curved;

                        // values for direction vector
                        mWork[ 0 ] = mX( 0, 0 ) + mX( 1, 0 ) - mX( 3, 0 ) - mX( 3, 0 );
                        mWork[ 1 ] = 0.5 * ( mX( 1, 0 ) - mX( 0, 0 ));
                        mWork[ 2 ] = mX( 0, 1 ) + mX( 1, 1 ) - mX( 3, 1 ) - mX( 3, 1 );
                        mWork[ 3 ] = 0.5 * ( mX( 1, 1 ) - mX( 0, 1 ));

                        break;
                    }
                    case ( 1 ) :
                    {
                        mComputeE = &EF_TRI6_TS::compute_E_edge1_curved;

                        // values for direction vector
                        mWork[ 0 ] = mX( 1, 0 ) + mX( 2, 0 ) - mX( 4, 0 ) - mX( 4, 0 );
                        mWork[ 1 ] = 0.5 * ( mX( 2, 0 ) - mX( 1, 0 ));
                        mWork[ 2 ] = mX( 1, 1 ) + mX( 2, 1 ) - mX( 4, 1 ) - mX( 4, 1 );
                        mWork[ 3 ] = 0.5 * ( mX( 2, 1 ) - mX( 1, 1 ));

                        break;
                    }
                    case ( 2 ) :
                    {
                        mComputeE = &EF_TRI6_TS::compute_E_edge2_curved;

                        // values for direction vector
                        mWork[ 0 ] = mX( 2, 0 ) + mX( 0, 0 ) - mX( 5, 0 ) - mX( 5, 0 );
                        mWork[ 1 ] = 0.5 * ( mX( 0, 0 ) - mX( 2, 0 ));
                        mWork[ 2 ] = mX( 2, 1 ) + mX( 0, 1 ) - mX( 5, 1 ) - mX( 5, 1 );
                        mWork[ 3 ] = 0.5 * ( mX( 0, 1 ) - mX( 2, 1 ));
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid Index on master" );
                    }
                }
            }
            else
            {
                // compute jacobian
                mJ( 0, 0 ) = mX( 0, 0 ) - mX( 2, 0 );
                mJ( 1, 0 ) = mX( 1, 0 ) - mX( 2, 0 );
                mJ( 0, 1 ) = mX( 0, 1 ) - mX( 2, 1 );
                mJ( 1, 1 ) = mX( 1, 1 ) - mX( 2, 1 );
                mJ *= 0.5 ;
                mDetJ = mJ( 0, 0 ) * mJ( 1, 1 ) - mJ( 0, 1 ) * mJ( 1, 0 );

                switch ( aElement->facet()->master_index())
                {
                    case ( 0 ) :
                    {
                        // link jacobian function
                        mComputeE = &EF_TRI6_TS::compute_E_edge_straight ;

                        // values for direction vector
                        mWork[ 4 ] = 0.5 * ( mX( 1, 0 ) - mX( 0, 0 ));
                        mWork[ 5 ] = 0.5 * ( mX( 1, 1 ) - mX( 0, 1 ));
                        mAbsDetJ = std::sqrt( mWork[ 4 ]*mWork[ 4 ] + mWork[ 5 ]*mWork[ 5 ] );

                        // values for constants k0 and k1
                        this->compute_nabla();
                        mWork[ 6 ] /= mDetJ * mAbsDetJ ;
                        mWork[ 7 ] /= mDetJ * mAbsDetJ ;
                        this->compute_edgepoly_edge0() ;

                        break;
                    }
                    case ( 1 ) :
                    {
                        mComputeE = &EF_TRI6_TS::compute_E_edge_straight;

                        // values for direction vector
                        mWork[ 4 ] = 0.5 * ( mX( 2, 0 ) - mX( 1, 0 ));
                        mWork[ 5 ] = 0.5 * ( mX( 2, 1 ) - mX( 1, 1 ));
                        mAbsDetJ = std::sqrt( mWork[ 4 ]*mWork[ 4 ] + mWork[ 5 ]*mWork[ 5 ] );
                        this->compute_nabla();
                        mWork[ 6 ] /= mDetJ * mAbsDetJ ;
                        mWork[ 7 ] /= mDetJ * mAbsDetJ ;
                        this->compute_edgepoly_edge1() ;
                        break;
                    }
                    case ( 2 ) :
                    {
                        mComputeE = &EF_TRI6_TS::compute_E_edge_straight;

                        // values for direction vector
                        mWork[ 4 ] = 0.5 * ( mX( 0, 0 ) - mX( 2, 0 ));
                        mWork[ 5 ] = 0.5 * ( mX( 0, 1 ) - mX( 2, 1 ));
                        mAbsDetJ = std::sqrt( mWork[ 4 ]*mWork[ 4 ] + mWork[ 5 ]*mWork[ 5 ] );
                        this->compute_nabla();
                        mWork[ 6 ] /= mDetJ * mAbsDetJ ;
                        mWork[ 7 ] /= mDetJ * mAbsDetJ ;
                        this->compute_edgepoly_edge2() ;
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid Index on master" );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        EF_TRI6_TS::compute_N_xi( const real aXi, const real aEta, Matrix< real > & aNxi )
        {
            aNxi( 0, 0 ) =  aXi + 0.5 ;
            aNxi( 1, 0 ) =  0.0 ;

            aNxi( 0, 1 ) =  0.0 ;
            aNxi( 1, 1 ) =  aEta + 0.5 ;

            aNxi( 0, 2 ) =  aEta + aXi + 0.5 ;
            aNxi( 1, 2 ) =  aEta + aXi + 0.5 ;

            aNxi( 0, 3 ) =  aEta + 1.0 ;
            aNxi( 1, 3 ) =  aXi + 1.0 ;

            aNxi( 0, 4 ) = -aEta - 1.0 ;
            aNxi( 1, 4 ) = -aEta - aEta - aXi - 1.0 ;

            aNxi( 0, 5 ) = -aEta - aXi - aXi - 1.0 ;
            aNxi( 1, 5 ) = -aXi - 1.0 ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_TRI6_TS::C( const uint aIndex )
        {
            BELFEM_ERROR( false, "function not implemented");
            return mC ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_TRI6_TS::B( const uint aIndex )
        {
            BELFEM_ERROR( false, "function not implemented");
            return mB ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        EF_TRI6_TS::CA( const uint aIndex )
        {
            BELFEM_ERROR( false, "function not implemented");
            return mC ;
        }

//------------------------------------------------------------------------------
    }
}