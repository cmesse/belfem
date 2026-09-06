/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_IF_ElementMapper.hpp"
#include "fn_inv2.hpp"
#include "fn_inv3.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
namespace belfem
{
    namespace fem
    {
        ElementMapper::ElementMapper()
        {
            mFactory = new InterpolationFunctionFactory() ;
        }

        ElementMapper::~ElementMapper()
        {
            for ( auto it : mFunctions )
            {
                delete it.second ;
            }
            delete mFactory ;
        }

        void
        ElementMapper::set_dimension( const uint aDim )
        {
            mDim = aDim ;
            mJ.set_size( aDim, aDim );
            mInvJ.set_size( aDim, aDim );
        }

        void
        ElementMapper::link( mesh::Element * aElement )
        {
            mElement = aElement ;

            uint n = aElement->number_of_nodes() ;
            uint tNatDim = aElement->dimension() ;

            if ( mDim == 0 )
            {
                this->set_dimension( tNatDim );
            }
            else
            {
                BELFEM_ERROR( tNatDim <= mDim, "Dimension mismatch" ) ;
            }

            mJ.set_size( tNatDim, tNatDim );
            mInvJ.set_size( tNatDim, tNatDim );
            mX.set_size( n, mDim ) ;

            for ( uint j=0; j<mDim; ++j )
            {
                for ( uint i=0; i<n; ++i )
                {
                    mX( i, j ) = aElement->node( i )->x( j );
                }
            }

            if ( mesh::geometry_type( aElement->type() ) == GeometryType::QUAD )
            {
                // mA: rows = {x, y}, cols = {1, xi, eta, xi*eta}
                mA.set_size( 2, 4 );

                mA( 0, 0 ) =  mX(0,0)+mX(1,0)+mX(2,0)+mX(3,0);
                mA( 0, 1 ) = -mX(0,0)+mX(1,0)+mX(2,0)-mX(3,0);
                mA( 0, 2 ) = -mX(0,0)-mX(1,0)+mX(2,0)+mX(3,0);
                mA( 0, 3 ) =  mX(0,0)-mX(1,0)+mX(2,0)-mX(3,0);

                mA( 1, 0 ) =  mX(0,1)+mX(1,1)+mX(2,1)+mX(3,1);
                mA( 1, 1 ) = -mX(0,1)+mX(1,1)+mX(2,1)-mX(3,1);
                mA( 1, 2 ) = -mX(0,1)-mX(1,1)+mX(2,1)+mX(3,1);
                mA( 1, 3 ) =  mX(0,1)-mX(1,1)+mX(2,1)-mX(3,1);

                mA *= 0.25 ;

                mIsAffin = std::abs( mA( 0, 3 ) ) < BELFEM_MESH_EPSILON &&
                           std::abs( mA( 1, 3 ) ) < BELFEM_MESH_EPSILON ;
            }

            if ( aElement->type() != mElementType )
            {
                mElementType = aElement->type() ;
                switch ( aElement->type() )
                {
                    case ElementType::TRI3 :
                    {
                        mFunEval = & ElementMapper::evaluate_tri3 ;
                        mFunCheck = & ElementMapper::inside_tri ;
                        mFunGuess = nullptr ;
                        break ;
                    }
                    case ElementType::QUAD4 :
                    {
                        mFunEval = & ElementMapper::evaluate_quad4 ;
                        mFunCheck = & ElementMapper::inside_quad ;
                        mFunGuess = nullptr ;
                        break ;
                    }
                    case ElementType::TET4 :
                    {
                        mFunEval = & ElementMapper::evaluate_tet4 ;
                        mFunCheck = & ElementMapper::inside_tet ;
                        mFunGuess = nullptr ;
                        break ;
                    }
                    default:
                    {
                        mFunEval = & ElementMapper::evaluate_general ;
                        switch ( mesh::geometry_type( aElement->type() ) )
                        {
                            case GeometryType::TRI :
                            {
                                mFunCheck = & ElementMapper::inside_tri ;
                                mFunGuess = & ElementMapper::evaluate_tri3 ;
                                mXi0 = { 1./3., 1./3. };
                                break ;
                            }
                            case GeometryType::QUAD :
                            {
                                mFunCheck = & ElementMapper::inside_quad ;
                                mFunGuess = & ElementMapper::evaluate_quad4 ;
                                mXi0   = { 0., 0. };
                                break ;
                            }
                            case GeometryType::TET :
                            {
                                mFunCheck = & ElementMapper::inside_tet ;
                                mFunGuess = & ElementMapper::guess_general ;
                                mXi0   = { 0.25, 0.25, 0.25 };
                                break ;
                            }
                            case GeometryType::PENTA :
                            {
                                mFunCheck = & ElementMapper::inside_penta ;
                                mFunGuess = & ElementMapper::guess_general ;
                                mXi0 = { 1./3., 1./3., 0.0 };
                                break;
                            }
                            case GeometryType::PYRA :
                            {
                                mFunCheck = & ElementMapper::inside_pyra ;
                                mFunGuess = & ElementMapper::guess_general ;
                                mXi0 = { 0., 0., 0.25 };
                                break;
                            }
                            case GeometryType::HEX :
                            {
                                mFunCheck = & ElementMapper::inside_hex ;
                                mFunGuess = & ElementMapper::guess_general ;
                                mXi0 = { 0., 0., 0. };
                                break;
                            }
                            default:
                            {
                                BELFEM_ERROR( false, "Unsupported geometry type" ) ;
                            }
                        }
                    }
                }
            }

            if ( ! mFunctions.key_exists( aElement->type() ) )
            {
                mFunction = mFactory->create_lagrange_function( aElement->type() ) ;
                mFunctions[ aElement->type() ] = mFunction ;
            }
            else
            {
                mFunction = mFunctions( aElement->type() ) ;
            }

            mRHS.set_size( aElement->dimension() );
            mN.set_size( 1, aElement->number_of_nodes() );
            mNxi.set_size( aElement->dimension(), aElement->number_of_nodes() );

        }

        bool
        ElementMapper::evaluate( const Vector< real > & aX, Vector< real > & aXi )
        {
            return ( this->*mFunEval )( aX, aXi ) ;
        }

        bool
        ElementMapper::evaluate_tri3( const Vector< real > & aX, Vector< real > & aXi )
        {
            mJ( 0, 0 ) = mX(0,0)-mX(2,0);
            mJ( 1, 0 ) = mX(0,1)-mX(2,1);
            mJ( 0, 1 ) = mX(1,0)-mX(2,0);
            mJ( 1, 1 ) = mX(1,1)-mX(2,1);

            mRHS( 0 ) = aX(0)-mX(2,0);
            mRHS( 1 ) = aX(1)-mX(2,1);

            inv2( mJ, mInvJ );
            aXi = mInvJ * mRHS ;

            if ( aXi( 0 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 0 ) + aXi( 1 )> 1.+BELFEM_MESH_EPSILON ) return false ;

            return true ;
        }

        bool
        ElementMapper::evaluate_quad4( const Vector< real > & aX, Vector< real > & aXi )
        {
            // mA layout: rows = {x, y}, cols = {1, xi, eta, xi*eta}
            if ( mIsAffin )
            {
                mJ( 0, 0 ) = mA( 0, 1 );
                mJ( 0, 1 ) = mA( 0, 2 );
                mJ( 1, 0 ) = mA( 1, 1 );
                mJ( 1, 1 ) = mA( 1, 2 );

                mRHS( 0 ) = aX( 0 ) - mA( 0, 0 );
                mRHS( 1 ) = aX( 1 ) - mA( 1, 0 );

                inv2( mJ, mInvJ );

                aXi = mInvJ * mRHS ;

                if ( std::abs( aXi( 0 ) ) > 1.0 + BELFEM_MESH_EPSILON ) return false ;
                if ( std::abs( aXi( 1 ) ) > 1.0 + BELFEM_MESH_EPSILON ) return false ;

                return true ;
            }
            else
            {
                // Eliminate xi from x = a0 + (a1 + a3*eta)*xi + a2*eta to get
                // A*eta^2 + B*eta + C = 0, with a_k = mA(0,k), b_k = mA(1,k).
                real a = mA(0,2)*mA(1,3) - mA(0,3)*mA(1,2) ;
                real b = mA(0,3)*( aX(1) - mA(1,0) )
                       - mA(1,3)*( aX(0) - mA(0,0) )
                       + mA(0,2)*mA(1,1) - mA(0,1)*mA(1,2) ;
                real c = mA(0,1)*( aX(1) - mA(1,0) )
                       - mA(1,1)*( aX(0) - mA(0,0) ) ;

                real & xi  = aXi( 0 ) ;
                real & eta = aXi( 1 ) ;

                if ( std::abs( a ) < BELFEM_MESH_EPSILON )
                {
                    eta = -c / b ;
                }
                else
                {
                    real d = b*b - 4*a*c ;

                    if ( d < 0. && d > -BELFEM_MESH_EPSILON )
                    {
                        d = 0. ;
                    }

                    eta = ( -b + std::sqrt( d ) ) / ( 2*a ) ;

                    if ( eta < -1.-BELFEM_MESH_EPSILON || eta > 1.+BELFEM_MESH_EPSILON )
                    {
                        eta = ( -b - std::sqrt( d ) ) / ( 2*a ) ;
                    }
                }

                xi = ( aX(0) - mA(0,0) - mA(0,2)*eta ) / ( mA(0,1) + mA(0,3)*eta ) ;

                if ( std::abs( aXi( 0 ) ) > 1.0 + BELFEM_MESH_EPSILON ) return false ;
                if ( std::abs( aXi( 1 ) ) > 1.0 + BELFEM_MESH_EPSILON ) return false ;

                return true ;
            }
        }

        bool
        ElementMapper::evaluate_tet4( const Vector< real > & aX, Vector< real > & aXi )
        {
            mJ( 0, 0 ) = mX(0,0)-mX(3,0);
            mJ( 1, 0 ) = mX(0,1)-mX(3,1);
            mJ( 2, 0 ) = mX(0,2)-mX(3,2);

            mJ( 0, 1 ) = mX(1,0)-mX(3,0);
            mJ( 1, 1 ) = mX(1,1)-mX(3,1);
            mJ( 2, 1 ) = mX(1,2)-mX(3,2);

            mJ( 0, 2 ) = mX(2,0)-mX(3,0);
            mJ( 1, 2 ) = mX(2,1)-mX(3,1);
            mJ( 2, 2 ) = mX(2,2)-mX(3,2);

            mRHS( 0 ) = aX(0)-mX(3,0);
            mRHS( 1 ) = aX(1)-mX(3,1);
            mRHS( 2 ) = aX(2)-mX(3,2);

            inv3( mJ, mInvJ );

            aXi = mInvJ * mRHS ;

            std::swap( aXi( 1 ), aXi( 2 ) );

            if ( aXi( 0 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 2 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 0 ) + aXi( 1 ) + aXi( 2 ) > 1.+BELFEM_MESH_EPSILON ) return false ;

            return true ;
        }

        bool
        ElementMapper::evaluate_general( const Vector< real > & aX, Vector< real > & aXi )
        {
            const uint tNumNodes = mN.n_cols() ;
            const uint tNatDim   = mJ.n_rows() ;

            // seed initial guess (return value is informational; aXi is always set)
            ( this->*mFunGuess )( aX, aXi );

            // residual r = N(xi) * X - aX at initial guess
            mFunction->N( aXi, mN );
            for ( uint k = 0; k < tNatDim; ++k )
            {
                real tSum = 0.0 ;
                for ( uint i = 0; i < tNumNodes; ++i )
                {
                    tSum += mN( 0, i ) * mX( i, k );
                }
                mRHS( k ) = tSum - aX( k );
            }
            real tEpsilon = norm( mRHS );

            index_t tCount = 0 ;
            while ( tEpsilon > BELFEM_EPSILON && tCount++ < 100 )
            {
                mFunction->dNdXi( aXi, mNxi );
                mJ = trans( mNxi * mX );

                if ( tNatDim == 2 )
                {
                    inv2( mJ, mInvJ );
                }
                else
                {
                    inv3( mJ, mInvJ );
                }

                aXi -= mInvJ * mRHS ;

                mFunction->N( aXi, mN );
                for ( uint k = 0; k < tNatDim; ++k )
                {
                    real tSum = 0.0 ;
                    for ( uint i = 0; i < tNumNodes; ++i )
                    {
                        tSum += mN( 0, i ) * mX( i, k );
                    }
                    mRHS( k ) = tSum - aX( k );
                }
                tEpsilon = norm( mRHS );
            }
            BELFEM_ERROR( tCount < 100, "ElementMapper::evaluate_general: failed to converge." ) ;

            return ( this->*mFunCheck )( aXi );
        }

        bool
        ElementMapper::guess_general(  const Vector< real > & aX, Vector< real > & aXi )
        {
            aXi = mXi0 ;
            return true ;
        }

        bool
        ElementMapper::inside_tri( const Vector< real > & aXi ) const
        {
            if ( aXi( 0 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 0 ) + aXi( 1 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            return true ;
        }

        bool
        ElementMapper::inside_quad( const Vector< real > & aXi ) const
        {
            if ( aXi( 0 ) < -1.-BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) < -1.-BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 0 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            return true ;
        }

        bool
        ElementMapper::inside_tet( const Vector< real > & aXi ) const
        {
            if ( aXi( 0 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 2 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 0 ) + aXi( 1 ) + aXi( 2 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            return true ;
        }

        bool
        ElementMapper::inside_penta( const Vector< real > & aXi ) const
        {
            if ( aXi( 0 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 0 ) + aXi( 1 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 2 ) < -1.-BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 2 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            return true ;
        }

        bool
        ElementMapper::inside_pyra( const Vector< real > & aXi ) const
        {
            if ( aXi( 0 ) < -1.-BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) < -1.-BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 0 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 2 ) < -BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 2 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            return true ;
        }

        bool
        ElementMapper::inside_hex( const Vector< real > & aXi ) const
        {
            if ( aXi( 0 ) < -1.-BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) < -1.-BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 2 ) < -1.-BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 0 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 1 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            if ( aXi( 2 ) > 1.+BELFEM_MESH_EPSILON ) return false ;
            return true ;
        }

        /**
        * computes the weights by evaluating the shape function
        * @param aXi
        * @param aWeights
        */
        void
        ElementMapper::weights( const Vector< real > & aXi, Vector< real > & aWeights )
        {
            mFunction->N( aXi, mN );

            uint n = mN.n_cols();
            for ( uint i = 0; i < n; ++i )
            {
                aWeights( i ) = mN( 0, i );
            }
        }

    }
}