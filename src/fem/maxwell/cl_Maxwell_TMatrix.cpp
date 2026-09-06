/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */
#include "cl_Maxwell_TMatrix.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            TMatrix::TMatrix( Mesh * aMesh ) :
                mGroup( new SideSet( ElementType::TRI6, ElementType::TET10, ElementType::TET10 ) ),
                mCalc( new Calculator( mGroup, aMesh ) )
            {
                mGroup->initialize_lookup_tables( 7 );

                mCalc->initialize_integration(  ElementType::TRI6, InterpolationType::LAGRANGE );
                mCalc->set_integration_order( 4 );
                mCalc->allocate();
                mNabla.set_size( 3, 4 );
                mResult6.set_size( 12 );
                mResult10.set_size( 2, 10 );
                mE1.set_size( 3 );
                mE2.set_size( 3 );
                mIndexLookup = { { 0, 2, 0, 0 }, { 2, 1, 3, 1 }, { 3, 3, 1, 2 } };
                mFacetLookup = { { 0, 1, 3, 4, 8, 7 }, {1, 2, 3, 5, 9, 8}, {0, 3, 2, 7, 9, 6}, {0, 2, 1, 6, 5, 4} };
            }

            TMatrix::~TMatrix()
            {
                delete mCalc;
                delete mGroup;
            }

            const Vector< real > &
            TMatrix::process( mesh::Facet * aFacet )
            {
                mCalc->link( aFacet );
                mResult10.fill( 0.0 );

                // get the indices for this facet
                uint i = mIndexLookup( 0, aFacet->index_on_master() );
                uint j = mIndexLookup( 1, aFacet->index_on_master() );
                uint k = mIndexLookup( 2, aFacet->index_on_master() );

                const Vector< real > & w = mCalc->integration()->weights() ;


                for ( uint p=0; p<mCalc->num_intpoints(); ++p )
                {
                    // compute the normal
                    const Vector< real > & n = mCalc->normal( p );

                    // gradient on slave
                    const Matrix< real > & B = mCalc->Bs( p );

                    // todo: these lines can be deleted
                    Matrix< real > P( mCalc->N( p ) * mCalc->X() );
                    Matrix< real > Q( mCalc->Nm( p ) * mCalc->Xm() );
                    Matrix< real > R( mCalc->Ns( p ) * mCalc->Xs() );

                    BELFEM_ASSERT(  std::sqrt(
                                      std::pow(( P(0,0)-Q(0,0)),2)
                                     +std::pow(( P(0,1)-Q(0,1)),2)
                                     +std::pow(( P(0,2)-Q(0,2)),2)) < 1e-10,
                                     "Master point does not match" );

                   BELFEM_ASSERT(  std::sqrt(
                         std::pow(( P(0,0)-R(0,0)),2)
                        +std::pow(( P(0,1)-R(0,1)),2)
                        +std::pow(( P(0,2)-R(0,2)),2)) < 1e-10, "Slave Point does not match");

                    // evaluate function for surface
                    this->compute_nedelec_function( i, j, k, p );

                    for ( uint q=0; q<10; ++q )
                    {
                       mResult10( 0, q ) -= w( p ) *
                            ( mE1( 0 ) * ( B(1,q) * n(2) - B(2,q) * n(1) )
                            + mE1( 1 ) * ( B(2,q) * n(0) - B(0,q) * n(2) )
                            + mE1( 2 ) * ( B(0,q) * n(1) - B(1,q) * n(0) ) ) * mCalc->dS( p ) ;

                       mResult10( 1, q ) -= w( p ) *
                                                  ( mE2( 0 ) * ( B(1,q) * n(2) - B(2,q) * n(1) )
                                                  + mE2( 1 ) * ( B(2,q) * n(0) - B(0,q) * n(2) )
                                                  + mE2( 2 ) * ( B(0,q) * n(1) - B(1,q) * n(0) ) ) * mCalc->dS( p ) ;

                    }
                }

                // extract values for facet
                k = 0 ;
                for ( i=0; i<2; ++i )
                {
                    for ( j=0; j<6; ++j )
                    {
                        real tValue =  mResult10( i, mFacetLookup( aFacet->index_on_master(), j ) );
                        mResult6( k++ ) = std::abs( tValue ) < 1e-12 ? 0.0 : tValue ;
                    }
                }

                // return the value
                return mResult6 ;
            }


            void
            TMatrix::compute_nabla( const uint aIndex )
            {
                const Matrix< real > & J = mCalc->Jm( aIndex );

                // the first three columns are simply the inverse of J
                mNabla(0,0) = J(1,1)*J(2,2)- J(1,2)*J(2,1) ;
                mNabla(1,0) = J(1,2)*J(2,0)-J(1,0)*J(2,2) ;
                mNabla(2,0) = J(1,0)*J(2,1)- J(1,1)*J(2,0) ;
                mNabla(0,1) = J(0,2)*J(2,1)-J(0,1)*J(2,2) ;
                mNabla(1,1) = J(0,0)*J(2,2)- J(0,2)*J(2,0) ;
                mNabla(2,1) = J(0,1)*J(2,0)-J(0,0)*J(2,1);
                mNabla(0,2) = J(0,1)*J(1,2)- J(0,2)*J(1,1) ;
                mNabla(1,2) = J(0,2)*J(1,0)-J(0,0)*J(1,2) ;
                mNabla(2,2) = J(0,0)*J(1,1)- J(0,1)*J(1,0) ;

                // compute the final column
                mNabla(0,3) = -(mNabla(0,0)+mNabla(0,1)+mNabla(0,2));
                mNabla(1,3) = -(mNabla(1,0)+mNabla(1,1)+mNabla(1,2));
                mNabla(2,3) = -(mNabla(2,0)+mNabla(2,1)+mNabla(2,2));

                // scaling
                mNabla /=  J(0,0)*(J(1,1)*J(2,2) - J(1,2)*J(2,1)) + J(0,1)*(J(1,2)*J(2,0)- J(1,0)*J(2,2))  + J(0,2)*(J(1,0)*J(2,1) - J(1,1)*J(2,0)) ;
            }

            void
            TMatrix::compute_nedelec_function( const uint aI, const uint aJ, const uint aK, const uint aIndex )
            {
                const Matrix< real > & points = mCalc->master_integration()->points();

                // compute the nablas
                this->compute_nabla( aIndex );

                real lambda_i = points(aI,aIndex );
                real lambda_j = points(aJ,aIndex );
                real lambda_k = points(aK,aIndex );

                // compute first function
                real u = 16. * lambda_j * lambda_k ;
                real v = -8. * lambda_i * lambda_k ;
                real w = -8. * lambda_i * lambda_j ;

                for ( uint l=0; l<3; ++l )
                {
                    mE1(l) = u * mNabla(l,aI )
                           + v * mNabla(l,aJ )
                           + w * mNabla(l,aK );
                }

                // compute second function
                u =  -8 * lambda_j * lambda_k ;
                v = 16. * lambda_i * lambda_k ;
                w = -8. * lambda_i * lambda_j ;

                for ( uint l=0; l<3; ++l )
                {
                    mE2(l) = u * mNabla(l,aI )
                           + v * mNabla(l,aJ )
                           + w * mNabla(l,aK );
                }
            }

        }
    }
}