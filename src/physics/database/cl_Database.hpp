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

#ifndef BELFEM_CL_DATABASE_HPP
#define BELFEM_CL_DATABASE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_TensorMeshFactory.hpp"
#include "cl_Mesh.hpp"
#include "hdf5_tools.hpp"

namespace belfem
{
    /**
     * @brief Precomputed lookup table on a tensor grid; evaluation, derivatives and HDF5 persistence.
     *
     * @ingroup grp_physics_database
     * @see @ref physics_database_database_usage_guide
     */
    class Database
    {
        const proc_t mCommRank ;
        const proc_t mCommSize ;
        const string mLabel ;

        Vector< real > mValues ;
        Matrix< index_t > mTopology ;
        TensorMeshConfig * mConfig = nullptr ;

        void ( Database::*mFunction2D ) ( const real xi, const real eta, real * N ) const = nullptr ;
        void ( Database::*mFunction3D ) ( const real xi, const real eta, const real zeta, real * N ) const = nullptr ;

        void ( Database::*mdFunction2Ddxi ) ( const real xi, const real eta, real * N ) const = nullptr ;
        void ( Database::*mdFunction3Ddxi ) ( const real xi, const real eta, const real zeta, real * N ) const = nullptr ;

        void ( Database::*mdFunction2Ddeta ) ( const real xi, const real eta, real * N ) const = nullptr ;
        void ( Database::*mdFunction3Ddeta ) ( const real xi, const real eta, const real zeta, real * N ) const = nullptr ;

        void ( Database::*mdFunction3Ddzeta ) ( const real xi, const real eta, const real zeta, real * N ) const = nullptr ;


        real * mN = nullptr ;
        uint  mNumNodesPerElement = 0 ;

    public:

        Database( const string & aFilePath, const string & aMaterial );

        Database( const hid_t aHDF5, const string & aLabel );

        Database( Mesh * aMesh, const string & aField, const bool aProject = true, const string aMaterial="" );

        ~Database();

        real
        evaluate( const real aX, const real aY ) const ;

        real
        evaluate( const real aX, const real aY, const real aZ ) const ;

        real
        evaluate_derivx( const real aX, const real aY ) const ;

        real
        evaluate_derivy( const real aX, const real aY ) const ;

        real
        evaluate_derivx( const real aX, const real aY, const real aZ ) const ;

        real
        evaluate_derivy( const real aX, const real aY, const real aZ ) const ;

        real
        evaluate_derivz( const real aX, const real aY, const real aZ ) const ;

        real
        min( const uint aDimension ) const ;

        real
        max( const uint aDimension ) const ;

        //! read access to the raw grid values ( in the stored scale, e.g.
        //! log10 for the jc/n tables ), for range checks at load time
        const Vector< real > &
        values() const
        {
            return mValues ;
        }

        void
        save( const hid_t aHDF5 );

    private:

        void
        load( const hid_t aHDF5 );

        void
        set_element_type();

        index_t
        nidx( const index_t i, const index_t j ) const;

        index_t
        nidx( const index_t i, const index_t j, const index_t k ) const;

        void
        eval_quad4( const real xi, const real eta, real * N ) const ;

        void
        eval_quad9( const real xi, const real eta, real * N ) const ;

        void
        eval_quad16( const real xi, const real eta, real * N ) const ;

        void
        eval_hex8( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        eval_hex27( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        eval_hex64( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        deval_quad4dxi( const real xi, const real eta, real * N ) const ;

        void
        deval_quad4deta( const real xi, const real eta, real * N ) const ;

        void
        deval_quad9dxi( const real xi, const real eta, real * N ) const ;

        void
        deval_quad9deta( const real xi, const real eta, real * N ) const ;


        void
        deval_quad16dxi( const real xi, const real eta, real * N ) const ;

        void
        deval_quad16deta( const real xi, const real eta, real * N ) const ;

        void
        deval_hex8dxi( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        deval_hex8deta( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        deval_hex8dzeta( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        deval_hex27dxi( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        deval_hex64dxi( const real xi, const real eta, const real zeta, real * N ) const ;


        void
        deval_hex27deta( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        deval_hex64deta( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        deval_hex27dzeta( const real xi, const real eta, const real zeta, real * N ) const ;

        void
        deval_hex64dzeta( const real xi, const real eta, const real zeta, real * N ) const ;
    };

    inline index_t
    Database::nidx( const index_t i, const index_t j ) const
    {
        return i * mConfig->num_nodes( 1 ) + j ;
    }

    index_t
    inline Database::nidx( const index_t i, const index_t j, const index_t k ) const
    {
        return mConfig->num_nodes( 2 ) * ( mConfig->num_nodes( 1 ) * i + j ) + k ;
    }

    inline
    real Database::evaluate( const real aX, const real aY ) const
    {
        // find element index
        index_t i = mConfig->element_ijk( 0, aX );
        index_t j = mConfig->element_ijk( 1, aY );
        index_t e = mConfig->element_index( i, j );

        // get center point of element
        real xc = mConfig->min( 0 ) + ( i + 0.5 ) * mConfig->element_step( 0 );
        real yc = mConfig->min( 1 ) + ( j + 0.5 ) * mConfig->element_step( 1 );

        // compute parameter coordinates
        real xi   = 2. * ( aX - xc ) * mConfig->inv_element_step( 0 ) ;
        real eta  = 2. * ( aY - yc ) * mConfig->inv_element_step( 1 ) ;

        // evaluate shape function
        ( this->*mFunction2D )( xi, eta, mN );

        // interpolate value
        real aValue = 0.0 ;
        for ( index_t n = 0; n < mNumNodesPerElement ; ++n )
        {
            aValue += mN[n] * mValues( mTopology( n, e ) );
        }
        return aValue;
    }

    inline
    real Database::evaluate_derivx( const real aX, const real aY ) const
    {
        // find element index
        index_t i = mConfig->element_ijk( 0, aX );
        index_t j = mConfig->element_ijk( 1, aY );
        index_t e = mConfig->element_index( i, j );

        // get center point of element
        real xc = mConfig->min( 0 ) + ( i + 0.5 ) * mConfig->element_step( 0 );
        real yc = mConfig->min( 1 ) + ( j + 0.5 ) * mConfig->element_step( 1 );

        // compute parameter coordinates
        real xi   = 2. * ( aX - xc ) * mConfig->inv_element_step( 0 ) ;
        real eta  = 2. * ( aY - yc ) * mConfig->inv_element_step( 1 ) ;

        // evaluate shape function
        ( this->*mdFunction2Ddxi )( xi, eta, mN );

        // interpolate value
        real aValue = 0.0 ;
        for ( index_t n = 0; n < mNumNodesPerElement ; ++n )
        {
            aValue += mN[n] * mValues( mTopology( n, e ) );
        }
        aValue *= 2. * mConfig->inv_element_step( 0 ) ;
        return aValue;
    }

    inline
    real Database::evaluate_derivy( const real aX, const real aY ) const
    {
        // find element index
        index_t i = mConfig->element_ijk( 0, aX );
        index_t j = mConfig->element_ijk( 1, aY );
        index_t e = mConfig->element_index( i, j );

        // get center point of element
        real xc = mConfig->min( 0 ) + ( i + 0.5 ) * mConfig->element_step( 0 );
        real yc = mConfig->min( 1 ) + ( j + 0.5 ) * mConfig->element_step( 1 );

        // compute parameter coordinates
        real xi   = 2. * ( aX - xc ) * mConfig->inv_element_step( 0 ) ;
        real eta  = 2. * ( aY - yc ) * mConfig->inv_element_step( 1 ) ;

        // evaluate shape function
        ( this->*mdFunction2Ddeta )( xi, eta, mN );

        // interpolate value
        real aValue = 0.0 ;
        for ( index_t n = 0; n < mNumNodesPerElement ; ++n )
        {
            aValue += mN[n] * mValues( mTopology( n, e ) );
        }
        aValue *= 2. * mConfig->inv_element_step( 1 ) ;
        return aValue;
    }

    inline
    real Database::evaluate( const real aX, const real aY, const real aZ ) const
    {
        // find element index
        index_t i = mConfig->element_ijk( 0, aX );
        index_t j = mConfig->element_ijk( 1, aY );
        index_t k = mConfig->element_ijk( 2, aZ );
        index_t e = mConfig->element_index( i, j, k );

        // get center point of element
        real xc = mConfig->min( 0 ) + ( i + 0.5 ) * mConfig->element_step( 0 );
        real yc = mConfig->min( 1 ) + ( j + 0.5 ) * mConfig->element_step( 1 );
        real zc = mConfig->min( 2 ) + ( k + 0.5 ) * mConfig->element_step( 2 );

        // compute parameter coordinates
        real xi   = 2. * ( aX - xc ) * mConfig->inv_element_step( 0 ) ;
        real eta  = 2. * ( aY - yc ) * mConfig->inv_element_step( 1 ) ;
        real zeta = 2. * ( aZ - zc ) * mConfig->inv_element_step( 2 ) ;

        // evaluate the shape function
        ( this->*mFunction3D )( xi, eta, zeta, mN );

        // interpolate value
        real aValue = 0.0 ;
        for ( index_t n = 0; n < mNumNodesPerElement ; ++n )
        {
            aValue += mN[n] * mValues( mTopology( n, e ) );
        }
        return aValue;
    }

    inline
    real Database::evaluate_derivx( const real aX, const real aY, const real aZ ) const
    {
        // find element index
        index_t i = mConfig->element_ijk( 0, aX );
        index_t j = mConfig->element_ijk( 1, aY );
        index_t k = mConfig->element_ijk( 2, aZ );
        index_t e = mConfig->element_index( i, j, k );

        // get center point of element
        real xc = mConfig->min( 0 ) + ( i + 0.5 ) * mConfig->element_step( 0 );
        real yc = mConfig->min( 1 ) + ( j + 0.5 ) * mConfig->element_step( 1 );
        real zc = mConfig->min( 2 ) + ( k + 0.5 ) * mConfig->element_step( 2 );

        // compute parameter coordinates
        real xi   = 2. * ( aX - xc ) * mConfig->inv_element_step( 0 ) ;
        real eta  = 2. * ( aY - yc ) * mConfig->inv_element_step( 1 ) ;
        real zeta = 2. * ( aZ - zc ) * mConfig->inv_element_step( 2 ) ;

        // evaluate the shape function
        ( this->*mdFunction3Ddxi )( xi, eta, zeta, mN );

        // interpolate value
        real aValue = 0.0 ;
        for ( index_t n = 0; n < mNumNodesPerElement ; ++n )
        {
            aValue += mN[n] * mValues( mTopology( n, e ) );
        }
        aValue *= 2. * mConfig->inv_element_step( 0 ) ;
        return aValue;
    }

    inline
    real Database::evaluate_derivy( const real aX, const real aY, const real aZ ) const
    {
        // find element index
        index_t i = mConfig->element_ijk( 0, aX );
        index_t j = mConfig->element_ijk( 1, aY );
        index_t k = mConfig->element_ijk( 2, aZ );
        index_t e = mConfig->element_index( i, j, k );

        // get center point of element
        real xc = mConfig->min( 0 ) + ( i + 0.5 ) * mConfig->element_step( 0 );
        real yc = mConfig->min( 1 ) + ( j + 0.5 ) * mConfig->element_step( 1 );
        real zc = mConfig->min( 2 ) + ( k + 0.5 ) * mConfig->element_step( 2 );

        // compute parameter coordinates
        real xi   = 2. * ( aX - xc ) * mConfig->inv_element_step( 0 ) ;
        real eta  = 2. * ( aY - yc ) * mConfig->inv_element_step( 1 ) ;
        real zeta = 2. * ( aZ - zc ) * mConfig->inv_element_step( 2 ) ;

        // evaluate the shape function
        ( this->*mdFunction3Ddeta )( xi, eta, zeta, mN );

        // interpolate value
        real aValue = 0.0 ;
        for ( index_t n = 0; n < mNumNodesPerElement ; ++n )
        {
            aValue += mN[n] * mValues( mTopology( n, e ) );
        }
        aValue *= 2. * mConfig->inv_element_step( 1 ) ;
        return aValue;
    }

    inline
    real Database::evaluate_derivz( const real aX, const real aY, const real aZ ) const
    {
        // find element index
        index_t i = mConfig->element_ijk( 0, aX );
        index_t j = mConfig->element_ijk( 1, aY );
        index_t k = mConfig->element_ijk( 2, aZ );
        index_t e = mConfig->element_index( i, j, k );

        // get center point of element
        real xc = mConfig->min( 0 ) + ( i + 0.5 ) * mConfig->element_step( 0 );
        real yc = mConfig->min( 1 ) + ( j + 0.5 ) * mConfig->element_step( 1 );
        real zc = mConfig->min( 2 ) + ( k + 0.5 ) * mConfig->element_step( 2 );

        // compute parameter coordinates
        real xi   = 2. * ( aX - xc ) * mConfig->inv_element_step( 0 ) ;
        real eta  = 2. * ( aY - yc ) * mConfig->inv_element_step( 1 ) ;
        real zeta = 2. * ( aZ - zc ) * mConfig->inv_element_step( 2 ) ;

        // evaluate the shape function
        ( this->*mdFunction3Ddzeta )( xi, eta, zeta, mN );

        // interpolate value
        real aValue = 0.0 ;
        for ( index_t n = 0; n < mNumNodesPerElement ; ++n )
        {
            aValue += mN[n] * mValues( mTopology( n, e ) );
        }
        aValue *= 2. * mConfig->inv_element_step( 2 ) ;
        return aValue;
    }


    inline void
    Database::eval_quad4( const real xi, const real eta, real * N ) const
    {
        N[0] = ( ( 1.0 - xi ) * ( 1.0 - eta ) ) * 0.25;
        N[1] = ( ( 1.0 + xi ) * ( 1.0 - eta ) ) * 0.25;
        N[2] = ( ( 1.0 + xi ) * ( 1.0 + eta ) ) * 0.25;
        N[3] = ( ( 1.0 - xi ) * ( 1.0 + eta ) ) * 0.25;
    }

    inline void
    Database::eval_quad9( const real xi, const real eta, real * N ) const
    {
        const real    c = xi * eta * 0.25;
        const real  xi2 = xi*xi;
        const real eta2 = eta*eta;

        N[0] = ( c * ( eta - 1.0 ) * (xi - 1.0) );
        N[1] = ( c * ( eta - 1.0 ) * (xi + 1.0) );
        N[2] = ( c * ( eta + 1.0 ) * (xi + 1.0) );
        N[3] = ( c * ( eta + 1.0 ) * (xi - 1.0) );
        N[4] = ( eta * ( 1.0 - xi2 ) * ( eta - 1.0 ) ) * 0.5;
        N[5] = ( xi * ( 1.0 - eta2)*( xi + 1.0 ) )*0.5;
        N[6] = ( eta * (1.0 - xi2)*( eta + 1.0 ) )*0.5;
        N[7] = ( xi*( 1.0 - eta2 )*( xi - 1.0 ) )*0.5;
        N[8] = ( eta2 - 1.0 )*( xi2 - 1.0 );
    }

    inline void
    Database::eval_quad16( const real xi, const real eta, real * N ) const
    {
        const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
        const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
        const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
        const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

        const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
        const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
        const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
        const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

        N[ 0 ] = a0*b0;
        N[ 1 ] = a3*b0;
        N[ 2 ] = a3*b3;
        N[ 3 ] = a0*b3;
        N[ 4 ] = a1*b0;
        N[ 5 ] = a2*b0;
        N[ 6 ] = a3*b1;
        N[ 7 ] = a3*b2;
        N[ 8 ] = a2*b3;
        N[ 9 ] = a1*b3;
        N[ 10 ] = a0*b2;
        N[ 11 ] = a0*b1;
        N[ 12 ] = a1*b1;
        N[ 13 ] = a2*b1;
        N[ 14 ] = a2*b2;
        N[ 15 ] = a1*b2;
    }

    inline void
    Database::eval_hex8( const real xi, const real eta, const real zeta, real * N ) const
    {
        N[0] =  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
        N[1] =    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
        N[2] =  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
        N[3] =    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
        N[4] =    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
        N[5] =  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
        N[6] =    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
        N[7] =  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
    }

    inline void
    Database::eval_hex27( const real xi, const real eta, const real zeta, real * N ) const
    {
        const real   xi2 = xi*xi;
        const real  eta2 = eta*eta;
        const real zeta2 = zeta*zeta;

        const real a = -0.25 * eta * zeta;
        const real b = -0.25 * xi * zeta;
        const real c = -0.25 * xi * eta;
        const real d = 0.125 * xi * eta * zeta;

        N[  0 ] = d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        N[  1 ] = d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        N[  2 ] = d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        N[  3 ] = d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        N[  4 ] = d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        N[  5 ] = d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        N[  6 ] = d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        N[  7 ] = d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        N[  8 ] = a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
        N[  9 ] = b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        N[ 10 ] = a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
        N[ 11 ] = b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        N[ 12 ] = c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );
        N[ 13 ] = c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );
        N[ 14 ] = c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );
        N[ 15 ] = c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );
        N[ 16 ] = a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
        N[ 17 ] = b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        N[ 18 ] = a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
        N[ 19 ] = b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        N[ 20 ] = -( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
        N[ 21 ] = ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
        N[ 22 ] = ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
        N[ 23 ] = ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.5;
        N[ 24 ] = ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.5;
        N[ 25 ] = ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.5;
        N[ 26 ] = ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.5;
    }

    inline void
    Database::eval_hex64( const real xi, const real eta, const real zeta, real * N ) const
    {
        const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
        const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
        const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
        const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

        const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
        const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
        const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
        const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

        const real c0 =  ( zeta*( 1.0 + 9.0 * zeta * ( 1.0 - zeta ) ) - 1.0 )*0.0625;
        const real c1 =  ( 9.0 - zeta * ( 27.0 + zeta*( 9.0 - 27.0*zeta ) ) )*0.0625;
        const real c2 =  ( 9.0 + zeta * ( 27.0 - zeta*( 9.0 + 27.0*zeta ) ) )*0.0625;
        const real c3 = ( -zeta*( 1.0 - 9.0 * zeta * ( 1.0 + zeta ) ) - 1.0 )*0.0625;

        N[  0 ] = a0 * b0 * c0;
        N[  1 ] = a3 * b0 * c0;
        N[  2 ] = a3 * b3 * c0;
        N[  3 ] = a0 * b3 * c0;
        N[  4 ] = a0 * b0 * c3;
        N[  5 ] = a3 * b0 * c3;
        N[  6 ] = a3 * b3 * c3;
        N[  7 ] = a0 * b3 * c3;
        N[  8 ] = a1 * b0 * c0;
        N[  9 ] = a2 * b0 * c0;
        N[ 10 ] = a0 * b1 * c0;
        N[ 11 ] = a0 * b2 * c0;
        N[ 12 ] = a0 * b0 * c1;
        N[ 13 ] = a0 * b0 * c2;
        N[ 14 ] = a3 * b1 * c0;
        N[ 15 ] = a3 * b2 * c0;
        N[ 16 ] = a3 * b0 * c1;
        N[ 17 ] = a3 * b0 * c2;
        N[ 18 ] = a2 * b3 * c0;
        N[ 19 ] = a1 * b3 * c0;
        N[ 20 ] = a3 * b3 * c1;
        N[ 21 ] = a3 * b3 * c2;
        N[ 22 ] = a0 * b3 * c1;
        N[ 23 ] = a0 * b3 * c2;
        N[ 24 ] = a1 * b0 * c3;
        N[ 25 ] = a2 * b0 * c3;
        N[ 26 ] = a0 * b1 * c3;
        N[ 27 ] = a0 * b2 * c3;
        N[ 28 ] = a3 * b1 * c3;
        N[ 29 ] = a3 * b2 * c3;
        N[ 30 ] = a2 * b3 * c3;
        N[ 31 ] = a1 * b3 * c3;
        N[ 32 ] = a1 * b1 * c0;
        N[ 33 ] = a1 * b2 * c0;
        N[ 34 ] = a2 * b2 * c0;
        N[ 35 ] = a2 * b1 * c0;
        N[ 36 ] = a1 * b0 * c1;
        N[ 37 ] = a2 * b0 * c1;
        N[ 38 ] = a2 * b0 * c2;
        N[ 39 ] = a1 * b0 * c2;
        N[ 40 ] = a0 * b1 * c1;
        N[ 41 ] = a0 * b1 * c2;
        N[ 42 ] = a0 * b2 * c2;
        N[ 43 ] = a0 * b2 * c1;
        N[ 44 ] = a3 * b1 * c1;
        N[ 45 ] = a3 * b2 * c1;
        N[ 46 ] = a3 * b2 * c2;
        N[ 47 ] = a3 * b1 * c2;
        N[ 48 ] = a2 * b3 * c1;
        N[ 49 ] = a1 * b3 * c1;
        N[ 50 ] = a1 * b3 * c2;
        N[ 51 ] = a2 * b3 * c2;
        N[ 52 ] = a1 * b1 * c3;
        N[ 53 ] = a2 * b1 * c3;
        N[ 54 ] = a2 * b2 * c3;
        N[ 55 ] = a1 * b2 * c3;
        N[ 56 ] = a1 * b1 * c1;
        N[ 57 ] = a2 * b1 * c1;
        N[ 58 ] = a2 * b2 * c1;
        N[ 59 ] = a1 * b2 * c1;
        N[ 60 ] = a1 * b1 * c2;
        N[ 61 ] = a2 * b1 * c2;
        N[ 62 ] = a2 * b2 * c2;
        N[ 63 ] = a1 * b2 * c2;
    }

    inline void
    Database::deval_quad4dxi( const real xi, const real eta, real * N ) const
    {
        N[0] = - 0.25 * ( 1.0 - eta ) ;
        N[1] =   0.25 * ( 1.0 - eta ) ;
        N[2] =   0.25 * ( 1.0 + eta ) ;
        N[3] = - 0.25 * ( 1.0 + eta ) ;
    }

    inline void
    Database::deval_quad4deta( const real xi, const real eta, real * N ) const
    {
        N[ 0 ] =  0.25 * ( xi - 1.0 );
        N[ 1 ] = -0.25 * ( xi + 1.0 );
        N[ 2 ] =  0.25 * ( xi + 1.0 );
        N[ 3 ] = -0.25 * ( xi - 1.0 );
    }

    inline void
    Database::deval_quad9dxi( const real xi, const real eta, real * N ) const
    {
        const real    c = xi * eta ;
        //const real  xi2 = xi*xi;
        const real eta2 = eta*eta;

        N[0] = ( eta * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
        N[1] = ( eta * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) ) * 0.25;
        N[2] = ( eta * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) ) * 0.25;
        N[3] = ( eta * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
        N[4] = - c * ( eta - 1.0 );
        N[5] = -( ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) ) * 0.5;
        N[6] = - c * ( eta + 1.0 );
        N[7] = -( ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) ) * 0.5;
        N[8] = 2.0 * xi * ( eta2 - 1.0 );
    }

    inline void
   Database::deval_quad9deta( const real xi, const real eta, real * N ) const
    {
        const real    c = xi * eta ;
        const real  xi2 = xi*xi;
        // const real eta2 = eta*eta;

        N[ 0 ] =  ( xi * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        N[ 1 ] =  ( xi * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
        N[ 2 ] =  ( xi * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;
        N[ 3 ] =  ( xi * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        N[ 4 ] =  -( ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) ) * 0.5;
        N[ 5 ] = - c * ( xi + 1.0 );
        N[ 6 ] =  -( ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) ) * 0.5;
        N[ 7 ] = - c * ( xi - 1.0 );
        N[ 8 ] = 2.0 * eta * ( xi2 - 1.0 );
    }

    inline void
    Database::deval_quad16dxi( const real xi, const real eta, real * N ) const
    {
        const real da0 = (   1.0 + xi*( 18.0 - 27.0*xi )) * 0.0625;
        const real da1 = ( -27.0 - xi*( 18.0 - 81.0*xi )) * 0.0625;
        const real da2 = (  27.0 - xi*( 18.0 + 81.0*xi )) * 0.0625;
        const real da3 = (  -1.0 + xi*( 18.0 + 27.0*xi )) * 0.0625;

        const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
        const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
        const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
        const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

        N[ 0 ] = da0*b0;
        N[ 1 ] = da3*b0;
        N[ 2 ] = da3*b3;
        N[ 3 ] = da0*b3;
        N[ 4 ] = da1*b0;
        N[ 5 ] = da2*b0;
        N[ 6 ] = da3*b1;
        N[ 7 ] = da3*b2;
        N[ 8 ] = da2*b3;
        N[ 9 ] = da1*b3;
        N[ 10 ] = da0*b2;
        N[ 11 ] = da0*b1;
        N[ 12 ] = da1*b1;
        N[ 13 ] = da2*b1;
        N[ 14 ] = da2*b2;
        N[ 15 ] = da1*b2;
    }




    inline void
    Database::deval_quad16deta( const real xi, const real eta, real * N ) const
    {
            // often used parameters
            const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
            const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
            const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
            const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

            const real db0 = (   1.0 + eta*( 18.0 - 27.0*eta )) * 0.0625;
            const real db1 = ( -27.0 - eta*( 18.0 - 81.0*eta )) * 0.0625;
            const real db2 = (  27.0 - eta*( 18.0 + 81.0*eta )) * 0.0625;
            const real db3 = (  -1.0 + eta*( 18.0 + 27.0*eta )) * 0.0625;

            // populate output matrix
            N[  0 ] = a0*db0;
            N[  1 ] = a3*db0;
            N[  2 ] = a3*db3;
            N[  3 ] = a0*db3;
            N[  4 ] = a1*db0;
            N[  5 ] = a2*db0;
            N[  6 ] = a3*db1;
            N[  7 ] = a3*db2;
            N[  8 ] = a2*db3;
            N[  9 ] = a1*db3;
            N[ 10 ] = a0*db2;
            N[ 11 ] = a0*db1;
            N[ 12 ] = a1*db1;
            N[ 13 ] = a2*db1;
            N[ 14 ] = a2*db2;
            N[ 15 ] = a1*db2;

    }


    inline void
    Database::deval_hex8dxi( const real xi, const real eta, const real zeta, real * N ) const
    {
        N[0] = -( eta - 1 ) * ( zeta - 1 ) * 0.125;
        N[1] =  ( eta - 1 ) * ( zeta - 1 ) * 0.125;
        N[2] = -( eta + 1 ) * ( zeta - 1 ) * 0.125;
        N[3] =  ( eta + 1 ) * ( zeta - 1 ) * 0.125;
        N[4] =  ( eta - 1 ) * ( zeta + 1 ) * 0.125;
        N[5] = -( eta - 1 ) * ( zeta + 1 ) * 0.125;
        N[6] =  ( eta + 1 ) * ( zeta + 1 ) * 0.125;
        N[7] = -( eta + 1 ) * ( zeta + 1 ) * 0.125;
    }

    inline void
   Database::deval_hex8deta( const real xi, const real eta, const real zeta, real * N ) const
    {
        N[ 0 ] = -(   xi - 1 ) * ( zeta - 1 ) * 0.125;
        N[ 1 ] =  (   xi + 1 ) * ( zeta - 1 ) * 0.125;
        N[ 2 ] = -(   xi + 1 ) * ( zeta - 1 ) * 0.125;
        N[ 3 ] =  (   xi - 1 ) * ( zeta - 1 ) * 0.125;
        N[ 4 ] =  (   xi - 1 ) * ( zeta + 1 ) * 0.125;
        N[ 5 ] = -(   xi + 1 ) * ( zeta + 1 ) * 0.125;
        N[ 6 ] =  (   xi + 1 ) * ( zeta + 1 ) * 0.125;
        N[ 7 ] = -(   xi - 1 ) * ( zeta + 1 ) * 0.125;
    }

    inline void
    Database::deval_hex8dzeta( const real xi, const real eta, const real zeta, real * N ) const
    {
        N[ 0 ] = -(  eta - 1 ) * (   xi - 1 ) * 0.125;
        N[ 1 ] =  (  eta - 1 ) * (   xi + 1 ) * 0.125;
        N[ 2 ] = -(  eta + 1 ) * (   xi + 1 ) * 0.125;
        N[ 3 ] =  (  eta + 1 ) * (   xi - 1 ) * 0.125;
        N[ 4 ] =  (  eta - 1 ) * (   xi - 1 ) * 0.125;
        N[ 5 ] = -(  eta - 1 ) * (   xi + 1 ) * 0.125;
        N[ 6 ] =  (  eta + 1 ) * (   xi + 1 ) * 0.125;
        N[ 7 ] = -(  eta + 1 ) * (   xi - 1 ) * 0.125;
    }

    inline void
    Database::deval_hex27dxi( const real xi, const real eta, const real zeta, real * N ) const
    {
        const real  eta2 = eta*eta;
        const real zeta2 = zeta*zeta;

        const real a =  0.125*eta*zeta;
        const real b =  0.125*xi*zeta;
        const real c =  0.125*xi*eta;
        const real d = -0.5*xi*eta*zeta;

        N[  0 ] = a * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
        N[  1 ] = a * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
        N[  2 ] = a * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
        N[  3 ] = a * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
        N[  4 ] = a * ( 2.0 * xi - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
        N[  5 ] = a * ( 2.0 * xi + 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
        N[  6 ] = a * ( 2.0 * xi + 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
        N[  7 ] = a * ( 2.0 * xi - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
        N[  8 ] = d * ( eta - 1.0 ) * ( zeta - 1.0 );
        N[  9 ] = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
        N[ 10 ] = d * ( eta + 1.0 ) * ( zeta - 1.0 );
        N[ 11 ] = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
        N[ 12 ] = - ( eta * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
        N[ 13 ] = - ( eta * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
        N[ 14 ] = - ( eta * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
        N[ 15 ] = - ( eta * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
        N[ 16 ] = d * ( eta - 1.0 ) * ( zeta + 1.0 );
        N[ 17 ] = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
        N[ 18 ] = d * ( eta + 1.0 ) * ( zeta + 1.0 );
        N[ 19 ] = - ( zeta * ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
        N[ 20 ] = - 2.0 * xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 );
        N[ 21 ] = 8.0 * b * ( eta2 - 1.0 ) * ( zeta - 1.0 );
        N[ 22 ] = 8.0 * b * ( eta2 - 1.0 ) * ( zeta + 1.0 );
        N[ 23 ] = ( ( eta2 - 1.0 ) * ( 2.0 * xi - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
        N[ 24 ] = ( ( eta2 - 1.0 ) * ( 2.0 * xi + 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
        N[ 25 ] = 8.0 * c * ( zeta2 - 1.0 ) * ( eta - 1.0 );
        N[ 26 ] = 8.0 * c * ( zeta2 - 1.0 ) * ( eta + 1.0 );
    }

    inline void
    Database::deval_hex27deta( const real xi, const real eta, const real zeta, real * N ) const
    {
        const real   xi2 = xi * xi ;
        const real zeta2 = zeta*zeta;

        const real a =  0.125*eta*zeta;
        const real b =  0.125*xi*zeta;
        const real c =  0.125*xi*eta;
        const real d = -0.5*xi*eta*zeta;

        N[  0 ] = b * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        N[  1 ] = b * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        N[  2 ] = b * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        N[  3 ] = b * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        N[  4 ] = b * ( 2.0 * eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        N[  5 ] = b * ( 2.0 * eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        N[  6 ] = b * ( 2.0 * eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        N[  7 ] = b * ( 2.0 * eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        N[  8 ] = - ( zeta * ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
        N[  9 ] = d * ( xi + 1.0 ) * ( zeta - 1.0 );
        N[ 10 ] =  - ( zeta * ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.25;
        N[ 11 ] = d * ( xi - 1.0 ) * ( zeta - 1.0 );
        N[ 12 ] = - ( xi * ( 2.0 * eta - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        N[ 13 ] = - ( xi * ( 2.0 * eta - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
        N[ 14 ] = - ( xi * ( 2.0 * eta + 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
        N[ 15 ] = - ( xi * ( 2.0 * eta + 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        N[ 16 ] = - ( zeta * ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
        N[ 17 ] = d * ( xi + 1.0 ) * ( zeta + 1.0 );
        N[ 18 ] = - ( zeta * ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.25;
        N[ 19 ] = d * ( xi - 1.0 ) * ( zeta + 1.0 );
        N[ 20 ] = - 2.0 * eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
        N[ 21 ] = 8.0 * a * ( xi2 - 1.0 ) * ( zeta - 1.0 );
        N[ 22 ] = 8.0 * a * ( xi2 - 1.0 ) * ( zeta + 1.0 );
        N[ 23 ] = 8.0 * c * ( zeta2 - 1.0 ) * ( xi - 1.0 );
        N[ 24 ] = 8.0 * c * ( zeta2 - 1.0 ) * ( xi + 1.0 );
        N[ 25 ] = ( ( 2.0 * eta - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
        N[ 26 ] = ( ( 2.0 * eta + 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) ) * 0.5;
    }

    inline void
    Database::deval_hex27dzeta( const real xi, const real eta, const real zeta, real * N ) const
    {
        const real   xi2 = xi * xi ;
        const real  eta2 = eta*eta;

        const real a =  0.125*eta*zeta;
        const real b =  0.125*xi*zeta;
        const real c =  0.125*xi*eta;
        const real d = -0.5*xi*eta*zeta;

        N[  0 ] = c * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );
        N[  1 ] = c * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );
        N[  2 ] = c * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );
        N[  3 ] = c * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );
        N[  4 ] = c * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );
        N[  5 ] = c * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );
        N[  6 ] = c * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );
        N[  7 ] = c * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );
        N[  8 ] = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( eta - 1.0 ) ) * 0.25;
        N[  9 ] = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( xi + 1.0 ) ) * 0.25;
        N[ 10 ] = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( eta + 1.0 ) ) * 0.25;
        N[ 11 ] = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        N[ 12 ] = d * ( eta - 1.0 ) * ( xi - 1.0 );
        N[ 13 ] = d * ( eta - 1.0 ) * ( xi + 1.0 );
        N[ 14 ] = d * ( eta + 1.0 ) * ( xi + 1.0 );
        N[ 15 ] = d * ( eta + 1.0 ) * ( xi - 1.0 );
        N[ 16 ] = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( eta - 1.0 ) ) * 0.25;
        N[ 17 ] = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( xi + 1.0 ) ) * 0.25;
        N[ 18 ] = - ( eta * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( eta + 1.0 ) ) * 0.25;
        N[ 19 ] = - ( xi * ( eta2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) * ( xi - 1.0 ) ) * 0.25;
        N[ 20 ] = - 2.0 * zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 );
        N[ 21 ] = ( ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( 2.0 * zeta - 1.0 ) ) * 0.5;
        N[ 22 ] = ( ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( 2.0 * zeta + 1.0 ) ) * 0.5;
        N[ 23 ] = 8.0 * b * ( eta2 - 1.0 ) * ( xi - 1.0 );
        N[ 24 ] = 8.0 * b * ( eta2 - 1.0 ) * ( xi + 1.0 );
        N[ 25 ] = 8.0 * a * ( xi2 - 1.0 ) * ( eta - 1.0 );
        N[ 26 ] = 8.0 * a * ( xi2 - 1.0 ) * ( eta + 1.0 );

    }

    inline void
    Database::deval_hex64dxi( const real xi, const real eta, const real zeta, real * N ) const
    {
        const real da0 = (   1.0 + xi*( 18.0 - 27.0*xi )) * 0.0625;
        const real da1 = ( -27.0 - xi*( 18.0 - 81.0*xi )) * 0.0625;
        const real da2 = (  27.0 - xi*( 18.0 + 81.0*xi )) * 0.0625;
        const real da3 = (  -1.0 + xi*( 18.0 + 27.0*xi )) * 0.0625;

        const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
        const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
        const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
        const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

        const real c0 =  ( zeta*( 1.0 + 9.0 * zeta * ( 1.0 - zeta ) ) - 1.0 )*0.0625;
        const real c1 =  ( 9.0 - zeta * ( 27.0 + zeta*( 9.0 - 27.0*zeta ) ) )*0.0625;
        const real c2 =  ( 9.0 + zeta * ( 27.0 - zeta*( 9.0 + 27.0*zeta ) ) )*0.0625;
        const real c3 = ( -zeta*( 1.0 - 9.0 * zeta * ( 1.0 + zeta ) ) - 1.0 )*0.0625;

        N[  0 ] = da0 * b0 * c0;
        N[  1 ] = da3 * b0 * c0;
        N[  2 ] = da3 * b3 * c0;
        N[  3 ] = da0 * b3 * c0;
        N[  4 ] = da0 * b0 * c3;
        N[  5 ] = da3 * b0 * c3;
        N[  6 ] = da3 * b3 * c3;
        N[  7 ] = da0 * b3 * c3;
        N[  8 ] = da1 * b0 * c0;
        N[  9 ] = da2 * b0 * c0;
        N[ 10 ] = da0 * b1 * c0;
        N[ 11 ] = da0 * b2 * c0;
        N[ 12 ] = da0 * b0 * c1;
        N[ 13 ] = da0 * b0 * c2;
        N[ 14 ] = da3 * b1 * c0;
        N[ 15 ] = da3 * b2 * c0;
        N[ 16 ] = da3 * b0 * c1;
        N[ 17 ] = da3 * b0 * c2;
        N[ 18 ] = da2 * b3 * c0;
        N[ 19 ] = da1 * b3 * c0;
        N[ 20 ] = da3 * b3 * c1;
        N[ 21 ] = da3 * b3 * c2;
        N[ 22 ] = da0 * b3 * c1;
        N[ 23 ] = da0 * b3 * c2;
        N[ 24 ] = da1 * b0 * c3;
        N[ 25 ] = da2 * b0 * c3;
        N[ 26 ] = da0 * b1 * c3;
        N[ 27 ] = da0 * b2 * c3;
        N[ 28 ] = da3 * b1 * c3;
        N[ 29 ] = da3 * b2 * c3;
        N[ 30 ] = da2 * b3 * c3;
        N[ 31 ] = da1 * b3 * c3;
        N[ 32 ] = da1 * b1 * c0;
        N[ 33 ] = da1 * b2 * c0;
        N[ 34 ] = da2 * b2 * c0;
        N[ 35 ] = da2 * b1 * c0;
        N[ 36 ] = da1 * b0 * c1;
        N[ 37 ] = da2 * b0 * c1;
        N[ 38 ] = da2 * b0 * c2;
        N[ 39 ] = da1 * b0 * c2;
        N[ 40 ] = da0 * b1 * c1;
        N[ 41 ] = da0 * b1 * c2;
        N[ 42 ] = da0 * b2 * c2;
        N[ 43 ] = da0 * b2 * c1;
        N[ 44 ] = da3 * b1 * c1;
        N[ 45 ] = da3 * b2 * c1;
        N[ 46 ] = da3 * b2 * c2;
        N[ 47 ] = da3 * b1 * c2;
        N[ 48 ] = da2 * b3 * c1;
        N[ 49 ] = da1 * b3 * c1;
        N[ 50 ] = da1 * b3 * c2;
        N[ 51 ] = da2 * b3 * c2;
        N[ 52 ] = da1 * b1 * c3;
        N[ 53 ] = da2 * b1 * c3;
        N[ 54 ] = da2 * b2 * c3;
        N[ 55 ] = da1 * b2 * c3;
        N[ 56 ] = da1 * b1 * c1;
        N[ 57 ] = da2 * b1 * c1;
        N[ 58 ] = da2 * b2 * c1;
        N[ 59 ] = da1 * b2 * c1;
        N[ 60 ] = da1 * b1 * c2;
        N[ 61 ] = da2 * b1 * c2;
        N[ 62 ] = da2 * b2 * c2;
        N[ 63 ] = da1 * b2 * c2;
    }

    inline void
    Database::deval_hex64deta( const real xi, const real eta, const real zeta, real * N ) const
    {
        // often used parameters
        const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
        const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
        const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
        const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

        const real c0 =  ( zeta*( 1.0 + 9.0 * zeta * ( 1.0 - zeta ) ) - 1.0 )*0.0625;
        const real c1 =  ( 9.0 - zeta * ( 27.0 + zeta*( 9.0 - 27.0*zeta ) ) )*0.0625;
        const real c2 =  ( 9.0 + zeta * ( 27.0 - zeta*( 9.0 + 27.0*zeta ) ) )*0.0625;
        const real c3 = ( -zeta*( 1.0 - 9.0 * zeta * ( 1.0 + zeta ) ) - 1.0 )*0.0625;

        const real db0 = (   1.0 + eta*( 18.0 - 27.0*eta )) * 0.0625;
        const real db1 = ( -27.0 - eta*( 18.0 - 81.0*eta )) * 0.0625;
        const real db2 = (  27.0 - eta*( 18.0 + 81.0*eta )) * 0.0625;
        const real db3 = (  -1.0 + eta*( 18.0 + 27.0*eta )) * 0.0625;

        N[ 0 ] = a0*c0*db0;
        N[ 1 ] = a3*c0*db0;
        N[ 2 ] = a3*c0*db3;
        N[ 3 ] = a0*c0*db3;
        N[ 4 ] = a0*c3*db0;
        N[ 5 ] = a3*c3*db0;
        N[ 6 ] = a3*c3*db3;
        N[ 7 ] = a0*c3*db3;
        N[ 8 ] = a1*c0*db0;
        N[ 9 ] = a2*c0*db0;
        N[ 10 ] = a0*c0*db1;
        N[ 11 ] = a0*c0*db2;
        N[ 12 ] = a0*c1*db0;
        N[ 13 ] = a0*c2*db0;
        N[ 14 ] = a3*c0*db1;
        N[ 15 ] = a3*c0*db2;
        N[ 16 ] = a3*c1*db0;
        N[ 17 ] = a3*c2*db0;
        N[ 18 ] = a2*c0*db3;
        N[ 19 ] = a1*c0*db3;
        N[ 20 ] = a3*c1*db3;
        N[ 21 ] = a3*c2*db3;
        N[ 22 ] = a0*c1*db3;
        N[ 23 ] = a0*c2*db3;
        N[ 24 ] = a1*c3*db0;
        N[ 25 ] = a2*c3*db0;
        N[ 26 ] = a0*c3*db1;
        N[ 27 ] = a0*c3*db2;
        N[ 28 ] = a3*c3*db1;
        N[ 29 ] = a3*c3*db2;
        N[ 30 ] = a2*c3*db3;
        N[ 31 ] = a1*c3*db3;
        N[ 32 ] = a1*c0*db1;
        N[ 33 ] = a1*c0*db2;
        N[ 34 ] = a2*c0*db2;
        N[ 35 ] = a2*c0*db1;
        N[ 36 ] = a1*c1*db0;
        N[ 37 ] = a2*c1*db0;
        N[ 38 ] = a2*c2*db0;
        N[ 39 ] = a1*c2*db0;
        N[ 40 ] = a0*c1*db1;
        N[ 41 ] = a0*c2*db1;
        N[ 42 ] = a0*c2*db2;
        N[ 43 ] = a0*c1*db2;
        N[ 44 ] = a3*c1*db1;
        N[ 45 ] = a3*c1*db2;
        N[ 46 ] = a3*c2*db2;
        N[ 47 ] = a3*c2*db1;
        N[ 48 ] = a2*c1*db3;
        N[ 49 ] = a1*c1*db3;
        N[ 50 ] = a1*c2*db3;
        N[ 51 ] = a2*c2*db3;
        N[ 52 ] = a1*c3*db1;
        N[ 53 ] = a2*c3*db1;
        N[ 54 ] = a2*c3*db2;
        N[ 55 ] = a1*c3*db2;
        N[ 56 ] = a1*c1*db1;
        N[ 57 ] = a2*c1*db1;
        N[ 58 ] = a2*c1*db2;
        N[ 59 ] = a1*c1*db2;
        N[ 60 ] = a1*c2*db1;
        N[ 61 ] = a2*c2*db1;
        N[ 62 ] = a2*c2*db2;
        N[ 63 ] = a1*c2*db2;
    }

    inline void
    Database::deval_hex64dzeta( const real xi, const real eta, const real zeta, real * N ) const
    {
        // often used parameters
        const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
        const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
        const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
        const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

        const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 ) * 0.0625;
        const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) ) * 0.0625;
        const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) ) * 0.0625;
        const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 ) * 0.0625;

        const real dc0 = (   1.0 + zeta*( 18.0 - 27.0*zeta )) * 0.0625;
        const real dc1 = ( -27.0 - zeta*( 18.0 - 81.0*zeta )) * 0.0625;
        const real dc2 = (  27.0 - zeta*( 18.0 + 81.0*zeta )) * 0.0625;
        const real dc3 = (  -1.0 + zeta*( 18.0 + 27.0*zeta )) * 0.0625;

        N[ 0 ] = a0*b0*dc0;
        N[ 1 ] = a3*b0*dc0;
        N[ 2 ] = a3*b3*dc0;
        N[ 3 ] = a0*b3*dc0;
        N[ 4 ] = a0*b0*dc3;
        N[ 5 ] = a3*b0*dc3;
        N[ 6 ] = a3*b3*dc3;
        N[ 7 ] = a0*b3*dc3;
        N[ 8 ] = a1*b0*dc0;
        N[ 9 ] = a2*b0*dc0;
        N[ 10 ] = a0*b1*dc0;
        N[ 11 ] = a0*b2*dc0;
        N[ 12 ] = a0*b0*dc1;
        N[ 13 ] = a0*b0*dc2;
        N[ 14 ] = a3*b1*dc0;
        N[ 15 ] = a3*b2*dc0;
        N[ 16 ] = a3*b0*dc1;
        N[ 17 ] = a3*b0*dc2;
        N[ 18 ] = a2*b3*dc0;
        N[ 19 ] = a1*b3*dc0;
        N[ 20 ] = a3*b3*dc1;
        N[ 21 ] = a3*b3*dc2;
        N[ 22 ] = a0*b3*dc1;
        N[ 23 ] = a0*b3*dc2;
        N[ 24 ] = a1*b0*dc3;
        N[ 25 ] = a2*b0*dc3;
        N[ 26 ] = a0*b1*dc3;
        N[ 27 ] = a0*b2*dc3;
        N[ 28 ] = a3*b1*dc3;
        N[ 29 ] = a3*b2*dc3;
        N[ 30 ] = a2*b3*dc3;
        N[ 31 ] = a1*b3*dc3;
        N[ 32 ] = a1*b1*dc0;
        N[ 33 ] = a1*b2*dc0;
        N[ 34 ] = a2*b2*dc0;
        N[ 35 ] = a2*b1*dc0;
        N[ 36 ] = a1*b0*dc1;
        N[ 37 ] = a2*b0*dc1;
        N[ 38 ] = a2*b0*dc2;
        N[ 39 ] = a1*b0*dc2;
        N[ 40 ] = a0*b1*dc1;
        N[ 41 ] = a0*b1*dc2;
        N[ 42 ] = a0*b2*dc2;
        N[ 43 ] = a0*b2*dc1;
        N[ 44 ] = a3*b1*dc1;
        N[ 45 ] = a3*b2*dc1;
        N[ 46 ] = a3*b2*dc2;
        N[ 47 ] = a3*b1*dc2;
        N[ 48 ] = a2*b3*dc1;
        N[ 49 ] = a1*b3*dc1;
        N[ 50 ] = a1*b3*dc2;
        N[ 51 ] = a2*b3*dc2;
        N[ 52 ] = a1*b1*dc3;
        N[ 53 ] = a2*b1*dc3;
        N[ 54 ] = a2*b2*dc3;
        N[ 55 ] = a1*b2*dc3;
        N[ 56 ] = a1*b1*dc1;
        N[ 57 ] = a2*b1*dc1;
        N[ 58 ] = a2*b2*dc1;
        N[ 59 ] = a1*b2*dc1;
        N[ 60 ] = a1*b1*dc2;
        N[ 61 ] = a2*b1*dc2;
        N[ 62 ] = a2*b2*dc2;
        N[ 63 ] = a1*b2*dc2;
    }


}
#endif //BELFEM_CL_DATABASE_HPP