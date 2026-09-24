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

#include <algorithm>
#include <cmath>
#include "commtools.hpp"
#include "cl_Database.hpp"
#include "cl_HDF5.hpp"
#include "assert.hpp"
#include "cl_Solver.hpp"
#include "cl_DatabaseProjector.hpp"

namespace belfem
{
    Database::Database( const string & aFilePath, const string & aMaterial ) :
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mLabel( aMaterial )
    {
        if ( comm_rank() == 0 )
        {
            HDF5 tFile( aFilePath, FileMode::OPEN_RDONLY );
            tFile.select_group( aMaterial );

            this->load( tFile.active_group() );
            tFile.close_active_group();
            tFile.close();
        }
        else
        {
            this->load( 0 );
        }

        this->set_element_type();
        TensorMeshFactory tFactory ;
        tFactory.create_topology( mConfig, mTopology );

    }

    Database::Database( const hid_t aHDF5, const string & aLabel ) :
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mLabel( aLabel )
    {
        this->load( aHDF5 );
        this->set_element_type();
        TensorMeshFactory tFactory ;
        tFactory.create_topology( mConfig, mTopology );
    }

    Database::Database( Mesh * aMesh, const string & aField, bool aProject, const string aMaterial ) :
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mLabel( aField )
    {
        uint tDimension ;
        uint tOrder ;
        Vector< real > tOrigin ;
        Vector< real > tNodeSteps ;
        Vector< index_t > tNumNodes ;

        if ( aProject )
        {
            database::Projector tProjector( aMesh );
            if ( aMaterial != "" ) tProjector.set_label( aMaterial );
            tProjector.project( aField, mValues );
        }
        else if ( mCommRank == 0 )
        {
            mValues = aMesh->field( aField )->data();
            share( mValues );
        }
        else
        {
            receive( mValues );
        }

        if ( mCommRank == 0 )
        {
            BELFEM_ERROR( aMesh->is_tensormesh(), "Database::Database( Mesh * aMesh, const string & aField ) expects a tensor mesh" );

            tDimension = aMesh->number_of_dimensions();
            tOrder     = aMesh->tensorconf()->order();
            tOrigin    = aMesh->tensorconf()->origin();
            tNodeSteps = aMesh->tensorconf()->step();
            tNumNodes  = aMesh->tensorconf()->num_nodes_vector();

            if ( mCommSize > 1 )
            {
                broadcast( tOrder );
                broadcast( tDimension );
                broadcast( tOrigin );
                broadcast( tNumNodes );
                broadcast( tNodeSteps );
            }
        }
        else
        {
            broadcast( tOrder );
            broadcast( tDimension );
            broadcast( tOrigin );
            broadcast( tNumNodes );
            broadcast( tNodeSteps );
        }

        mConfig = new TensorMeshConfig( tOrder, tNumNodes, tNodeSteps, tOrigin  );
        this->set_element_type();
        TensorMeshFactory tFactory ;
        tFactory.create_topology( mConfig, mTopology );
    }

    Database::~Database()
    {
        if ( mN != nullptr )
        {
            free( mN );
        }
        if ( mConfig != nullptr )
        {
            delete mConfig ;
        }
    }

    void Database::save( const hid_t aHDF5 )
    {
        if ( mCommRank != 0 ) return;

        hid_t tID = aHDF5 ;
        herr_t tStatus ;
        hdf5::save_scalar_to_file( tID, "order", mConfig->order(), tStatus );
        hdf5::save_vector_to_file( tID, "origin", mConfig->origin(), tStatus );
        hdf5::save_vector_to_file( tID, "points", mConfig->num_nodes_vector(), tStatus );
        hdf5::save_vector_to_file( tID, "step", mConfig->step(), tStatus );
        hdf5::save_vector_to_file( tID, "values", mValues, tStatus );
    }

    void
    Database::load( const hid_t aHDF5 )
    {
        uint tDimension ;
        uint32_t tOrder ;

        Vector< uint32_t > tNumPoints ;
        Vector< real >     tNodeSteps ;
        Vector< real >     tOrigin ;

        if ( mCommRank == 0 )
        {
            hid_t tID = aHDF5 ;
            herr_t tStatus ;
            hdf5::load_scalar_from_file( tID, "order", tOrder, tStatus );

            hdf5::load_vector_from_file( tID, "origin", tOrigin, tStatus );
            tDimension = tOrigin.length() ;

            hdf5::load_vector_from_file( tID, "points", tNumPoints, tStatus );
            hdf5::load_vector_from_file( tID, "step", tNodeSteps, tStatus );
            hdf5::load_vector_from_file( tID, "values", mValues, tStatus );

            if ( mCommSize > 1 )
            {
                broadcast( tOrder );
                broadcast( tDimension );
                broadcast( tOrigin );
                broadcast( tNumPoints );
                broadcast( tNodeSteps );
                share( mValues );
            }
        }
        else
        {
            broadcast( tOrder );
            broadcast( tDimension );
            broadcast( tOrigin );
            broadcast( tNumPoints );
            broadcast( tNodeSteps );
            receive( mValues );
        }

        mConfig = new TensorMeshConfig( tOrder, tNumPoints, tNodeSteps, tOrigin  );
    }

    void
    Database::set_element_type()
    {
        switch ( mConfig->element_type() )
        {
            case ElementType::QUAD4:
            {
                mFunction2D      = & Database::eval_quad4 ;
                mdFunction2Ddxi  = & Database::deval_quad4dxi ;
                mdFunction2Ddeta = & Database::deval_quad4deta ;
                mFunction3D      = nullptr ;
                mdFunction3Ddxi  = nullptr ;
                break ;
            }
            case ElementType::QUAD9:
            {
                mFunction2D      = & Database::eval_quad9 ;
                mdFunction2Ddxi  = & Database::deval_quad9dxi ;
                mdFunction2Ddeta = & Database::deval_quad9deta ;
                mFunction3D      = nullptr ;
                mdFunction3Ddxi  = nullptr ;
                break ;
            }
            case ElementType::QUAD16:
            {
                mFunction2D      = & Database::eval_quad16 ;
                mdFunction2Ddxi  = & Database::deval_quad16dxi ;
                mdFunction2Ddeta = & Database::deval_quad16deta ;
                mFunction3D      = nullptr ;
                mdFunction3Ddxi  = nullptr ;
                break ;
            }
            case ElementType::HEX8:
            {
                mFunction2D       = nullptr ;
                mdFunction2Ddxi   = nullptr ;
                mFunction3D       = & Database::eval_hex8 ;
                mdFunction3Ddxi   = & Database::deval_hex8dxi ;
                mdFunction3Ddeta  = & Database::deval_hex8deta ;
                mdFunction3Ddzeta = & Database::deval_hex8dzeta ;
                break ;
            }
            case ElementType::HEX27:
            {
                mFunction2D       = nullptr ;
                mdFunction2Ddxi   = nullptr ;
                mFunction3D       = & Database::eval_hex27 ;
                mdFunction3Ddxi   = & Database::deval_hex27dxi ;
                mdFunction3Ddeta  = & Database::deval_hex27deta ;
                mdFunction3Ddzeta = & Database::deval_hex27dzeta ;
                break ;
            }
            case ElementType::HEX64:
            {
                mFunction2D       = nullptr ;
                mdFunction2Ddxi   = nullptr ;
                mFunction3D       = & Database::eval_hex64 ;
                mdFunction3Ddxi   = & Database::deval_hex64dxi ;
                mdFunction3Ddeta  = & Database::deval_hex64deta ;
                mdFunction3Ddzeta = & Database::deval_hex64dzeta ;
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "unsupported element type in database %s", mLabel.c_str() );
            }
        }

        mNumNodesPerElement = mesh::number_of_nodes( mConfig->element_type() );

        if ( mN != nullptr ) free( mN );
        mN = ( real * ) malloc( mNumNodesPerElement * sizeof( real ) );

        if ( mConfig->dimension() == 3 )
        {
            this->build_tensor_map();
        }
    }

    uint
    Database::dimension() const
    {
        return mConfig->dimension() ;
    }

    real
    Database::lagrange_1d( const uint aOrder, const uint a, const real t )
    {
        const real tA = -1.0 + 2.0 * a / aOrder ;
        real tValue = 1.0 ;
        for ( uint b = 0; b <= aOrder; ++b )
        {
            if ( b == a ) continue ;
            const real tB = -1.0 + 2.0 * b / aOrder ;
            tValue *= ( t - tB ) / ( tA - tB ) ;
        }
        return tValue ;
    }

    real
    Database::dlagrange_1d( const uint aOrder, const uint a, const real t )
    {
        const real tA = -1.0 + 2.0 * a / aOrder ;
        real tSum = 0.0 ;
        for ( uint c = 0; c <= aOrder; ++c )
        {
            if ( c == a ) continue ;
            const real tC = -1.0 + 2.0 * c / aOrder ;
            real tProduct = 1.0 / ( tA - tC ) ;
            for ( uint b = 0; b <= aOrder; ++b )
            {
                if ( b == a || b == c ) continue ;
                const real tB = -1.0 + 2.0 * b / aOrder ;
                tProduct *= ( t - tB ) / ( tA - tB ) ;
            }
            tSum += tProduct ;
        }
        return tSum ;
    }

    void
    Database::build_tensor_map()
    {
        const uint p = mConfig->order() ;
        BELFEM_ERROR( p >= 1 && p <= 3, "unsupported order %u in database %s",
                      ( unsigned int ) p, mLabel.c_str() );

        mTensorIndexX.set_size( mNumNodesPerElement, BELFEM_UINT_MAX );
        mTensorIndexY.set_size( mNumNodesPerElement, BELFEM_UINT_MAX );
        mTensorIndexZ.set_size( mNumNodesPerElement, BELFEM_UINT_MAX );

        // the node of a tensor grid point is the one whose shape function is
        // 1 there; every other one is 0
        for ( uint c = 0; c <= p; ++c )
        {
            for ( uint b = 0; b <= p; ++b )
            {
                for ( uint a = 0; a <= p; ++a )
                {
                    ( this->*mFunction3D )( -1.0 + 2.0 * a / p,
                                            -1.0 + 2.0 * b / p,
                                            -1.0 + 2.0 * c / p, mN );
                    for ( uint n = 0; n < mNumNodesPerElement; ++n )
                    {
                        if ( std::abs( mN[ n ] - 1.0 ) < 1e-9 )
                        {
                            mTensorIndexX( n ) = a ;
                            mTensorIndexY( n ) = b ;
                            mTensorIndexZ( n ) = c ;
                        }
                    }
                }
            }
        }
        for ( uint n = 0; n < mNumNodesPerElement; ++n )
        {
            BELFEM_ERROR( mTensorIndexX( n ) != BELFEM_UINT_MAX,
                          "node %u of the database element has no tensor position",
                          ( unsigned int ) n );
        }

        // Bernstein coefficients from nodal values: invert A( i, k ) = B_k( s_i ),
        // s_i = i / p, by Gauss-Jordan elimination ( at most 4 x 4, once )
        real tA[ 4 ][ 8 ] ;
        for ( uint i = 0; i <= p; ++i )
        {
            const real s = static_cast< real >( i ) / p ;
            real tBinomial = 1.0 ;
            for ( uint k = 0; k <= p; ++k )
            {
                if ( k > 0 ) tBinomial = tBinomial * ( p - k + 1 ) / k ;
                tA[ i ][ k ] = tBinomial * std::pow( s, k ) * std::pow( 1.0 - s, p - k );
                tA[ i ][ p + 1 + k ] = ( i == k ) ? 1.0 : 0.0 ;
            }
        }
        for ( uint c = 0; c <= p; ++c )
        {
            uint tPivot = c ;
            for ( uint r = c + 1; r <= p; ++r )
            {
                if ( std::abs( tA[ r ][ c ] ) > std::abs( tA[ tPivot ][ c ] ) ) tPivot = r ;
            }
            for ( uint k = 0; k < 2 * ( p + 1 ); ++k ) std::swap( tA[ c ][ k ], tA[ tPivot ][ k ] );
            const real tDiag = tA[ c ][ c ] ;
            for ( uint k = 0; k < 2 * ( p + 1 ); ++k ) tA[ c ][ k ] /= tDiag ;
            for ( uint r = 0; r <= p; ++r )
            {
                if ( r == c ) continue ;
                const real tFactor = tA[ r ][ c ] ;
                for ( uint k = 0; k < 2 * ( p + 1 ); ++k ) tA[ r ][ k ] -= tFactor * tA[ c ][ k ] ;
            }
        }
        for ( uint i = 0; i < 4; ++i )
        {
            for ( uint k = 0; k < 4; ++k )
            {
                mToBernstein[ i ][ k ] = ( i <= p && k <= p ) ? tA[ i ][ p + 1 + k ] : 0.0 ;
            }
        }
    }

    index_t
    Database::locate( const real aX, const real aY, const real aZ,
                      real & aXi, real & aEta, real & aZeta ) const
    {
        const index_t i = mConfig->element_ijk( 0, aX );
        const index_t j = mConfig->element_ijk( 1, aY );
        const index_t k = mConfig->element_ijk( 2, aZ );

        const real xc = mConfig->min( 0 ) + ( i + 0.5 ) * mConfig->element_step( 0 );
        const real yc = mConfig->min( 1 ) + ( j + 0.5 ) * mConfig->element_step( 1 );
        const real zc = mConfig->min( 2 ) + ( k + 0.5 ) * mConfig->element_step( 2 );

        aXi   = 2. * ( aX - xc ) * mConfig->inv_element_step( 0 ) ;
        aEta  = 2. * ( aY - yc ) * mConfig->inv_element_step( 1 ) ;
        aZeta = 2. * ( aZ - zc ) * mConfig->inv_element_step( 2 ) ;

        return mConfig->element_index( i, j, k );
    }

    real
    Database::evaluate_derivxy( const real aX, const real aY, const real aZ ) const
    {
        BELFEM_ASSERT( mConfig->dimension() == 3, "evaluate_derivxy needs a 3-D database" );
        real xi, eta, zeta ;
        const index_t e = this->locate( aX, aY, aZ, xi, eta, zeta );
        const uint p = mConfig->order() ;

        // the 1-D factors per axis, once per call
        real tX[ 4 ], tY[ 4 ], tZ[ 4 ] ;
        for ( uint a = 0; a <= p; ++a )
        {
            tX[ a ] = dlagrange_1d( p, a, xi );
            tY[ a ] = dlagrange_1d( p, a, eta );
            tZ[ a ] = lagrange_1d( p, a, zeta );
        }

        real tValue = 0.0 ;
        for ( uint n = 0; n < mNumNodesPerElement; ++n )
        {
            tValue += mValues( mTopology( n, e ) )
                    * tX[ mTensorIndexX( n ) ] * tY[ mTensorIndexY( n ) ] * tZ[ mTensorIndexZ( n ) ] ;
        }
        return tValue * 4. * mConfig->inv_element_step( 0 ) * mConfig->inv_element_step( 1 ) ;
    }

    real
    Database::evaluate_derivyz( const real aX, const real aY, const real aZ ) const
    {
        BELFEM_ASSERT( mConfig->dimension() == 3, "evaluate_derivyz needs a 3-D database" );
        real xi, eta, zeta ;
        const index_t e = this->locate( aX, aY, aZ, xi, eta, zeta );
        const uint p = mConfig->order() ;

        // the 1-D factors per axis, once per call
        real tX[ 4 ], tY[ 4 ], tZ[ 4 ] ;
        for ( uint a = 0; a <= p; ++a )
        {
            tX[ a ] = lagrange_1d( p, a, xi );
            tY[ a ] = dlagrange_1d( p, a, eta );
            tZ[ a ] = dlagrange_1d( p, a, zeta );
        }

        real tValue = 0.0 ;
        for ( uint n = 0; n < mNumNodesPerElement; ++n )
        {
            tValue += mValues( mTopology( n, e ) )
                    * tX[ mTensorIndexX( n ) ] * tY[ mTensorIndexY( n ) ] * tZ[ mTensorIndexZ( n ) ] ;
        }
        return tValue * 4. * mConfig->inv_element_step( 1 ) * mConfig->inv_element_step( 2 ) ;
    }

    real
    Database::upper_bound_derivy_on_lower_y_face() const
    {
        BELFEM_ERROR( mConfig->dimension() == 3, "the face bound needs a 3-D database" );
        const uint p = mConfig->order() ;
        real tBound = -BELFEM_REAL_MAX ;

        for ( uint k = 0; k < mConfig->num_elements( 2 ); ++k )
        {
            for ( uint i = 0; i < mConfig->num_elements( 0 ); ++i )
            {
                const index_t e = mConfig->element_index( i, 0, k );

                // nodal values over ( x, z ) of the face polynomial of df/deta at
                // eta = -1: the element is Lagrange in eta, so this is exact
                real tC[ 4 ][ 4 ] = { { 0 } } ;
                for ( uint n = 0; n < mNumNodesPerElement; ++n )
                {
                    tC[ mTensorIndexX( n ) ][ mTensorIndexZ( n ) ] +=
                        mValues( mTopology( n, e ) ) * dlagrange_1d( p, mTensorIndexY( n ), -1.0 );
                }
                for ( uint a = 0; a <= p; ++a )
                {
                    for ( uint c = 0; c <= p; ++c )
                    {
                        real tB = 0.0 ;
                        for ( uint i2 = 0; i2 <= p; ++i2 )
                        {
                            for ( uint k2 = 0; k2 <= p; ++k2 )
                            {
                                tB += mToBernstein[ a ][ i2 ] * tC[ i2 ][ k2 ] * mToBernstein[ c ][ k2 ] ;
                            }
                        }
                        // a non-finite table has no certificate; NaN fails
                        // every comparison the caller makes with it
                        if ( ! std::isfinite( tB ) ) return BELFEM_QUIET_NAN ;
                        tBound = std::max( tBound, tB );
                    }
                }
            }
        }
        return tBound * 2. * mConfig->inv_element_step( 1 ) ;
    }

    real
    Database::lower_bound_on_lower_y_face() const
    {
        BELFEM_ERROR( mConfig->dimension() == 3, "the face bound needs a 3-D database" );
        const uint p = mConfig->order() ;
        real tBound = BELFEM_REAL_MAX ;

        for ( uint k = 0; k < mConfig->num_elements( 2 ); ++k )
        {
            for ( uint i = 0; i < mConfig->num_elements( 0 ); ++i )
            {
                const index_t e = mConfig->element_index( i, 0, k );

                real tC[ 4 ][ 4 ] = { { 0 } } ;
                for ( uint n = 0; n < mNumNodesPerElement; ++n )
                {
                    if ( mTensorIndexY( n ) != 0 ) continue ;
                    tC[ mTensorIndexX( n ) ][ mTensorIndexZ( n ) ] = mValues( mTopology( n, e ) );
                }
                for ( uint a = 0; a <= p; ++a )
                {
                    for ( uint c = 0; c <= p; ++c )
                    {
                        real tB = 0.0 ;
                        for ( uint i2 = 0; i2 <= p; ++i2 )
                        {
                            for ( uint k2 = 0; k2 <= p; ++k2 )
                            {
                                tB += mToBernstein[ a ][ i2 ] * tC[ i2 ][ k2 ] * mToBernstein[ c ][ k2 ] ;
                            }
                        }
                        if ( ! std::isfinite( tB ) ) return BELFEM_QUIET_NAN ;
                        tBound = std::min( tBound, tB );
                    }
                }
            }
        }
        return tBound ;
    }

    real
    Database::min( const uint aDimension ) const
    {
        return mConfig->min( aDimension );
    }

    real
    Database::max( const uint aDimension ) const
    {
        return mConfig->max( aDimension );
    }
}
