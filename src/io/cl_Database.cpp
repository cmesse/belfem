//
// Created by christian on 10/8/22.
//
#include "commtools.hpp"
#include "cl_Database.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    Database::Database( const std::string & aFilename,
                        const std::string & aTablename )
    {
        if( comm_rank() == 0 )
        {
            this->load_from_file( aFilename, aTablename );
        }

        // send data from master to others
        this->distribute_data() ;

        // link the shape function
        this->link_interpolation_function();
    }

//----------------------------------------------------------------------------

    // the database constructor (passing database)
    Database::Database( HDF5 & aFile ,
                        const std::string & aTablename )
    {
        if( comm_rank() == 0 )
        {
            this->load_from_file(aFile, aTablename);
        }

        // send data from master to others
        this->distribute_data() ;

        // link the shape function
        this->link_interpolation_function();
    }

//----------------------------------------------------------------------------

    Database::~Database()
    {
        this->free();
    }

//----------------------------------------------------------------------------

    void
    Database::free()
    {
        if( mNumPoints != nullptr )
        {
            ::free( mNumPoints );
            mNumPoints = nullptr ;
        }
        if( mXmax != nullptr )
        {
            ::free( mXmax );
            mXmax = nullptr ;
        }
        if( mXmin != nullptr )
        {
            ::free( mXmin );
            mXmin = nullptr ;
        }
        if( mX != nullptr )
        {
            ::free( mX );
            mX = nullptr ;
        }
        if ( mStep != nullptr )
        {
            ::free( mStep );
            mStep = nullptr ;
        }
        if ( mInvStep != nullptr )
        {
            ::free( mInvStep );
            mInvStep = nullptr ;
        }
        if ( mValues != nullptr )
        {
            ::free( mValues );
            mValues = nullptr ;
        }
        if( mWork != nullptr )
        {
            ::free( mWork );
            mWork = nullptr ;
        }
        if( mCoeffs != nullptr )
        {
            ::free( mCoeffs );
            mCoeffs = nullptr ;
        }

        mInterpolationOrder = 0 ;
        mNumberOfDimensions = 0 ;
        mMemorySize = 0 ;
        mNumberOfPointsPerCell = 0 ;
    }

//----------------------------------------------------------------------------

    void
    Database::load_from_file(
            const std::string & aFilename,
            const std::string & aTablename )
    {
        // open the Database
        HDF5 tFile( aFilename, FileMode::OPEN_RDONLY );

        this->load_from_file( tFile, aTablename );
        tFile.close();
    }

//----------------------------------------------------------------------------

    void
    Database::load_from_file(
            HDF5 & aFile,
            const std::string & aTablename )
    {
        // reset data
        this->free();

        // select the group
        aFile.select_group( aTablename );

        // get the number of dimensions
        aFile.load_data( "dimension", mNumberOfDimensions );

        // get the interpolation order
        aFile.load_data( "order", mInterpolationOrder );

        // order, but as real
        mOrder = real( mInterpolationOrder );

        // compute the inverse of the order
        mInvOrder = 1.0 / mOrder ;

        // read the number of samples
        mNumPoints = ( uint * ) malloc( mNumberOfDimensions * sizeof ( uint  ) ) ;
        aFile.load_data( "numpoints", mNumPoints, mNumberOfDimensions );

        // read the stepsize
        mStep =  ( real * ) malloc( mNumberOfDimensions * sizeof ( real ) ) ;
        aFile.load_data( "step", mStep, mNumberOfDimensions );

        // inverse stepsize
        mInvStep = ( real * ) malloc( mNumberOfDimensions * sizeof ( real ) ) ;
        for( uint k = 0; k<mNumberOfDimensions; ++k )
        {
            mInvStep[ k ] = 1.0 / mStep[ k ];
        }

        // read the offset
        mXmin = ( real * ) malloc( mNumberOfDimensions * sizeof ( real ) ) ;
        aFile.load_data( "offset", mXmin, mNumberOfDimensions );

        // compute the maximum values
        mXmax = ( real * ) malloc( mNumberOfDimensions * sizeof ( real ) ) ;
        for( uint d=0; d<mNumberOfDimensions; ++d )
        {
            mXmax[ d ] = mXmin[ d ] + mStep[ d ] *  real( mNumPoints[ d ] - 1 ) * mInvOrder ;
        }

        // allocate the work vector for the coordinates
        mX = ( real * ) malloc( mNumberOfDimensions * sizeof ( real ) ) ;

        for( uint d=0; d<mNumberOfDimensions; ++d )
        {
            mMemorySize *= mNumPoints[ d ];
        }

        // compute the memory size
        mMemorySize = 1 ;
        for( uint d=0; d<mNumberOfDimensions; ++d )
        {
            mMemorySize *= mNumPoints[ d ];
        }

        if( mOrder == 0 )
        {
            // read the precomputed coefficients
            std::size_t tSize = 4 * mMemorySize ;
            mCoeffs = ( real * ) malloc( tSize * sizeof ( real  ) ) ;
            aFile.load_data( "coeffs", mCoeffs, tSize );
        }
        else
        {
            // read the data
            mValues = ( real * ) malloc(mMemorySize * sizeof( real ) );
            aFile.load_data("values", mValues, mMemorySize );
        }
        aFile.close_active_group();
    }

//----------------------------------------------------------------------------

    void
    Database::distribute_data()
    {
        if( comm_size() == 1 )
        {
            return;
        }

        // wait for other procs
        comm_barrier() ;

        // prepare info tag
        Vector< uint > tIntInfo( 8, 0 );
        Vector< real > tRealInfo;
        if( comm_rank() == 0 )
        {
            uint tCount = 0 ;

            tIntInfo( tCount++ ) = mInterpolationOrder ;
            tIntInfo( tCount++ ) = mNumberOfDimensions ;
            tIntInfo( tCount++ ) = mMemorySize ;

            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                tIntInfo( tCount++ ) = mNumPoints[ k ];
            }

            distribute( tIntInfo );

            tRealInfo.set_size( 3 * mNumberOfDimensions );

            // reset counter
            tCount = 0 ;

            // send min and max data
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                tRealInfo( tCount++ ) = mXmin[ k ];
            }
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                tRealInfo( tCount++ ) = mXmax[ k ];
            }
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                tRealInfo( tCount++ ) = mStep[ k ];
            }


            distribute( tRealInfo );

            // catch special case
            if( mOrder == 0 )
            {
                tRealInfo.set_size( 4 * mMemorySize );
                std::copy( mCoeffs, mCoeffs + 4 * mMemorySize, tRealInfo.begin() );
            }
            else
            {
                // prepare data
                tRealInfo.set_size( mMemorySize );
                std::copy(mValues, mValues + mMemorySize, tRealInfo.begin() );
            }

            distribute( tRealInfo );
        }
        else
        {
            uint tCount = 0 ;


            distribute( tIntInfo );

            // unpack info
            mInterpolationOrder    = tIntInfo( tCount++ );
            mNumberOfDimensions    = tIntInfo( tCount++ );
            mMemorySize            = tIntInfo( tCount++ ) ;

            mNumPoints = ( uint * ) malloc( mNumberOfDimensions * sizeof ( uint ) ) ;
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mNumPoints[ k ] = tIntInfo( tCount++ ) ;
            }

            // order, but as real
            mOrder = real( mInterpolationOrder );

            // compute the inverse of the order
            mInvOrder = 1.0 / mOrder ;

            // receive real values
            distribute( tRealInfo );

            // allocate data
            mXmin = ( real * ) malloc( mNumberOfDimensions * sizeof ( real  ) ) ;
            mXmax = ( real * ) malloc( mNumberOfDimensions * sizeof ( real  ) ) ;
            mStep = ( real * ) malloc( mNumberOfDimensions * sizeof ( real  ) ) ;
            mInvStep = ( real * ) malloc( mNumberOfDimensions * sizeof ( real  ) ) ;

            // allocate work array
            mX = ( real * ) malloc( mNumberOfDimensions * sizeof ( real ) ) ;

            // populate data
            tCount = 0 ;
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mXmin[ k ] = tRealInfo( tCount++ );
            }
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mXmax[ k ] = tRealInfo( tCount++ );
            }
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mStep[ k ] = tRealInfo( tCount++ );
                mInvStep[ k ] = 1.0 / mStep[ k ];
            }

            distribute( tRealInfo );
            if( mOrder == 0 )
            {
                // allocate memory for coefficients
                mCoeffs = ( real * ) malloc(4 * mMemorySize * sizeof( real ) );

                // copy data from vector
                std::copy(tRealInfo.begin(), tRealInfo.end(), mCoeffs );
            }
            else
            {
                // allocate value memory
                mValues = ( real * ) malloc(mMemorySize * sizeof( real ) );

                // copy data from vector
                std::copy( tRealInfo.begin(), tRealInfo.end(), mValues );
            }
        }

        // wait for other procs
        comm_barrier() ;
    }


//----------------------------------------------------------------------------

    double
    Database::compute( const real aX )
    {
        mX[ 0 ] = aX ;
        return ( this->*mInterpolate )( mX );
    }

//----------------------------------------------------------------------------

    double
    Database::compute( const double aX, const double aY )
    {
        mX[ 0 ] = aX ;
        mX[ 1 ] = aY ;
        return ( this->*mInterpolate )( mX );
    }

//----------------------------------------------------------------------------

    double
    Database::compute( const double aX, const double aY, const double aZ )
    {
        mX[ 0 ] = aX ;
        mX[ 1 ] = aY ;
        mX[ 2 ] = aZ ;
        return ( this->*mInterpolate )( mX );
    }

//----------------------------------------------------------------------------

    void
    Database::link_interpolation_function()
    {
        switch( mNumberOfDimensions )
        {
            case( 1 ) :
            {
                switch ( mInterpolationOrder )
                {
                    case ( 0 ) :
                    {
                        mInterpolate = & Database::interpolate_1dspline ;
                        mNumberOfPointsPerCell = 0 ;
                        break ;
                    }
                    case ( 1 ) :
                    {
                        mInterpolate = &Database::interpolate_line2;
                        mNumberOfPointsPerCell = 2;
                        break;
                    }
                    case ( 2 ) :
                    {
                        mInterpolate = &Database::interpolate_line3;
                        mNumberOfPointsPerCell = 3;
                        break;
                    }
                    case ( 3 ) :
                    {
                        mInterpolate = &Database::interpolate_line4;
                        mNumberOfPointsPerCell = 4;
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false , "Invalid interpolation order: %u",
                                          (unsigned int) mInterpolationOrder);
                    }
                }
                break ;
            }
            case( 2 ) :
            {
                switch ( mInterpolationOrder )
                {
                    case( 1 ) :
                    {
                        mInterpolate = & Database::interpolate_quad4 ;
                        mNumberOfPointsPerCell = 4 ;
                        break ;
                    }
                    case( 2 ) :
                    {
                        mInterpolate = & Database::interpolate_quad9 ;
                        mNumberOfPointsPerCell = 9 ;
                        break ;
                    }
                    case( 3 ) :
                    {
                        mInterpolate = & Database::interpolate_quad16 ;
                        mNumberOfPointsPerCell = 16 ;
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false , "Invalid interpolation order: %u",
                                      (unsigned int) mInterpolationOrder);
                    }
                }
                break ;
            }
            case( 3 ) :
            {
                switch ( mInterpolationOrder )
                {
                    case( 1 ) :
                    {
                        mInterpolate = &  Database::interpolate_hex8 ;
                        mNumberOfPointsPerCell = 8 ;
                        break ;
                    }
                    case( 2 ) :
                    {
                        mInterpolate = & Database::interpolate_hex27 ;
                        mNumberOfPointsPerCell = 27 ;
                        break ;
                    }
                    case( 3 ) :
                    {
                        mInterpolate = & Database::interpolate_hex64 ;
                        mNumberOfPointsPerCell = 64 ;
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false , "Invalid interpolation order: %u",
                                      (unsigned int) mInterpolationOrder );
                    }
                }
                break ;
            }
            default :
            {
                BELFEM_ERROR( false , "Invalid interpolation order: %u",
                              (unsigned int) mInterpolationOrder);
            }
        }

        if( mNumberOfPointsPerCell > 0 )
        {
            // allocate work vector
            mWork = (double *)
                    malloc(mNumberOfPointsPerCell * sizeof(double));
        }
    }
    
//----------------------------------------------------------------------------

    real
    Database::interpolate_line2( const real * x )
    {
        // identify the indices of the first point of the cell
        const uint i = this->find_cell_linear( x, 0 );

        // compute parameter coordinates
        const real xi  = this->xtoxi_linear( x, i, 0 );

        // compute the interpolation
        return 0.5 * (  mValues[ i ] * ( 1.0 - xi )
                        + mValues[ i + 1 ] * ( 1.0 + xi ) );
    }

//----------------------------------------------------------------------------

    real
    Database::interpolate_line3( const real * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = this->find_cell_higher_order( x, 0 );

        // step 2: collect the node data
        mWork[ 0 ] = mValues[ i ];
        mWork[ 1 ] = mValues[ i+2 ];
        mWork[ 2 ] = mValues[ i+1 ];

        // step 3: compute parameter coordinate
        const real xi  = this->xtoxi_higher_order( x, i, 0 );
        real xi2 = xi*xi ;

        // step 4: compute the interpolation
        mWork[ 0 ] *=  xi2 - xi ;
        mWork[ 1 ] *= xi2 + xi  ;
        mWork[ 2 ] *= 1.0 - xi2;

        return 0.5 * ( mWork[ 0 ] +  mWork[ 1 ] ) +  mWork[ 2 ] ;
    }

//----------------------------------------------------------------------------

    real
    Database::interpolate_line4( const real * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = this->find_cell_higher_order( x, 0 );

        // step : compute parameter coordinate
        const real xi  = this->xtoxi_higher_order( x, i, 0 );

        // compute helpers
        mWork[ 0 ] = 1.0 - xi ;
        mWork[ 1 ] = 1.0 + xi ;
        mWork[ 2 ] = 3.0 * xi - 1.0 ;
        mWork[ 3 ] = 3.0 * xi + 1.0 ;

        // compute interpolation
        return  (   mWork[ 2 ] * mWork[ 3 ] *
                    ( mWork[ 0 ] * mValues[ i ] + mWork[ 1 ] * mValues[ i+3 ] )
                    + 9.0 * mWork[ 0 ] * mWork[ 1 ] * ( mWork[ 3 ] * mValues[ i+2 ]
                                                        - mWork[ 2 ] * mValues[ i+1 ] ) ) *  0.06250 ;
    }

//----------------------------------------------------------------------------

    real
    Database::interpolate_quad4( const real * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = this->find_cell_linear( x, 0 );
        const uint j = this->find_cell_linear( x, 1 );

        // step 2: collect the node data
        mWork[ 0 ] = mValues[ this->index( i, j ) ];
        mWork[ 1 ] = mValues[ this->index( i+1, j ) ];
        mWork[ 2 ] = mValues[ this->index( i+1, j+1 ) ];
        mWork[ 3 ] = mValues[ this->index( i, j+1 ) ];

        // step 3: compute parameter coordinates
        const real xi  = this->xtoxi_linear( x, i, 0 );
        const real eta = this->xtoxi_linear( x, j, 1 );

        // step 4: compute the interpolation
        mWork[ 0 ] *= ( ( 1.0 - xi ) * ( 1.0 - eta ) ) ;
        mWork[ 1 ] *= ( ( 1.0 + xi ) * ( 1.0 - eta ) ) ;
        mWork[ 2 ] *= ( ( 1.0 + xi ) * ( 1.0 + eta ) ) ;
        mWork[ 3 ] *= ( ( 1.0 - xi ) * ( 1.0 + eta ) ) ;

        // compute and return the value
        return 0.25 * std::accumulate( mWork, mWork + 4, 0.0 ) ;
    }

//----------------------------------------------------------------------------

    real
    Database::interpolate_quad9( const real * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = this->find_cell_higher_order( x, 0 );
        const uint j = this->find_cell_higher_order( x, 1 );

        // step 2: collect the node data
        mWork[ 0 ] = mValues[ this->index( i  , j   ) ];
        mWork[ 1 ] = mValues[ this->index( i+2, j   ) ];
        mWork[ 2 ] = mValues[ this->index( i+2, j+2 ) ];
        mWork[ 3 ] = mValues[ this->index( i  , j+2 ) ];
        mWork[ 4 ] = mValues[ this->index( i+1, j   ) ];
        mWork[ 5 ] = mValues[ this->index( i+2, j+1 ) ];
        mWork[ 6 ] = mValues[ this->index( i+1, j+2 ) ];
        mWork[ 7 ] = mValues[ this->index( i  , j+1 ) ];
        mWork[ 8 ] = mValues[ this->index( i+1, j+1 ) ];

        // step 3: compute parameter coordinates
        const real xi  = this->xtoxi_quadratic( x, i, 0 );
        const real eta = this->xtoxi_quadratic( x, j, 1 );

        // step 4: compute the interpolation
        const real    c = xi * eta * 0.25;
        const real  xi2 = xi*xi;
        const real eta2 = eta*eta;
        mWork[ 0 ] *= ( c * ( eta - 1.0 ) * (xi - 1.0) );
        mWork[ 1 ] *= ( c * ( eta - 1.0 ) * (xi + 1.0) );
        mWork[ 2 ] *= ( c * ( eta + 1.0 ) * (xi + 1.0) );
        mWork[ 3 ] *= ( c * ( eta + 1.0 ) * (xi - 1.0) );
        mWork[ 4 ] *= ( eta * ( 1.0 - xi2 ) * ( eta - 1.0 ) ) * 0.5;
        mWork[ 5 ] *= ( xi * ( 1.0 - eta2)*( xi + 1.0 ) )*0.5;
        mWork[ 6 ] *= ( eta * (1.0 - xi2)*( eta + 1.0 ) )*0.5;
        mWork[ 7 ] *= ( xi*( 1.0 - eta2 )*( xi - 1.0 ) )*0.5;
        mWork[ 8 ] *= ( eta2 - 1.0 )*( xi2 - 1.0 );

        return std::accumulate( mWork, mWork + 9, 0.0 );
    }

//----------------------------------------------------------------------------

    real
    Database::interpolate_quad16( const real * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = this->find_cell_higher_order( x, 0 );
        const uint j = this->find_cell_higher_order( x, 1 );

        // step 2: collect the node data
        mWork[  0 ] = mValues[ this->index( i  , j   ) ];
        mWork[  1 ] = mValues[ this->index( i+3, j   ) ];
        mWork[  2 ] = mValues[ this->index( i+3, j+3 ) ];
        mWork[  3 ] = mValues[ this->index( i  , j+3 ) ];
        mWork[  4 ] = mValues[ this->index( i+1, j   ) ];
        mWork[  5 ] = mValues[ this->index( i+2, j   ) ];
        mWork[  6 ] = mValues[ this->index( i+3, j+1 ) ];
        mWork[  7 ] = mValues[ this->index( i+3, j+2 ) ];
        mWork[  8 ] = mValues[ this->index( i+2, j+3 ) ];
        mWork[  9 ] = mValues[ this->index( i+1, j+3 ) ];
        mWork[ 10 ] = mValues[ this->index( i  , j+2 ) ];
        mWork[ 11 ] = mValues[ this->index( i  , j+1 ) ];
        mWork[ 12 ] = mValues[ this->index( i+1, j+1 ) ];
        mWork[ 13 ] = mValues[ this->index( i+2, j+1 ) ];
        mWork[ 14 ] = mValues[ this->index( i+2, j+2 ) ];
        mWork[ 15 ] = mValues[ this->index( i+1, j+2 ) ];

        // step 3: compute parameter coordinates
        const real xi  = this->xtoxi_higher_order( x, i, 0 );
        const real eta = this->xtoxi_higher_order( x, j, 1 );

        // step 4: compute the interpolation
        const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
        const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
        const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
        const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

        const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
        const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
        const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
        const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

        mWork[ 0 ] *= a0*b0;
        mWork[ 1 ] *= a3*b0;
        mWork[ 2 ] *= a3*b3;
        mWork[ 3 ] *= a0*b3;
        mWork[ 4 ] *= a1*b0;
        mWork[ 5 ] *= a2*b0;
        mWork[ 6 ] *= a3*b1;
        mWork[ 7 ] *= a3*b2;
        mWork[ 8 ] *= a2*b3;
        mWork[ 9 ] *= a1*b3;
        mWork[ 10 ] *= a0*b2;
        mWork[ 11 ] *= a0*b1;
        mWork[ 12 ] *= a1*b1;
        mWork[ 13 ] *= a2*b1;
        mWork[ 14 ] *= a2*b2;
        mWork[ 15 ] *= a1*b2;

        return std::accumulate( mWork, mWork + 16, 0.0 );
    }


//----------------------------------------------------------------------------

    real
    Database::interpolate_hex8( const real * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = this->find_cell_linear( x, 0 );
        const uint j = this->find_cell_linear( x, 1 );
        const uint k = this->find_cell_linear( x, 2 );

        // step 2: collect the node data
        mWork[ 0 ] = mValues[ this->index( i, j, k ) ];
        mWork[ 1 ] = mValues[ this->index( i+1, j, k ) ];
        mWork[ 2 ] = mValues[ this->index( i+1, j+1, k ) ];
        mWork[ 3 ] = mValues[ this->index( i, j+1, k ) ];
        mWork[ 4 ] = mValues[ this->index( i, j, k+1 ) ];
        mWork[ 5 ] = mValues[ this->index( i+1, j, k+1 ) ];
        mWork[ 6 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
        mWork[ 7 ] = mValues[ this->index( i, j+1, k+1 ) ];

        // step 3: compute parameter coordinates
        const real  xi  = this->xtoxi_linear( x, i, 0 );
        const real  eta = this->xtoxi_linear( x, j, 1 );
        const real zeta = this->xtoxi_linear( x, k, 2 );

        // step 4: compute the interpolation
        mWork[ 0 ] *=  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[ 1 ] *=    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[ 2 ] *=  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[ 3 ] *=    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[ 4 ] *=    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        mWork[ 5 ] *=  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[ 6 ] *=    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[ 7 ] *=  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );

        // compute and return the value
        return 0.125 * std::accumulate( mWork, mWork + 8, 0.0 );
    }

//----------------------------------------------------------------------------

    real
    Database::interpolate_hex27( const real * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = this->find_cell_higher_order( x, 0 );
        const uint j = this->find_cell_higher_order( x, 1 );
        const uint k = this->find_cell_higher_order( x, 2 );

        // step 2: collect the node data
        mWork[  0 ] = mValues[ this->index( i  , j  , k   ) ];
        mWork[  1 ] = mValues[ this->index( i+2, j  , k   ) ];
        mWork[  2 ] = mValues[ this->index( i+2, j+2, k+2 ) ];
        mWork[  3 ] = mValues[ this->index( i  , j+2, k+2 ) ];
        mWork[  4 ] = mValues[ this->index( i  , j  , k   ) ];
        mWork[  5 ] = mValues[ this->index( i+2, j  , k   ) ];
        mWork[  6 ] = mValues[ this->index( i+2, j+2, k+2 ) ];
        mWork[  7 ] = mValues[ this->index( i  , j+2, k+2 ) ];
        mWork[  8 ] = mValues[ this->index( i+1, j  , k   ) ];
        mWork[  9 ] = mValues[ this->index( i+2, j+1, k+1 ) ];
        mWork[ 10 ] = mValues[ this->index( i+1, j+2, k+2 ) ];
        mWork[ 11 ] = mValues[ this->index( i  , j+1, k+1 ) ];
        mWork[ 12 ] = mValues[ this->index( i  , j  , k   ) ];
        mWork[ 13 ] = mValues[ this->index( i+2, j  , k   ) ];
        mWork[ 14 ] = mValues[ this->index( i+2, j+2, k+2 ) ];
        mWork[ 15 ] = mValues[ this->index( i  , j+2, k+2 ) ];
        mWork[ 16 ] = mValues[ this->index( i+1, j  , k   ) ];
        mWork[ 17 ] = mValues[ this->index( i+2, j+1, k+1 ) ];
        mWork[ 18 ] = mValues[ this->index( i+1, j+2, k+2 ) ];
        mWork[ 19 ] = mValues[ this->index( i  , j+1, k+1 ) ];
        mWork[ 20 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
        mWork[ 21 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
        mWork[ 22 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
        mWork[ 23 ] = mValues[ this->index( i  , j+1, k+1 ) ];
        mWork[ 24 ] = mValues[ this->index( i+2, j+1, k+1 ) ];
        mWork[ 25 ] = mValues[ this->index( i+1, j  , k   ) ];
        mWork[ 26 ] = mValues[ this->index( i+1, j+2, k+2 ) ];

        // step 3: compute parameter coordinates
        const real   xi = this->xtoxi_quadratic( x, i, 0 );
        const real  eta = this->xtoxi_quadratic( x, j, 1 );
        const real zeta = this->xtoxi_quadratic( x, k, 2 );

        // step 4: compute the interpolation

        const real   xi2 = xi*xi;
        const real  eta2 = eta*eta;
        const real zeta2 = zeta*zeta;

        const real a = -0.25 * eta * zeta;
        const real b = -0.25 * xi * zeta;
        const real c = -0.25 * xi * eta;
        const real d = 0.125 * xi * eta * zeta;

        mWork[  0 ] *=  d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[  1 ] *=  d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[  2 ] *=  d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[  3 ] *=  d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[  4 ] *=  d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        mWork[  5 ] *=  d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[  6 ] *=  d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[  7 ] *=  d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        mWork[  8 ] *=  a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
        mWork[  9 ] *=  b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[ 10 ] *=  a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
        mWork[ 11 ] *=  b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[ 12 ] *=  c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );
        mWork[ 13 ] *=  c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );
        mWork[ 14 ] *=  c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );
        mWork[ 15 ] *=  c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );
        mWork[ 16 ] *=  a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
        mWork[ 17 ] *=  b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[ 18 ] *=  a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
        mWork[ 19 ] *=  b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        mWork[ 20 ] *=  -( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
        mWork[ 21 ] *=  ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
        mWork[ 22 ] *=  ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
        mWork[ 23 ] *=  ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.5;
        mWork[ 24 ] *=  ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.5;
        mWork[ 25 ] *=  ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.5;
        mWork[ 26 ] *=  ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.5;

        return std::accumulate( mWork, mWork + 27, 0.0 );
    }
//----------------------------------------------------------------------------

    real
    Database::interpolate_hex64( const real * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = this->find_cell_higher_order( x, 0 );
        const uint j = this->find_cell_higher_order( x, 1 );
        const uint k = this->find_cell_higher_order( x, 2 );

        // step 2: collect the node data
        mWork[  0 ] = mValues[ this->index( i  , j  , k   ) ];
        mWork[  1 ] = mValues[ this->index( i+3, j  , k   ) ];
        mWork[  2 ] = mValues[ this->index( i+3, j+3, k   ) ];
        mWork[  3 ] = mValues[ this->index( i  , j+3, k   ) ];
        mWork[  4 ] = mValues[ this->index( i  , j  , k+3 ) ];
        mWork[  5 ] = mValues[ this->index( i+3, j  , k+3 ) ];
        mWork[  6 ] = mValues[ this->index( i+3, j+3, k+3 ) ];
        mWork[  7 ] = mValues[ this->index( i  , j+3, k+3 ) ];
        mWork[  8 ] = mValues[ this->index( i+1, j  , k   ) ];
        mWork[  9 ] = mValues[ this->index( i+2, j  , k   ) ];
        mWork[ 10 ] = mValues[ this->index( i  , j+1, k   ) ];
        mWork[ 11 ] = mValues[ this->index( i  , j+2, k   ) ];
        mWork[ 12 ] = mValues[ this->index( i  , j  , k+1 ) ];
        mWork[ 13 ] = mValues[ this->index( i  , j  , k+2 ) ];
        mWork[ 14 ] = mValues[ this->index( i+3, j+1, k   ) ];
        mWork[ 15 ] = mValues[ this->index( i+3, j+2, k   ) ];
        mWork[ 16 ] = mValues[ this->index( i+3, j  , k+1 ) ];
        mWork[ 17 ] = mValues[ this->index( i+3, j  , k+2 ) ];
        mWork[ 18 ] = mValues[ this->index( i+2, j+3, k   ) ];
        mWork[ 19 ] = mValues[ this->index( i+1, j+3, k   ) ];
        mWork[ 20 ] = mValues[ this->index( i+3, j+3, k+1 ) ];
        mWork[ 21 ] = mValues[ this->index( i+3, j+3, k+2 ) ];
        mWork[ 22 ] = mValues[ this->index( i  , j+3, k+1 ) ];
        mWork[ 23 ] = mValues[ this->index( i  , j+3, k+2 ) ];
        mWork[ 24 ] = mValues[ this->index( i+1, j  , k+3 ) ];
        mWork[ 25 ] = mValues[ this->index( i+2, j  , k+3 ) ];
        mWork[ 26 ] = mValues[ this->index( i  , j+1, k+3 ) ];
        mWork[ 27 ] = mValues[ this->index( i  , j+2, k+3 ) ];
        mWork[ 28 ] = mValues[ this->index( i+3, j+1, k+3 ) ];
        mWork[ 29 ] = mValues[ this->index( i+3, j+2, k+3 ) ];
        mWork[ 30 ] = mValues[ this->index( i+2, j+3, k+3 ) ];
        mWork[ 31 ] = mValues[ this->index( i+1, j+3, k+3 ) ];
        mWork[ 32 ] = mValues[ this->index( i+1, j+1, k   ) ];
        mWork[ 33 ] = mValues[ this->index( i+1, j+2, k   ) ];
        mWork[ 34 ] = mValues[ this->index( i+2, j+2, k   ) ];
        mWork[ 35 ] = mValues[ this->index( i+2, j+1, k   ) ];
        mWork[ 36 ] = mValues[ this->index( i+1, j  , k+1 ) ];
        mWork[ 37 ] = mValues[ this->index( i+2, j  , k+1 ) ];
        mWork[ 38 ] = mValues[ this->index( i+2, j  , k+2 ) ];
        mWork[ 39 ] = mValues[ this->index( i+1, j  , k+2 ) ];
        mWork[ 40 ] = mValues[ this->index( i  , j+1, k+1 ) ];
        mWork[ 41 ] = mValues[ this->index( i  , j+1, k+2 ) ];
        mWork[ 42 ] = mValues[ this->index( i  , j+2, k+2 ) ];
        mWork[ 43 ] = mValues[ this->index( i  , j+2, k+1 ) ];
        mWork[ 44 ] = mValues[ this->index( i+3, j+1, k+1 ) ];
        mWork[ 45 ] = mValues[ this->index( i+3, j+2, k+1 ) ];
        mWork[ 46 ] = mValues[ this->index( i+3, j+2, k+2 ) ];
        mWork[ 47 ] = mValues[ this->index( i+3, j+1, k+2 ) ];
        mWork[ 48 ] = mValues[ this->index( i+2, j+3, k+1 ) ];
        mWork[ 49 ] = mValues[ this->index( i+1, j+3, k+1 ) ];
        mWork[ 50 ] = mValues[ this->index( i+1, j+3, k+2 ) ];
        mWork[ 51 ] = mValues[ this->index( i+2, j+3, k+2 ) ];
        mWork[ 52 ] = mValues[ this->index( i+1, j+1, k+3 ) ];
        mWork[ 53 ] = mValues[ this->index( i+2, j+1, k+3 ) ];
        mWork[ 54 ] = mValues[ this->index( i+2, j+2, k+3 ) ];
        mWork[ 55 ] = mValues[ this->index( i+1, j+2, k+3 ) ];
        mWork[ 56 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
        mWork[ 57 ] = mValues[ this->index( i+2, j+1, k+1 ) ];
        mWork[ 58 ] = mValues[ this->index( i+2, j+2, k+1 ) ];
        mWork[ 59 ] = mValues[ this->index( i+1, j+2, k+1 ) ];
        mWork[ 60 ] = mValues[ this->index( i+1, j+1, k+2 ) ];
        mWork[ 61 ] = mValues[ this->index( i+2, j+1, k+2 ) ];
        mWork[ 62 ] = mValues[ this->index( i+2, j+2, k+2 ) ];
        mWork[ 63 ] = mValues[ this->index( i+1, j+2, k+2 ) ];

        // step 3: compute parameter coordinates
        const real  xi  = this->xtoxi_quadratic( x, i, 0 );
        const real  eta = this->xtoxi_quadratic( x, j, 1 );
        const real zeta = this->xtoxi_quadratic( x, k, 2 );

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

        mWork[  0 ] *= a0 * b0 * c0;
        mWork[  1 ] *= a3 * b0 * c0;
        mWork[  2 ] *= a3 * b3 * c0;
        mWork[  3 ] *= a0 * b3 * c0;
        mWork[  4 ] *= a0 * b0 * c3;
        mWork[  5 ] *= a3 * b0 * c3;
        mWork[  6 ] *= a3 * b3 * c3;
        mWork[  7 ] *= a0 * b3 * c3;
        mWork[  8 ] *= a1 * b0 * c0;
        mWork[  9 ] *= a2 * b0 * c0;
        mWork[ 10 ] *= a0 * b1 * c0;
        mWork[ 11 ] *= a0 * b2 * c0;
        mWork[ 12 ] *= a0 * b0 * c1;
        mWork[ 13 ] *= a0 * b0 * c2;
        mWork[ 14 ] *= a3 * b1 * c0;
        mWork[ 15 ] *= a3 * b2 * c0;
        mWork[ 16 ] *= a3 * b0 * c1;
        mWork[ 17 ] *= a3 * b0 * c2;
        mWork[ 18 ] *= a2 * b3 * c0;
        mWork[ 19 ] *= a1 * b3 * c0;
        mWork[ 20 ] *= a3 * b3 * c1;
        mWork[ 21 ] *= a3 * b3 * c2;
        mWork[ 22 ] *= a0 * b3 * c1;
        mWork[ 23 ] *= a0 * b3 * c2;
        mWork[ 24 ] *= a1 * b0 * c3;
        mWork[ 25 ] *= a2 * b0 * c3;
        mWork[ 26 ] *= a0 * b1 * c3;
        mWork[ 27 ] *= a0 * b2 * c3;
        mWork[ 28 ] *= a3 * b1 * c3;
        mWork[ 29 ] *= a3 * b2 * c3;
        mWork[ 30 ] *= a2 * b3 * c3;
        mWork[ 31 ] *= a1 * b3 * c3;
        mWork[ 32 ] *= a1 * b1 * c0;
        mWork[ 33 ] *= a1 * b2 * c0;
        mWork[ 34 ] *= a2 * b2 * c0;
        mWork[ 35 ] *= a2 * b1 * c0;
        mWork[ 36 ] *= a1 * b0 * c1;
        mWork[ 37 ] *= a2 * b0 * c1;
        mWork[ 38 ] *= a2 * b0 * c2;
        mWork[ 39 ] *= a1 * b0 * c2;
        mWork[ 40 ] *= a0 * b1 * c1;
        mWork[ 41 ] *= a0 * b1 * c2;
        mWork[ 42 ] *= a0 * b2 * c2;
        mWork[ 43 ] *= a0 * b2 * c1;
        mWork[ 44 ] *= a3 * b1 * c1;
        mWork[ 45 ] *= a3 * b2 * c1;
        mWork[ 46 ] *= a3 * b2 * c2;
        mWork[ 47 ] *= a3 * b1 * c2;
        mWork[ 48 ] *= a2 * b3 * c1;
        mWork[ 49 ] *= a1 * b3 * c1;
        mWork[ 50 ] *= a1 * b3 * c2;
        mWork[ 51 ] *= a2 * b3 * c2;
        mWork[ 52 ] *= a1 * b1 * c3;
        mWork[ 53 ] *= a2 * b1 * c3;
        mWork[ 54 ] *= a2 * b2 * c3;
        mWork[ 55 ] *= a1 * b2 * c3;
        mWork[ 56 ] *= a1 * b1 * c1;
        mWork[ 57 ] *= a2 * b1 * c1;
        mWork[ 58 ] *= a2 * b2 * c1;
        mWork[ 59 ] *= a1 * b2 * c1;
        mWork[ 60 ] *= a1 * b1 * c2;
        mWork[ 61 ] *= a2 * b1 * c2;
        mWork[ 62 ] *= a2 * b2 * c2;
        mWork[ 63 ] *= a1 * b2 * c2;

        return std::accumulate( mWork, mWork + 64, 0.0 );
    }
    
//----------------------------------------------------------------------------

}