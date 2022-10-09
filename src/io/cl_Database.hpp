//
// Created by christian on 10/8/22.
//

#ifndef BELFEM_CL_DATABASE_HPP
#define BELFEM_CL_DATABASE_HPP

#include "typedefs.hpp"
#include "cl_HDF5.hpp"

namespace belfem
{
    class Database
    {
        uint mInterpolationOrder = 0 ;
        uint mNumberOfDimensions = 0 ;
        uint mNumberOfPointsPerCell = 0 ;

        uint mMemorySize = 0 ;

        uint * mNumPoints = nullptr ;

        // inverse of order
        real mOrder = 0 ;
        real mInvOrder = 0 ;

        // miniumum value in dataset
        real * mXmin = nullptr ;

        // maximum value in dataset
        real * mXmax = nullptr ;

        // stepwidth
        real * mStep = nullptr ;

        // inverse stepwidth
        real * mInvStep = nullptr ;

        // the actual data of the set
        real * mValues = nullptr ;

        // coefficients, only for 1D cubic
        real * mCoeffs = nullptr ;

        // work vector for the interpolation
        real * mWork = nullptr ;

        // work vector for the coordinates
        real  * mX = nullptr ;

        // pointer to the function that collects the data
        real
        ( Database::*mInterpolate )( const real * x ) ;

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        // the database constructor (standalone via filename)
        Database( const std::string & aFilename,
                  const std::string & aTablename );

//----------------------------------------------------------------------------

        // the database constructor (passing database)
        Database( HDF5 & aHdf5File ,
        const std::string & aTablename );

//----------------------------------------------------------------------------

        // the destructor
        ~Database() ;

//------------------------------------------------------------------------

        void
        load_from_file( const std::string & aFilename,
                        const std::string & aTablename );

//------------------------------------------------------------------------

        void
        load_from_file( HDF5          & aFile,
                        const std::string & aTablename );

//------------------------------------------------------------------------

        real
        compute( const real aX ) ;

//----------------------------------------------------------------------------

        real
        compute( const real aX, const real aY ) ;

//----------------------------------------------------------------------------

        real
        compute( const real aX, const real aY, const real aZ ) ;
        
//----------------------------------------------------------------------------
    private:
//----------------------------------------------------------------------------

        /**
         * clear the memory
         */
        void
        free();

//----------------------------------------------------------------------------

        void
        distribute_data();

//----------------------------------------------------------------------------

        void
        link_interpolation_function();

//----------------------------------------------------------------------------

        real
        interpolate_line2( const real * x );

//----------------------------------------------------------------------------

        real
        interpolate_line3( const real * x );

//----------------------------------------------------------------------------

        real
        interpolate_line4( const real * x );

//----------------------------------------------------------------------------

        real
        interpolate_quad4( const real * x );

//----------------------------------------------------------------------------

        real
        interpolate_quad9( const real * x );

//----------------------------------------------------------------------------

        real
        interpolate_quad16( const real * x );

//----------------------------------------------------------------------------

        real
        interpolate_hex8( const real * x );

//----------------------------------------------------------------------------

        real
        interpolate_hex27( const real * x );

//----------------------------------------------------------------------------

        real
        interpolate_hex64( const real * x );

//----------------------------------------------------------------------------

        /**
         * the index function for the 2d case
         */
        inline uint
        index( const uint i, const uint j ) const
        {
            return mNumPoints[ 0 ] * j + i ;
        }

//----------------------------------------------------------------------------

        /**
         * the index function for the 3d case
         */
        inline uint
        index( const uint i, const uint j, const uint k ) const
        {
            return mNumPoints[ 0 ] *  ( k * mNumPoints[ 1 ] +  j ) + i ;
        }

//----------------------------------------------------------------------------

        // only for 1d cubic spline
        inline uint
        find_cell_1dspline( const real x )
        {
            long int pivot = std::floor( ( x - mXmin[ 0 ] ) * mInvStep[ 0 ] ) ;

            if( pivot < 0 )
            {
                return 0 ;
            }
            else if ( pivot >= mNumPoints[ 0 ] )
            {
                return 4 * ( mNumPoints[ 0 ] - 1 );
            }
            else
            {
                return 4 * pivot ;
            }
        }

//------------------------------------------------------------------------

        inline uint
        find_cell_linear( const real * x, const uint dimension ) const
        {
            if( x[ dimension ] < mXmin[ dimension ] )
            {
                return 0 ;
            }
            else if ( mXmax[ dimension ] < x[ dimension ] )
            {
                return mNumPoints[ dimension ] - 2 ;
            }
            else
            {
                return uint( (x[ dimension ] - mXmin[ dimension ] ) * mInvStep[ dimension ] );
            }
        }

//------------------------------------------------------------------------

        inline uint
        find_cell_higher_order( const real * x, const unsigned int dimension ) const
        {
            if( x[ dimension ] < mXmin[ dimension ] )
            {
                return 0 ;
            }
            else if ( mXmax[ dimension ] < x[ dimension ] )
            {
                return mNumPoints[ dimension ] - mInterpolationOrder - 1 ;
            }
            else
            {
                return mInterpolationOrder * uint( ( x[ dimension ] - mXmin[ dimension ] ) * mInvStep[ dimension ] );
            }
        }

//------------------------------------------------------------------------

        inline real
        xtoxi_linear( const real * x, const uint index, const unsigned int dimension ) const
        {
            return ( x[ dimension ] + x[ dimension ]
                     - mXmin[ dimension ] - mXmin[ dimension ] ) * mInvStep[ dimension ]
                   - index - index - 1 ;
        }

//------------------------------------------------------------------------

        inline real
        xtoxi_quadratic( const real * x, const uint index, const unsigned int dimension ) const
        {
            return ( x[ dimension ] + x[ dimension ]
                     - mXmin[ dimension ] - mXmin[ dimension ] ) * mInvStep[ dimension ]
                   - index - 1 ;
        }

//------------------------------------------------------------------------

        inline real
        xtoxi_higher_order( const real * x, const uint index, const unsigned int dimension ) const
        {
            return 2 * ( x[ dimension ] - mXmin[ dimension ]
                         - mInvOrder * index * mStep[ dimension ] )
                   * mInvStep[ dimension ] - 1 ;
        }


//------------------------------------------------------------------------


        inline real
        interpolate_1dspline( const real * x )
        {
            real X = *x ;

            // get index
            std::size_t k = this->find_cell_1dspline( X );

            return ( (   mCoeffs[ k ]     * X
                         + mCoeffs[ k+1 ] ) * X
                     + mCoeffs[ k+2 ] ) * X
                   + mCoeffs[ k+3 ]  ;
        }
        
//----------------------------------------------------------------------------
    };
}
#endif //BELFEM_CL_DATABASE_HPP
