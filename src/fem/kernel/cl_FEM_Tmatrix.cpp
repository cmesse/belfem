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
#include "cl_FEM_Tmatrix.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Tmatrix::Tmatrix( const Matrix< belfem::real > & aMatrix ) :
            mNumRows( aMatrix.n_rows() ),
            mNumCols( aMatrix.n_cols() ),
            mNumNonZeros( this->count_nnz( aMatrix ) )
        {

            mValues = ( real * ) malloc( mNumNonZeros * sizeof( real ));

            // allocate pointer array
            mPointers = ( uint * ) malloc(( mNumRows + 1 ) * sizeof( uint ));

            // first entry
            mPointers[ 0 ] = 0;

            // allocate column array
            mIndices = ( uint * ) malloc(( mNumNonZeros ) * sizeof( uint ));

            // position in columns
            uint tStep = 0;

            uint tCount ;

            // populate pointers array
            for ( uint i = 0; i < mNumRows; ++i )
            {
                // reset counter
                tCount = 0;

                // loop over all columns
                for ( uint j = 0; j < mNumCols; ++j )
                {
                    if ( abs( aMatrix( i, j )) > BELFEM_EPSILON )
                    {
                        // write column
                        mIndices[ tStep ] = j;

                        // write value
                        mValues[ tStep++ ] = aMatrix( i, j );

                        // increment counter
                        ++tCount;
                    }
                }

                // count entries
                mPointers[ i + 1 ] = mPointers[ i ] + tCount;
            }
        }

//------------------------------------------------------------------------------

        Tmatrix::~Tmatrix()
        {
            free( mPointers );
            free( mValues );
            free( mIndices );
        }

//------------------------------------------------------------------------------

        uint
        Tmatrix::count_nnz( const Matrix< real > & aMatrix )
        {
            uint aNNZ = 0 ;

            for( uint i=0; i<aMatrix.n_rows(); ++i )
            {
                for( uint j=0; j<aMatrix.n_cols(); ++j )
                {
                    if( abs( aMatrix( i, j )) > BELFEM_EPSILON )
                    {
                        ++aNNZ ;
                    }
                }
            }

            return aNNZ ;
        }

//------------------------------------------------------------------------------

        void
        Tmatrix::project( const Matrix< real > & aA, Matrix< real > & aB ) const
        {
            BELFEM_ASSERT( aA.n_rows() == mNumRows && aA.n_cols() == mNumRows, "Invalid dimesion of input matrix" );

            // Set size of aB and initialize all elements to 0.0
            aB.set_size( mNumCols, mNumCols, 0.0 );

            // Compute B = T' * A * T

            // loop over all rows
            for ( uint i = 0; i < mNumRows; ++i)
            {
                // Iterate over non-zero elements in row i of T
                for ( uint p = mPointers[ i ]; p < mPointers[ i + 1 ]; ++p )
                {
                    uint j = mIndices[ p ];
                    real T_ji = mValues[ p ];

                    for( uint k=0; k<mNumRows; ++k )
                    {
                        for( uint q=mPointers[ k ]; q<mPointers[ k+1 ]; ++q )
                        {
                            uint l = mIndices[ q ];
                            real T_kl = mValues[ q ];

                            aB( j, l ) += T_ji * aA( i, k ) * T_kl ;
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Tmatrix::project( const Vector< real > & aA, Vector< real > & aB ) const
        {
            BELFEM_ASSERT( aA.length() == mNumRows, "Invalid dimension of input matrix" );

            // performs the operation B = T'*A
            // Ensure the output vector is of the correct size and initialize it to zero
            aB.set_size( mNumCols, 0.0 );

            // Perform the matrix-vector multiplication
            for ( uint i = 0; i < mNumRows; ++i)
            {
                // Iterate over non-zero elements in row i of T
                for ( uint j = mPointers[ i ]; j < mPointers[ i + 1 ]; ++j)
                {
                    aB( mIndices[ j ] ) += mValues[ j ] * aA( i );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Tmatrix::get_matrix( Matrix< real > & aMatrix ) const
        {
            aMatrix.set_size( mNumRows, mNumCols, 0.0 );

            // populate data
            for ( uint i = 0; i < mNumRows; ++i)
            {
                // Iterate over non-zero elements in row i of T
                for ( uint j = mPointers[ i ]; j < mPointers[ i + 1 ]; ++j)
                {
                    aMatrix( i, mIndices[ j ] ) = mValues[ j ] ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Tmatrix::print() const
        {
            /*std::cout << "n   : " << mNumRows << std::endl ;
            std::cout << "m   : " << mNumCols << std::endl ;
            std::cout << "nnz : " << mNumNonZeros << std::endl ;

            std::cout << "pointers :" << std::endl ;
            for( uint k=0; k<=mNumRows; ++k )
            {
                std::cout << "  " << k << ", " << mPointers[ k ] << std::endl ;
            }
            std::cout <<  std::endl << "indices & values :" << std::endl ;
            for( uint k=0; k<mNumNonZeros; ++k )
            {
                std::cout << "  " << k << ", " << mIndices[ k ] << ", " << mValues[ k ] << std::endl ;
            }*/

            Matrix< real > T( mNumRows, mNumCols, 0.0 );

            this->get_matrix( T );

            T.print("T");
        }

//------------------------------------------------------------------------------

    }
}