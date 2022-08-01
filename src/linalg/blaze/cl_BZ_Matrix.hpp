//
// Created by Christian Messe on 03.09.19.
//

#ifndef BELFEM_CL_BZ_MATRIX_HPP
#define BELFEM_CL_BZ_MATRIX_HPP

#include "typedefs.hpp"
#include "assert.hpp"
namespace belfem
{
    template< typename T >
    class Matrix
    {
//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        typedef blaze::DynamicMatrix< T, BLAZE_DEFAULT_STORAGE_ORDER > MatrixType;

//------------------------------------------------------------------------------
    private :
//------------------------------------------------------------------------------

        // member class of underlying vector implementation
        MatrixType mMatrix;

//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        /**
         * empty constructor
         */
        Matrix() = default;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
          * Constructor without fill value
          */
        Matrix( const size_t aNumRows,
                   const size_t aNumCols ) :
                mMatrix( aNumRows, aNumCols ) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
          * Constructor with fill value
          */
        Matrix( const size_t aNumRows,
                   const size_t aNumCols,
                   const real aValue ) :
            mMatrix( aNumRows, aNumCols, aValue )
        {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * Constructor with initializer list
         */
        Matrix( const std::initializer_list<std::initializer_list<T> > & aInitList )
                : mMatrix( aInitList )
        {

        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
        * Constructor from expression
        */
        Matrix( const MatrixType & aExpression )
                : mMatrix( aExpression ) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
          * Constructor from expression
          */
        template< typename ET >
        Matrix( const blaze::Expression<ET> & aExpression )
                : mMatrix( aExpression ) {}

//------------------------------------------------------------------------------

        /**
         * empty destructor
         */
        virtual ~Matrix() = default;

//------------------------------------------------------------------------------
// MEMORY
//------------------------------------------------------------------------------

        /*
         * expose the underlying raw pointer
         */

        inline T *
        data()
        {
            return mMatrix.data();
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /*
         * expose the underlying raw pointer ( const version )
         */
        inline const T *
        data() const
        {
            return mMatrix.data();
        }

//------------------------------------------------------------------------------

        /*
         * expose the underlying matrix implementation ( writable version )
         */
        inline MatrixType &
        matrix_data()
        {
            return mMatrix;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /*
         * expose the underlying matrix implementation ( const version )
         */
        inline const MatrixType &
        matrix_data() const
        {
            return mMatrix;
        }

//------------------------------------------------------------------------------
// UTILITIES
//------------------------------------------------------------------------------

        inline void
        fill( const T & aValue )
        {
            mMatrix = aValue;
        }

//------------------------------------------------------------------------------

        inline void
        set_size( const size_t aNumRows, const size_t aNumCols )
        {
            mMatrix.resize( aNumRows, aNumCols, false );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline void
        set_size( const size_t aNumRows, const size_t aNumCols, const T & aValue )
        {
            this->set_size( aNumRows, aNumCols );
            this->fill( aValue );
        }

//------------------------------------------------------------------------------
// SIZE CHECKS
//------------------------------------------------------------------------------

        inline size_t
        n_rows() const
        {
            return mMatrix.rows();
        }

//------------------------------------------------------------------------------

        inline size_t
        n_cols() const
        {
            return mMatrix.columns();
        }

//------------------------------------------------------------------------------

        /**
         * length of data container
         */
        inline size_t
        capacity() const
        {
            return mMatrix.capacity();
        }

//------------------------------------------------------------------------------
// OPERATORS
//------------------------------------------------------------------------------

        /**
         * access operator ( writable version )
         */
        inline T &
        operator()( const size_t aRowIndex, const size_t aColIndex )
        {
            BELFEM_ASSERT( aRowIndex< this->n_rows(),
                          "Row index %lu out of bounds, which must be smaller than %lu.",
                          ( long unsigned int ) aRowIndex,
                          ( long unsigned int ) this->n_rows() );

            BELFEM_ASSERT( aColIndex< this->n_cols(),
                          "Col index %lu out of bounds, which must be smaller than %lu.",
                          ( long unsigned int ) aColIndex,
                          ( long unsigned int ) this->n_cols() );

            return mMatrix( aRowIndex, aColIndex );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * access operator ( const version )
         */
        virtual inline const T &
        operator()( const size_t aRowIndex, const size_t aColIndex ) const
        {
            BELFEM_ASSERT( aRowIndex< this->n_rows(),
                          "Row index %lu out of bounds, which must be smaller than %lu.",
                          ( long unsigned int ) aRowIndex,
                          ( long unsigned int ) this->n_rows() );

            BELFEM_ASSERT( aColIndex< this->n_cols(),
                          "Col index %lu out of bounds, which must be smaller than %lu.",
                          ( long unsigned int ) aColIndex,
                          ( long unsigned int ) this->n_cols() );

            return mMatrix( aRowIndex, aColIndex );
        }

//------------------------------------------------------------------------------
// Cols and Rows
//------------------------------------------------------------------------------

        inline auto
        row( const size_t aRowIndex )
            ->decltype( blaze::rows( mMatrix, { 0 } ) )
        {
            BELFEM_ASSERT( aRowIndex < this->n_rows(),
                          "Row index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aRowIndex,
                          ( long unsigned int ) this->n_rows() );

            std::array < size_t, 1 > tArray;
            tArray[ 0 ] = aRowIndex;
            return blaze::rows( mMatrix, tArray );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline auto
        row( const size_t aRowIndex ) const
            ->decltype( blaze::row( mMatrix, aRowIndex ) )
        {
            BELFEM_ASSERT( aRowIndex < this->n_rows(),
                          "Row index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aRowIndex,
                          ( long unsigned int ) this->n_rows() );

            return blaze::row( mMatrix, aRowIndex );
        }

//------------------------------------------------------------------------------

        inline auto
        col( const size_t aColIndex )
            ->decltype( blaze::column( mMatrix, aColIndex )  )
        {
            BELFEM_ASSERT( aColIndex < this->n_cols(),
                          "Col index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aColIndex,
                          ( long unsigned int ) this->n_cols() );

            return blaze::column( mMatrix, aColIndex );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline auto
        col( const size_t aColIndex ) const
            ->decltype( blaze::column( mMatrix, aColIndex )  )
        {
            BELFEM_ASSERT( aColIndex < this->n_cols(),
                          "Col index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aColIndex,
                          ( long unsigned int ) this->n_cols() );

            return blaze::column( mMatrix, aColIndex );
        }

//------------------------------------------------------------------------------

        inline auto
        submat( const size_t aFirstRow,
                const size_t aFirstCol,
                const size_t aLastRow,
                const size_t aLastCol )
            -> decltype( blaze::submatrix( mMatrix,  aFirstRow, aFirstCol, aLastRow, aLastCol ) )
        {
            BELFEM_ASSERT( aFirstRow < this->n_cols(),
                          "First row index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aFirstRow,
                          ( long unsigned int ) this->n_rows() );

            BELFEM_ASSERT( aLastRow < this->n_cols(),
                          "Last row index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aLastRow,
                          ( long unsigned int ) this->n_rows() );

            BELFEM_ASSERT( aFirstRow < this->n_cols(),
                          "First col index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aFirstCol,
                          ( long unsigned int ) this->n_cols() );

            BELFEM_ASSERT( aLastRow < this->n_cols(),
                          "Last col index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aLastCol,
                          ( long unsigned int ) this->n_cols() );

            return blaze::submatrix( mMatrix,  aFirstRow, aFirstCol, aLastRow, aLastCol );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline auto
        submat( const size_t aFirstRow,
                const size_t aFirstCol,
                const size_t aLastRow,
                const size_t aLastCol ) const
            -> decltype( blaze::submatrix( mMatrix,  aFirstRow, aFirstCol, aLastRow, aLastCol ) )
        {
            BELFEM_ASSERT( aFirstRow < this->n_cols(),
                          "First row index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aFirstRow,
                          ( long unsigned int ) this->n_rows() );

            BELFEM_ASSERT( aLastRow < this->n_cols(),
                          "Last row index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aLastRow,
                          ( long unsigned int ) this->n_rows() );

            BELFEM_ASSERT( aFirstRow < this->n_cols(),
                          "First col index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aFirstCol,
                          ( long unsigned int ) this->n_cols() );

            BELFEM_ASSERT( aLastRow < this->n_cols(),
                          "Last col index %lu out of bounds, which must be less than %lu.",
                          ( long unsigned int ) aLastCol,
                          ( long unsigned int ) this->n_cols() );

            return blaze::submatrix( mMatrix,  aFirstRow, aFirstCol, aLastRow, aLastCol );
        }

//------------------------------------------------------------------------------

        inline void
        set_row( const size_t aRowIndex, const Vector <T> & aVector )
        {

            // make sure that length is correct
            BELFEM_ASSERT( aVector.length() == this->n_cols(),
                          "Wrong length of vector : %lu ( expect %lu )",
                          ( long unsigned int ) aVector.length(),
                          ( long unsigned int ) this->n_cols() );

            // get pointer to data
            const T * tSource = aVector.data();

            // manual copy of data
            for( std::size_t k=0; k<aVector.length(); ++k )
            {
                mMatrix( aRowIndex, k ) = tSource[ k ];
            }
        }

//------------------------------------------------------------------------------

        inline void
        set_col( const size_t aRowIndex, const Vector <T> & aVector )
        {
            // make sure that length is correct
            BELFEM_ASSERT( aVector.length() == this->n_rows(),
                          "Wrong length of vector : %lu ( expect %lu )",
                          ( long unsigned int ) aVector.length(),
                          ( long unsigned int ) this->n_rows() );

            std::copy(  aVector.data(), aVector.data() + aVector.length(),
                    this->col( aRowIndex ).data() );
        }

//------------------------------------------------------------------------------
// internal operators
//------------------------------------------------------------------------------

        /**
         * Fill with value
         */
        Matrix< T > &
        operator=( const T & aValue )
        {
            this->fill( aValue );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * Fill with expression
         */
        template < typename ET >
        Matrix< T > &
        operator=( const ET & aExpression )
        {
            mMatrix = aExpression;
            return *this;
        }

//------------------------------------------------------------------------------
// ADD OPERATORS
//------------------------------------------------------------------------------

        Matrix< T > &
        operator+=( const T & aValue )
        {
            mMatrix += aValue;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Matrix< T > &
        operator+=( const Matrix< T > & aMatrix )
        {
            mMatrix += aMatrix.matrix_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef BELFEM_GCC
        template < typename ET >
        Matrix< T > &
        operator+=( const ET & aExpression )
        {
            mMatrix += aExpression;
            return *this;
        }
#endif

//------------------------------------------------------------------------------
// SUBTRACT OPERATORS
//------------------------------------------------------------------------------

        Matrix< T > &
        operator-=( const T & aValue )
        {
            mMatrix -= aValue;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Matrix< T > &
        operator-=( const Matrix< T > & aMatrix )
        {
            mMatrix -= aMatrix.matrix_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef BELFEM_GCC
       Matrix< T > &
       operator-=( const MatrixType & aExpression )
       {
           mMatrix -= aExpression;
           return *this;
       }
#endif

//------------------------------------------------------------------------------
// MULTIPLY OPERATORS
//------------------------------------------------------------------------------

        Matrix< T > &
        operator*=( const T & aValue )
        {
            mMatrix *= aValue;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Matrix< T > &
        operator*=( const Matrix< T > & aMatrix )
        {
            mMatrix *= aMatrix.matrix_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef BELFEM_GCC
        Matrix< T > &
        operator*=( const MatrixType & aExpression )
        {
            mMatrix *= aExpression;
            return *this;
        }
#endif

//------------------------------------------------------------------------------
// DIVIDE OPERATORS
//------------------------------------------------------------------------------

        Matrix< T > &
        operator/=( const T & aValue )
        {
            mMatrix /= aValue;
            return *this;
        }

//------------------------------------------------------------------------------
// PRINT OPERATION
//------------------------------------------------------------------------------

        inline void
        print( const std::string aLabel="Matrix" ) const;

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------
// turn off annoying waning
#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#endif
    template < typename T > void
    Matrix< T >::print( const std::string aLabel ) const
    {
        FILE * tOutFile = stdout;

        fprintf( tOutFile, "%s = [ ... \n", aLabel.c_str() );

        // get a ref of this
        const Matrix< T > & tThis = *this;

        uint tRows = tThis.n_rows();
        uint tCols = tThis.n_cols();

        for( uint i=0; i< tRows; ++i )
        {
            for( uint j=0; j< tCols; ++j )
            {
                if( j < tCols-1 )
                {
                    fprintf( tOutFile, "%d, ", ( int ) tThis( i, j ) );
                }
                else
                {
                    fprintf( tOutFile, "%d; ", ( int ) tThis( i, j ) );
                }
            }

            if( i < tRows-1 )
            {
                fprintf( tOutFile, "...\n" );
            }
            else
            {
                fprintf( tOutFile, "];\n" );
            }
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline void
    Matrix<real>::print( const std::string aLabel ) const
    {
        FILE * tOutFile = stdout;

        fprintf( tOutFile, "%s = [ ... \n", aLabel.c_str() );

        const Matrix< real > & tThis = *this;

        uint tRows = tThis.n_rows();
        uint tCols = tThis.n_cols();

        for( uint i=0; i< tRows; ++i )
        {
            for( uint j=0; j< tCols; ++j )
            {
                if( j < tCols-1 )
                {
                    fprintf( tOutFile, "%+.15e, ", ( double ) tThis( i, j ) );
                }
                else
                {
                    fprintf( tOutFile, "%+.15e; ", ( double ) tThis( i, j ) );
                }
            }

            if( i < tRows-1 )
            {
                fprintf( tOutFile, "...\n" );
            }
            else
            {
                fprintf( tOutFile, "];\n" );
            }
        }
    }

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
}
#endif //BELFEM_CL_BZ_MATRIX_HPP
