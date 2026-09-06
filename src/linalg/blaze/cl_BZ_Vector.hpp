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

#ifndef BELFEM_CL_BZ_VECTOR_HPP
#define BELFEM_CL_BZ_VECTOR_HPP

#ifdef BELFEM_INTEL
#pragma warning push
#pragma warning disable 3058
#endif

#include <memory>
#include <algorithm>
#include "blaze_config.hpp"
#include <blaze/Blaze.h>

#include "typedefs.hpp"
#include "assert.hpp"
#include <blaze/math/Column.h>
#include "cl_Cell.hpp"

namespace belfem
{
    /**
     * @brief Column vector. The backend (Armadillo or Blaze) is chosen at build time;
     * both provide the same belfem::Vector interface.
     *
     * @ingroup grp_linalg
     * @see @ref linalg_lapack_usage_guide
     */
    template< typename T >
    class Vector
    {
//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        typedef blaze::DynamicVector<T, blaze::columnVector> VectorType;

//------------------------------------------------------------------------------
    private :
//------------------------------------------------------------------------------

        // member class of underlying vector implementation
        VectorType mVector;

//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        /**
         * empty constructor
         */
        Vector() = default;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * Constructor without fill value
         */
        Vector( const size_t aNumRows ) :
                mVector( aNumRows ) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * Constructor with fill value
         */
        Vector( const size_t aNumRows, const T & aValue ) :
                mVector( aNumRows, aValue ) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * Constructor with initializer list
         */
        Vector( std::initializer_list<T> aInitList )
        {
            if( aInitList.size() == 1 )
            {
                // Single element: interpret as Vector(1, value)
                mVector.resize( 1, false );
                mVector[ 0 ] = *aInitList.begin();
            }
            else
            {
                // Multiple elements: interpret as list of values
                mVector.resize( aInitList.size(), false );
                size_t i = 0;
                for( const T & val : aInitList )
                {
                    mVector[ i++ ] = val;
                }
            }
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
        * Constructor from expression
        */
        Vector( const VectorType & aExpression ) :
                mVector( aExpression ) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector( const Cell< T > & aCell )
        {
            mVector.resize( aCell.size(), false );
            std::copy( aCell.begin(), aCell.end(), mVector.begin() );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector( const std::vector< T > & aVector )
        {
            mVector.resize( aVector.size(), false );
            std::copy( aVector.begin(), aVector.end(), mVector.begin() );
        }


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * Constructor from column vector expression (TF=false)
         */
        template< typename VT >
        Vector( const blaze::DenseVector<VT, false> & aExpression )
                : mVector( static_cast< const VT& >( aExpression ) ) {}

        /**
         * Constructor from row vector expression (TF=true), transposes to column
         */
        template< typename VT >
        Vector( const blaze::DenseVector<VT, true> & aExpression )
        {
            const VT & tRow = static_cast< const VT& >( aExpression );
            mVector.resize( tRow.size(), false );
            for( size_t k = 0; k < tRow.size(); ++k )
            {
                mVector[ k ] = tRow[ k ];
            }
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
        * Constructor from column
        */
        Vector( const blaze::Columns<blaze::DynamicMatrix<T, true>, true, true, false> & aColumn )
        {
            // get length of vector
            std::size_t tLength = aColumn.rows();

            // make sure that this is a column matrix
            BELFEM_ASSERT( aColumn.columns() == 1,
                          "aColumn must be a single column, but has %lu columns.",
                          ( long unsigned int ) aColumn.columns() );

            // create data
            mVector.resize( tLength , false );

            std::copy(  aColumn.data(), aColumn.data() + tLength, mVector.data() );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
        * Constructor from row
        */
        Vector( const blaze::Rows<blaze::DynamicMatrix<T, true>, false, true, false> & aRow )
        {
            // get length of vector
            std::size_t tLength = aRow.columns();

            // make sure that this is a row matrix
            BELFEM_ASSERT( aRow.rows() == 1,
                          "aColumn must be a single row, but has %lu rows.",
                          ( long unsigned int ) aRow.rows());

            // create data
            mVector.resize( tLength, false );

            // get pointer to data
            T * tTarget = mVector.data();

            // manual copy of data
            for ( std::size_t k = 0; k < tLength; ++k )
            {
                tTarget[ k ] = aRow.at( 0, k );
            }
        }

//------------------------------------------------------------------------------

        /**
        *copy constructor
        */
        Vector( const Vector< T > & aVector ) :
            mVector( aVector.mVector )
        {}

        /**
         *move constructor
         */
        Vector( Vector< T > && aVector ) noexcept :
                mVector( std::move( aVector.mVector ) )
        {}
        
//------------------------------------------------------------------------------

        /**
         * empty destructor
         */
        virtual ~Vector() = default;

//------------------------------------------------------------------------------
// MEMORY
//------------------------------------------------------------------------------

        /**
         * expose the underlying raw pointer ( writable version )
         */
        T *
        data()
        {
            return mVector.data();
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * expose the underlying raw pointer ( const version )
         */
        const T *
        data() const
        {
            return mVector.data();
        }

//------------------------------------------------------------------------------

        /**
        * expose the underlying matrix implementation ( writable version )
        */
        VectorType &
        vector_data()
        {
            return mVector;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * expose the underlying matrix implementation ( const version )
         */
        const VectorType &
        vector_data() const
        {
            return mVector;
        }

//------------------------------------------------------------------------------
// UTILITIES
//------------------------------------------------------------------------------

        /**
         * write value into all entries of the vector
         */
        void
        fill( const T & aValue )
        {
            mVector = aValue;
        }

//------------------------------------------------------------------------------

        /**
         * change the size of the vector
         */
        void
        set_size( const size_t aNumRows )
        {
            mVector.resize( aNumRows, false );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        set_size( const size_t aNumRows, const T & aValue )
        {
            this->set_size( aNumRows );
            this->fill( aValue );
        }

//------------------------------------------------------------------------------

        /**
         * get the length of the vector
         */
        size_t
        length() const
        {
            return mVector.size();
        }

//------------------------------------------------------------------------------
// ACCESS OPERATORS
//------------------------------------------------------------------------------

        /**
        * move assignment operator
        */
        Vector< T > &
        operator=( Vector< T > && aVector ) noexcept
        {
            if ( this != & aVector )
            {
                mVector = std::move( aVector.mVector );
            }
            return *this;
        }

        /**
         * copy assignment operator
         */
        Vector< T > &
        operator=( const Vector< T > & aVector )
        {
            if ( this != &aVector )
            {
                mVector = aVector.mVector;
            }
            return *this;
        }

//------------------------------------------------------------------------------
        T &
        operator()( const size_t aIndex )
        {
            BELFEM_ASSERT( aIndex < this->length(),
                          "Index %lu out of bounds, which must be smaller than %lu.",
                          ( long unsigned int ) aIndex,
                          ( long unsigned int ) this->length());

            return mVector[ aIndex ];
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const T &
        operator()( const size_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < this->length(),
                          "Index %lu out of bounds, which must be smaller than %lu.",
                          ( long unsigned int ) aIndex,
                          ( long unsigned int ) this->length());

            return mVector[ aIndex ];
        }

//------------------------------------------------------------------------------
// ITERATORS
//------------------------------------------------------------------------------

        T *
        begin()
        {
            return mVector.data();
        }

        T *
        end()
        {
            return mVector.data() + this->length();
        }

        const T *
        begin() const
        {
            return mVector.data();
        }

        auto
        end() const
        {
            return mVector.data() + this->length();
        }


//------------------------------------------------------------------------------
// EQUAL OPERATORS
//------------------------------------------------------------------------------

        Vector<T> &
        operator=( const T & aValue )
        {
            this->fill( aValue );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector<T> &
        operator=( std::initializer_list<T> aInitList )
        {
            if( aInitList.size() == 1 )
            {
                // Single element: interpret as set_size(1, value)
                this->set_size( 1, *aInitList.begin() );
            }
            else
            {
                // Multiple elements: interpret as list of values
                this->set_size( aInitList.size() );
                size_t i = 0;
                for( const T & val : aInitList )
                {
                    mVector[ i++ ] = val;
                }
            }
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET >
        Vector<T> &
        operator=( const ET & aExpression )
        {
            mVector = aExpression;
            return *this;
        }

//------------------------------------------------------------------------------
// ADD OPERATORS
//------------------------------------------------------------------------------

        Vector<T> &
        operator+=( const T & aValue )
        {
            std::for_each( mVector.data(), mVector.data() + this->length(),
                           [ aValue ]( T & tVal )
                           { tVal += aValue; } );
            return *this;
        }


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        //Vector<T> &
        //operator+=( const VectorType & aExpression )
       // {
       //     mVector += aExpression;
       //     return *this;
       // }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector<T> &
        operator+=( const Vector<T> & aVector )
        {
            mVector += aVector.vector_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        //template < typename ET >
        //Vector<T> &
        //operator+=( const ET & aExpression )
        //{
        //    mVector += aExpression;
        //    return *this;
        //}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET >
        Vector<T> &
        operator+=( const blaze::DMatScalarMultExpr< ET, T, false> & aExpression )
        {
            mVector += blaze::column( aExpression, 0UL );
            return *this;
        }

//------------------------------------------------------------------------------
// SUBTRACT OPERATORS
//------------------------------------------------------------------------------

        Vector<T> &
        operator-=( const T & aValue )
        {
            std::for_each( mVector.data(), mVector.data() + this->length(),
                           [ aValue ]( T & tVal )
                           { tVal -= aValue; } );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector<T> &
        operator-=( const Vector<T> & aVector )
        {
            mVector -= aVector.vector_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        //Vector<T> &
        //operator-=( const VectorType & aExpression )
       // {
       //     mVector -= aExpression;
       //     return *this;
       // }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        //template < typename ET >
        //Vector<T> &
        //operator-=( const ET & aExpression )
        //{
        //    mVector -= aExpression;
        //    return *this;
        //}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET >
        Vector<T> &
        operator-=( const blaze::DMatScalarMultExpr< ET, T, false> & aExpression )
        {
            mVector -= blaze::column( aExpression, 0UL );
            return *this;
        }

//------------------------------------------------------------------------------
// MULTIPLY OPERATORS
//------------------------------------------------------------------------------

        Vector<T> &
        operator*=( const T & aValue )
        {
            mVector *= aValue;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector<T> &
        operator*=( const Vector<T> & aVector )
        {
            mVector *= aVector.vector_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector<T> &
        operator*=( const VectorType & aExpression )
        {
            mVector *= aExpression;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator%=( const Vector< T > & aVector )
        {
            size_t tN = this->length();

            BELFEM_ASSERT( aVector.length() == tN,
                          "Vector sizes do not match ( %lu and  %lu ).",
                          ( long unsigned int ) tN,
                          ( long unsigned int ) aVector.length() );

            for( size_t k=0; k<tN; ++k )
            {
                mVector[ k ] *= aVector( k );
            }
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET >
        Vector<T> &
        operator*=( const ET & aExpression )
        {
            mVector *= aExpression;
            return *this;
        }

//------------------------------------------------------------------------------
// DEVIDE OPERATORS
//------------------------------------------------------------------------------

        Vector<T> &
        operator/=( const T & aValue )
        {
            mVector /= aValue;
            return *this;
        }

//------------------------------------------------------------------------------
// PRINT OPERATION
//------------------------------------------------------------------------------

        void
        print( const std::string aLabel="Vector" ) const;

//------------------------------------------------------------------------------
    };

// turn off annoying waning
#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#endif

    template < typename T > void
    Vector< T >::print( const std::string aLabel ) const
    {
        FILE * tOutFile = stdout;

        fprintf( tOutFile, "\n%s = [ ... \n", aLabel.c_str() );

        // get a ref of this
        const Vector< T > & tThis = *this;

        uint tLength = tThis.length();

        for( uint i=0; i< tLength; ++i )
        {

            fprintf( tOutFile, "%d; ", ( int ) tThis( i ) );


            if( i < tLength-1 )
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
    Vector< real >::print( const std::string aLabel ) const
    {
        FILE * tOutFile = stdout;

        fprintf( tOutFile, "%s = [ ... \n", aLabel.c_str() );

        const Vector< real > & tThis = *this;

        uint tLength = tThis.length();

        for( uint i=0; i< tLength; ++i )
        {

            fprintf( tOutFile, "%+.15e; ", ( double ) tThis( i ) );


            if( i < tLength-1 )
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

#ifdef BELFEM_INTEL
#pragma warning pop
#endif

#endif //BELFEM_CL_BZ_VECTOR_HPP
