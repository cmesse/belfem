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

#ifndef BELFEM_CL_AR_VECTOR_HPP
#define BELFEM_CL_AR_VECTOR_HPP

#include <memory>
#include "armadillo.hpp"

#include "typedefs.hpp"
#include "assert.hpp"
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
    template < typename T >
    class Vector
    {
//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        typedef arma::Mat <T> VectorType;

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
                mVector( aNumRows, 1 ) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * Constructor with fill value
         */
        Vector( const size_t aNumRows, const T & aValue ) :
                mVector( aNumRows, 1 )
        {
            this->fill( aValue );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * Constructor with initializer list
         */
        Vector( std::initializer_list<T> aInitList )
        {
            if( aInitList.size() == 1 )
            {
                // Single element: interpret as Vector(1, value)
                mVector.set_size( 1, 1 );
                mVector( 0 ) = *aInitList.begin();
            }
            else
            {
                // Multiple elements: interpret as list of values
                mVector.set_size( aInitList.size(), 1 );
                size_t i = 0;
                for( const T & val : aInitList )
                {
                    mVector( i++ ) = val;
                }
            }
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector( const Cell< T > & aCell )
        {
            mVector.set_size( aCell.size(), 1 );
            std::copy( aCell.begin(), aCell.end(), mVector.begin() );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector( const std::vector< T > & aVector )
        {
            mVector.set_size( aVector.size(), 1 );
            std::copy( aVector.begin(), aVector.end(), mVector.begin() );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
        * Constructor from expression
        */
        Vector( const VectorType & aExpression ) :
        mVector( aExpression )
        {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
        * Constructor from expression
        */
        template < typename ET, typename OP>
        Vector( const arma::Op<ET,OP> & aExpression )
                : mVector( aExpression ) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
        * Constructor from column
        */
        Vector( const arma::subview_col< T > & aCol ) :
                mVector( aCol )
        {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
        * Constructor from row
        */
        Vector( const arma::subview_row< T > & aRow ) :
                mVector( arma::trans( aRow ) )
        {}

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
            return mVector.memptr();
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * expose the underlying raw pointer ( const version )
         */
        const T *
        data() const
        {
            return mVector.memptr();
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
            mVector.fill( aValue );
        }

//------------------------------------------------------------------------------

        /**
         * change the size of the vector
         */
        void
        set_size( const size_t aNumRows )
        {
            mVector.set_size( aNumRows, 1 );
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
            return mVector.n_rows;
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

            return mVector( aIndex );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const T &
        operator()( const size_t aIndex ) const
        {
            BELFEM_ASSERT( aIndex < this->length(),
                          "Index %lu out of bounds, which must be smaller than %lu.",
                          ( long unsigned int ) aIndex,
                          ( long unsigned int ) this->length());

            return mVector( aIndex );
        }

//------------------------------------------------------------------------------
// ITERATORS
//------------------------------------------------------------------------------

        auto
        begin() -> decltype( mVector.begin() )
        {
            return mVector.begin();
        }

        auto
        end() -> decltype( mVector.end() )
        {
            return mVector.end();
        }

        auto
        begin() const -> decltype( mVector.begin() )
        {
            return mVector.begin();
        }

        auto
        end() const -> decltype( mVector.end() )
        {
            return mVector.end();
        }

//------------------------------------------------------------------------------
// EQUAL OPERATORS
//------------------------------------------------------------------------------

        Vector< T > &
        operator=( const T & aValue )
        {
            this->fill( aValue );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
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
                    mVector( i++ ) = val;
                }
            }
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET >
        Vector< T > &
        operator=( const ET & aExpression )
        {
            mVector = aExpression;
            return *this;
        }

//------------------------------------------------------------------------------
// ADD OPERATORS
//------------------------------------------------------------------------------

        Vector< T > &
        operator+=( const T & aValue )
        {
            mVector += aValue;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator+=( const Vector< T > & aVector )
        {
            mVector += aVector.vector_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator+=( const VectorType & aExpression )
        {
            mVector += aExpression;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET >
        Vector<T> &
        operator+=( const ET & aExpression )
        {
            mVector += aExpression;
            return *this;
        }

//------------------------------------------------------------------------------
// SUBTRACT OPERATORS
//------------------------------------------------------------------------------

        Vector< T > &
        operator-=( const T & aValue )
        {
            mVector -= aValue;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator-=( const Vector< T > & aVector )
        {
            mVector -= aVector.vector_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator-=( const VectorType & aExpression )
        {
            mVector -= aExpression;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET >
        Vector<T> &
        operator-=( const ET & aExpression )
        {
            mVector -= aExpression;
            return *this;
        }

//------------------------------------------------------------------------------
// MULTIPLY OPERATORS
//------------------------------------------------------------------------------

        Vector< T > &
        operator*=( const T & aValue )
        {
            mVector *= aValue;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator*=( const Vector< T > & aVector )
        {
            mVector *= aVector.vector_data();
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator*=( const VectorType & aExpression )
        {
            mVector *= aExpression;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator%=( const Vector< T > & aVector )
        {
            BELFEM_ASSERT( aVector.length() == this->length(),
                          "Vector sizes do not match ( %lu and  %lu ).",
                          ( long unsigned int ) this->length(),
                          ( long unsigned int ) aVector.length() );

            mVector %= aVector.vector_data();
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

        Vector< T > &
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

            fprintf( tOutFile, "%d; ",
                ( int ) tThis( i ) );


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
#endif //BELFEM_CL_AR_VECTOR_HPP
