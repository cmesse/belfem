//
// Created by Christian Messe on 03.09.19.
//

#ifndef BELFEM_CL_AR_VECTOR_HPP
#define BELFEM_CL_AR_VECTOR_HPP

#include <memory>
#include "armadillo.hpp"

#include "typedefs.hpp"
#include "assert.hpp"

namespace belfem
{
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
        Vector( const std::initializer_list<T> & aInitList ) :
                mVector( arma::strans( arma::Mat<T>( { aInitList } )))
        {}

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
         * empty destructor
         */
        virtual ~Vector() = default;

//------------------------------------------------------------------------------
// MEMORY
//------------------------------------------------------------------------------

        /**
         * expose the underlying raw pointer ( writable version )
         */
        inline T *
        data()
        {
            return mVector.memptr();
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * expose the underlying raw pointer ( const version )
         */
        inline const T *
        data() const
        {
            return mVector.memptr();
        }

//------------------------------------------------------------------------------

        /**
        * expose the underlying matrix implementation ( writable version )
        */
        inline VectorType &
        vector_data()
        {
            return mVector;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * expose the underlying matrix implementation ( const version )
         */
        inline const VectorType &
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
        inline void
        fill( const T & aValue )
        {
            mVector.fill( aValue );
        }

//------------------------------------------------------------------------------

        /**
         * change the size of the vector
         */
        inline void
        set_size( const size_t aNumRows )
        {
            mVector.set_size( aNumRows, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline void
        set_size( const size_t aNumRows, const T & aValue )
        {
            this->set_size( aNumRows );
            this->fill( aValue );
        }

//------------------------------------------------------------------------------

        /**
         * get the length of the vector
         */
        inline size_t
        length() const
        {
            return mVector.n_rows;
        }

//------------------------------------------------------------------------------
// ACCESS OPERATORS
//------------------------------------------------------------------------------

        inline T &
        operator()( const size_t aIndex )
        {
            BELFEM_ASSERT( aIndex < this->length(),
                          "Index %lu out of bounds, which must be smaller than %lu.",
                          ( long unsigned int ) aIndex,
                          ( long unsigned int ) this->length());

            return mVector( aIndex );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline const T &
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

        inline auto
        begin() -> decltype( mVector.begin() )
        {
            return mVector.begin();
        }

        inline auto
        end() -> decltype( mVector.end() )
        {
            return mVector.end();
        }

        inline auto
        begin() const -> decltype( mVector.begin() )
        {
            return mVector.begin();
        }

        inline auto
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

        inline void
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
#endif //BELFEM_CL_AR_VECTOR_HPP
