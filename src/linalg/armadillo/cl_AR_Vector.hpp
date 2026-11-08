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
#include <type_traits>
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
    class Vector final
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

        /**
         * the orientation contract of a column vector: a row view, an
         * element-wise expression or a transpose whose evaluated shape is
         * 1 x n ( n > 1 ) is transposed into the column. Only operands whose
         * Proxy is lazy are inspected ( eOp, eGlue, mtOp, mtGlue, the two
         * transpose Ops ); a product, a Mat, an arma::Row or any other view
         * is handed to the backend as it is, and in a debug build an
         * assignment that leaves a row behind is asserted.
         */
        template < typename ET >
        struct IsTransposeOp { static constexpr bool value = false; };

        template < typename T1 >
        struct IsTransposeOp< arma::Op< T1, arma::op_htrans > > { static constexpr bool value = true; };

        template < typename T1 >
        struct IsTransposeOp< arma::Op< T1, arma::op_strans > > { static constexpr bool value = true; };

        template < typename ET >
        struct IsOrientedOperand
        {
            static constexpr bool value =
                arma::is_eOp< ET >::value || arma::is_eGlue< ET >::value ||
                arma::is_mtOp< ET >::value || arma::is_mtGlue< ET >::value ||
                IsTransposeOp< ET >::value ;
        };

        template < typename ET >
        typename std::enable_if< ! IsOrientedOperand< ET >::value, void >::type
        assign_oriented( const ET & aExpression )
        {
            mVector = aExpression;
            BELFEM_ASSERT( mVector.n_cols == 1 || mVector.n_elem == 0,
                           "Vector: a %lu x %lu result was assigned; products and raw views are not transposed",
                           ( unsigned long ) mVector.n_rows, ( unsigned long ) mVector.n_cols );
        }

        template < typename ET >
        typename std::enable_if< IsOrientedOperand< ET >::value, void >::type
        assign_oriented( const ET & aExpression )
        {
            const arma::Proxy< ET > tP( aExpression );
            if ( tP.get_n_rows() == 1 && tP.get_n_cols() > 1 )
            {
                mVector = arma::strans( aExpression );
            }
            else
            {
                mVector = aExpression;
            }
        }

        template < typename ET >
        typename std::enable_if< ! IsOrientedOperand< ET >::value, void >::type
        plus_oriented( const ET & aExpression ) { mVector += aExpression; }

        template < typename ET >
        typename std::enable_if< IsOrientedOperand< ET >::value, void >::type
        plus_oriented( const ET & aExpression )
        {
            const arma::Proxy< ET > tP( aExpression );
            if ( tP.get_n_rows() == 1 && tP.get_n_cols() > 1 )
            {
                mVector += arma::strans( aExpression );
            }
            else
            {
                mVector += aExpression;
            }
        }

        template < typename ET >
        typename std::enable_if< ! IsOrientedOperand< ET >::value, void >::type
        minus_oriented( const ET & aExpression ) { mVector -= aExpression; }

        template < typename ET >
        typename std::enable_if< IsOrientedOperand< ET >::value, void >::type
        minus_oriented( const ET & aExpression )
        {
            const arma::Proxy< ET > tP( aExpression );
            if ( tP.get_n_rows() == 1 && tP.get_n_cols() > 1 )
            {
                mVector -= arma::strans( aExpression );
            }
            else
            {
                mVector -= aExpression;
            }
        }

//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        Vector() = default;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector( const size_t aNumRows ) :
                mVector( aNumRows, 1 ) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector( const size_t aNumRows, const T & aValue ) :
                mVector( aNumRows, 1 )
        {
            this->fill( aValue );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

        Vector( const VectorType & aExpression ) :
        mVector( aExpression )
        {}

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
                mVector( arma::strans( aRow ) )
        {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * constructor from any Armadillo expression; a 1 x n result is a row
         * and is transposed into the column, as the row-view constructor does
         */
        template< typename ET >
        Vector( const arma::Base< T, ET > & aExpression )
        {
            this->assign_oriented( aExpression.get_ref() );
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
        ~Vector() = default;

//------------------------------------------------------------------------------
// MEMORY
//------------------------------------------------------------------------------

        /** Returns a borrowed raw buffer under the layout contract stated
         *  in cl_Vector.hpp; may be nullptr when empty. */
        T *
        data()
        {
            return mVector.memptr();
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

        template < typename ET,
                   typename = typename std::enable_if< arma::is_arma_type< ET >::value >::type >
        Vector< T > &
        operator=( const ET & aExpression )
        {
            this->assign_oriented( aExpression );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< T > &
        operator=( const arma::subview_row< T > & aRow )
        {
            mVector = arma::strans( aRow );
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

        /**
         * a row view is transposed into the column; strans, so complex
         * entries are not conjugated ( the same orientation Blaze gives )
         */
        Vector< T > &
        operator+=( const arma::subview_row< T > & aRow )
        {
            mVector += arma::strans( aRow );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET,
                   typename = typename std::enable_if< arma::is_arma_type< ET >::value >::type >
        Vector<T> &
        operator+=( const ET & aExpression )
        {
            this->plus_oriented( aExpression );
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

        Vector< T > &
        operator-=( const arma::subview_row< T > & aRow )
        {
            mVector -= arma::strans( aRow );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < typename ET,
                   typename = typename std::enable_if< arma::is_arma_type< ET >::value >::type >
        Vector<T> &
        operator-=( const ET & aExpression )
        {
            this->minus_oriented( aExpression );
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
