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

#ifndef BELFEM_CL_TENSOR_HPP
#define BELFEM_CL_TENSOR_HPP

#include <cstring>    // for std::memcpy
#include <algorithm>  // for std::sort std::fill_n
#include <functional> // for std::plus std::minus

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "fn_TR_mat_to_ten.hpp"
#include "fn_TR_ten_to_mat.hpp"
#include "fn_TR_contract42.hpp"
#include "fn_TR_contract44.hpp"
#include "fn_TR_equal_equal.hpp"

#include "fn_TR_fill.hpp"
#include "fn_compliance_matrix.hpp"

#include "fn_inv.hpp"
namespace belfem
{
    /**
     * @brief Third- or fourth-order tensor container; the constitutive helpers
     *        (contraction, rotation, inversion, Voigt conversion) are specialised to 3x3x3x3.
     *
     * @ingroup grp_math_tensor
     * @see @ref math_tensor_tensor_usage_guide
     */
    template< typename T >
    class Tensor
    {
//----------------------------------------------------------------------------

        // indices
        const index_t mSizeI;
        const index_t mSizeJ;
        const index_t mSizeK;
        const index_t mSizeL;

        // memory jumps
        const index_t mOffsetJ;
        const index_t mOffsetK;
        const index_t mOffsetL;

        const index_t mOrder;
        const index_t mCapacity;

        T * mData = nullptr;

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        // create an empty tensor
        Tensor( const index_t aSizeI,
                const index_t aSizeJ,
                const index_t aSizeK,
                const index_t aSizeL ) :
            mSizeI( aSizeI ),
                    mSizeJ( aSizeJ ),
                    mSizeK( aSizeK ),
                    mSizeL( aSizeL ),
                    mOffsetJ( aSizeI ),
                    mOffsetK( aSizeI * aSizeJ ),
                    mOffsetL( aSizeI * aSizeJ * aSizeK ),
            mOrder( 4 ),
            mCapacity( aSizeI * aSizeJ * aSizeK * aSizeL )
        {
            BELFEM_ASSERT( mCapacity > 0, "Tensor sizes must not be zero" );

            mData = ( T * ) malloc( ( mCapacity ) * sizeof( T ) );
        }

        // create an empty tensor
        Tensor( const index_t aSizeI,
                const index_t aSizeJ,
                const index_t aSizeK ) :
            mSizeI( aSizeI ),
                    mSizeJ( aSizeJ ),
                    mSizeK( aSizeK ),
                    mSizeL( 1 ),
                    mOffsetJ( aSizeI ),
                    mOffsetK( aSizeI * aSizeJ ),
                    mOffsetL( aSizeI * aSizeJ * aSizeK ),
            mOrder( 3 ),
            mCapacity( aSizeI * aSizeJ * aSizeK )
        {
            BELFEM_ASSERT( mCapacity > 0, "Tensor sizes must not be zero" );

            mData = ( T * ) malloc( ( mCapacity ) * sizeof( T ) );
        }

        // create an empty tensor and initialize it with values
        Tensor( const index_t aSizeI,
                const index_t aSizeJ,
                const index_t aSizeK,
                const index_t aSizeL,
                const real aValue ) :
                Tensor( aSizeI, aSizeJ, aSizeK, aSizeL )
        {
            this->fill( aValue );
        }

        // create a tensor form an elasticity matrix
        Tensor( const Matrix< real > & aElasticityMatrix ) :
            Tensor( 3, 3, 3, 3 )
        {
            BELFEM_ASSERT(    aElasticityMatrix.n_rows() == 6
                          && aElasticityMatrix.n_cols() == 6,
                          "Matrix must be allocated as 6x6" );

            tensor::mat_to_ten( aElasticityMatrix, mData );
        }

//----------------------------------------------------------------------------

        // copy constructor
        Tensor( const Tensor< T > & aOther ) :
                mSizeI( aOther.size_i() ),
                mSizeJ( aOther.size_j() ),
                mSizeK( aOther.size_k() ),
                mSizeL( aOther.size_l() ),
                mOffsetJ( aOther.mOffsetJ ),
                mOffsetK( aOther.mOffsetK ),
                mOffsetL( aOther.mOffsetL ),
                mOrder( aOther.mOrder ),
                mCapacity( aOther.mCapacity )
        {
            mData = ( T * ) malloc( mCapacity * sizeof( T ) );
            std::memcpy( mData, aOther.mData, mCapacity * sizeof( T ) );
        }

//----------------------------------------------------------------------------

        // move constructor
        Tensor( Tensor< T > && aOther ) :
                mSizeI( aOther.size_i() ),
                mSizeJ( aOther.size_j() ),
                mSizeK( aOther.size_k() ),
                mSizeL( aOther.size_l() ),
                mOffsetJ( aOther.mOffsetJ ),
                mOffsetK( aOther.mOffsetK ),
                mOffsetL( aOther.mOffsetL ),
                mOrder( aOther.mOrder ),
                mCapacity( aOther.mCapacity )
        {
            mData = aOther.mData ;
            aOther.mData = nullptr ;
        }

//----------------------------------------------------------------------------

        ~Tensor()
        {
            free( mData );
        }
        
//------------------------------------------------------------------------------
// SIZES
//------------------------------------------------------------------------------
        
        index_t
        size_i() const
        {
            return mSizeI ;
        }
        
        index_t
        size_j() const
        {
            return mSizeJ ;
        }

        index_t
        size_k() const
        {
            return mSizeK ;
        }

        index_t
        size_l() const
        {
            return mSizeL ;
        }
        index_t
        order() const
        {
            return mOrder ;
        }

//------------------------------------------------------------------------------
// MEMORY
//------------------------------------------------------------------------------

        /**
         * expose the underlying raw pointer
         */
        T *
        data()
        {
            return mData ;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * expose the underlying raw pointer ( const version )
         */
        const T *
        data() const
        {
            return mData ;
        }

//------------------------------------------------------------------------------
// UTILITIES
//------------------------------------------------------------------------------

        /**
         * fill all values
         */
        void
        fill( const T aValue )
        {
            // populate pointer array
            std::fill_n( mData, mCapacity, aValue );
        }

//------------------------------------------------------------------------------

        /**
         * fill tensor in an isotropic way
         */
        void
        fill( const T aA, const T aB )
        {
            BELFEM_ASSERT( this->is_3333(), "Tensor must be 3x3x3x3" );
            tensor::fill( mData, aA, aB );
        }
        
//------------------------------------------------------------------------------

        /**
         * special funcition to create an isotropic elasticity tensor
         */
        void
        fill_isotropic_elasticity( const T E, const T nu )
        {
            BELFEM_ASSERT( this->is_3333(), "Tensor must be 3x3x3x3" );

            // bulk modulus
            T K  = E / ( 3. * ( 1.0 - 2.0 * nu ) );

            // shear modulus
            T G  = E / ( 2.0 * ( 1.0 + nu ) );

            tensor::fill( mData, K, G );
        }

//------------------------------------------------------------------------------

        void
        fill_orthotropic_elasticity(
                const T aYoung1,
                const T aYoung2,
                const T aYoung3,
                const T aPoisson23,
                const T aPoisson13,
                const T aPoisson12,
                const T aShear23,
                const T aShear13,
                const T aShear12
                )
        {
            BELFEM_ASSERT( this->is_3333(), "Tensor must be 3x3x3x3" );

            // compliance matrix
            Matrix< T > tS( 6, 6 );
            compliance_matrix( aYoung1, aYoung2, aYoung3,
                               aPoisson23, aPoisson13, aPoisson12,
                               aShear23, aShear13, aShear12, tS );

            // elasticity matrix
            Matrix< T > tC = inv( tS );

            // polulate data
            tensor::mat_to_ten( tC, mData );
        }

//------------------------------------------------------------------------------

        /**
         * memory size
         */
        index_t
        capacity() const
        {
            return mCapacity ;
        }

//------------------------------------------------------------------------------

        /**
         * returns true if this is a 3x3x3x3 tensor
         */
         bool
         is_3333() const
         {
             return mSizeI == 3 && mSizeJ == 3 && mSizeK == 3 && mSizeL == 3 ;
         }

//------------------------------------------------------------------------------

        void
        print( const string aLabel="Tensor")
        {
             index_t tCount = 0 ;
             fprintf( stdout, "    %s:\n", aLabel.c_str() );

             for( index_t l=0; l<mSizeL; ++l )
             {
                 for( index_t k=0; k<mSizeK; ++k )
                 {
                     for( index_t j=0; j<mSizeJ; ++j )
                     {
                         for( index_t i=0; i<mSizeI; ++i )
                         {
                            fprintf( stdout, "    %u : ( %u, %u, %u, %u ) = %12.3f\n",
                                     ( unsigned int ) tCount,
                                     ( unsigned int ) i,
                                     ( unsigned int ) j,
                                     ( unsigned int ) k,
                                     ( unsigned int ) l,
                                     ( double ) mData[ tCount ] );
                            ++tCount ;
                         }
                     }
                 }
             }
        }

//------------------------------------------------------------------------------
// ACCESS OPERATORS
//------------------------------------------------------------------------------

        /**
        * access operator ( writable version )
        */
        T &
        operator()(
                const index_t I,
                const index_t J,
                const index_t K )
        {
            BELFEM_ASSERT( mOrder == 3, "Tensor order must be of order 3" );

            BELFEM_ASSERT( I < mSizeI, "Index i out of bounds ( %u vs %u )",
                          ( unsigned int ) I,
                          ( unsigned int ) mSizeI );

            BELFEM_ASSERT( J < mSizeJ, "Index j out of bounds ( %u vs %u )",
                          ( unsigned int ) J,
                          ( unsigned int ) mSizeJ );

            BELFEM_ASSERT( K < mSizeK, "Index k out of bounds ( %u vs %u )",
                          ( unsigned int ) K,
                          ( unsigned int ) mSizeK );

            return mData[ K * mOffsetK + J * mOffsetJ + I ] ;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * access operator ( const version )
         */
        const T &
        operator()(
                const index_t I,
                const index_t J,
                const index_t K ) const
        {
            BELFEM_ASSERT( mOrder == 3, "Tensor order must be of order 3" );

            BELFEM_ASSERT( I < mSizeI, "Index i out of bounds ( %u vs %u )",
                          ( unsigned int ) I,
                          ( unsigned int ) mSizeI );

            BELFEM_ASSERT( J < mSizeJ, "Index j out of bounds ( %u vs %u )",
                          ( unsigned int ) J,
                          ( unsigned int ) mSizeJ );

            BELFEM_ASSERT( K < mSizeK, "Index k out of bounds ( %u vs %u )",
                          ( unsigned int ) K,
                          ( unsigned int ) mSizeK );

            return mData[ K * mOffsetK + J * mOffsetJ + I ] ;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * access operator ( writable version )
         */
        T &
        operator()(
                const index_t I,
                const index_t J,
                const index_t K,
                const index_t L )
        {
            BELFEM_ASSERT( mOrder == 4, "Tensor order must be of order 4" );

            BELFEM_ASSERT( I < mSizeI, "Index i out of bounds ( %u vs %u )",
                          ( unsigned int ) I,
                          ( unsigned int ) mSizeI );

            BELFEM_ASSERT( J < mSizeJ, "Index j out of bounds ( %u vs %u )",
                          ( unsigned int ) J,
                          ( unsigned int ) mSizeJ );

            BELFEM_ASSERT( K < mSizeK, "Index k out of bounds ( %u vs %u )",
                          ( unsigned int ) K,
                          ( unsigned int ) mSizeK );

            BELFEM_ASSERT( L < mSizeL, "Index l out of bounds ( %u vs %u )",
                          ( unsigned int ) L,
                          ( unsigned int ) mSizeL );

            return mData[ L * mOffsetL + K * mOffsetK + J * mOffsetJ + I ] ;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * access operator ( const version )
         */
        const T &
                operator()(
                        const index_t I,
                        const index_t J,
                        const index_t K,
                        const index_t L ) const
        {
            BELFEM_ASSERT( mOrder == 4, "Tensor order must be of order 4" );

            BELFEM_ASSERT( I < mSizeI, "Index i out of bounds ( %u vs %u )",
                          ( unsigned int ) I,
                          ( unsigned int ) mSizeI );

            BELFEM_ASSERT( J < mSizeJ, "Index j out of bounds ( %u vs %u )",
                          ( unsigned int ) J,
                          ( unsigned int ) mSizeJ );

            BELFEM_ASSERT( K < mSizeK, "Index k out of bounds ( %u vs %u )",
                          ( unsigned int ) K,
                          ( unsigned int ) mSizeK );

            BELFEM_ASSERT( L < mSizeL, "Index l out of bounds ( %u vs %u )",
                          ( unsigned int ) L,
                          ( unsigned int ) mSizeL );

            return mData[ L * mOffsetL + K * mOffsetK + J * mOffsetJ + I ] ;
        }

//------------------------------------------------------------------------------
// Equal Operators
//------------------------------------------------------------------------------

        /**
         * copy assignment
         */
        Tensor< T > &
        operator=( const Tensor< T > & aTensor )
        {
            if( this == &aTensor ) return *this;

            BELFEM_ASSERT( mOrder == aTensor.order(),
                          "Tensor orders must match for copy assignment" );

            BELFEM_ASSERT( mSizeI == aTensor.size_i(),
                       "Tensor dimensions must match for copy assignment" );

            BELFEM_ASSERT( mSizeJ == aTensor.size_j(),
                         "Tensor dimensions must match for copy assignment" );

            BELFEM_ASSERT( mSizeK == aTensor.size_k(),
                         "Tensor dimensions must match for copy assignment" );

            BELFEM_ASSERT( mSizeL == aTensor.size_l(),
                         "Tensor dimensions must match for copy assignment" );

            std::memcpy( mData, aTensor.data(), mCapacity * sizeof( T ) );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * move assignment
         */
        Tensor< T > &
        operator=( Tensor< T > && aTensor )
        {
            if( this == &aTensor ) return *this;

            BELFEM_ASSERT( mOrder == aTensor.order(),
                          "Tensor orders must match for move assignment" );

            BELFEM_ASSERT( mSizeI == aTensor.size_i(),
                          "Tensor dimensions must match for move assignment" );

            BELFEM_ASSERT( mSizeJ == aTensor.size_j(),
                         "Tensor dimensions must match for move assignment" );

            BELFEM_ASSERT( mSizeK == aTensor.size_k(),
                         "Tensor dimensions must match for move assignment" );

            BELFEM_ASSERT( mSizeL == aTensor.size_l(),
                         "Tensor dimensions must match for move assignment" );

            free( mData );
            mData = aTensor.mData ;
            aTensor.mData = nullptr ;
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Tensor< T > &
        operator=( const T & aScalar )
        {
            this->fill( aScalar );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Tensor< T > &
        operator=( const Matrix< T > & aElasticityMatrix )
        {
            BELFEM_ASSERT( this->is_3333(), "The tensor must be a 3x3x3x3 tensor" );

            tensor::mat_to_ten( aElasticityMatrix, mData );

            return *this;
        }

//------------------------------------------------------------------------------
// Addition operators
//------------------------------------------------------------------------------

        Tensor< T > &
        operator+=( const T & aScalar )
        {
            std::for_each( mData, mData + mCapacity,
                           [ aScalar ]( T & tVal )
                           { tVal += aScalar; } );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Tensor< T > &
        operator+=( const Tensor< T > & aTensor )
        {
            BELFEM_ASSERT( mOrder == aTensor.order(),
                       "Tensor orders must match for addition operator" );

            BELFEM_ASSERT( mSizeI == aTensor.size_i(),
                         "Tensor dimensions must match for addition operator" );

            BELFEM_ASSERT( mSizeJ == aTensor.size_j(),
                         "Tensor dimensions must match for addition operator" );

            BELFEM_ASSERT( mSizeK == aTensor.size_k(),
                         "Tensor dimensions must match for addition operator" );

            BELFEM_ASSERT( mSizeL == aTensor.size_l(),
                         "Tensor dimensions must match for addition operator" );

            std::transform( mData, mData + mCapacity, aTensor.data(),
                            mData, std::plus< T >() );
            return *this;
        }

//------------------------------------------------------------------------------
// Subtraction operators
//------------------------------------------------------------------------------

        Tensor< T > &
        operator-=( const T & aScalar )
        {
            std::for_each( mData, mData + mCapacity,
                           [ aScalar ]( T & tVal )
                           { tVal -= aScalar; } );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Tensor< T > &
        operator-=( const Tensor< T > & aTensor )
        {
            BELFEM_ASSERT( mOrder == aTensor.order(),
                       "Tensor orders must match for subtraction operator" );

            BELFEM_ASSERT( mSizeI == aTensor.size_i(),
                           "Tensor dimensions must match for subtraction operator" );

            BELFEM_ASSERT( mSizeJ == aTensor.size_j(),
                         "Tensor dimensions must match for subtraction operator" );

            BELFEM_ASSERT( mSizeK == aTensor.size_k(),
                         "Tensor dimensions must match for subtraction operator" );

            BELFEM_ASSERT( mSizeL == aTensor.size_l(),
                         "Tensor dimensions must match for subtraction operator" );

            std::transform( mData, mData + mCapacity, aTensor.data(),
                            mData, std::minus< T >() );
            return *this;
        }

//------------------------------------------------------------------------------
// Multiplication
//------------------------------------------------------------------------------

        Tensor< T > &
        operator*=( const T & aScalar )
        {
            std::for_each( mData, mData + mCapacity,
                           [ aScalar ]( T & tVal )
                           { tVal *= aScalar; } );
            return *this;
        }

//------------------------------------------------------------------------------
// Division
//------------------------------------------------------------------------------

        Tensor< T > &
        operator/=( const T & aScalar )
        {
            std::for_each( mData, mData + mCapacity,
                           [ aScalar ]( T & tVal )
                           { tVal /= aScalar; } );
            return *this;
        }

//------------------------------------------------------------------------------
//   Contraction
//------------------------------------------------------------------------------

        /**
         * contract with  other tensor
         * @param aB
         * @param aC
         */
        inline
        void
        ddot( const Tensor< T > & aB, Tensor< T > & aC )
        {
            BELFEM_ASSERT( this->is_3333(), "operating tensor must be 3x3x3x3" );
            BELFEM_ASSERT( aB.is_3333(), "argument tensor must be 3x3x3x3" );
            BELFEM_ASSERT( aC.is_3333(), "target tensor must be 3x3x3x3" );

            tensor::contract44( mData, aB.data(), aC.data() );
        }

//------------------------------------------------------------------------------

        /**
         * contract with 3x3 matrox
         * @param aB
         * @param aC
         */
        inline
        void
        ddot( const Matrix< T > & aB, Matrix< T > & aC )
        {
            BELFEM_ASSERT( this->is_3333(),
                          "operating tensor must be 3x3x3x3" );

            BELFEM_ASSERT(      aB.n_rows() == 3
                            && aB.n_cols() == 3,
                            "argument matrix must be allocated as 3x3" );

            BELFEM_ASSERT(      aC.n_rows() == 3
                               && aC.n_cols() == 3,
                               "target matrix must be allocated as 3x3" );

            tensor::contract42( mData, aB, aC );
        }

//------------------------------------------------------------------------------
// Conversion
//------------------------------------------------------------------------------

        /**
         * converts a tensor to the elastitity matrix in Voigt notation
         * @param aMatrix
         */
        void
        to_matrix( Matrix< real > & aMatrix )
        {
            BELFEM_ASSERT(   aMatrix.n_rows() == 6
                         && aMatrix.n_cols() == 6,
                 "Matrix must be allocated as 6x6" );

            BELFEM_ASSERT( this->is_3333(), "The tensor must be a 3x3x3x3 tensor" );

            tensor::ten_to_mat( mData, aMatrix );
        }
    };
//------------------------------------------------------------------------------

    template< typename T >
    Tensor< T >
    operator+( const Tensor< T > & aB,
               const Tensor< T > & aC )
    {
        BELFEM_ASSERT( aB.order() == aC.order(),
                      "Tensor orders must match for addition operator" );

        BELFEM_ASSERT( aB.size_i() == aC.size_i(),
                       "Tensor dimensions must match for addition operator" );

        BELFEM_ASSERT( aB.size_j() == aC.size_j(),
                       "Tensor dimensions must match for addition operator" );

        BELFEM_ASSERT( aB.size_k() == aC.size_k(),
                       "Tensor dimensions must match for addition operator" );

        BELFEM_ASSERT( aB.size_l() == aC.size_l(),
                       "Tensor dimensions must match for addition operator" );

        Tensor< T > aA( aB );
        aA += aC;

        return aA ;
    }

//------------------------------------------------------------------------------

    template< typename T >
    Tensor< T >
    operator-( const Tensor< T > & aB,
               const Tensor< T > & aC )
    {
        BELFEM_ASSERT( aB.order() == aC.order(),
                             "Tensor orders must match for subtraction operator" );

        BELFEM_ASSERT( aB.size_i() == aC.size_i(),
                       "Tensor dimensions must match for subtraction operator" );

        BELFEM_ASSERT( aB.size_j() == aC.size_j(),
                       "Tensor dimensions must match for subtraction operator" );

        BELFEM_ASSERT( aB.size_k() == aC.size_k(),
                       "Tensor dimensions must match for subtraction operator" );

        BELFEM_ASSERT( aB.size_l() == aC.size_l(),
                       "Tensor dimensions must match for subtraction operator" );

        Tensor< T > aA( aB );
        aA -= aC;
        return aA ;
    }

//------------------------------------------------------------------------------
// contraction
//------------------------------------------------------------------------------

    template< typename T >
    Tensor< T >
    operator%( const Tensor< T > & aA,
               const Tensor< T > & aB )
    {
        BELFEM_ASSERT( aA.is_3333() && aB.is_3333(),
                      "Both tensors must be of size 3x3x3x3" );

        Tensor< T > aC( 3, 3, 3, 3 );
        tensor::contract44( aA.data(), aB.data(), aC.data() );
        return aC ;
    }

//------------------------------------------------------------------------------

    template< typename T >
    Matrix< T >
    operator%( const Tensor< T > & aA,
               const Matrix< T > & aB )
    {
        BELFEM_ASSERT( aA.is_3333(),
                      "Tensor A must be of size 3x3x3x3" );

        BELFEM_ASSERT( aB.n_rows() == 3 && aB.n_cols() == 3,
            "when contracting a 3x3x3x3 tensor with a matrix, latter one must be 3x3" );

        Matrix< T > aC ( 3, 3 );
        tensor::contract42( aA.data(), aB, aC );
        return aC ;
    }

//------------------------------------------------------------------------------

    template< typename T >
    bool
    operator==( const Tensor< T > & aA,
                const Tensor< T > & aB )
    {
        if(    aA.order()  != aB.order()
            || aA.size_i() != aB.size_i()
            || aA.size_j() != aB.size_j()
            || aA.size_k() != aB.size_k()
            || aA.size_l() != aB.size_l() ) return false;

        return tensor::equal_equal( aA.data(), aB.data(), aA.capacity() );
    }

//------------------------------------------------------------------------------
} /* namespace belfem */
#endif //BELFEM_CL_TENSOR_HPP
