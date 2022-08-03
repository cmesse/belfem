//
// Created by Christian Messe on 17.12.20.
//

#ifndef BELFEM_CL_TENSOR_HPP
#define BELFEM_CL_TENSOR_HPP

#include <cstring>    // for std::memcpy
#include <algorithm>  // for std::sort std::fill_n

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "fn_TR_mat_to_ten.hpp"
#include "fn_TR_ten_to_mat.hpp"
#include "fn_TR_add44.hpp"
#include "fn_TR_subtract44.hpp"
#include "fn_TR_contract42.hpp"
#include "fn_TR_contract44.hpp"
#include "fn_TR_equal_equal.hpp"
#include "fn_TR_fill_isotropic.hpp"
#include "fn_compliance_matrix.hpp"
#include "fn_TR_kelvin_christoffel.hpp"
#include "fn_inv.hpp"
namespace belfem
{
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

            tensor::mat_to_ten( mData, aElasticityMatrix.data() );
        }

//----------------------------------------------------------------------------

        ~Tensor()
        {
            free( mData );
        }

//------------------------------------------------------------------------------
// MEMORY
//------------------------------------------------------------------------------

        /**
         * expose the underlying raw pointer
         */
        inline T *
        data()
        {
            return mData ;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * expose the underlying raw pointer ( const version )
         */
        inline const T *
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
        inline void
        fill( const T aValue )
        {
            // populate pointer array
            std::fill_n( mData, 81, aValue );
        }

//------------------------------------------------------------------------------

        /**
         * fill tensor in an isotropic way
         */
        inline void
        fill( const T aKappa, const T aGamma )
        {
            BELFEM_ASSERT( this->is_3333(), "Tensor must be 3x3x3x3" );
            tensor::fill_isotropic( mData, aKappa, aGamma / 3.0 );
        }

//------------------------------------------------------------------------------

        /**
         * special funcition to create an isotropic elasticity tensor
         */
        inline void
        fill_isotropic_elasticity( const T aYoung, const T aPoisson )
        {
            BELFEM_ASSERT( this->is_3333(), "Tensor must be 3x3x3x3" );

            // fill with K=Bulk Modulus, Gamma = Shear Modulus
            tensor::fill_isotropic(
                    mData,
                    aYoung / ( 3. - 6. * aPoisson ) ,
                    aYoung / ( 6. + 6. * aPoisson ) );
        }

//------------------------------------------------------------------------------

        inline void
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
            tensor::mat_to_ten( mData, tC.data() );
        }

//------------------------------------------------------------------------------

        /**
         * memory size
         */
        inline index_t
        capacity() const
        {
            return mCapacity ;
        }

//------------------------------------------------------------------------------

        /**
         * returns true if this is a 3x3x3x3 tensor
         */
         inline bool
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
                     for( index_t j=0; j<mSizeL; ++j )
                     {
                         for( index_t i=0; i<mSizeL; ++i )
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
        inline T &
        operator()(
                const index_t I,
                const index_t J,
                const index_t K,
                const index_t L )
        {
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
        inline const T &
                operator()(
                        const index_t I,
                        const index_t J,
                        const index_t K,
                        const index_t L ) const
        {
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
         * Copy constructor
         */
        inline Tensor< T > &
        operator=( const Tensor< T > & aTensor )
        {
            std::memcpy( mData, aTensor.data(), ( mCapacity ) * sizeof( T ) );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline Tensor< T > &
        operator=( const T & aScalar )
        {
            this->fill( aScalar );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline Tensor< T > &
        operator=( const Matrix< T > & aElasticityMatrix )
        {
            BELFEM_ASSERT( this->is_3333(), "The tensor must be a 3x3x3x3 tensor" );

            tensor::mat_to_ten( mData, aElasticityMatrix.data() );
        }

//------------------------------------------------------------------------------
// Addition operators
//------------------------------------------------------------------------------

        inline Tensor< T > &
        operator+=( const T & aScalar )
        {
            std::for_each( mData, mData + mCapacity,
                           [ aScalar ]( T & tVal )
                           { tVal += aScalar; } );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline Tensor< T > &
        operator+=( const Tensor< T > & aTensor )
        {
            BELFEM_ASSERT( this->is_3333() && aTensor.is_3333(),
                          "Both tensors must be of size 3x3x3x3" );

            tensor::add( mData, aTensor.data() );
            return *this;
        }

//------------------------------------------------------------------------------
// Subtraction operators
//------------------------------------------------------------------------------

        inline Tensor< T > &
        operator-=( const T & aScalar )
        {
            std::for_each( mData, mData + mCapacity,
                           [ aScalar ]( T & tVal )
                           { tVal -= aScalar; } );
            return *this;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline Tensor< T > &
        operator-=( const Tensor< T > & aTensor )
        {
            BELFEM_ASSERT( this->is_3333() && aTensor.is_3333(),
                          "Both tensors must be of size 3x3x3x3" );

            tensor::subtract( mData, aTensor.data() );
            return *this;
        }

//------------------------------------------------------------------------------
// Multiplication
//------------------------------------------------------------------------------

        inline Tensor< T > &
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

        inline Tensor< T > &
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

            tensor::contract44( mData, aB.data(), aC.data() );
        }

//------------------------------------------------------------------------------
// Conversion
//------------------------------------------------------------------------------

        /**
         * converts a tensor to the elastitity matrix in Voigt notation
         * @param aMatrix
         */
        inline void
        to_matrix( Matrix< real > & aMatrix )
        {
            BELFEM_ASSERT(   aMatrix.n_rows() == 6
                         && aMatrix.n_cols() == 6,
                 "Matrix must be allocated as 6x6" );

            BELFEM_ASSERT( this->is_3333(), "The tensor must be a 3x3x3x3 tensor" );

            tensor::ten_to_mat( aMatrix.data(), mData );
        }
    };
//------------------------------------------------------------------------------

    template< typename T >
    inline Tensor< T >
    operator+( const Tensor< T > & aB,
               const Tensor< T > & aC )
    {
        BELFEM_ASSERT( aB.is_3333() && aC.is_3333(),
                      "Both tensors must be of size 3x3x3x3" );

        Tensor< T > aA ;
        tensor::add( aA.data(), aB.data(), aC.data() );
        return aA ;
    }

//------------------------------------------------------------------------------

    template< typename T >
    inline Tensor< T >
    operator-( const Tensor< T > & aB,
               const Tensor< T > & aC )
    {
        BELFEM_ASSERT( aB.is_3333() && aC.is_3333(),
                      "Both tensors must be of size 3x3x3x3" );
        Tensor< T > aA( 3, 3, 3, 3 );
        tensor::subtract( aA.data(), aB.data(), aC.data() );
        return aA ;
    }

//------------------------------------------------------------------------------
// contraction
//------------------------------------------------------------------------------

    template< typename T >
    inline Tensor< T >
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
    inline Matrix< T >
    operator%( const Tensor< T > & aA,
               const Matrix< T > & aB )
    {
        BELFEM_ASSERT( aA.is_3333(),
                      "Tensor A must be of size 3x3x3x3" );

        BELFEM_ASSERT( aB.n_rows() == 3 && aB.n_cols() == 3,
            "when contracting a 3x3x3x3 tensor with a matrix, latter one must be 3x3" );

        Matrix< T > aC ( 3, 3 );
        tensor::contract42( aA.data(), aB.data(), aC.data() );
        return aC ;
    }

//------------------------------------------------------------------------------

    template< typename T >
    inline bool
    operator==( const Tensor< T > & aA,
                const Tensor< T > & aB )
    {
        return tensor::equal_equal( aA.data(), aB.data(), aA.capacity() );
    }

//------------------------------------------------------------------------------
} /* namespace belfem */
#endif //BELFEM_CL_TENSOR_HPP
