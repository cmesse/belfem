/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/**
 * @file
 * @brief Shared helpers for the LAPACK wrappers, such as leading_dimension().
 * @ingroup grp_linalg
 *
 * Not a solver itself.
 */

#ifndef BELFEM_LAPACKTOOLS_HPP
#define BELFEM_LAPACKTOOLS_HPP

#include <algorithm>
#include <complex>
#include <limits>
#include <type_traits>

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"

#ifdef BELFEM_MKL
#include <mkl_types.h>
#endif

#ifdef BELFEM_BLAZE
#include <blaze/math/blas/Types.h>
#endif


namespace belfem
{
    namespace lapack
    {
//------------------------------------------------------------------------------

        // helper for compile-time errors in unsupported template instantiations
        template< typename T >
        constexpr bool dependent_false = false ;

//------------------------------------------------------------------------------

        // underlying real type of a scalar: T for real T, T::value_type for
        // std::complex< T >. LAPACK keeps some arrays real valued even in
        // the complex flavors ( e.g. the singular values and rwork of gesvd )
        template< typename T >
        struct real_type { typedef T type ; };

        template< typename T >
        struct real_type< std::complex< T > > { typedef T type ; };

        template< typename T >
        using real_t = typename real_type< T >::type ;

        // complementary complex type of a scalar: std::complex<T> for real
        // T, T itself for std::complex<T> — e.g. the eigenvalues of geev
        // are complex valued in all flavors
        template< typename T >
        struct cplx_type { typedef std::complex< T > type ; };

        template< typename T >
        struct cplx_type< std::complex< T > > { typedef std::complex< T > type ; };

        template< typename T >
        using cplx_t = typename cplx_type< T >::type ;

//------------------------------------------------------------------------------

        // complex datatype for the extern "C" prototypes: MKL declares its
        // LAPACK interface with MKL_Complex8/16, everything else uses the
        // Fortran real-pair convention ( two consecutive floats/doubles ).
        // Under Blaze, Blaze's own clapack headers put real-pair prototypes
        // into every TU ( mkl_cblas.h does not define INTEL_MKL_VERSION, so
        // Blaze never suppresses them ) — the glue must match Blaze there,
        // also when MKL is the linked library
#if defined( BELFEM_MKL ) && ! defined( BELFEM_BLAZE )
        typedef MKL_Complex8  cplx_float_t ;
        typedef MKL_Complex16 cplx_double_t ;
#else
        typedef float  cplx_float_t ;
        typedef double cplx_double_t ;
#endif

        // hidden Fortran argument carrying the length of char* parameters
        // ( gfortran/ifort append it after the regular argument list;
        //   Blaze declares it explicitly, so we must match )
#ifdef BELFEM_BLAZE
        typedef blaze::fortran_charlen_t fortran_charlen_t ;
#else
        typedef size_t fortran_charlen_t ;
#endif

        // std::complex is guaranteed layout-compatible ( C++11, [complex.numbers] )
#if defined( BELFEM_MKL ) && ! defined( BELFEM_BLAZE )
        static_assert( sizeof( std::complex< float > ) == sizeof( cplx_float_t ),
            "sizes of std::complex< float > and MKL_Complex8 do not match" );

        static_assert( sizeof( std::complex< double > ) == sizeof( cplx_double_t ),
            "sizes of std::complex< double > and MKL_Complex16 do not match" );
#else
        static_assert( sizeof( std::complex< float > ) == 2UL * sizeof( cplx_float_t ),
            "std::complex< float > is not layout-compatible with LAPACK" );

        static_assert( sizeof( std::complex< double > ) == 2UL * sizeof( cplx_double_t ),
            "std::complex< double > is not layout-compatible with LAPACK" );
#endif

        // SCLS builds all third-party libraries with one consistent integer
        // width, so int_t is also the LAPACK integer. These checks back up
        // that invariant.
#ifdef BELFEM_BLAZE
        static_assert( sizeof( int_t ) == sizeof( blaze::blas_int_t ),
            "sizes of belfem::int_t and blaze::blas_int_t do not match" );
#endif

#if defined( BELFEM_MKL ) && ( defined( MKL_ILP64 ) != defined( BELFEM_INT64 ) )
#error "MKL_ILP64 and BELFEM_INT64 must be set together, check USE_MKL_64BIT_API"
#endif

#ifdef BELFEM_MKL
        static_assert( sizeof( int_t ) == sizeof( MKL_INT ),
            "sizes of belfem::int_t and MKL_INT do not match" );
#endif

//------------------------------------------------------------------------------

        /**
         * logical length of a vector operand, as passed to LAPACK as LDB;
         * vector storage is contiguous under both backends
         */
        template< typename T >
        int_t
        leading_dimension( const Vector< T > & A )
        {
            BELFEM_ASSERT( A.length() <= ( size_t ) std::numeric_limits< int_t >::max(),
                "vector length exceeds LAPACK integer range ( %lu )",
                ( long unsigned int ) A.length() );

            return std::max< int_t >( 1, ( int_t ) A.length() );
        }

//------------------------------------------------------------------------------

        /**
         * stride between the columns of an owning, column-major matrix,
         * as passed to LAPACK as LDA/LDB/LDC. This is a property of the
         * stored array only — it does NOT change when the operand is used
         * transposed ( trans only changes which logical dimension LAPACK
         * validates the stride against ). Under Blaze, columns are padded
         * for SIMD alignment, so the stride may exceed the row count —
         * passing n_rows() instead would silently corrupt the result.
         */
        template< typename T >
        int_t
        leading_dimension( const Matrix< T > & A )
        {
#ifdef BELFEM_ARMADILLO
            BELFEM_ASSERT( A.n_rows() <= ( size_t ) std::numeric_limits< int_t >::max(),
                "matrix row count exceeds LAPACK integer range ( %lu )",
                ( long unsigned int ) A.n_rows() );

            return std::max< int_t >( 1, ( int_t ) A.n_rows() );
#elif BELFEM_BLAZE
            BELFEM_ASSERT( A.matrix_data().spacing() <= ( size_t ) std::numeric_limits< int_t >::max(),
                "matrix column stride exceeds LAPACK integer range ( %lu )",
                ( long unsigned int ) A.matrix_data().spacing() );

            return std::max< int_t >( 1, ( int_t ) A.matrix_data().spacing() );
#else
#error "leading_dimension() : no matrix backend selected"
#endif
        }

//------------------------------------------------------------------------------

        /**
         * lapack reports the optimal work size in the first entry of the work
         * array, which stays real valued for the complex flavors
         */
        inline int_t
        work_size( const float & aValue )
        {
            return static_cast< int_t >( aValue );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline int_t
        work_size( const double & aValue )
        {
            return static_cast< int_t >( aValue );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline int_t
        work_size( const std::complex< float > & aValue )
        {
            return static_cast< int_t >( aValue.real() );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline int_t
        work_size( const std::complex< double > & aValue )
        {
            return static_cast< int_t >( aValue.real() );
        }
        
//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_LAPACKTOOLS_HPP
