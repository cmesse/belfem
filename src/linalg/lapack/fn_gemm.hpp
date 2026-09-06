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

/**
 * @file
 * @brief Matrix-matrix product C := alpha*op(A)*op(B) + beta*C (BLAS ?gemm).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_GEMM_HPP
#define BELFEM_FN_GEMM_HPP

#include "assert.hpp"
#include "lapacktools.hpp"

//------------------------------------------------------------------------------
namespace belfem
{
    namespace lapack
    {

#ifdef __cplusplus
        extern "C"
        {
#endif
//------------------------------------------------------------------------------

            void
            sgemm_(  char   * transa,
                     char   * transb,
                     int_t  * m,
                     int_t  * n,
                     int_t  * k,
                     float  * alpha,
                     float  * a,
                     int_t  * lda,
                     float  * b,
                     int_t  * ldb,
                     float  * beta,
                     float  * c,
                     int_t  * ldc,
                 fortran_charlen_t lta,
                 fortran_charlen_t ltb
                      );

//------------------------------------------------------------------------------

            void
            dgemm_(  char   * transa,
                     char   * transb,
                     int_t  * m,
                     int_t  * n,
                     int_t  * k,
                     double * alpha,
                     double * a,
                     int_t    * lda,
                     double * b,
                     int_t    * ldb,
                     double * beta,
                     double * c,
                     int_t  * ldc,
                 fortran_charlen_t lta,
                 fortran_charlen_t ltb );
//------------------------------------------------------------------------------

        void
        cgemm_(  char         * transa,
                 char         * transb,
                 int_t        * m,
                 int_t        * n,
                 int_t        * k,
                 cplx_float_t * alpha,
                 cplx_float_t * a,
                 int_t        * lda,
                 cplx_float_t * b,
                 int_t        * ldb,
                 cplx_float_t * beta,
                 cplx_float_t * c,
                 int_t        * ldc,
                 fortran_charlen_t lta,
                 fortran_charlen_t ltb );

//------------------------------------------------------------------------------

        void
        zgemm_(  char          * transa,
                 char          * transb,
                 int_t         * m,
                 int_t         * n,
                 int_t         * k,
                 cplx_double_t * alpha,
                 cplx_double_t * a,
                 int_t         * lda,
                 cplx_double_t * b,
                 int_t         * ldb,
                 cplx_double_t * beta,
                 cplx_double_t * c,
                 int_t         * ldc,
                 fortran_charlen_t lta,
                 fortran_charlen_t ltb );

//------------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

//------------------------------------------------------------------------------

        template< typename T >
        void
        gemm(    const char     * transa,
                 const char     * transb,
                 const int_t    * m,
                 const int_t    * n,
                 const int_t    * k,
                 const T        * alpha,
                 const T        * a,
                 const int_t    * lda,
                 const T        * b,
                 const int_t    * ldb,
                 const T        * beta,
                       T        * c,
                 const int_t    * ldc )
        {
            static_assert( dependent_false< T >,
                "gemm not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gemm(    const char     * transa,
                 const char     * transb,
                 const int_t    * m,
                 const int_t    * n,
                 const int_t    * k,
                 const float    * alpha,
                 const float    * a,
                 const int_t    * lda,
                 const float    * b,
                 const int_t    * ldb,
                 const float    * beta,
                       float    * c,
                 const int_t    * ldc )
        {
            sgemm_( const_cast< char * >  ( transa ),
                    const_cast< char * >  ( transb ),
                    const_cast< int_t * > ( m ),
                    const_cast< int_t * > ( n ),
                    const_cast< int_t * > ( k ),
                    const_cast< float * > ( alpha ),
                    const_cast< float * > ( a ),
                    const_cast< int_t * > ( lda ),
                    const_cast< float * > ( b ),
                    const_cast< int_t * > ( ldb ),
                    const_cast< float * > ( beta ),
                                            c,
                    const_cast< int_t * >  ( ldc ),
                    1, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gemm(    const char     * transa,
                 const char     * transb,
                 const int_t    * m,
                 const int_t    * n,
                 const int_t    * k,
                 const double   * alpha,
                 const double   * a,
                 const int_t    * lda,
                 const double   * b,
                 const int_t    * ldb,
                 const double   * beta,
                       double   * c,
                 const int_t    * ldc )
        {
            dgemm_( const_cast< char * >  ( transa ),
                    const_cast< char * >  ( transb ),
                    const_cast< int_t * > ( m ),
                    const_cast< int_t * > ( n ),
                    const_cast< int_t * > ( k ),
                    const_cast< double * >( alpha ),
                    const_cast< double * >( a ),
                    const_cast< int_t * > ( lda ),
                    const_cast< double * >( b ),
                    const_cast< int_t * > ( ldb ),
                    const_cast< double * >( beta ),
                                            c,
                    const_cast< int_t * > ( ldc ),
                    1, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gemm(    const char                       * transa,
                 const char                       * transb,
                 const int_t                      * m,
                 const int_t                      * n,
                 const int_t                      * k,
                 const std::complex< float >    * alpha,
                 const std::complex< float >    * a,
                 const int_t                      * lda,
                 const std::complex< float >    * b,
                 const int_t                      * ldb,
                 const std::complex< float >    * beta,
                       std::complex< float >    * c,
                 const int_t    * ldc )
        {
            cgemm_( const_cast< char * >  ( transa ),
                    const_cast< char * >  ( transb ),
                    const_cast< int_t * > ( m ),
                    const_cast< int_t * > ( n ),
                    const_cast< int_t * > ( k ),
                    const_cast< cplx_float_t * >( reinterpret_cast< const cplx_float_t * >( alpha ) ),
                    const_cast< cplx_float_t * >( reinterpret_cast< const cplx_float_t * >( a ) ),
                    const_cast< int_t * > ( lda ),
                    const_cast< cplx_float_t * >( reinterpret_cast< const cplx_float_t * >( b ) ),
                    const_cast< int_t * > ( ldb ),
                    const_cast< cplx_float_t * >( reinterpret_cast< const cplx_float_t * >( beta ) ),
                    reinterpret_cast< cplx_float_t * >( c ),
                    const_cast< int_t * > ( ldc ),
                    1, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gemm(    const char                       * transa,
                 const char                       * transb,
                 const int_t                      * m,
                 const int_t                      * n,
                 const int_t                      * k,
                 const std::complex< double >   * alpha,
                 const std::complex< double >   * a,
                 const int_t                      * lda,
                 const std::complex< double >   * b,
                 const int_t                      * ldb,
                 const std::complex< double >   * beta,
                       std::complex< double >   * c,
                 const int_t    * ldc )
        {
            zgemm_( const_cast< char * >  ( transa ),
                    const_cast< char * >  ( transb ),
                    const_cast< int_t * > ( m ),
                    const_cast< int_t * > ( n ),
                    const_cast< int_t * > ( k ),
                    const_cast< cplx_double_t * >( reinterpret_cast< const cplx_double_t * >( alpha ) ),
                    const_cast< cplx_double_t * >( reinterpret_cast< const cplx_double_t * >( a ) ),
                    const_cast< int_t * > ( lda ),
                    const_cast< cplx_double_t * >( reinterpret_cast< const cplx_double_t * >( b ) ),
                    const_cast< int_t * > ( ldb ),
                    const_cast< cplx_double_t * >( reinterpret_cast< const cplx_double_t * >( beta ) ),
                    reinterpret_cast< cplx_double_t * >( c ),
                    const_cast< int_t * > ( ldc ),
                    1, 1 );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */
//------------------------------------------------------------------------------

    /**
     * @brief matrix-matrix product C := alpha * op(A) * op(B) + beta * C
     *        via BLAS ?gemm, where op is the identity ( trans = 'N' ) or
     *        the ( conjugate ) transpose ( 'T' / 'C' )
     *
     * @param[in]     A      left operand, stored form
     * @param[in]     B      right operand, stored form
     * @param[in,out] C      result; resized to op(A) * op(B) when beta is
     *                       zero, otherwise it must already have that shape
     * @param[in]     alpha  scaling of the product ( default 1 )
     * @param[in]     beta   scaling of the accumulator ( default 0 )
     * @param[in]     transa 'N', 'T' or 'C' for op(A)
     * @param[in]     transb 'N', 'T' or 'C' for op(B)
     */
    template< typename T >
    void
    gemm(
        const Matrix< T > & A,
        const Matrix< T > & B,
        Matrix< T > & C,
        const T alpha = 1.0,
        const T beta = 0.0,
        const char transa = 'N',
        const char transb = 'N' )
    {
        // dimensions of op(A) ( m x k ) and op(B) ( k x n )
        const bool tNoTransA = ( transa == 'N' || transa == 'n' );
        const bool tNoTransB = ( transb == 'N' || transb == 'n' );

        BELFEM_ASSERT( transa == 'N' || transa == 'T' || transa == 'C' ||
                       transa == 'n' || transa == 't' || transa == 'c',
                        "unsupported transa flag '%c'", transa );

        BELFEM_ASSERT( transb == 'N' || transb == 'T' || transb == 'C' ||
                       transb == 'n' || transb == 't' || transb == 'c',
                        "unsupported transb flag '%c'", transb );

        int_t m = ( int_t ) ( tNoTransA ? A.n_rows() : A.n_cols() );
        int_t k = ( int_t ) ( tNoTransA ? A.n_cols() : A.n_rows() );
        int_t n = ( int_t ) ( tNoTransB ? B.n_cols() : B.n_rows() );

        BELFEM_ASSERT( ( int_t ) ( tNoTransB ? B.n_rows() : B.n_cols() ) == k,
            "Inner dimensions of op(A) and op(B) do not match ( %u vs %u )",
            ( unsigned int ) k,
            ( unsigned int ) ( tNoTransB ? B.n_rows() : B.n_cols() ) );

        if ( beta == static_cast< T >( 0.0 ) )
        {
            C.set_size( m, n, 0.0 );
        }

        BELFEM_ASSERT( ( int_t ) C.n_rows() == m,
            "Number of rows of op(A) and C do not match ( %u vs %u )",
            ( unsigned int ) m, ( unsigned int ) C.n_rows() );

        BELFEM_ASSERT( ( int_t ) C.n_cols() == n,
            "Number of cols of op(B) and C do not match ( %u vs %u )",
            ( unsigned int ) n, ( unsigned int ) C.n_cols() );

        // leading dimensions describe the stored arrays and are
        // independent of the transposition flags
        int_t lda = lapack::leading_dimension( A );
        int_t ldb = lapack::leading_dimension( B );
        int_t ldc = lapack::leading_dimension( C );

        lapack::gemm( &transa,
                 &transb,
                 &m,
                 &n,
                 &k,
                 &alpha,
                 A.data(),
                 &lda,
                 B.data(),
                 &ldb,
                 &beta,
                 C.data(),
                 &ldc );
    }

//------------------------------------------------------------------------------

} /* end namespace belfem */
//------------------------------------------------------------------------------
#endif //BELFEM_FN_GEMM_HPP
