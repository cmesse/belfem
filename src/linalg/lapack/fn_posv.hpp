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
 * @brief Solves a real symmetric or complex Hermitian positive-definite system by
 *        Cholesky factorization (LAPACK ?posv).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_POSV_HPP
#define BELFEM_FN_POSV_HPP


#include "lapacktools.hpp"

namespace belfem
{
    namespace lapack
    {
//------------------------------------------------------------------------------
#ifdef __cplusplus
        extern "C"
        {
#endif
        void
        sposv_( char   * uplo,
                int_t  * n,
                int_t  * nrhs,
                float  * a,
                int_t  * lda,
                float  * b,
                int_t  * ldb,
                int_t  * info,
                fortran_charlen_t l );

//------------------------------------------------------------------------------

        void
        dposv_( char   * uplo,
                int_t  * n,
                int_t  * nrhs,
                double * a,
                int_t  * lda,
                double * b,
                int_t  * ldb,
                int_t  * info,
                fortran_charlen_t l );

//------------------------------------------------------------------------------

        // the complex flavors expect a hermitian positive definite matrix
        void
        cposv_( char         * uplo,
                int_t        * n,
                int_t        * nrhs,
                cplx_float_t * a,
                int_t        * lda,
                cplx_float_t * b,
                int_t        * ldb,
                int_t        * info,
                fortran_charlen_t l );

//------------------------------------------------------------------------------

        void
        zposv_( char          * uplo,
                int_t         * n,
                int_t         * nrhs,
                cplx_double_t * a,
                int_t         * lda,
                cplx_double_t * b,
                int_t         * ldb,
                int_t         * info,
                fortran_charlen_t l );

//------------------------------------------------------------------------------

#ifdef __cplusplus
        }
#endif

        template< typename T >
        void
        posv( const char   * uplo,
              const int_t  * n,
              const int_t  * nrhs,
                    T      * a,
              const int_t  * lda,
                    T      * b,
              const int_t  * ldb,
                    int_t  * info )
        {
            static_assert( dependent_false< T >,
                "posv not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        posv( const char   * uplo,
              const int_t  * n,
              const int_t  * nrhs,
                    float  * a,
              const int_t  * lda,
                    float  * b,
              const int_t  * ldb,
                    int_t  * info )
        {
            sposv_(
                const_cast< char * > ( uplo ),
                const_cast< int_t * >( n ),
                const_cast< int_t * >( nrhs ),
                a,
                const_cast< int_t * >( lda ),
                b,
                const_cast< int_t * >( ldb ),
                info,
                1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        posv( const char   * uplo,
              const int_t  * n,
              const int_t  * nrhs,
                    double * a,
              const int_t  * lda,
                    double * b,
              const int_t  * ldb,
                    int_t  * info )
        {
            dposv_(
                const_cast< char * > ( uplo ),
                const_cast< int_t * >( n ),
                const_cast< int_t * >( nrhs ),
                a,
                const_cast< int_t * >( lda ),
                b,
                const_cast< int_t * >( ldb ),
                info,
                1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        posv( const char             * uplo,
              const int_t            * n,
              const int_t            * nrhs,
              std::complex< float >  * a,
              const int_t            * lda,
              std::complex< float >  * b,
              const int_t            * ldb,
                    int_t            * info )
        {
            cposv_(
                const_cast< char * > ( uplo ),
                const_cast< int_t * >( n ),
                const_cast< int_t * >( nrhs ),
                reinterpret_cast< cplx_float_t * > ( a ),
                const_cast< int_t * >( lda ),
                reinterpret_cast< cplx_float_t * > ( b ),
                const_cast< int_t * >( ldb ),
                info,
                1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        posv( const char             * uplo,
              const int_t            * n,
              const int_t            * nrhs,
              std::complex< double > * a,
              const int_t            * lda,
              std::complex< double > * b,
              const int_t            * ldb,
                    int_t            * info )
        {
            zposv_(
                 const_cast< char * > ( uplo ),
                 const_cast< int_t * >( n ),
                 const_cast< int_t * >( nrhs ),
                 reinterpret_cast< cplx_double_t * > ( a ),
                 const_cast< int_t * >( lda ),
                 reinterpret_cast< cplx_double_t * > ( b ),
                 const_cast< int_t * >( ldb ),
                 info,
                 1 );
        }
    }

//------------------------------------------------------------------------------

    /**
     * @brief solve A * x = b for a symmetric ( real ) or Hermitian
     *        ( complex ) positive definite matrix via LAPACK ?posv
     *        ( Cholesky factorization, no pivoting )
     *
     * @param[in,out] A      positive definite matrix; overwritten with
     *                       the Cholesky factor
     * @param[in,out] B      right hand side on input, solution on output
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: the leading minor of order info
     *                is not positive definite
     */
    template< typename T >
    int_t
    posv( Matrix< T > & A, Vector< T > & B, const bool AbortOnError = true )
    {
        BELFEM_ASSERT( A.n_rows() == B.length(),
                      "Number of rows of matrix does not match." );
        BELFEM_ASSERT( A.n_cols() == B.length(),
                      "Number of cols of matrix does not match." );

        // size of matrix
        int_t n = ( int_t ) A.n_rows();

        // leading dimensions
        int_t lda = lapack::leading_dimension( A );
        int_t ldb = lapack::leading_dimension( B );

        char uplo = 'L';

        int_t nrhs = 1;

        // error code
        int_t info = 0;

        // call lapack
        lapack::posv(
                &uplo,
                &n,
                &nrhs,
                A.data(),
                &lda,
                B.data(),
                &ldb,
                &info );

        BELFEM_ERROR( info == 0 || ! AbortOnError,
                   "LAPACK posv has thrown an error: %i", ( int ) info );

        return info ;
    }

//------------------------------------------------------------------------------

    /**
     * @brief solve A * X = B for multiple right hand sides of a symmetric
     *        ( real ) or Hermitian ( complex ) positive definite matrix,
     *        see the vector version above
     *
     * @param[in,out] A      positive definite matrix; overwritten with
     *                       the Cholesky factor
     * @param[in,out] B      right hand sides in columns on input,
     *                       solutions on output
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: the leading minor of order info
     *                is not positive definite
     */
    template< typename T >
    int_t
    posv( Matrix< T > & A, Matrix< T > & B, const bool AbortOnError = true )
    {
        BELFEM_ASSERT( A.n_rows() == B.n_rows(),
                      "Number of rows of matrix does not match." );
        BELFEM_ASSERT( A.n_cols() == B.n_rows(),
                      "Number of cols of matrix does not match." );
        // size of matrix
        int_t n    = ( int_t ) B.n_rows();

        int_t nrhs = ( int_t ) B.n_cols();

        // leading dimensions
        int_t lda = lapack::leading_dimension( A );
        int_t ldb = lapack::leading_dimension( B );

        // error code
        int_t info = 0;

        char uplo = 'L';

        // call lapack
        lapack::posv(
                &uplo,
                &n,
                &nrhs,
                A.data(),
                &lda,
                B.data(),
                &ldb,
                &info );

        BELFEM_ERROR( info == 0 || ! AbortOnError,
                   "LAPACK posv has thrown an error: %i", ( int ) info );

        return info;
    }
    
//------------------------------------------------------------------------------
}


#endif //BELFEM_FN_POSV_HPP
