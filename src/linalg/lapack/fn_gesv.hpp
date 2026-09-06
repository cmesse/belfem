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
 * @brief Solves a square linear system by LU factorization (LAPACK ?gesv).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_GESV_HPP
#define BELFEM_FN_GESV_HPP

#include "lapacktools.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace lapack
    {
//------------------------------------------------------------------------------
#ifdef __cplusplus
        extern "C"
        {
#endif
//------------------------------------------------------------------------------

            void
            sgesv_( int_t   * n,
                    int_t   * nrhs,
                    float   * a,
                    int_t   * lda,
                    int_t   * ipiv,
                    float   * b,
                    int_t   * ldb,
                    int_t   * info );

//------------------------------------------------------------------------------

            void
            dgesv_( int_t   * n,
                    int_t   * nrhs,
                    double  * a,
                    int_t   * lda,
                    int_t   * ipiv,
                    double  * b,
                    int_t   * ldb,
                    int_t   * info );

//------------------------------------------------------------------------------

            void
            cgesv_( int_t          * n,
                    int_t          * nrhs,
                    cplx_float_t   * a,
                    int_t          * lda,
                    int_t          * ipiv,
                    cplx_float_t   * b,
                    int_t          * ldb,
                    int_t          * info );

//------------------------------------------------------------------------------

            void
            zgesv_( int_t          * n,
                    int_t          * nrhs,
                    cplx_double_t  * a,
                    int_t          * lda,
                    int_t          * ipiv,
                    cplx_double_t  * b,
                    int_t          * ldb,
                    int_t          * info );

//------------------------------------------------------------------------------
#ifdef __cplusplus
        }
#endif
//------------------------------------------------------------------------------

        template< typename T >
        void
        gesv( const int_t * n,
              const int_t * nrhs,
                    T     * a,
              const int_t * lda,
                    int_t * ipiv,
                    T     * b,
              const int_t * ldb,
                    int_t * info )
        {
            static_assert( dependent_false< T >,
                "gesv not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gesv( const int_t * n,
              const int_t * nrhs,
                    float * a,
              const int_t * lda,
                    int_t * ipiv,
                    float * b,
              const int_t * ldb,
                    int_t * info )
        {
            sgesv_(
                const_cast< int_t * > ( n ),
                const_cast< int_t * > ( nrhs ),
                a,
                const_cast< int_t * > ( lda ),
                ipiv,
                b,
                const_cast< int_t * > ( ldb ),
                info );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gesv( const int_t  * n,
              const int_t  * nrhs,
                    double * a,
              const int_t  * lda,
                    int_t  * ipiv,
                    double * b,
              const int_t  * ldb,
                    int_t  * info )
        {
            dgesv_(
                const_cast< int_t * > ( n ),
                const_cast< int_t * > ( nrhs ),
                a,
                const_cast< int_t * > ( lda ),
                ipiv,
                b,
                const_cast< int_t * > ( ldb ),
                info  );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gesv( const int_t                 * n,
              const int_t                 * nrhs,
                    std::complex< float > * a,
              const int_t                 * lda,
                    int_t                 * ipiv,
                    std::complex< float > * b,
              const int_t                 * ldb,
                    int_t                 * info )
        {
            cgesv_(
                const_cast< int_t * > ( n ),
                const_cast< int_t * > ( nrhs ),
                reinterpret_cast< cplx_float_t * >( a ),
                const_cast< int_t * > ( lda ),
                ipiv,
                reinterpret_cast< cplx_float_t * >( b ),
                const_cast< int_t * > ( ldb ),
                info );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gesv( const int_t                  * n,
              const int_t                  * nrhs,
                    std::complex< double > * a,
              const int_t                  * lda,
                    int_t                  * ipiv,
                    std::complex< double > * b,
              const int_t                  * ldb,
                    int_t                  * info )
        {
            zgesv_(
               const_cast< int_t * > ( n ),
               const_cast< int_t * > ( nrhs ),
               reinterpret_cast< cplx_double_t * >( a ),
               const_cast< int_t * > ( lda ),
               ipiv,
               reinterpret_cast< cplx_double_t * >( b ),
               const_cast< int_t * > ( ldb ) ,
               info );
        }

//------------------------------------------------------------------------------
    }
//------------------------------------------------------------------------------

    /**
     * @brief solve the square linear system A * x = b via LAPACK ?gesv
     *        ( LU factorization with partial pivoting )
     *
     * @param[in,out] A      square system matrix; overwritten with the
     *                       LU factors
     * @param[in,out] B      right hand side on input, solution on output
     * @param[out]    Pivot  pivot indices, at least n entries
     * @param[in]     AbortOnError  abort on info != 0 ( default ); pass
     *                       false to receive info instead, e.g. inside
     *                       iterative schemes that recover from a
     *                       singular matrix
     * @return info   0 on success; > 0: U(info,info) is exactly zero
     */
    template< typename T >
    int_t
    gesv( Matrix< T > & A, Vector< T > & B, Vector< int_t > & Pivot,
          const bool AbortOnError = true )
    {
        BELFEM_ASSERT( A.n_rows() == A.n_cols(),
            "Matrix A must be square ( is %lu x %lu )",
            ( long unsigned int ) A.n_rows(),
            ( long unsigned int ) A.n_cols() );

        BELFEM_ASSERT( B.length() == A.n_rows(),
            "Length of right hand side does not match ( %lu vs %lu )",
            ( long unsigned int ) B.length(),
            ( long unsigned int ) A.n_rows() );

        BELFEM_ASSERT( Pivot.length() >= A.n_rows(),
            "Pivot vector is too short ( %lu, need %lu )",
            ( long unsigned int ) Pivot.length(),
            ( long unsigned int ) A.n_rows() );

        // size of matrix
        int_t n = ( int_t ) A.n_rows();

        int_t nrhs = 1;

        int_t lda = lapack::leading_dimension( A );
        int_t ldb = lapack::leading_dimension( B );

        // error code
        int_t info = 0;

        // call lapack
        lapack::gesv(   &n,
                &nrhs,
                A.data(),
                &lda,
                Pivot.data(),
                B.data(),
                &ldb,
                &info );

        BELFEM_ERROR( info == 0 || ! AbortOnError,
                   "LAPACK gesv has thrown an error: %i", ( int ) info );

        return info ;
    }

//------------------------------------------------------------------------------

    /**
     * @brief solve A * X = B for multiple right hand sides via LAPACK
     *        ?gesv, see the vector version above
     *
     * @param[in,out] A      square system matrix; overwritten with the
     *                       LU factors
     * @param[in,out] B      right hand sides in columns on input,
     *                       solutions on output
     * @param[out]    Pivot  pivot indices, at least n entries
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: U(info,info) is exactly zero
     */
    template< typename T >
    int_t
    gesv( Matrix< T > & A, Matrix< T > & B, Vector< int_t > & Pivot, const bool AbortOnError = true )
    {
        BELFEM_ASSERT( A.n_rows() == A.n_cols(),
            "Matrix A must be square ( is %lu x %lu )",
            ( long unsigned int ) A.n_rows(),
            ( long unsigned int ) A.n_cols() );

        BELFEM_ASSERT( B.n_rows() == A.n_rows(),
            "Number of rows of right hand side does not match ( %lu vs %lu )",
            ( long unsigned int ) B.n_rows(),
            ( long unsigned int ) A.n_rows() );

        BELFEM_ASSERT( Pivot.length() >= A.n_rows(),
            "Pivot vector is too short ( %lu, need %lu )",
            ( long unsigned int ) Pivot.length(),
            ( long unsigned int ) A.n_rows() );

        // size of matrix
        int_t n = ( int_t ) A.n_rows();

        int_t nrhs = ( int_t ) B.n_cols();

        int_t lda = lapack::leading_dimension( A );
        int_t ldb = lapack::leading_dimension( B );

        // error code
        int_t info = 0;

        // call lapack
        lapack::gesv(   &n,
                & nrhs,
                A.data(),
                & lda,
                Pivot.data(),
                B.data(),
                & ldb,
                & info );

        BELFEM_ERROR( info == 0 || ! AbortOnError,
                   "LAPACK gesv has thrown an error: %i", ( int ) info );

        return info ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_GESV_HPP
