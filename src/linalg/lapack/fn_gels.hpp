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
 * @brief Least-squares or minimum-norm solution of a full-rank system (LAPACK ?gels).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_GELS_HPP
#define BELFEM_FN_GELS_HPP

#include <algorithm>
#include <complex>

#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "lapacktools.hpp"

//------------------------------------------------------------------------------
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
        sgels_( char   * trans,
                int_t  * m,
                int_t  * n,
                int_t  * nrhs,
                float  * a,
                int_t  * lda,
                float  * b,
                int_t  * ldb,
                float  * work,
                int_t  * lwork,
                int_t  * info,
                fortran_charlen_t lt );

//------------------------------------------------------------------------------

        void
        dgels_( char   * trans,
                int_t  * m,
                int_t  * n,
                int_t  * nrhs,
                double * a,
                int_t  * lda,
                double * b,
                int_t  * ldb,
                double * work,
                int_t  * lwork,
                int_t  * info,
                fortran_charlen_t lt );

//------------------------------------------------------------------------------

        void
        cgels_( char         * trans,
                int_t        * m,
                int_t        * n,
                int_t        * nrhs,
                cplx_float_t * a,
                int_t        * lda,
                cplx_float_t * b,
                int_t        * ldb,
                cplx_float_t * work,
                int_t        * lwork,
                int_t        * info,
                fortran_charlen_t lt );

//------------------------------------------------------------------------------

        void
        zgels_( char          * trans,
                int_t         * m,
                int_t         * n,
                int_t         * nrhs,
                cplx_double_t * a,
                int_t         * lda,
                cplx_double_t * b,
                int_t         * ldb,
                cplx_double_t * work,
                int_t         * lwork,
                int_t         * info,
                fortran_charlen_t lt );

//------------------------------------------------------------------------------

#ifdef __cplusplus
        }
#endif

//------------------------------------------------------------------------------

        template< typename T >
        void
        gels( const char     * trans,
              const int_t    * m,
              const int_t    * n,
              const int_t    * nrhs,
                    T        * a,
              const int_t    * lda,
                    T        * b,
              const int_t    * ldb,
                    T        * work,
              const int_t    * lwork,
                    int_t    * info )
        {
            static_assert( dependent_false< T >,
                "gels not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gels( const char     * trans,
              const int_t    * m,
              const int_t    * n,
              const int_t    * nrhs,
                    float    * a,
              const int_t    * lda,
                    float    * b,
              const int_t    * ldb,
                    float    * work,
              const int_t    * lwork,
                    int_t    * info )
        {
            sgels_(
                const_cast< char * >( trans ),
                const_cast< int_t * >( m ),
                const_cast< int_t * >( n ),
                const_cast< int_t * >( nrhs ),
                a,
                const_cast< int_t * >( lda ),
                b,
                const_cast< int_t * >( ldb ),
                work,
                const_cast< int_t * >( lwork ),
                info,
                1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gels( const char     * trans,
              const int_t    * m,
              const int_t    * n,
              const int_t    * nrhs,
                    double   * a,
              const int_t    * lda,
                    double   * b,
              const int_t    * ldb,
                    double   * work,
              const int_t    * lwork,
                    int_t    * info )
        {
            dgels_(
                const_cast< char * >( trans ),
                const_cast< int_t * >( m ),
                const_cast< int_t * >( n ),
                const_cast< int_t * >( nrhs ),
                a,
                const_cast< int_t * >( lda ),
                b,
                const_cast< int_t * >( ldb ),
                work,
                const_cast< int_t * >( lwork ),
                info,
                1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gels( const char            * trans,
              const int_t           * m,
              const int_t           * n,
              const int_t           * nrhs,
              std::complex< float > * a,
              const int_t           * lda,
              std::complex< float > * b,
              const int_t           * ldb,
              std::complex< float > * work,
              const int_t           * lwork,
                    int_t           * info )
        {
            cgels_(
                const_cast< char * >( trans ),
                const_cast< int_t * >( m ),
                const_cast< int_t * >( n ),
                const_cast< int_t * >( nrhs ),
                reinterpret_cast< cplx_float_t * > ( a ),
                const_cast< int_t * >( lda ),
                reinterpret_cast< cplx_float_t * > ( b ),
                const_cast< int_t * >( ldb ),
                reinterpret_cast< cplx_float_t * > ( work ),
                const_cast< int_t * >( lwork ),
                info,
                1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        gels( const char             * trans,
              const int_t            * m,
              const int_t            * n,
              const int_t            * nrhs,
              std::complex< double > * a,
              const int_t            * lda,
              std::complex< double > * b,
              const int_t            * ldb,
              std::complex< double > * work,
              const int_t            * lwork,
                    int_t            * info )
        {
            zgels_(
                const_cast< char * >( trans ),
                const_cast< int_t * >( m ),
                const_cast< int_t * >( n ),
                const_cast< int_t * >( nrhs ),
                reinterpret_cast< cplx_double_t * > ( a ),
                const_cast< int_t * >( lda ),
                reinterpret_cast< cplx_double_t * > ( b ),
                const_cast< int_t * >( ldb ),
                reinterpret_cast< cplx_double_t * > ( work ),
                const_cast< int_t * >( lwork ),
                info,
                1 );
        }

//------------------------------------------------------------------------------
    }
//------------------------------------------------------------------------------

    /**
     * @brief solve the least squares problem min || A * x - b || via
     *        LAPACK ?gels ( QR or LQ factorization, A must have full rank )
     *
     * @param[in,out] A     m x n system matrix; overwritten with the QR
     *                      ( m >= n ) or LQ ( m < n ) factorization
     * @param[in,out] B     right hand side on input ( first m entries ),
     *                      solution on output ( first n entries ); must be
     *                      allocated with at least max( m, n ) entries
     * @param[in,out] Work  scratch, grown to the optimal size on first use
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: a diagonal entry of the triangular
     *                factor is zero, A is rank deficient and no solution
     *                was computed
     */
    template< typename T >
    int_t
    gels( Matrix< T > & A, Vector< T > & B, Vector< T > & Work, const bool AbortOnError = true )
    {
        int_t m = A.n_rows();
        int_t n = A.n_cols();

        // B carries the right hand side ( m ) as well as the solution ( n )
        BELFEM_ASSERT( static_cast< int_t >( B.length() ) >= std::max( m, n ),
              "Length of rhs vector is %i, but must be at least %i.",
              ( int ) B.length(), ( int ) std::max( m, n ) );

        char trans = 'N';

        int_t nrhs = 1 ;
        int_t lda = lapack::leading_dimension( A );

        // leading dimension of B, taken from the actual allocation
        int_t ldb = lapack::leading_dimension( B );

        int_t mn = std::min( m, n );
        int_t lwork = std::max< int_t >( 1, mn + std::max( mn, nrhs ) );
        int_t info = 0 ;

        // check length of work array
        if ( static_cast< int_t >( Work.length() ) < lwork )
        {
            // ask lapack for the optimal size, this leaves A and B untouched
            Work.set_size( 1 );

            int_t query = -1 ;

            lapack::gels( &trans, &m, &n, &nrhs, A.data(), &lda, B.data(), &ldb, Work.data(), &query, &info );

            BELFEM_ERROR( info == 0 || ! AbortOnError,
                 "LAPACK gels workspace query has thrown an error: %i", ( int ) info );

            if ( info != 0 ) return info ;

            lwork = lapack::work_size( Work( 0 ) );

            Work.set_size( lwork );
        }
        else
        {
            // use the full buffer the caller has provided
            lwork = ( int_t ) Work.length();
        }

        lapack::gels(
            &trans,
            &m,
            &n,
            &nrhs,
            A.data(),
            &lda,
            B.data(),
            &ldb,
            Work.data(),
            &lwork, &info );

        // info > 0 : a diagonal entry of the triangular factor is zero, so A
        //            does not have full rank and no solution was computed
        BELFEM_ERROR( info == 0 || ! AbortOnError,
              "LAPACK gels has thrown an error: %i", ( int ) info );

        return info;
    }

//------------------------------------------------------------------------------

    /**
     * @brief solve min || A * X - B || for multiple right hand sides, see
     *        the single right hand side version above
     *
     * @param[in,out] A     m x n system matrix; overwritten with the QR
     *                      or LQ factorization
     * @param[in,out] B     right hand sides in columns; must have at
     *                      least max( m, n ) rows
     * @param[in,out] Work  scratch, grown to the optimal size on first use
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: A is rank deficient
     */
    template< typename T >
    int_t
    gels( Matrix< T > & A, Matrix< T > & B, Vector< T > & Work, const bool AbortOnError = true )
    {
        int_t m = A.n_rows();
        int_t n = A.n_cols();

        // B carries the right hand sides ( m ) as well as the solutions ( n )
        BELFEM_ASSERT( static_cast< int_t >( B.n_rows() ) >= std::max( m, n ),
             "Number of rows of rhs matrix is %i, but must be at least %i.",
             ( int ) B.n_rows(), ( int ) std::max( m, n ) );

        char trans = 'N';

        int_t nrhs = B.n_cols() ;

        int_t lda = lapack::leading_dimension( A );
        int_t ldb = lapack::leading_dimension( B );

        int_t mn = std::min( m, n );
        int_t lwork = std::max< int_t >( 1, mn + std::max( mn, nrhs ) );
        int_t info = 0 ;

        // check length of work array
        if ( static_cast< int_t >( Work.length() ) < lwork )
        {
            // ask lapack for the optimal size, this leaves A and B untouched
            Work.set_size( 1 );

            int_t query = -1 ;

            lapack::gels( &trans, &m, &n, &nrhs, A.data(), &lda, B.data(), &ldb, Work.data(), &query, &info );

            BELFEM_ERROR( info == 0 || ! AbortOnError,
                 "LAPACK gels workspace query has thrown an error: %i", ( int ) info );

            if ( info != 0 ) return info ;

            lwork = lapack::work_size( Work( 0 ) );

            Work.set_size( lwork );
        }
        else
        {
            // use the full buffer the caller has provided
            lwork = ( int_t ) Work.length();
        }

        lapack::gels( &trans, &m, &n, &nrhs, A.data(), &lda, B.data(), &ldb, Work.data(), &lwork, &info );

        // info > 0 : a diagonal entry of the triangular factor is zero, so A
        //            does not have full rank and no solution was computed
        BELFEM_ERROR( info == 0 || ! AbortOnError,
              "LAPACK gels has thrown an error: %i", ( int ) info );

        return info;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_GELS_HPP
