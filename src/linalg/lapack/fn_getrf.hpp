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
 * @brief LU factorization with partial pivoting, A = P*L*U (LAPACK ?getrf).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_LAPACK_GETRF_HPP
#define BELFEM_FN_LAPACK_GETRF_HPP

#include "assert.hpp"
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
            sgetrf_( int_t    * m,
                     int_t    * n,
                     float    * a,
                     int_t    * lda,
                     int_t    * ipiv,
                     int_t    * info );

//------------------------------------------------------------------------------

            void
            dgetrf_( int_t    * m,
                     int_t    * n,
                     double   * a,
                     int_t    * lda,
                     int_t    * ipiv,
                     int_t    * info );

//------------------------------------------------------------------------------

            void
            cgetrf_( int_t        * m,
                     int_t        * n,
                     cplx_float_t * a,
                     int_t        * lda,
                     int_t        * ipiv,
                     int_t        * info );

//------------------------------------------------------------------------------

            void
            zgetrf_( int_t         * m,
                     int_t         * n,
                     cplx_double_t * a,
                     int_t         * lda,
                     int_t         * ipiv,
                     int_t         * info );

//------------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

//------------------------------------------------------------------------------

        template< typename T >
        void
        getrf( const int_t    * m,
               const int_t    * n,
                         T    * a,
               const int_t    * lda,
                     int_t    * ipiv,
                     int_t    * info )
        {
            static_assert( dependent_false< T >,
                "getrf not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getrf( const int_t    * m,
               const int_t    * n,
                     float    * a,
               const int_t    * lda,
                     int_t    * ipiv,
                     int_t    * info )
        {
            sgetrf_( const_cast< int_t * >( m ),
                     const_cast< int_t * >( n ),
                     a,
                     const_cast< int_t * >( lda ),
                     ipiv,
                     info  );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getrf( const int_t    * m,
               const int_t    * n,
                     double   * a,
               const int_t    * lda,
                     int_t    * ipiv,
                     int_t    * info )
        {
            dgetrf_( const_cast< int_t * >( m ),
                     const_cast< int_t * >( n ),
                     a,
                     const_cast< int_t * >( lda ),
                     ipiv,
                     info  );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getrf( const int_t                   * m,
               const int_t                   * n,
                     std::complex< float > * a,
               const int_t                   * lda,
                     int_t                   * ipiv,
                     int_t                   * info )
        {
            cgetrf_( const_cast< int_t * >( m ),
                     const_cast< int_t * >( n ),
                     reinterpret_cast< cplx_float_t * >( a ),
                     const_cast< int_t * >( lda ),
                     ipiv,
                     info  );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getrf( const int_t                    * m,
               const int_t                    * n,
                     std::complex< double > * a,
               const int_t                    * lda,
                     int_t                    * ipiv,
                     int_t                    * info )
        {
            zgetrf_( const_cast< int_t * >( m ),
                     const_cast< int_t * >( n ),
                     reinterpret_cast< cplx_double_t * >( a ),
                     const_cast< int_t * >( lda ),
                     ipiv,
                     info  );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */
//------------------------------------------------------------------------------

    /**
     * @brief LU factorization with partial pivoting via LAPACK ?getrf,
     *        A = P * L * U; use together with getri() to invert a matrix
     *
     * @param[in,out] A      matrix to factorize; overwritten with L and U
     * @param[in,out] Pivot  pivot indices; grown to min( m, n ) if too
     *                       short. Feed the result unchanged into getri()
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: U(info,info) is exactly zero, the
     *                factorization is complete but U is singular
     */
    template< typename T >
    int_t
    getrf( Matrix< T > & A , Vector< int_t > & Pivot, const bool AbortOnError = true )
    {
        int_t m    = A.n_rows() ;
        int_t n    = A.n_cols() ;
        int_t lda  = lapack::leading_dimension( A );
        int_t info = 0 ;

        if ( ( int_t ) Pivot.length() < std::min( m, n ) )
        {
            Pivot.set_size( std::min( m, n ), 0 );
        }

        lapack::getrf( &m, &n, A.data(), &lda, Pivot.data(), & info );

        BELFEM_ERROR( info == 0 || ! AbortOnError,
                  "LAPACK getrf has thrown an error: %i", ( int ) info );

        return info ;
    }

//------------------------------------------------------------------------------

} /* end namespace belfem */

#endif //BELFEM_FN_LAPACK_GETRF_HPP
