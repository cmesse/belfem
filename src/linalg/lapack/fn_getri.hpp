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
 * @brief Inverts a square matrix in place from its LU factors (LAPACK ?getri).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_GETRI_HPP
#define BELFEM_FN_GETRI_HPP

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
            sgetri_( int_t  * n,
                     float  * a,
                     int_t  * lda,
                     int_t  * ipiv,
                     float  * work,
                     int_t  * lwork,
                     int_t  * info );

//------------------------------------------------------------------------------

            void
            dgetri_( int_t  * n,
                 double  * a,
                 int_t  * lda,
                 int_t  * ipiv,
                 double  * work,
                 int_t  * lwork,
                 int_t  * info );
//------------------------------------------------------------------------------

            void
            cgetri_( int_t         * n,
                     cplx_float_t  * a,
                     int_t         * lda,
                     int_t         * ipiv,
                     cplx_float_t  * work,
                     int_t         * lwork,
                     int_t         * info );

//------------------------------------------------------------------------------

            void
            zgetri_( int_t     * n,
                 cplx_double_t * a,
                 int_t         * lda,
                 int_t         * ipiv,
                 cplx_double_t * work,
                 int_t         * lwork,
                 int_t         * info );


//------------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

//------------------------------------------------------------------------------

        template< typename T >
        void
        getri(  const int_t   * n,
                      T       * a,
                const int_t   * lda,
                      int_t   * ipiv,
                      T       * work,
                const int_t   * lwork,
                      int_t   * info )
        {
            static_assert( dependent_false< T >,
                "getri not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getri(  const int_t   * n,
                      float   * a,
                const int_t   * lda,
                      int_t   * ipiv,
                      float   * work,
                const int_t   * lwork,
                      int_t   * info )
        {
            sgetri_( const_cast< int_t * >( n ),
                     a,
                     const_cast< int_t * >( lda ),
                     ipiv,
                     work,
                     const_cast< int_t * >( lwork ),
                     info );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getri( const int_t    * n,
                     double  * a,
               const int_t   * lda,
                     int_t   * ipiv,
                     double  * work,
               const int_t   * lwork,
                     int_t   * info )
        {
            dgetri_( const_cast< int_t * >( n ),
                     a,
                     const_cast< int_t * >( lda ),
                     ipiv,
                     work,
                     const_cast< int_t * >( lwork ),
                     info );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getri( const int_t   * n,
               std::complex< float > * a,
               const int_t   * lda,
                     int_t   * ipiv,
               std::complex< float > * work,
               const int_t   * lwork,
                     int_t   * info )
        {
            cgetri_( const_cast< int_t * >( n ),
                     reinterpret_cast< cplx_float_t * >( a ),
                     const_cast< int_t * >( lda ),
                     ipiv,
                     reinterpret_cast< cplx_float_t * >( work ),
                     const_cast< int_t * >( lwork ),
                     info );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void
        getri( const int_t   * n,
               std::complex< double > * a,
               const int_t   * lda,
                     int_t   * ipiv,
               std::complex< double > * work,
               const int_t   * lwork,
                     int_t   * info )
        {
            zgetri_( const_cast< int_t * >( n ),
                     reinterpret_cast< cplx_double_t * >( a ),
                     const_cast< int_t * >( lda ),
                     ipiv,
                     reinterpret_cast< cplx_double_t * >( work ),
                     const_cast< int_t * >( lwork ),
                     info );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */
//------------------------------------------------------------------------------

    /**
     * @brief invert a square matrix in place via LAPACK ?getri, using the
     *        LU factorization computed by getrf()
     *
     * @param[in,out] A      must hold the getrf() factorization on entry;
     *                       overwritten with the inverse
     * @param[in]     Pivot  pivot indices exactly as produced by getrf()
     * @param[in,out] Work   scratch, grown to the optimal size on first
     *                       use; reused unchanged when already large enough
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: U(info,info) is exactly zero, the
     *                matrix is singular and no inverse was computed
     */
    template< typename T >
    int_t
    getri( Matrix< T > & A , Vector< int_t > & Pivot, Vector< T > & Work, const bool AbortOnError = true )
    {
        BELFEM_ASSERT( A.n_rows() == A.n_cols(),
            "Matrix A must be square ( is %lu x %lu )",
            ( long unsigned int ) A.n_rows(),
            ( long unsigned int ) A.n_cols() );

        BELFEM_ASSERT( Pivot.length() >= A.n_rows(),
            "Pivot vector is too short ( %lu, need %lu )",
            ( long unsigned int ) Pivot.length(),
            ( long unsigned int ) A.n_rows() );

        int_t n    = ( int_t ) A.n_cols() ;
        int_t lda  = lapack::leading_dimension( A );
        int_t info = 0 ;

        int_t lwork = std::max< int_t >( 1, n );

        if ( static_cast< int_t >( Work.length() ) < lwork )
        {
            // ask lapack for the optimal size, this leaves A untouched
            Work.set_size( 1 );

            int_t query = -1 ;

            lapack::getri( &n, A.data(), &lda, Pivot.data(), Work.data(), &query, &info );

            BELFEM_ERROR( info == 0 || ! AbortOnError,
                 "LAPACK getri workspace query has thrown an error: %i", ( int ) info );

            if ( info != 0 ) return info ;

            lwork = lapack::work_size( Work( 0 ) );

            Work.set_size( lwork );
        }
        else
        {
            // use the full buffer the caller has provided
            lwork = ( int_t ) Work.length();
        }

        lapack::getri( &n, A.data(), &lda, Pivot.data(), Work.data(), &lwork, &info );

        BELFEM_ERROR( info == 0 || ! AbortOnError,
                  "LAPACK getri has thrown an error: %i", ( int ) info );

        return info ;
    }
} /* end namespace belfem */
#endif //BELFEM_FN_GETRI_HPP
