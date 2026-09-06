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
 * @brief Singular value decomposition A = U * diag(S) * VT (LAPACK ?gesvd).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_GESVD_HPP
#define BELFEM_FN_GESVD_HPP

#include "assert.hpp"
#include "lapacktools.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace lapack
    {
#ifdef __cplusplus
        extern "C"
        {
#endif
// -----------------------------------------------------------------------------

        // two gesvd particularities: the singular values s ( and the
        // complex-only scratch rwork ) stay REAL valued in all flavors,
        // and the two job chars carry two hidden Fortran length arguments

        void
        sgesvd_(
            char          * jobu,
            char          * jobvt,
            int_t         * m,
            int_t         * n,
            float         * a,
            int_t         * lda,
            float         * s,
            float         * u,
            int_t         * ldu,
            float         * vt,
            int_t         * ldvt,
            float         * work,
            int_t         * lwork,
            int_t         * info,
            fortran_charlen_t lu,
            fortran_charlen_t lvt );

// -----------------------------------------------------------------------------

        void
        dgesvd_(
            char          * jobu,
            char          * jobvt,
            int_t         * m,
            int_t         * n,
            double        * a,
            int_t         * lda,
            double        * s,
            double        * u,
            int_t         * ldu,
            double        * vt,
            int_t         * ldvt,
            double        * work,
            int_t         * lwork,
            int_t         * info,
            fortran_charlen_t lu,
            fortran_charlen_t lvt );

// -----------------------------------------------------------------------------

        void
        cgesvd_(
            char          * jobu,
            char          * jobvt,
            int_t         * m,
            int_t         * n,
            cplx_float_t  * a,
            int_t         * lda,
            float         * s,
            cplx_float_t  * u,
            int_t         * ldu,
            cplx_float_t  * vt,
            int_t         * ldvt,
            cplx_float_t  * work,
            int_t         * lwork,
            float         * rwork,
            int_t         * info,
            fortran_charlen_t lu,
            fortran_charlen_t lvt );

// -----------------------------------------------------------------------------

        void
        zgesvd_(
            char          * jobu,
            char          * jobvt,
            int_t         * m,
            int_t         * n,
            cplx_double_t * a,
            int_t         * lda,
            double        * s,
            cplx_double_t * u,
            int_t         * ldu,
            cplx_double_t * vt,
            int_t         * ldvt,
            cplx_double_t * work,
            int_t         * lwork,
            double        * rwork,
            int_t         * info,
            fortran_charlen_t lu,
            fortran_charlen_t lvt );

#ifdef __cplusplus
        }
#endif
// -----------------------------------------------------------------------------

        // unified dispatch: rwork is referenced by the complex flavors only,
        // the real specializations ignore it

        template< typename T >
        void gesvd(
            const char        * jobu,
            const char        * jobvt,
            const int_t       * m,
            const int_t       * n,
                  T           * a,
            const int_t       * lda,
                  real_t< T > * s,
                  T           * u,
            const int_t       * ldu,
                  T           * vt,
            const int_t       * ldvt,
                  T           * work,
            const int_t       * lwork,
                  real_t< T > * rwork,
                  int_t       * info )
        {
            static_assert( dependent_false< T >,
                "gesvd not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void gesvd(
            const char  * jobu,
            const char  * jobvt,
            const int_t * m,
            const int_t * n,
                  float * a,
            const int_t * lda,
                  float * s,
                  float * u,
            const int_t * ldu,
                  float * vt,
            const int_t * ldvt,
                  float * work,
            const int_t * lwork,
                  float * /* rwork */,
                  int_t * info )
        {
            sgesvd_(
                const_cast< char  * >( jobu ),
                const_cast< char  * >( jobvt ),
                const_cast< int_t * >( m ),
                const_cast< int_t * >( n ),
                a,
                const_cast< int_t * >( lda ),
                s,
                u,
                const_cast< int_t * >( ldu ),
                vt,
                const_cast< int_t * >( ldvt ),
                work,
                const_cast< int_t * >( lwork ),
                info,
                1, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void gesvd(
            const char   * jobu,
            const char   * jobvt,
            const int_t  * m,
            const int_t  * n,
                  double * a,
            const int_t  * lda,
                  double * s,
                  double * u,
            const int_t  * ldu,
                  double * vt,
            const int_t  * ldvt,
                  double * work,
            const int_t  * lwork,
                  double * /* rwork */,
                  int_t  * info )
        {
            dgesvd_(
                const_cast< char  * >( jobu ),
                const_cast< char  * >( jobvt ),
                const_cast< int_t * >( m ),
                const_cast< int_t * >( n ),
                a,
                const_cast< int_t * >( lda ),
                s,
                u,
                const_cast< int_t * >( ldu ),
                vt,
                const_cast< int_t * >( ldvt ),
                work,
                const_cast< int_t * >( lwork ),
                info,
                1, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void gesvd(
            const char                  * jobu,
            const char                  * jobvt,
            const int_t                 * m,
            const int_t                 * n,
                  std::complex< float > * a,
            const int_t                 * lda,
                  float                 * s,
                  std::complex< float > * u,
            const int_t                 * ldu,
                  std::complex< float > * vt,
            const int_t                 * ldvt,
                  std::complex< float > * work,
            const int_t                 * lwork,
                  float                 * rwork,
                  int_t                 * info )
        {
            cgesvd_(
                const_cast< char  * >( jobu ),
                const_cast< char  * >( jobvt ),
                const_cast< int_t * >( m ),
                const_cast< int_t * >( n ),
                reinterpret_cast< cplx_float_t * >( a ),
                const_cast< int_t * >( lda ),
                s,
                reinterpret_cast< cplx_float_t * >( u ),
                const_cast< int_t * >( ldu ),
                reinterpret_cast< cplx_float_t * >( vt ),
                const_cast< int_t * >( ldvt ),
                reinterpret_cast< cplx_float_t * >( work ),
                const_cast< int_t * >( lwork ),
                rwork,
                info,
                1, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void gesvd(
            const char                   * jobu,
            const char                   * jobvt,
            const int_t                  * m,
            const int_t                  * n,
                  std::complex< double > * a,
            const int_t                  * lda,
                  double                 * s,
                  std::complex< double > * u,
            const int_t                  * ldu,
                  std::complex< double > * vt,
            const int_t                  * ldvt,
                  std::complex< double > * work,
            const int_t                  * lwork,
                  double                 * rwork,
                  int_t                  * info )
        {
            zgesvd_(
                const_cast< char  * >( jobu ),
                const_cast< char  * >( jobvt ),
                const_cast< int_t * >( m ),
                const_cast< int_t * >( n ),
                reinterpret_cast< cplx_double_t * >( a ),
                const_cast< int_t * >( lda ),
                s,
                reinterpret_cast< cplx_double_t * >( u ),
                const_cast< int_t * >( ldu ),
                reinterpret_cast< cplx_double_t * >( vt ),
                const_cast< int_t * >( ldvt ),
                reinterpret_cast< cplx_double_t * >( work ),
                const_cast< int_t * >( lwork ),
                rwork,
                info,
                1, 1 );
        }

// -----------------------------------------------------------------------------
    } /* end namespace lapack */

    /**
     * @brief singular value decomposition A = U * diag( S ) * VT via
     *        LAPACK ?gesvd
     *
     * @param[in,out] A     m x n input matrix; destroyed by the call
     * @param[out]    S     the min( m, n ) singular values, descending —
     *                      always real valued, also for complex A
     * @param[out]    U     left singular vectors; resized to m x m for
     *                      jobu = 'A', m x min( m, n ) for 'S', untouched
     *                      for 'N'
     * @param[out]    VT    transposed right singular vectors; resized to
     *                      n x n for jobvt = 'A', min( m, n ) x n for 'S',
     *                      untouched for 'N'
     * @param[in,out] Work  single real-valued scratch, grown to the
     *                      optimal size on first use: LAPACK's work array
     *                      in its head ( for complex T reinterpreted in
     *                      place, two reals per entry ), the complex-only
     *                      5*min( m, n ) real scratch rwork in its tail —
     *                      no internal allocation
     * @param[in]     jobu  'A', 'S' or 'N' ( 'O' is not supported )
     * @param[in]     jobvt 'A', 'S' or 'N' ( 'O' is not supported )
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: the QR iteration did not converge
     */
    template< typename T >
    int_t
    gesvd(
        Matrix< T > & A,
        Vector< lapack::real_t< T > > & S,
        Matrix< T > & U,
        Matrix< T > & VT,
        Vector< lapack::real_t< T > > & Work,
        const char jobu  = 'A',
        const char jobvt = 'A',
        const bool AbortOnError = true )
    {
        BELFEM_ASSERT( jobu == 'A' || jobu == 'S' || jobu == 'N',
            "unsupported jobu flag '%c'", jobu );
        BELFEM_ASSERT( jobvt == 'A' || jobvt == 'S' || jobvt == 'N',
            "unsupported jobvt flag '%c'", jobvt );

        // reals per LAPACK work entry
        constexpr int_t tRealsPerT =
            std::is_same< T, lapack::real_t< T > >::value ? 1 : 2 ;

        int_t m  = ( int_t ) A.n_rows();
        int_t n  = ( int_t ) A.n_cols();
        int_t mn = std::min( m, n );

        int_t lda = lapack::leading_dimension( A );

        S.set_size( mn );

        // U is m x m for 'A' and m x mn for 'S'; lapack requires
        // ldu >= 1 even when U is not referenced ( jobu = 'N' )
        int_t ldu = 1 ;
        if ( jobu != 'N' )
        {
            U.set_size( m, jobu == 'A' ? m : mn );
            ldu = lapack::leading_dimension( U );
        }

        // VT is n x n for 'A' and mn x n for 'S'
        int_t ldvt = 1 ;
        if ( jobvt != 'N' )
        {
            VT.set_size( jobvt == 'A' ? n : mn, n );
            ldvt = lapack::leading_dimension( VT );
        }

        // real scratch in the tail, referenced by the complex flavors only
        int_t tRWorkSize = tRealsPerT == 1 ? 0 : 5 * mn ;

        // minimum work sizes: real flavors need
        // max( 3*mn + max( m, n ), 5*mn ), complex ones 2*mn + max( m, n )
        int_t lwork = std::max< int_t >( 1, tRealsPerT == 1 ?
            std::max( 3 * mn + std::max( m, n ), 5 * mn ) :
            2 * mn + std::max( m, n ) );

        int_t info = 0 ;

        // required buffer: lwork entries of T in the head plus the rwork tail
        if ( static_cast< int_t >( Work.length() ) < tRealsPerT * lwork + tRWorkSize )
        {
            // ask lapack for the optimal size ( one T entry plus the
            // tail ); this leaves A untouched. Grow-only: never shrink a
            // buffer that already covers the query call
            if ( static_cast< int_t >( Work.length() ) < tRealsPerT + tRWorkSize )
            {
                Work.set_size( tRealsPerT + tRWorkSize );
            }

            int_t query = -1 ;

            lapack::gesvd( &jobu, &jobvt, &m, &n, A.data(), &lda, S.data(),
                U.data(), &ldu, VT.data(), &ldvt,
                reinterpret_cast< T * >( Work.data() ), &query,
                Work.data() + tRealsPerT, &info );

            BELFEM_ERROR( info == 0 || ! AbortOnError,
                "LAPACK gesvd workspace query has thrown an error: %i", ( int ) info );

            if ( info != 0 ) return info ;

            // the optimal size sits in the real part of the first entry.
            // Keep the reference minimum as floor: a vendor query may
            // undercut it ( MKL: 7 vs 15 for a 4x3 'A','A' ), and a buffer
            // sized below the length test above would re-enter this branch
            // on every call instead of being reused
            lwork = std::max( lwork, lapack::work_size( Work( 0 ) ) );

            Work.set_size( tRealsPerT * lwork + tRWorkSize );
        }
        else
        {
            // use the full buffer the caller has provided
            lwork = ( ( int_t ) Work.length() - tRWorkSize ) / tRealsPerT ;
        }

        // work segment in the head, rwork in the tail
        lapack::gesvd( &jobu, &jobvt, &m, &n, A.data(), &lda, S.data(),
            U.data(), &ldu, VT.data(), &ldvt,
            reinterpret_cast< T * >( Work.data() ), &lwork,
            Work.data() + tRealsPerT * lwork, &info );

        // info > 0 : the QR iteration did not converge
        BELFEM_ERROR( info == 0 || ! AbortOnError,
            "LAPACK gesvd has thrown an error: %i", ( int ) info );

        return info;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_GESVD_HPP
