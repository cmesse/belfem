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
 * @brief Eigenvalues and eigenvectors of a general square matrix (LAPACK ?geev).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_GEEV_HPP
#define BELFEM_FN_GEEV_HPP

#include "assert.hpp"
#include "lapacktools.hpp"
#include "cl_Vector.hpp"

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

        // geev particularities: the real flavors return the eigenvalues as
        // separate real and imaginary arrays ( wr, wi ) and have no rwork;
        // the complex flavors take one complex w plus a REAL rwork of size
        // 2*n. The two job chars carry two hidden Fortran length arguments.

        void
        sgeev_( char     * jobvl,
                char     * jobvr,
                int_t    * n,
                float    * a,
                int_t    * lda,
                float    * wr,
                float    * wi,
                float    * vl,
                int_t    * ldvl,
                float    * vr,
                int_t    * ldvr,
                float    * work,
                int_t    * lwork,
                int_t    * info,
                fortran_charlen_t lvl,
                fortran_charlen_t lvr );

//------------------------------------------------------------------------------

        void
        dgeev_( char     * jobvl,
                char     * jobvr,
                int_t    * n,
                double   * a,
                int_t    * lda,
                double   * wr,
                double   * wi,
                double   * vl,
                int_t    * ldvl,
                double   * vr,
                int_t    * ldvr,
                double   * work,
                int_t    * lwork,
                int_t    * info,
                fortran_charlen_t lvl,
                fortran_charlen_t lvr );

//------------------------------------------------------------------------------

        void
        cgeev_( char         * jobvl,
                char         * jobvr,
                int_t        * n,
                cplx_float_t * a,
                int_t        * lda,
                cplx_float_t * w,
                cplx_float_t * vl,
                int_t        * ldvl,
                cplx_float_t * vr,
                int_t        * ldvr,
                cplx_float_t * work,
                int_t        * lwork,
                float        * rwork,
                int_t        * info,
                fortran_charlen_t lvl,
                fortran_charlen_t lvr );

//------------------------------------------------------------------------------

        void
        zgeev_( char          * jobvl,
                char          * jobvr,
                int_t         * n,
                cplx_double_t * a,
                int_t         * lda,
                cplx_double_t * w,
                cplx_double_t * vl,
                int_t         * ldvl,
                cplx_double_t * vr,
                int_t         * ldvr,
                cplx_double_t * work,
                int_t         * lwork,
                double        * rwork,
                int_t         * info,
                fortran_charlen_t lvl,
                fortran_charlen_t lvr );

//------------------------------------------------------------------------------
#ifdef __cplusplus
        }
#endif
//------------------------------------------------------------------------------

        // unified dispatch. The real scratch rwork must hold 2*n entries in
        // ALL flavors: the real specializations use it as wr/wi and pack the
        // eigenvalues into the complex w afterwards, the complex ones pass
        // it through as LAPACK's actual rwork.

        template< typename T >
        void geev(
            const char        * jobvl,
            const char        * jobvr,
            const int_t       * n,
                  T           * a,
            const int_t       * lda,
                  cplx_t< T > * w,
                  T           * vl,
            const int_t       * ldvl,
                  T           * vr,
            const int_t       * ldvr,
                  T           * work,
            const int_t       * lwork,
                  real_t< T > * rwork,
                  int_t       * info )
        {
            static_assert( dependent_false< T >,
                "geev not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void geev(
           const char                   * jobvl,
           const char                   * jobvr,
           const int_t                  * n,
                 float                  * a,
           const int_t                  * lda,
                 std::complex< float >  * w,
                 float                  * vl,
           const int_t                  * ldvl,
                 float                  * vr,
           const int_t                  * ldvr,
                 float                  * work,
           const int_t                  * lwork,
                 float                  * rwork,
                 int_t                  * info )
        {
            sgeev_(
                const_cast< char * >( jobvl ),
                const_cast< char * >( jobvr ),
                const_cast< int_t * >( n ),
                a,
                const_cast< int_t * >( lda ),
                rwork,
                rwork + *n,
                vl,
                const_cast< int_t * >( ldvl ),
                vr,
                const_cast< int_t * >( ldvr ),
                work,
                const_cast< int_t * >( lwork ),
                info,
                1, 1 );

            // pack wr/wi into the complex eigenvalue vector; a workspace
            // query does not touch wr/wi, so there is nothing to pack
            if ( *lwork != -1 )
            {
                for ( int_t k = 0; k < *n; ++k )
                {
                    w[ k ] = std::complex< float >( rwork[ k ], rwork[ k + *n ] );
                }
            }
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void geev(
           const char                   * jobvl,
           const char                   * jobvr,
           const int_t                  * n,
                 double                 * a,
           const int_t                  * lda,
                 std::complex< double > * w,
                 double                 * vl,
           const int_t                  * ldvl,
                 double                 * vr,
           const int_t                  * ldvr,
                 double                 * work,
           const int_t                  * lwork,
                 double                 * rwork,
                 int_t                  * info )
        {
            dgeev_(
                const_cast< char * >( jobvl ),
                const_cast< char * >( jobvr ),
                const_cast< int_t * >( n ),
                a,
                const_cast< int_t * >( lda ),
                rwork,
                rwork + *n,
                vl,
                const_cast< int_t * >( ldvl ),
                vr,
                const_cast< int_t * >( ldvr ),
                work,
                const_cast< int_t * >( lwork ),
                info,
                1, 1 );

            if ( *lwork != -1 )
            {
                for ( int_t k = 0; k < *n; ++k )
                {
                    w[ k ] = std::complex< double >( rwork[ k ], rwork[ k + *n ] );
                }
            }
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void geev(
           const char                   * jobvl,
           const char                   * jobvr,
           const int_t                  * n,
                 std::complex< float >  * a,
           const int_t                  * lda,
                 std::complex< float >  * w,
                 std::complex< float >  * vl,
           const int_t                  * ldvl,
                 std::complex< float >  * vr,
           const int_t                  * ldvr,
                 std::complex< float >  * work,
           const int_t                  * lwork,
                 float                  * rwork,
                 int_t                  * info )
        {
            cgeev_(
                const_cast< char * >( jobvl ),
                const_cast< char * >( jobvr ),
                const_cast< int_t * >( n ),
                reinterpret_cast< cplx_float_t * >( a ),
                const_cast< int_t * >( lda ),
                reinterpret_cast< cplx_float_t * >( w ),
                reinterpret_cast< cplx_float_t * >( vl ),
                const_cast< int_t * >( ldvl ),
                reinterpret_cast< cplx_float_t * >( vr ),
                const_cast< int_t * >( ldvr ),
                reinterpret_cast< cplx_float_t * >( work ),
                const_cast< int_t * >( lwork ),
                rwork,
                info,
                1, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void geev(
           const char                   * jobvl,
           const char                   * jobvr,
           const int_t                  * n,
                 std::complex< double > * a,
           const int_t                  * lda,
                 std::complex< double > * w,
                 std::complex< double > * vl,
           const int_t                  * ldvl,
                 std::complex< double > * vr,
           const int_t                  * ldvr,
                 std::complex< double > * work,
           const int_t                  * lwork,
                 double                 * rwork,
                 int_t                  * info )
        {
            zgeev_(
                const_cast< char * >( jobvl ),
                const_cast< char * >( jobvr ),
                const_cast< int_t * >( n ),
                reinterpret_cast< cplx_double_t * >( a ),
                const_cast< int_t * >( lda ),
                reinterpret_cast< cplx_double_t * >( w ),
                reinterpret_cast< cplx_double_t * >( vl ),
                const_cast< int_t * >( ldvl ),
                reinterpret_cast< cplx_double_t * >( vr ),
                const_cast< int_t * >( ldvr ),
                reinterpret_cast< cplx_double_t * >( work ),
                const_cast< int_t * >( lwork ),
                rwork,
                info,
                1, 1 );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */

    /**
     * @brief eigenvalues and eigenvectors of a general square matrix,
     *        A * v = lambda * v, via LAPACK ?geev
     *
     * @param[in,out] A     square input matrix; destroyed by the call
     * @param[out]    W     the n eigenvalues — always complex valued
     * @param[out]    VL    left eigenvectors, resized to n x n for
     *                      jobvl = 'V'; untouched for 'N'
     * @param[out]    VR    right eigenvectors, likewise for jobvr.
     *                      For real T, a complex conjugate eigenvalue pair
     *                      ( W(j), W(j+1) ) stores its eigenvectors
     *                      LAPACK-packed across two consecutive columns:
     *                      v_j = VR(:,j) + i * VR(:,j+1) and
     *                      v_{j+1} = conj( v_j )
     * @param[in,out] Work  single real-valued scratch, grown to the
     *                      optimal size on first use: LAPACK's work array
     *                      in its head ( for complex T reinterpreted in
     *                      place, two reals per entry ), the 2*n real
     *                      scratch ( wr/wi for real T, rwork for complex
     *                      T ) in its tail — no internal allocation
     * @param[in]     jobvl 'V' or 'N'
     * @param[in]     jobvr 'V' or 'N'
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: the QR iteration failed —
     *                eigenvalues info+1 .. n have converged, no
     *                eigenvectors were computed
     */
    template< typename T >
    int_t
    geev(
        Matrix< T > & A,
        Vector< lapack::cplx_t< T > > & W,
        Matrix< T > & VL,
        Matrix< T > & VR,
        Vector< lapack::real_t< T > > & Work,
        const char jobvl = 'V',
        const char jobvr = 'V',
        const bool AbortOnError = true )
    {
        BELFEM_ASSERT( A.n_rows() == A.n_cols(),
            "Matrix A must be square ( is %lu x %lu )",
            ( long unsigned int ) A.n_rows(),
            ( long unsigned int ) A.n_cols() );
        BELFEM_ASSERT( jobvl == 'V' || jobvl == 'N',
            "unsupported jobvl flag '%c'", jobvl );
        BELFEM_ASSERT( jobvr == 'V' || jobvr == 'N',
            "unsupported jobvr flag '%c'", jobvr );

        // reals per LAPACK work entry
        constexpr int_t tRealsPerT =
            std::is_same< T, lapack::real_t< T > >::value ? 1 : 2 ;

        int_t n   = ( int_t ) A.n_rows();
        int_t lda = lapack::leading_dimension( A );

        W.set_size( n );

        // lapack requires ldvl/ldvr >= 1 even when not referenced
        int_t ldvl = 1 ;
        if ( jobvl == 'V' )
        {
            VL.set_size( n, n );
            ldvl = lapack::leading_dimension( VL );
        }

        int_t ldvr = 1 ;
        if ( jobvr == 'V' )
        {
            VR.set_size( n, n );
            ldvr = lapack::leading_dimension( VR );
        }

        // minimum work sizes: real flavors need 4n with eigenvectors and
        // 3n without, the complex ones need 2n
        int_t lwork = std::max< int_t >( 1, tRealsPerT == 1 ?
            ( ( jobvl == 'V' || jobvr == 'V' ) ? 4 * n : 3 * n ) : 2 * n );

        int_t info = 0 ;

        // required buffer: lwork entries of T in the head plus 2*n reals
        if ( static_cast< int_t >( Work.length() ) < tRealsPerT * lwork + 2 * n )
        {
            // ask lapack for the optimal size; one T entry plus the tail
            Work.set_size( tRealsPerT + 2 * n );

            int_t query = -1 ;

            lapack::geev( &jobvl, &jobvr, &n, A.data(), &lda, W.data(),
                VL.data(), &ldvl, VR.data(), &ldvr,
                reinterpret_cast< T * >( Work.data() ), &query,
                Work.data() + tRealsPerT, &info );

            BELFEM_ERROR( info == 0 || ! AbortOnError,
                "LAPACK geev workspace query has thrown an error: %i", ( int ) info );

            if ( info != 0 ) return info ;

            // the optimal size sits in the real part of the first entry
            lwork = lapack::work_size( Work( 0 ) );

            Work.set_size( tRealsPerT * lwork + 2 * n );
        }
        else
        {
            // use the full buffer the caller has provided
            lwork = ( ( int_t ) Work.length() - 2 * n ) / tRealsPerT ;
        }

        // work segment in the head, 2*n real scratch in the tail
        lapack::geev( &jobvl, &jobvr, &n, A.data(), &lda, W.data(),
            VL.data(), &ldvl, VR.data(), &ldvr,
            reinterpret_cast< T * >( Work.data() ), &lwork,
            Work.data() + tRealsPerT * lwork, &info );

        // info > 0 : the QR iteration failed; eigenvalues info+1 .. n
        //            have converged, no eigenvectors were computed
        BELFEM_ERROR( info == 0 || ! AbortOnError,
            "LAPACK geev has thrown an error: %i", ( int ) info );

        return info ;
    }

//------------------------------------------------------------------------------

    /**
     * @brief eigenvalues only of a general square matrix, see the full
     *        version above
     *
     * @param[in,out] A     square input matrix; destroyed by the call
     * @param[out]    W     the n eigenvalues — always complex valued
     * @param[in,out] Work  single real-valued scratch as in the full
     *                      version; no internal allocation
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: the QR iteration failed
     */
    template< typename T >
    int_t
    geev(
        Matrix< T > & A,
        Vector< lapack::cplx_t< T > > & W,
        Vector< lapack::real_t< T > > & Work,
        const bool AbortOnError = true )
    {
        BELFEM_ASSERT( A.n_rows() == A.n_cols(),
            "Matrix A must be square ( is %lu x %lu )",
            ( long unsigned int ) A.n_rows(),
            ( long unsigned int ) A.n_cols() );

        constexpr int_t tRealsPerT =
            std::is_same< T, lapack::real_t< T > >::value ? 1 : 2 ;

        char jobv = 'N' ;

        int_t n   = ( int_t ) A.n_rows();
        int_t lda = lapack::leading_dimension( A );

        W.set_size( n );

        // not referenced under jobvl = jobvr = 'N', but ld >= 1 is required
        int_t ldv = 1 ;

        int_t lwork = std::max< int_t >( 1, tRealsPerT == 1 ? 3 * n : 2 * n );

        int_t info = 0 ;

        if ( static_cast< int_t >( Work.length() ) < tRealsPerT * lwork + 2 * n )
        {
            Work.set_size( tRealsPerT + 2 * n );

            int_t query = -1 ;

            lapack::geev( &jobv, &jobv, &n, A.data(), &lda, W.data(),
                static_cast< T * >( nullptr ), &ldv,
                static_cast< T * >( nullptr ), &ldv,
                reinterpret_cast< T * >( Work.data() ), &query,
                Work.data() + tRealsPerT, &info );

            BELFEM_ERROR( info == 0 || ! AbortOnError,
                "LAPACK geev workspace query has thrown an error: %i", ( int ) info );

            if ( info != 0 ) return info ;

            lwork = lapack::work_size( Work( 0 ) );

            Work.set_size( tRealsPerT * lwork + 2 * n );
        }
        else
        {
            lwork = ( ( int_t ) Work.length() - 2 * n ) / tRealsPerT ;
        }

        lapack::geev( &jobv, &jobv, &n, A.data(), &lda, W.data(),
            static_cast< T * >( nullptr ), &ldv,
            static_cast< T * >( nullptr ), &ldv,
            reinterpret_cast< T * >( Work.data() ), &lwork,
            Work.data() + tRealsPerT * lwork, &info );

        BELFEM_ERROR( info == 0 || ! AbortOnError,
            "LAPACK geev has thrown an error: %i", ( int ) info );

        return info ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_GEEV_HPP
