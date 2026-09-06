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
 * @brief Schur decomposition of a general square matrix (LAPACK ?gees).
 * @ingroup grp_linalg
 *
 * Thin wrapper over the Fortran routine: BELFEM containers are passed straight through
 * after their leading dimensions are worked out. Arguments are not copied -- see each
 * parameter for what is overwritten in place.
 */

#ifndef BELFEM_FN_GEES_HPP
#define BELFEM_FN_GEES_HPP

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

            // gees particularities: like geev, the real flavors return the
            // eigenvalues as wr/wi and have no rwork, while the complex ones
            // take one complex w plus a REAL rwork ( size n, not 2*n ). The
            // sort callback SELECT is a Fortran LOGICAL function that takes
            // its arguments by reference: two reals for s/d, one complex for
            // c/z. bwork is a LOGICAL array, referenced only when sorting.

            typedef int_t ( * sgees_select_t )( const float *, const float * );
            typedef int_t ( * dgees_select_t )( const double *, const double * );
            typedef int_t ( * cgees_select_t )( const cplx_float_t * );
            typedef int_t ( * zgees_select_t )( const cplx_double_t * );

//------------------------------------------------------------------------------

            void sgees_(
                char * jobvs,
                char * sort,
                sgees_select_t select,
                int_t * n,
                float * a,
                int_t * lda,
                int_t * sdim,
                float * wr,
                float * wi,
                float * vs,
                int_t * ldvs,
                float * work,
                int_t * lwork,
                int_t * bwork,
                int_t * info,
                fortran_charlen_t lj,
                fortran_charlen_t ls );

//------------------------------------------------------------------------------

            void dgees_(
                char * jobvs,
                char * sort,
                dgees_select_t select,
                int_t  * n,
                double * a,
                int_t  * lda,
                int_t  * sdim,
                double * wr,
                double * wi,
                double * vs,
                int_t  * ldvs,
                double * work,
                int_t  * lwork,
                int_t  * bwork,
                int_t  * info,
                fortran_charlen_t lj,
                fortran_charlen_t ls );

//------------------------------------------------------------------------------

            void cgees_(
                char         * jobvs,
                char         * sort,
                cgees_select_t select,
                int_t        * n,
                cplx_float_t * a,
                int_t        * lda,
                int_t        * sdim,
                cplx_float_t * w,
                cplx_float_t * vs,
                int_t        * ldvs,
                cplx_float_t * work,
                int_t        * lwork,
                float        * rwork,
                int_t        * bwork,
                int_t        * info,
                fortran_charlen_t lj,
                fortran_charlen_t ls );

//------------------------------------------------------------------------------

            void zgees_(
                char * jobvs,
                char * sort,
                zgees_select_t select,
                int_t * n,
                cplx_double_t * a,
                int_t * lda,
                int_t * sdim,
                cplx_double_t * w,
                cplx_double_t * vs,
                int_t * ldvs,
                cplx_double_t * work,
                int_t * lwork,
                double * rwork,
                int_t * bwork,
                int_t * info,
                fortran_charlen_t lj,
                fortran_charlen_t ls );

//------------------------------------------------------------------------------
#ifdef __cplusplus
        }
#endif
//------------------------------------------------------------------------------

        // user-facing select callback type: two real arguments for the real
        // flavors, one complex argument for the complex ones ( Fortran
        // passes by reference in both cases )
        template< typename T >
        struct gees_select
        {
            typedef int_t ( * type )( const T *, const T * ) ;
        };

        template< typename R >
        struct gees_select< std::complex< R > >
        {
            typedef int_t ( * type )( const std::complex< R > * ) ;
        };

        template< typename T >
        using gees_select_t = typename gees_select< T >::type ;

//------------------------------------------------------------------------------

        // unified dispatch. As in geev, rwork is the real scratch: wr/wi
        // ( 2*n ) for the real flavors, LAPACK's rwork ( n ) for the
        // complex ones.

        template< typename T >
        void gees(
            const char            * jobvs,
            const char            * sort,
                  gees_select_t< T > select,
            const int_t           * n,
                  T               * a,
            const int_t           * lda,
                  int_t           * sdim,
                  cplx_t< T >     * w,
                  T               * vs,
            const int_t           * ldvs,
                  T               * work,
            const int_t           * lwork,
                  real_t< T >     * rwork,
                  int_t           * bwork,
                  int_t           * info )
        {
            static_assert( dependent_false< T >,
                "gees not implemented for selected data type" );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void gees(
            const char            * jobvs,
            const char            * sort,
                  sgees_select_t    select,
            const int_t           * n,
                  float           * a,
            const int_t           * lda,
                  int_t           * sdim,
                  std::complex< float > * w,
                  float           * vs,
            const int_t           * ldvs,
                  float           * work,
            const int_t           * lwork,
                  float           * rwork,
                  int_t           * bwork,
                  int_t           * info )
        {
            sgees_(
                const_cast< char * >( jobvs ),
                const_cast< char * >( sort ),
                select,
                const_cast< int_t * >( n ),
                a,
                const_cast< int_t * >( lda ),
                sdim,
                rwork,
                rwork + *n,
                vs,
                const_cast< int_t * >( ldvs ),
                work,
                const_cast< int_t * >( lwork ),
                bwork,
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
        inline void gees(
            const char            * jobvs,
            const char            * sort,
                  dgees_select_t    select,
            const int_t           * n,
                  double          * a,
            const int_t           * lda,
                  int_t           * sdim,
                  std::complex< double > * w,
                  double          * vs,
            const int_t           * ldvs,
                  double          * work,
            const int_t           * lwork,
                  double          * rwork,
                  int_t           * bwork,
                  int_t           * info )
        {
            dgees_(
                const_cast< char * >( jobvs ),
                const_cast< char * >( sort ),
                select,
                const_cast< int_t * >( n ),
                a,
                const_cast< int_t * >( lda ),
                sdim,
                rwork,
                rwork + *n,
                vs,
                const_cast< int_t * >( ldvs ),
                work,
                const_cast< int_t * >( lwork ),
                bwork,
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
        inline void gees(
            const char            * jobvs,
            const char            * sort,
                  gees_select_t< std::complex< float > > select,
            const int_t           * n,
                  std::complex< float > * a,
            const int_t           * lda,
                  int_t           * sdim,
                  std::complex< float > * w,
                  std::complex< float > * vs,
            const int_t           * ldvs,
                  std::complex< float > * work,
            const int_t           * lwork,
                  float           * rwork,
                  int_t           * bwork,
                  int_t           * info )
        {
            cgees_(
                const_cast< char * >( jobvs ),
                const_cast< char * >( sort ),
                // function-pointer type pun: the callback parameter is
                // layout-compatible ( std::complex< float > vs
                // cplx_float_t ), same cast direction as the data arrays;
                // strictly UB in ISO C++, universally sound on the SysV
                // ABI ( one pointer argument, integer-width return )
                reinterpret_cast< cgees_select_t >( select ),
                const_cast< int_t * >( n ),
                reinterpret_cast< cplx_float_t * >( a ),
                const_cast< int_t * >( lda ),
                sdim,
                reinterpret_cast< cplx_float_t * >( w ),
                reinterpret_cast< cplx_float_t * >( vs ),
                const_cast< int_t * >( ldvs ),
                reinterpret_cast< cplx_float_t * >( work ),
                const_cast< int_t * >( lwork ),
                rwork,
                bwork,
                info,
                1, 1 );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline void gees(
            const char            * jobvs,
            const char            * sort,
                  gees_select_t< std::complex< double > > select,
            const int_t           * n,
                  std::complex< double > * a,
            const int_t           * lda,
                  int_t           * sdim,
                  std::complex< double > * w,
                  std::complex< double > * vs,
            const int_t           * ldvs,
                  std::complex< double > * work,
            const int_t           * lwork,
                  double          * rwork,
                  int_t           * bwork,
                  int_t           * info )
        {
            zgees_(
                const_cast< char * >( jobvs ),
                const_cast< char * >( sort ),
                reinterpret_cast< zgees_select_t >( select ),
                const_cast< int_t * >( n ),
                reinterpret_cast< cplx_double_t * >( a ),
                const_cast< int_t * >( lda ),
                sdim,
                reinterpret_cast< cplx_double_t * >( w ),
                reinterpret_cast< cplx_double_t * >( vs ),
                const_cast< int_t * >( ldvs ),
                reinterpret_cast< cplx_double_t * >( work ),
                const_cast< int_t * >( lwork ),
                rwork,
                bwork,
                info,
                1, 1 );
        }

//------------------------------------------------------------------------------
    } /* end namespace lapack */

    /**
     * @brief Schur decomposition of a general square matrix via LAPACK
     *        ?gees: A = Z * T * Z^T for the real flavors,
     *        A = Z * T * Z^H for the complex ones
     *
     * @param[in,out] A      square input matrix; overwritten with the
     *                       Schur form T — upper triangular for complex T,
     *                       quasi upper triangular with 2x2 blocks for the
     *                       complex conjugate eigenvalue pairs of the real
     *                       flavors
     * @param[out]    W      the n eigenvalues — always complex valued
     * @param[out]    VS     Schur vectors Z, resized to n x n for
     *                       jobvs = 'V'; untouched for 'N'
     * @param[in,out] Work   single real-valued scratch, grown to the
     *                       optimal size on first use: LAPACK's work array
     *                       in its head ( for complex T reinterpreted in
     *                       place ), the real scratch in its tail — wr/wi
     *                       ( 2*n ) for real T, rwork ( n ) for complex T.
     *                       No internal allocation.
     * @param[in,out] BWork  LOGICAL scratch, referenced only when sorting;
     *                       grown to n entries for sort = 'S'
     * @param[in]     select moves the selected eigenvalues to the top left
     *                       of the Schur form when sort = 'S'; pass
     *                       nullptr with sort = 'N'.
     *                       NOTE for real T: a complex conjugate pair is
     *                       selected as a whole if select is true for
     *                       EITHER member, and it counts as TWO towards
     *                       sdim — so sdim can differ from the number of
     *                       true returns, and can come out odd when reals
     *                       and pairs mix
     * @param[in]     jobvs  'V' or 'N'
     * @param[in]     sort   'S' or 'N'
     * @param[out]    aSdim  if given, receives the number of selected
     *                       eigenvalues
     * @param[in]     AbortOnError  abort on info != 0, or return info
     * @return info   0 on success; > 0: <= n the QR iteration failed,
     *                n+1 reordering failed because some eigenvalues are
     *                too close to separate, n+2 roundoff changed the
     *                leading eigenvalues after reordering
     */
    template< typename T >
    int_t
    gees(
        Matrix< T > & A,
        Vector< lapack::cplx_t< T > > & W,
        Matrix< T > & VS,
        Vector< lapack::real_t< T > > & Work,
        Vector< int_t > & BWork,
        lapack::gees_select_t< T > select = nullptr,
        const char jobvs = 'V',
        const char sort  = 'N',
        int_t * aSdim = nullptr,
        const bool AbortOnError = true )
    {
        BELFEM_ASSERT( A.n_rows() == A.n_cols(),
            "Matrix A must be square ( is %lu x %lu )",
            ( long unsigned int ) A.n_rows(),
            ( long unsigned int ) A.n_cols() );
        BELFEM_ASSERT( jobvs == 'V' || jobvs == 'N',
            "unsupported jobvs flag '%c'", jobvs );
        BELFEM_ASSERT( sort == 'S' || sort == 'N',
            "unsupported sort flag '%c'", sort );
        BELFEM_ASSERT( ( sort == 'S' ) == ( select != nullptr ),
            "sort = 'S' requires a select callback, sort = 'N' forbids it" );

        // reals per LAPACK work entry
        constexpr int_t tRealsPerT =
            std::is_same< T, lapack::real_t< T > >::value ? 1 : 2 ;

        int_t n   = ( int_t ) A.n_rows();
        int_t lda = lapack::leading_dimension( A );

        W.set_size( n );

        // lapack requires ldvs >= 1 even when not referenced
        int_t ldvs = 1 ;
        if ( jobvs == 'V' )
        {
            VS.set_size( n, n );
            ldvs = lapack::leading_dimension( VS );
        }

        // the LOGICAL scratch is referenced only when sorting
        if ( sort == 'S' && static_cast< int_t >( BWork.length() ) < n )
        {
            BWork.set_size( n );
        }

        // real tail: wr/wi for the real flavors, rwork for the complex ones
        const int_t tTail = ( tRealsPerT == 1 ) ? 2 * n : n ;

        // minimum work sizes: real flavors 3n, complex 2n
        int_t lwork = std::max< int_t >( 1, tRealsPerT == 1 ? 3 * n : 2 * n );

        int_t sdim = 0 ;
        int_t info = 0 ;

        // check length of work array: lwork entries of T plus the tail
        if ( static_cast< int_t >( Work.length() ) < tRealsPerT * lwork + tTail )
        {
            // ask lapack for the optimal size; one T entry plus the tail
            Work.set_size( tRealsPerT + tTail );

            int_t query = -1 ;

            lapack::gees( &jobvs, &sort, select, &n, A.data(), &lda, &sdim,
                W.data(), VS.data(), &ldvs,
                reinterpret_cast< T * >( Work.data() ), &query,
                Work.data() + tRealsPerT, BWork.data(), &info );

            BELFEM_ERROR( info == 0 || ! AbortOnError,
                "LAPACK gees workspace query has thrown an error: %i", ( int ) info );

            if ( info != 0 ) return info ;

            // the optimal size sits in the real part of the first entry
            lwork = lapack::work_size( Work( 0 ) );

            Work.set_size( tRealsPerT * lwork + tTail );
        }
        else
        {
            // use the full buffer the caller has provided
            lwork = ( ( int_t ) Work.length() - tTail ) / tRealsPerT ;
        }

        // work segment in the head, real scratch in the tail
        lapack::gees( &jobvs, &sort, select, &n, A.data(), &lda, &sdim,
            W.data(), VS.data(), &ldvs,
            reinterpret_cast< T * >( Work.data() ), &lwork,
            Work.data() + tRealsPerT * lwork, BWork.data(), &info );

        BELFEM_ERROR( info == 0 || ! AbortOnError,
            "LAPACK gees has thrown an error: %i", ( int ) info );

        if ( aSdim != nullptr )
        {
            *aSdim = sdim ;
        }

        return info ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_FN_GEES_HPP
