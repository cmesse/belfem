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
#ifndef BELFEM_FN_FEM_ANDERSON_MIXING_HPP
#define BELFEM_FN_FEM_ANDERSON_MIXING_HPP

#include <cmath>

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_ShiftRegister.hpp"
#include "fn_gels.hpp"

namespace belfem
{
    namespace fem
    {
        //! reject a mixing solution whose largest coefficient exceeds this
        //! bound: the least-squares system has gone near-collinear and the
        //! extrapolation would be dominated by noise. Engineering default,
        //! meaningful only together with the column normalization below;
        //! recalibrate from the gamma norms logged in validation runs.
        constexpr real gAndersonGammaMax = 1.0e2 ;

//------------------------------------------------------------------------------

        /**
         * one type-II Anderson mixing step (Walker & Ni 2011) for the
         * fixed-point iteration x <- G( x ) with residual r = G( x ) - x :
         *
         *   gamma   = argmin || r_k - DR * gamma ||_2
         *   x_{k+1} = x_k + beta r_k - ( DX + beta DR ) gamma
         *
         * where the difference columns DR / DX are formed from the committed
         * history (newest aXHistory.size() snapshots) chained with the current
         * pair ( aX, aR ). The least-squares solve runs on column-normalized
         * DR via QR ( gels ); gamma is un-scaled afterwards.
         *
         * Failure policy: a rank-deficient solve ( info > 0 ), a non-finite
         * or oversized gamma ( gAndersonGammaMax ), or a vanishing difference
         * column drops the OLDEST column and retries once; if the retry fails
         * too, the plain relaxed Picard step  x + beta r  is written instead.
         * An illegal-argument error ( info < 0 ) is a programming error and
         * always aborts.
         *
         * With an empty history ( or on fallback ) the result is exactly the
         * relaxed Picard step: x_{k+1} = ( 1 - beta ) x + beta G( x ).
         *
         * All scratch is caller-owned and must be presized: aDeltaR to
         * n x depth, aRhs to n, aGamma and aColNorm to depth. aWork is grown
         * to the optimal gels size on first use and reused thereafter.
         *
         * @return number of history columns used ( 0 = plain step )
         */
        inline uint
        anderson_mixing_step(
            const Vector< real >                    & aX,        // x_k
            const Vector< real >                    & aR,        // r_k = G(x_k) - x_k
            const ShiftRegister< Vector< real > >   & aXHistory, // committed x, (0) = newest
            const ShiftRegister< Vector< real > >   & aRHistory, // committed r, (0) = newest
            const real                                aBeta,     // mixing = live relaxation
                  Matrix< real >                    & aDeltaR,   // scratch, destroyed by gels
                  Vector< real >                    & aRhs,      // scratch, destroyed by gels
                  Vector< real >                    & aWork,     // gels workspace
                  Vector< real >                    & aGamma,    // scratch, un-scaled solution
                  Vector< real >                    & aColNorm,  // scratch, column norms
                  Vector< real >                    & aXNew )    // out
        {
            const index_t tN = aX.length() ;

            BELFEM_ASSERT( aR.length() == tN, "size mismatch of x and r" );
            BELFEM_ASSERT( aXHistory.size() == aRHistory.size(),
                "x and r histories out of sync ( %u vs %u )",
                ( unsigned int ) aXHistory.size(),
                ( unsigned int ) aRHistory.size() );

            // the plain relaxed Picard step, also the fallback
            auto tPlainStep = [ & ]()
            {
                for ( index_t k = 0; k < tN; ++k )
                {
                    aXNew( k ) = aX( k ) + aBeta * aR( k );
                }
            };

            // one attempt with the newest tUse difference columns; returns
            // true if aGamma holds a usable un-scaled solution
            auto tTrySolve = [ & ]( const uint tUse ) -> bool
            {
                // chain point i ( 0 = oldest used ): history( tUse - 1 - i )
                // for i < tUse, the current pair for i = tUse
                for ( uint j = 0; j < tUse; ++j )
                {
                    const Vector< real > & tR1 = ( j + 1 == tUse ) ?
                        aR : aRHistory( tUse - 2 - j );
                    const Vector< real > & tR0 = aRHistory( tUse - 1 - j );

                    real tNorm = 0.0 ;
                    for ( index_t k = 0; k < tN; ++k )
                    {
                        const real tValue = tR1( k ) - tR0( k );
                        aDeltaR( k, j ) = tValue ;
                        tNorm += tValue * tValue ;
                    }
                    tNorm = std::sqrt( tNorm );

                    // a vanishing difference column can not be normalized;
                    // treat like rank deficiency
                    if ( tNorm < BELFEM_EPSILON )
                    {
                        return false ;
                    }
                    aColNorm( j ) = tNorm ;

                    const real tScale = 1.0 / tNorm ;
                    for ( index_t k = 0; k < tN; ++k )
                    {
                        aDeltaR( k, j ) *= tScale ;
                    }
                }

                // gels destroys the rhs; the solution lands in its head
                for ( index_t k = 0; k < tN; ++k )
                {
                    aRhs( k ) = aR( k );
                }

                // call the typed backend directly with the logical column
                // count: the scratch matrix keeps its full allocation and the
                // leading dimension makes the size mismatch legal
                char  tTrans = 'N' ;
                int_t tM     = tN ;
                int_t tCols  = tUse ;
                int_t tNrhs  = 1 ;
                int_t tLda   = lapack::leading_dimension( aDeltaR );
                int_t tLdb   = lapack::leading_dimension( aRhs );
                int_t tInfo  = 0 ;

                int_t tLwork = ( int_t ) aWork.length() ;
                if ( tLwork < 2 * tCols )
                {
                    // workspace query, then grow the buffer once
                    int_t tQuery = -1 ;
                    aWork.set_size( 1 );
                    lapack::gels( &tTrans, &tM, &tCols, &tNrhs,
                        aDeltaR.data(), &tLda, aRhs.data(), &tLdb,
                        aWork.data(), &tQuery, &tInfo );
                    BELFEM_ERROR( tInfo == 0,
                        "gels workspace query failed: %i", ( int ) tInfo );
                    aWork.set_size( lapack::work_size( aWork( 0 ) ) );
                    tLwork = ( int_t ) aWork.length() ;
                }

                lapack::gels( &tTrans, &tM, &tCols, &tNrhs,
                    aDeltaR.data(), &tLda, aRhs.data(), &tLdb,
                    aWork.data(), &tLwork, &tInfo );

                // an illegal argument is a programming error, always fatal
                BELFEM_ERROR( tInfo >= 0,
                    "gels reports an illegal argument: %i", ( int ) tInfo );

                // rank deficient: caller shrinks the window
                if ( tInfo > 0 )
                {
                    return false ;
                }

                // un-scale and guard the solution
                real tGammaMax = 0.0 ;
                for ( uint j = 0; j < tUse; ++j )
                {
                    const real tValue = aRhs( j ) / aColNorm( j );
                    if ( ! std::isfinite( tValue ) )
                    {
                        return false ;
                    }
                    aGamma( j ) = tValue ;
                    tGammaMax = std::max( tGammaMax, std::abs( tValue ) );
                }
                return tGammaMax <= gAndersonGammaMax ;
            };

            // window: newest tUse columns; drop the oldest and retry once.
            // A window wider than the system would make the least-squares
            // problem underdetermined and overrun the rhs buffer ( gels needs
            // max( m, n ) rows ), so clamp — this must hold in release too
            uint tUse = aXHistory.size() ;
            if ( ( index_t ) tUse > tN )
            {
                tUse = ( uint ) tN ;
            }

            bool tOK = tUse > 0 ? tTrySolve( tUse ) : false ;
            if ( ! tOK && tUse > 1 )
            {
                --tUse ;
                tOK = tTrySolve( tUse );
            }
            if ( ! tOK )
            {
                tPlainStep();
                return 0 ;
            }

            // x_{k+1} = x + beta r - sum_j gamma_j ( dx_j + beta dr_j ),
            // differences re-formed from the histories ( gels destroyed the
            // normalized copy )
            tPlainStep();
            for ( uint j = 0; j < tUse; ++j )
            {
                const bool tHead = ( j + 1 == tUse );
                const Vector< real > & tX1 = tHead ? aX : aXHistory( tUse - 2 - j );
                const Vector< real > & tX0 = aXHistory( tUse - 1 - j );
                const Vector< real > & tR1 = tHead ? aR : aRHistory( tUse - 2 - j );
                const Vector< real > & tR0 = aRHistory( tUse - 1 - j );

                const real tGamma = aGamma( j );
                for ( index_t k = 0; k < tN; ++k )
                {
                    aXNew( k ) -= tGamma * (   ( tX1( k ) - tX0( k ) )
                                    + aBeta * ( tR1( k ) - tR0( k ) ) );
                }
            }
            return tUse ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_FEM_ANDERSON_MIXING_HPP
