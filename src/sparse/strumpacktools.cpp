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
#include "commtools.hpp"
#include "strumpacktools.hpp"
namespace belfem
{
    namespace sparse
    {
        string
        strumpack_message( StrumpackReturnCode aCode )
        {
#ifdef BELFEM_STRUMPACK
            switch( aCode )
            {
                case strumpack::ReturnCode::SUCCESS             : return "Operation completed successfully." ;
                case strumpack::ReturnCode::MATRIX_NOT_SET      : return "The input matrix was not set." ;
                case strumpack::ReturnCode::REORDERING_ERROR    : return "The matrix reordering failed." ;
                case strumpack::ReturnCode::ZERO_PIVOT          : return "A zero pivot was encountered. If matrix matching was disabled ( matching : off ), re-enable it." ;
                case strumpack::ReturnCode::NO_CONVERGENCE      : return "The iterative solver did not converge." ;
                case strumpack::ReturnCode::INACCURATE_INERTIA  : return "Inertia could not be computed." ;
                default : return "Unknown error.";
            }
#else
            return "we are not linked to STRUMPACK" ;
#endif
        }

        void
        set_strumpack_options( const SolverParameters & aParams,
                               StrumpackOptions & aOpts,
                               const index_t aNumRows )
        {
#ifdef BELFEM_STRUMPACK

            // Replace near-zero pivots with a threshold value during factorization.
            // This improves robustness for ill-conditioned matrices (e.g. after
            // edge deduplication) at slight accuracy cost, which the outer
            // iterative solver (REFINE or GMRES) will correct.
            aOpts.enable_replace_tiny_pivots();

            // MC64 matching is on by default: its permutation prevents exact
            // zero pivots in the coupled h-phi Jacobian, which neither
            // equilibration nor replace_tiny_pivots can absorb. In MPI mode
            // STRUMPACK gathers the full matrix on rank 0 for it ( once per
            // initialize ); ( matching : off ) opts out for verified cases.
            if ( ! aParams.use_matrix_matching() )
            {
                aOpts.set_matching( strumpack::MatchingJob::NONE );
            }

            // METIS_NodeNDP returns the separator tree; METIS_NodeND does not,
            // so STRUMPACK falls back to rebuilding a supernodal tree from the
            // elimination tree ( "supernodal tree was built from etree" in its
            // log ) and that reconstruction is what makes the tree deep. On
            // tapestack3d it reached 56 levels, at which point STRUMPACK warned
            // that it "does not handle this safely, which could lead to
            // segmentation faults due to stack overflows" and recommended this
            // exact flag. A deep tree is also the likely reason the factor
            // memory concentrated on a single rank ( 23 GiB vs 2.8 GiB ).
            //
            // Applies to SERIAL METIS only -- PARMETIS and SCOTCH ignore it, so
            // it is safe to set unconditionally. METIS_NodeNDP is undocumented
            // in METIS and declared by STRUMPACK itself, hence the opt-out.
            if ( aParams.use_metis_nodendp() )
            {
                aOpts.enable_METIS_NodeNDP();
            }
            else
            {
                aOpts.enable_METIS_NodeND();
            }

            switch ( aParams.reordering_method() )
            {
                case ReorderingMethod::NATURAL :
                {
                    aOpts.set_reordering_method( strumpack::ReorderingStrategy::NATURAL );
                    break;
                }
                // parmetis / ptscotch select what metis / scotch already
                // mean here: the parallel library above the serial one
                case ReorderingMethod::METIS :
                case ReorderingMethod::PARMETIS :
                {
                    if ( gComm.size() > 1 )
                    {
                        aOpts.set_reordering_method( strumpack::ReorderingStrategy::PARMETIS );
                    }
                    else
                    {
                        aOpts.set_reordering_method( strumpack::ReorderingStrategy::METIS );
                    }

                    break;
                }
                case ReorderingMethod::SCOTCH :
                case ReorderingMethod::PTSCOTCH :
                {
                    if ( gComm.size() > 1 )
                    {
                        aOpts.set_reordering_method( strumpack::ReorderingStrategy::PTSCOTCH );
                    }
                    else
                    {
                        aOpts.set_reordering_method( strumpack::ReorderingStrategy::SCOTCH );
                    }
                    break;
                }
                default:
                {
                    // pass
                }
            }

            switch ( aParams.compression_method() )
            {
                case CompressionMethod::BLR :
                {
                    aOpts.set_compression(strumpack::CompressionType::BLR);

                    // larger systems tolerate larger low-rank leaf blocks
                    aOpts.BLR_options().set_leaf_size(
                        aNumRows < 100000 ? 128 :
                        aNumRows < 250000 ? 256 : 512 );

                    // STRUMPACK's default BLR tolerance ( 1e-4 ) makes the
                    // factorization a lossy preconditioner whose GMRES can
                    // stagnate silently; that noise stalls the Newton loop
                    // ( observed 2026-07-06, N = 72k, all rank counts ).
                    aOpts.BLR_options().set_rel_tol(
                        aParams.compression_cutoff() );
                    break;
                }
                default:
                {
                    // OFF and AUTOMATIC: compression disabled. BLR pays off
                    // only for very large systems and must be requested
                    // explicitly via ( compression scheme : blr ).
                    aOpts.set_compression( strumpack::CompressionType::NONE );
                }
            }

            switch ( aParams.krylov_method() )
            {
                case KrylovMethod::AUTO :
                {
                    // GMRES preconditioned by the factorization, NOT
                    // STRUMPACK's own AUTO ( which picks Richardson-style
                    // REFINE on the exact path ). Ruled 2026-08-16 for
                    // uniform krylov-method semantics across libraries,
                    // and because GMRES is strictly more robust here: on
                    // a well-conditioned matrix the first preconditioned
                    // iterate lands near machine precision ( same cost ),
                    // while on an ill-conditioned one ( cond ~ 5e11 on
                    // tapestack3d: the raw factor gives ~4 digits )
                    // GMRES minimizes the residual where plain refinement
                    // provably stalled ( the refinement floor trace ). The
                    // maxit cap below applies to this path too.
                    // ( krylov method : preonly ) skips the outer loop
                    // entirely — an expert setting for measured cases
                    // only; on an ill-conditioned matrix it returns the
                    // four-digit factor solution, see §4.1
                    aOpts.set_Krylov_solver(
                        strumpack::KrylovSolver::PREC_GMRES );
                    break;
                }
                case KrylovMethod::PREONLY :
                {
                    aOpts.set_Krylov_solver( strumpack::KrylovSolver::DIRECT );
                    break;
                }
                case KrylovMethod::BCGS :
                {
                    if ( aParams.preconditioner() == Preconditioner::NONE )
                    {
                        aOpts.set_Krylov_solver( strumpack::KrylovSolver::BICGSTAB );
                    }
                    else
                    {
                        aOpts.set_Krylov_solver( strumpack::KrylovSolver::PREC_BICGSTAB );
                    }
                    break;
                }
                case KrylovMethod::GMRES :
                {
                    if ( aParams.preconditioner() == Preconditioner::NONE )
                    {
                        aOpts.set_Krylov_solver( strumpack::KrylovSolver::GMRES );
                    }
                    else
                    {
                        aOpts.set_Krylov_solver( strumpack::KrylovSolver::PREC_GMRES );
                    }
                    break ;
                }
                default:
                {
                    // CG/CGS/IBCGS/TFQMR are PETSc methods with no STRUMPACK
                    // mapping. Falling through silently would leave the
                    // library's own KrylovSolver::AUTO, which picks
                    // Richardson-style REFINE on the exact path — the deck
                    // would run a solver nobody asked for
                    BELFEM_ERROR( false,
                        "STRUMPACK has no mapping for krylov method '%s' "
                        "( supported: auto, preonly, gmres, bcgs )",
                        to_string( aParams.krylov_method() ).c_str() );
                }
            }
            // ALWAYS apply the shared relative tolerance ( class default
            // 1e-10 ), mirroring the abs_tol treatment below. Until
            // 2026-08-18 this was gated on have_relative_tolerance(): a
            // REFINE-era stall fix ( an unmeetable relative target
            // on a numerically zero RHS ground IterativeRefinementMPI
            // against the library's 5000-iteration cap ). That geometry is
            // gone: AUTO maps to PREC_GMRES since 2026-08-16, maxit is
            // capped at 50 below, and abs_tol provides the small-RHS
            // escape. The gate's cost was real: an unstated deck inherited
            // the library rel_tol 1e-6, and the tapestack3d A/B
            // ( 2026-08-18, 1e-8 vs 1e-10 from the same restart ) proved
            // the loose exit test WAS the printed nonlinear residual —
            // Picard "stalled" at the linear exit ( -85 dB ), Newton was
            // invoked on exit-test noise, and the timestep collapsed. At
            // 1e-10 the same solves finish their dive ( ~1e-15 ), steps
            // converge in two Picard iterates, and the collapse never
            // happens. Do not re-introduce the guard.
            aOpts.set_rel_tol( aParams.relative_tolerance() );

            // ALWAYS override the library's absolute tolerance: STRUMPACK
            // defaults abs_tol to 1e-10, and on a small right-hand side
            // ( transient startup, ||b|| ~ 1e-2 ) that exit test sits
            // exactly at a 1e-11 RELATIVE nonlinear target — the outer
            // GMRES then stops a factor ~300 short of what the Newton
            // loop needs and convergence degenerates into an overshoot
            // lottery ( jury 2026-08-17: measured linear finals
            // clustered at 6e-11..1.4e-10 with the stall entry at
            // floor/||b|| = -85 dB exactly as predicted ). 1e-14 restores
            // the headroom; the maxit cap below bounds what it may cost
            aOpts.set_abs_tol( aParams.absolute_tolerance() );

            // bound the outer loop either way: an unmeetable deck-stated
            // target must cost minutes, not grind toward the library cap
            // of 5000 ( the hour-long refinement stall ). NOTE the cap's real semantics,
            // source-verified 2026-08-16: STRUMPACK returns SUCCESS at
            // maxit with the best-effort solution ( the GMRES kernels
            // return only the residual and the callers discard it ) — so
            // the cap is a cost bound, NOT a named failure; the nonlinear
            // loop above judges the true residual
            aOpts.set_maxit( 50 );
#endif
        }

    }
}