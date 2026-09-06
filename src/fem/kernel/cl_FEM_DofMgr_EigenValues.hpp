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
#ifndef CL_FEM_DOFMGR_EIGENVALUES_HPP
#define CL_FEM_DOFMGR_EIGENVALUES_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Bitset.hpp"
#include "cl_Solver.hpp"
#include "arpacktools.hpp"

namespace belfem
{
    class SpMatrix ;

#ifdef BELFEM_PARPACK
    namespace sparse
    {
        template< typename T > class DistMatrixCSR ;
    }
#endif

    //! how an eigenvalue diagnostic ended. The MASTER decides this and
    //! BROADCASTS it: every branch that may call run() again has to be taken
    //! from the same value on every rank, or the ranks part company inside a
    //! collective. Deriving it rank-locally from Ritz data is NOT safe -- in a
    //! multi-rank build without PARPACK the workers' Ritz buffers are stale by
    //! construction ( see lambda_real() below )
    enum class EigenOutcome
    {
        Ok                  = 0,   //!< a usable value came back
        NotConverged        = 1,   //!< ARPACK ran out of restarts
        NotPositiveDefinite = 2,   //!< shift-invert returned lambda_min <= 0: the ratio is not a condition number
        Complex             = 3,   //!< the winning Ritz value is not real
        Inconsistent        = 4,   //!< shift-invert returned lambda_min > rho, which cannot be
        Infeasible          = 5,   //!< no legal ncv fits the basis budget, or no usable solver
        SolverFailed        = 6    //!< the factorization or a frozen solve failed. NOT
                                   //!< NotConverged ( that is ARPACK's ) and NOT Infeasible
                                   //!< ( that is size/build ) -- conflating them would send
                                   //!< the reader after the wrong cause
    };

    enum class ArpackWhat
    {
        LambaMinMagn = 0,
        LambaMaxMagn = 1,
        LambaMinReal = 2,
        LambaMaxReal = 3,
        LambaMinImag = 4,
        LambaMaxImag = 5
    };

    namespace fem
    {
        class DofManager ;

        namespace dofmgr
        {
            class EigenValues
            {
                const proc_t mCommRank ;
                const proc_t mCommSize ;
                DofManager * mParent ;

                // matrix data
                SpMatrix   * mK = nullptr ;
                SpMatrix   * mM = nullptr ;
                int_t  mN   = 0 ;
                int_t  mNNZ = 0 ;


                // options
                int_t          mNumMinVals   = 1 ;
                int_t          mNumMaxVals   = 1 ;
                int_t          mNumMaxIter   = 300 ;

                // the three tuning knobs used to be hardwired here and were
                // reachable from nowhere -- no caller in the tree ever invoked
                // the setters. They are sized automatically now ( configure() ),
                // and these flags are what keeps an explicit caller in charge:
                // once a setter has been used, the heuristic leaves that knob
                // alone for the rest of the run
                bool           mSubspaceSizeExplicit = false ;
                bool           mToleranceExplicit    = false ;
                bool           mNumMaxIterExplicit   = false ;

                // ARPACK accepts a Ritz value when its error bound falls below
                // mEpsilon * |lambda|, and that bound cannot go below the
                // backward error of the matvec, roughly eps_mach * |lambda_max|.
                // The small end is therefore only reachable while
                //
                //     mEpsilon  >  eps_mach * kappa      ( ~ 1e-16 * kappa )
                //
                // so 1e-7 caps the usable condition number at about 1e9 -- below
                // what an h-phi Jacobian reaches on a refined mesh, where SM then
                // grinds to mNumMaxIter and converges nothing. 1e-4 lifts that
                // ceiling to ~1e12 and still pins lambda to four digits, which is
                // far more than a conditioning DIAGNOSTIC needs.
                // Tighten only with that inequality in hand
                real           mEpsilon      = 1e-4 ;

                // floor for the Krylov subspace handed to ( p )dnaupd. Remark 4
                // of dnaupd: raising ncv at fixed nev usually REDUCES the total
                // OP*x count, at the price of an n x ncv basis
                int_t          mSubspaceSize = 20 ;

                // Tolerance for the FOLDED run, which is a different job from
                // the plain one. lambda_min falls out of a cancellation,
                // lambda_min = sigma - mu, and mu is of order sigma, so an
                // error of tol*sigma in mu lands undivided on lambda_min:
                //
                //     rel_err( lambda_min )  ~  max( tol, eps_mach ) * kappa
                //
                // 1e-4 would therefore be worthless at the small end ( it caps
                // kappa at ~1e4 before the answer is all noise ), while the
                // FOLDED problem is an exterior one and can afford a tight
                // tolerance -- it converges in a handful of restarts either way.
                // The floor is eps_mach*kappa and cannot be bought off
                real           mFoldEpsilon  = 1e-10 ;

                // NOTE: the spectral FOLD that this class used to run for the
                // small end has been REPLACED by shift-invert
                // ( run_shift_invert ), because an affine shift cannot improve
                // eigenvalue SEPARATION and separation is what sets the
                // convergence rate -- measured, see compute_conditioning().
                // The sigma argument on the Fortran drivers is deliberately
                // KEPT: it costs nothing, the folded operator is still correct
                // where it converges, and the shift-invert shims share those
                // drivers' argument shape

                // ceiling on the Arnoldi basis, which is n x ncv doubles. When
                // not even the smallest legal ncv fits, the diagnostic reports
                // Infeasible rather than allocating past it
                size_t         mBasisBudget  = 512UL * 1024UL * 1024UL ;

                // output
                Vector< int_t > mInfo ;

                // the SIGNED real and imaginary part of the Ritz value that won
                // the last reduction. The function value is a magnitude, which
                // is all the caller used to need ; the fold needs the sign to
                // tell a positive definite spectrum from one straddling zero,
                // and the imaginary part to notice a complex pair. Copied from
                // the winning SLOT, so the two always describe the same value
                real            mExtremalReal = BELFEM_QUIET_NAN ;
                real            mExtremalImag = BELFEM_QUIET_NAN ;

                // set once an end of the spectrum is proven out of reach, and
                // never cleared: retrying it every timestep cost 4.7 s per step
                // on tapestack3d and produced nothing. Assigned ONLY from a
                // broadcast outcome -- see EigenOutcome
                bool            mDiagnosticUnavailable = false ;

                // what ended the last compute_conditioning(), for the caller
                // and for the one-time message
                EigenOutcome    mOutcome = EigenOutcome::Ok ;

                // which of the two runs is in flight, for the messages. A
                // failure that does not say whether the plain end or the fold
                // died cannot be acted on -- measured 2026-08-28, when exactly
                // that ambiguity cost a diagnostic cycle
                const char *    mRunLabel = "eigenvalue" ;

                // which ARPACK family drives this field. TOLD, not detected:
                // the thermal Jacobian is symmetric by construction and the
                // magnetic h-phi one is not ( Christian, 2026-08-28 ). Detecting
                // it would mean comparing A against A^T across a distribution
                // that may be transposed column blocks already ( see the CSC
                // note in run_parpack ), which is a bigger job than the caller
                // simply knowing. Defaults to the SAFE answer: the nonsymmetric
                // driver is correct for a symmetric matrix, merely slower
                bool            mSymmetric = false ;

                // ---- shift-invert ( mode 3 ) state -------------------------
                // ARPACK's reverse-communication arrays. Members, not locals:
                // the loop hands them to Fortran across many calls and dseupd
                // requires them untouched between the last step and the
                // extraction, so nothing transient may own them
                Vector< real >  mSiResid ;      // [ n ]
                Vector< real >  mSiBasis ;      // [ n * ncv ], the Lanczos basis
                Vector< real >  mSiWorkD ;      // [ 3 * n ]
                Vector< real >  mSiWorkL ;      // [ ncv * ( ncv + 8 ) ]
                Vector< int_t > mSiIparam ;     // [ 11 ]
                Vector< int_t > mSiIpntr ;      // [ 11 ]
                Vector< real >  mSiLambda ;     // [ nev ], eigenvalues of A

                // right-hand side and solution of one inverse application.
                // Sized once ; the solve writes into them every iteration
                Vector< real >  mSiRhs ;
                Vector< real >  mSiLhs ;

                // dedicated solver for the inverse operator. NOT mSolver,
                // which belongs to compute_lambda_max() and is built from the
                // parent's solver type -- this one is MUMPS specifically,
                // because it is the wrapper that can hold a factorization
                // still and solve against it repeatedly
                Solver *        mShiftInvertSolver = nullptr ;

                // the tolerance the last run was actually given. configure()
                // picks between mEpsilon and mFoldEpsilon and an explicit
                // set_tolerance() overrides both, so recomputing that choice at
                // the point of use duplicates the rule and lets the two drift.
                // Recorded once, read where it is reported
                real            mEffectiveTolerance = BELFEM_QUIET_NAN ;

                // set by a driver when configure() found no legal ncv. It is
                // the difference between "this iterate did not converge", which
                // gets a few strikes, and "no run of this size can ever be
                // attempted", which latches at once -- and a bare NaN cannot
                // tell the two apart. Cleared at the top of every run
                bool            mSubspaceInfeasible = false ;

                // consecutive failures. Every outcome but Infeasible gets a few
                // strikes ( a non-positive iterate may be transient: the Jacobian
                // is rebuilt every step ); Infeasible latches at once because it
                // depends on size and budget, not on the numbers. Every rank
                // counts the same broadcast outcome
                int_t           mFailureCount = 0 ;
                int_t           mMaxFailures  = 3 ;

                // which driver produced the last result. It depends on the
                // rank count AND on the build ( multi-rank without PARPACK
                // still runs ARPACK on the master ), so a caller cannot
                // derive it and has to be told
                string          mBackendLabel = "ARPACK" ;

                // work data
                string          mWhich;
                real mLambdaMax = BELFEM_QUIET_NAN ;

                Vector< real >  mLambdaReal ;
                Vector< real >  mLambdaImag ;

                Vector< real >  mXreal ;
                Vector< real >  mYreal ;
                Vector< real >  mZreal ;

                Vector< real >  mXimag ;
                Vector< real >  mYimag ;
                Vector< real >  mZimag ;

                Vector< cplx >  mX ;
                Vector< cplx >  mY ;
                Vector< cplx >  mZ ;

                bool mMatrixFlag = false ;
                bool mFirstRun = true ;

                Solver * mSolver = nullptr ;

                int_t mOriginalBase = 0 ;

#ifdef BELFEM_PARPACK
                // row-distributed view of the master's Jacobian. The master
                // owns the full matrix ; every other rank holds only its
                // assembly submatrix, so the row blocks PARPACK needs have to
                // be scattered from rank 0
                sparse::DistMatrixCSR< int_t > * mDistMatrix = nullptr ;

                // one-based copies of the distributed pattern. DistMatrixCSR
                // is zero-based by construction ( distribute_values() forces
                // SpMatrixIndexingBase::Cpp on the source ), and its arrays
                // are handed out const, so the shift lives here
                Vector< int_t > mParpackPointers ;
                Vector< int_t > mParpackIndices ;

                // structural fingerprint of the matrix the pattern was built
                // from. NOT reset() -- that runs after every Jacobian
                // computation, so rebuilding there would redistribute the
                // pattern once per Newton iteration
                const SpMatrix * mDistSource = nullptr ;
                int_t mDistNumRows = 0 ;
                int_t mDistNumNonzeros = 0 ;

#endif

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                EigenValues( DofManager * aParent ) ;

                ~EigenValues();

                // NON-COPYABLE, NON-MOVABLE. This class owns four raw pointers
                // and deletes all four in its destructor -- mSolver,
                // mShiftInvertSolver, mM and ( under PARPACK ) mDistMatrix --
                // so the implicit copy would be a shallow pointer copy and the
                // second destructor a double free. Nothing copies it today
                // ( DofManager holds it by pointer and news it once ), which is
                // exactly why the hazard was latent rather than loud.
                //
                // Same shape and same remedy as SpMatrix ( cl_SpMatrix.hpp )
                // and Solver ( cl_Solver.hpp ). Deleting them makes any future attempt
                // a compile error instead of a run-time double free
                EigenValues( const EigenValues & ) = delete ;
                EigenValues( EigenValues && ) = delete ;
                EigenValues & operator=( const EigenValues & ) = delete ;
                EigenValues & operator=( EigenValues && ) = delete ;

//------------------------------------------------------------------------------

                void
                set_num_minvals( const index_t aNumMinVals );

//------------------------------------------------------------------------------

                void
                set_num_maxvals( const index_t aNumMaxVals );

//------------------------------------------------------------------------------

                /**
                 * relative accuracy demanded of each Ritz value. The small end
                 * of the spectrum is only reachable while this stays above
                 * eps_mach * kappa -- see the note at mEpsilon before lowering
                 * it. Non-positive values are rejected: ARPACK reads <= 0 as
                 * "use machine precision", which is exactly the unreachable
                 * setting on an ill conditioned matrix
                 */
                void
                set_tolerance( const real aTolerance );

                /**
                 * floor for the Krylov subspace. Larger usually means fewer
                 * matrix-vector products but a bigger basis ( n x ncv )
                 */
                void
                set_subspace_size( const index_t aSubspaceSize );

                /**
                 * restart budget handed to ( p )dnaupd. There is no matvec
                 * budget: one restart does NOT cost ncv - nev products
                 * ( ARPACK boosts nev internally between restarts, so np
                 * shrinks; measured 10.1 and 12.2 products per restart at
                 * ncv - nev = 19 ), so this bounds the restart count and the
                 * achieved product count is REPORTED rather than modeled
                 */
                void
                set_max_iterations( const index_t aNumMaxIter );

                /**
                 * declare the matrix symmetric, selecting dsaupd / dseupd over
                 * dnaupd / dneupd. Lanczos is cheaper and its workspace is
                 * smaller, and for a symmetric matrix the spectral ratio IS
                 * kappa_2 -- which the nonsymmetric driver cannot promise.
                 *
                 * This is an ASSERTION BY THE CALLER, not a measurement. A
                 * matrix wrongly declared symmetric gives a wrong answer with
                 * no diagnostic: dsaupd only ever references one triangle's
                 * worth of information through the matvec, so an unsymmetric
                 * part is silently ignored rather than detected
                 */
                void
                set_symmetric( const bool aSymmetric );

                /**
                 * whether the spectral ratio this object returns may be called
                 * kappa_2. That needs symmetry AND a positive definite
                 * spectrum -- real eigenvalues alone are not enough
                 */
                bool
                is_symmetric() const ;

//------------------------------------------------------------------------------

                real
                compute_conditioning() ;

//------------------------------------------------------------------------------

                real
                compute_lambda_max() ;

//------------------------------------------------------------------------------

                /**
                 * compute the eigenvalues at one end of the spectrum and
                 * return the extremal magnitude found there. how many are
                 * computed is set by set_num_minvals / set_num_maxvals ;
                 * the full set is exposed through lambda_real / lambda_imag
                 */
                real
                compute_smallest_eigenvalues();

                real
                compute_largest_eigenvalues();

//------------------------------------------------------------------------------

                /**
                 * real and imaginary parts of the values found by the last
                 * call above. Which ranks hold them depends on the build:
                 *
                 *  - serial, or PARPACK: valid on every rank. pdneupd returns
                 *    the same spectrum everywhere
                 *  - multi-rank without PARPACK: valid on the MASTER ONLY.
                 *    The master computes alone and only the scalar extremum
                 *    is broadcast, so a worker's buffers stay stale
                 *
                 * The same split applies to number_of_converged_values().
                 */
                const Vector< real > &
                lambda_real() const ;

                const Vector< real > &
                lambda_imag() const ;

                /**
                 * how many of the requested values actually converged
                 */
                int_t
                number_of_converged_values() const ;

                /**
                 * "ARPACK" or "PARPACK", whichever drove the last call
                 */
                const string &
                backend_label() const ;

                /**
                 * how the last compute_conditioning() ended. Identical on every
                 * rank -- it is broadcast before it is stored
                 */
                EigenOutcome
                outcome() const ;

//------------------------------------------------------------------------------

                void
                reset();

//------------------------------------------------------------------------------
            private:
//------------------------------------------------------------------------------

                void
                compute_matrices() ;

                /**
                 * size ncv, maxit and tol for one run from ( n, nev, which end ).
                 * aJob picks the product budget: an INTERIOR request ( job 0,
                 * 'SM' unfolded ) is the expensive one and keeps the larger
                 * budget the public compute_smallest_eigenvalues() has always
                 * had. aFolded selects the tight tolerance the cancellation needs.
                 * Returns false when no legal ncv fits mBasisBudget, in which
                 * case the diagnostic is Infeasible and nothing is attempted
                 */
                bool
                configure( const int_t aJob,
                           const bool  aFolded,
                           const int_t aNumEigenValues,
                           int_t     & aSubspaceSize,
                           int_t     & aNumMaxIter,
                           real      & aTolerance ) ;

                /**
                 * aSigma = 0 runs the plain operator, anything else runs the
                 * fold OP = aSigma*I - A. The same value must reach every rank
                 */
                real
                run( const int_t aJob, const real aSigma );

                void
                link_matrix();

                real
                run_arpack( const int_t aJob, const real aSigma );

                real
                run_parpack( const int_t aJob, const real aSigma );

                /**
                 * pick the extremum out of the converged set and remember the
                 * winning slot's signed parts in mExtremalReal / mExtremalImag.
                 * Shared by both drivers so the two cannot drift apart
                 */
                real
                reduce_extremum( const int_t aJob, const int_t aNumConverged );

                /**
                 * store the outcome, count the strike, latch when the end is
                 * proven unreachable, say so once, and hand back a NaN.
                 * aOutcome must already be the SAME on every rank
                 */
                real
                report_unavailable( const int_t aOutcome );

                /**
                 * smallest eigenvalue by SHIFT-INVERT: ARPACK mode 3 with
                 * sigma = 0, so OP = A^-1 and the wanted end of A is the
                 * dominant end of OP.
                 *
                 * This is the only Krylov route to the small end of a
                 * discretized diffusion operator. A spectral FOLD cannot get
                 * there: a shift is affine and preserves ratios of
                 * differences, so it moves the wanted eigenvalue to the
                 * exterior without improving its separation. Measured on the
                 * tapestack3d thermal Jacobian, the small-end relative gap is
                 * 9.5e-10 folded against 2.8e-2 inverted.
                 *
                 * COLLECTIVE. The master drives ARPACK ; every rank enters the
                 * solve for each inverse application, steered by a broadcast
                 * command word. Returns lambda_min, or NaN with aOutcome set
                 */
                real
                run_shift_invert( int_t & aOutcome );

                void
                restore_indexing_base();
//------------------------------------------------------------------------------

            };

//------------------------------------------------------------------------------

            inline const Vector< real > &
            EigenValues::lambda_real() const
            {
                return mLambdaReal ;
            }

//------------------------------------------------------------------------------

            inline const Vector< real > &
            EigenValues::lambda_imag() const
            {
                return mLambdaImag ;
            }

//------------------------------------------------------------------------------

            inline int_t
            EigenValues::number_of_converged_values() const
            {
                return mInfo( arpack::gInfoNumConverged ) ;
            }

//------------------------------------------------------------------------------

            inline const string &
            EigenValues::backend_label() const
            {
                return mBackendLabel ;
            }

//------------------------------------------------------------------------------

            inline bool
            EigenValues::is_symmetric() const
            {
                return mSymmetric ;
            }

//------------------------------------------------------------------------------

            inline EigenOutcome
            EigenValues::outcome() const
            {
                return mOutcome ;
            }

//------------------------------------------------------------------------------
        }
    }
}
#endif //CL_FEM_DOFMGR_EIGENVALUES_HPP
