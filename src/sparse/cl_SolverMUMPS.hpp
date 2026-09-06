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

#ifndef BELFEM_CL_SOLVERMUMPS_HPP
#define BELFEM_CL_SOLVERMUMPS_HPP

#include "cl_Cell.hpp"
#include "cl_SolverWrapper.hpp"
#include "cl_SolverParameters.hpp"

namespace belfem
{
    namespace solver
    {

        class MUMPS : public Wrapper
        {
            // parameter object
            const SolverParameters * mParams ;

            // Rank of HOST
            const proc_t mMasterRank ;

            // id of solver
            int_t mSolverID = 0 ;

            // vector containing user parameters
            Vector< int_t >  mIParameters ;
            Vector< real >   mRParameters ;

            // vector containing debug information ( rank-local INFO: a
            // failing rank holds the true code, the others the propagated -1;
            // still read by check_warnings(), which is rank-local by design )
            Vector< int_t >  mInfo ;

            // rank-uniform INFOG: same code and supplementary value on every
            // rank, so the -9 retry ladder and the error arms key on it --
            // a branch taken from mInfo would strand ranks outside the
            // collective call. Vector< int_t > to match its sibling mInfo,
            // not a linear-algebra object
            Vector< int_t >  mInfoG ;

            Vector< real > mRInfoG ;

            // pointer to matrix
            SpMatrix * mMatrix = nullptr ;

            //! frozen-factorization scope. While armed, a solve runs JOB 3
            //! ( solve only ) against the factors the last successful JOB 5/6
            //! built, instead of rebuilding them.
            //!
            //! This is its OWN struct rather than extra loose members on
            //! purpose: MUMPS re-declares mMatrix above, shadowing
            //! Wrapper::mMatrix, and hanging the record off either of those
            //! two would make it ambiguous which one it describes.
            //!
            //! Only the master can fill mIdentity -- on a worker the matrix
            //! argument is a local submatrix that the distributed solve never
            //! reads, so a worker records mArmed alone
            struct FrozenFactorization
            {
                bool             mArmed    = false ;

                //! true only where mIdentity below is meaningful
                bool             mHaveIdentity = false ;

                const SpMatrix * mSource   = nullptr ;
                index_t          mNumRows  = 0 ;
                index_t          mNumCols  = 0 ;
                index_t          mNumNonZeros = 0 ;
                const real     * mValues   = nullptr ;
                const int_t    * mPointers = nullptr ;
            };

            FrozenFactorization mFrozen ;

            //! per-process budget in MB the machine afforded this instance,
            //! measured ONCE in initialize() -- before MUMPS has allocated
            //! anything, so the number is the total the instance may take,
            //! which is what ICNTL(23) bounds. 0 = could not be measured on
            //! some rank. Rank-uniform ( MIN-reduced ); read by
            //! escalate_workspace() and print_soft_fail()
            int_t mMeasuredBudgetMB = 0 ;

            //! a Fortran instance exists in the mumpstools registry that
            //! free() owes a JOB = -2. Set by a SUCCESSFUL initialize(),
            //! cleared by free(), and read in exactly one place -- free()'s
            //! teardown gate. It is NOT the wrapper's usability flag: that is
            //! Wrapper::mIsInitialized, reached through is_initialized(), and
            //! the two are separate members with separate meanings
            bool mInitialized = false ;

            //! UNSYMMETRIC. It used to default to GeneralSymmetric, which is
            //! MUMPS SYM = 2 -- see the guard in initialize() for why that is
            //! not something BELFEM can honour
            SymmetryMode mSymmetryMode = SymmetryMode::Unsymmetric ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            MUMPS( const SolverParameters * aParams, const proc_t aMasterRank = 0 );

//------------------------------------------------------------------------------

            ~MUMPS() override ;

//------------------------------------------------------------------------------

            //! MUMPS keeps its factors in the saved DMUMPS_STRUC until JOB -2,
            //! so solving against them is a matter of asking for JOB 3
            bool
            supports_factorization_reuse() const override ;

            void
            freeze_factorization( const SpMatrix & aMatrix ) override ;

            void
            unfreeze_factorization() override ;

            bool
            factorization_is_frozen() const override ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * pick the MUMPS JOB for a solve against aMatrix: 3 while frozen,
             * else 6 for a new matrix and 5 for a repeat. Shared by BOTH
             * solve overloads so their policies cannot drift apart.
             * Validates the frozen record with always-active errors
             */
            int_t
            select_job( SpMatrix & aMatrix );

            /**
             * drop the frozen scope because the factors are no longer
             * trustworthy. Called from both soft-failure paths and free()
             */
            void
            invalidate_factorization();

//------------------------------------------------------------------------------

            /**
             * the workspace-failure retry's single decision point, shared by
             * BOTH solve overloads so their policies cannot drift apart
             * ( same argument as select_job ). Reads the rank-uniform mInfoG
             * and applies mumps::next_workspace_action(): first the machine
             * budget ( mMeasuredBudgetMB ) as ICNTL(23), then the ICNTL(14)
             * ladder ( bounded by gMaxMemoryRelaxation ), give-up on -19.
             * Returns true when the caller should re-enter the collective
             * solve, with aJob rewritten to 5. Every input is rank-uniform,
             * so every rank returns the same answer; nothing collective
             * happens inside
             */
            bool
            escalate_workspace( int_t & aJob );

            /**
             * per-process memory budget in MB the machine affords this
             * instance: available memory / ranks on this node / live MUMPS
             * instances * safety, reduced to the MIN over all ranks so every
             * rank writes the same ICNTL(23). 0 when any rank could not
             * measure; never 0 for a machine that WAS measured, however full
             * ( 1 MB floor ), so "over the limit" stays distinguishable from
             * "unknown". COLLECTIVE; called from initialize()
             */
            int_t
            memory_budget_mb();

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

            void
            solve(
                    SpMatrix & aMatrix,
                    Vector< real > & aLHS,
                    Vector< real > & aRHS ) override ;

//------------------------------------------------------------------------------

            void
            solve(
                    SpMatrix & aMatrix,
                    Matrix< real > & aLHS,
                    Matrix< real > & aRHS ) override ;

//------------------------------------------------------------------------------

            void
            free() override ;

//------------------------------------------------------------------------------

            void
            set_reordering(
                    const MumpsSerialReodrdering   aSerial,
                    const MumpsParallelReodrdering aParallel
            );

//------------------------------------------------------------------------------

            void
            set_block_low_ranking(
                    const MumpsBlockLowRanking aBLK,
                    const real            aEpsilon
            );

//------------------------------------------------------------------------------

            void
            set_error_analysis(
                    const MumpsErrorAnalysis aSetting
            );

//------------------------------------------------------------------------------

            /**
             * returns the determinant, if supported by the solver
             * and computation was requested
             */
            real
            get_determinant() const override ;

//------------------------------------------------------------------------------

            /**
            * returns the conditioning number, if supported by the solver
            * and computation was requested
            */
            real
            get_cond1() const override ;

//------------------------------------------------------------------------------

            /**
             * returns the conditioning number, if supported by the solver
             * and computation was requested
            */
            real
            get_cond2() const override ;

//------------------------------------------------------------------------------

            /**
             * RINFOG( 9 ): estimated forward error of the last solve
             */
            real
            get_forward_error() const override ;

//------------------------------------------------------------------------------

            /**
             * RINFOG( 7 ) + RINFOG( 8 ): the componentwise backward error.
             * MUMPS reports it in two halves, split by whether a row's
             * denominator is trustworthy, and they sum
             */
            real
            get_backward_error() const override ;

//------------------------------------------------------------------------------

            /**
             * RINFOG( 8 ) alone: omega2. Zero means COND2's term drops out of
             * the forward error. It does NOT prove MUMPS skipped the COND2
             * estimator -- see the implementation for why the two differ
             */
            real
            get_omega2() const override ;

//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            initialize(
                    SpMatrix & aMatrix,
                    const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                    const int_t aNumRhsColumns=1 ) override ;


//------------------------------------------------------------------------------

            string
            error_message(
                    const int_t   * aInfo,
                    const int_t   & aN,
                    const int_t   & aNNZ ) ;

//------------------------------------------------------------------------------

            // decodes a positive INFO(1) warning ( a usable solution is still
            // returned ) into one human-readable line per active warning bit
            void
            warning_message(
                    const int_t    * aInfo,
                    Cell< string > & aWarnings ) ;

//------------------------------------------------------------------------------
        private :
//------------------------------------------------------------------------------

            void
            init_defaults();

//------------------------------------------------------------------------------

            // handles a positive INFO(1) after a solve: raises a hard error on
            // the out-of-range-index bit ( +1 ) and reports the remaining
            // warnings per rank. NOT used by free(), which is teardown
            void
            check_warnings();

//------------------------------------------------------------------------------

            // rank-0 log line for the soft-fail path, box-styled to match the
            // controller's timestep box, with a compact INFO decode ( the full
            // decode lives in error_message() on the hard path )
            void
            print_soft_fail();

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_SOLVERMUMPS_HPP
