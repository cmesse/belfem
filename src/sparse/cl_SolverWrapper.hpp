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

#ifndef BELFEM_CL_SOLVERWRAPPER_HPP
#define BELFEM_CL_SOLVERWRAPPER_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SpMatrix.hpp"
#include "en_SolverEnums.hpp"
#include "fn_create_graph_from_matrix.hpp"
#include "cl_Communicator.hpp"

namespace belfem
{
    namespace solver
    {
        /**
         * parent class for solver specific data
         */
        class Wrapper : public CommunicationObject
        {
//------------------------------------------------------------------------------

            // rank of this proc
            const proc_t mCommRank ;

            const proc_t mCommSize ;

            const bool mUsesMPI ;

            // label of this solver type
            const string mLabel ;

            // flag telling if we have been initialized
            bool mIsInitialized = false ;

            //! soft-fail contract: when armed, a failed factorization/solve
            //! is RECORDED instead of raising BELFEM_ERROR, and the caller
            //! ( the nonlinear Controller ) treats the event as a failed
            //! trial and cuts the timestep. Default off — standalone users
            //! keep the loud abort. Rank-uniformity: MUMPS propagates errors
            //! to every rank ( non-failing ranks see INFO(1) = -1 ), and the
            //! STRUMPACK collective solve returns the same ReturnCode on all
            //! ranks, so the recorded failure is uniform.
            bool mSoftFail = false ;
            bool mFailed   = false ;

        protected:

            //! subclass hooks for the soft-fail contract
            bool
            soft_fail() const
            {
                return mSoftFail ;
            }

            void
            flag_failure()
            {
                mFailed = true ;
            }

            SpMatrix * mMatrix = nullptr ;
            Vector< real > * mX = nullptr ;
            Vector< real > * mY = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Wrapper( const string & aLabel, const bool aUsesMPI );

//------------------------------------------------------------------------------

            virtual ~Wrapper();

//------------------------------------------------------------------------------

            /**
             * returns the name of the solver as string
             * @return
             */
            const string &
            label() const ;

//------------------------------------------------------------------------------

            /**
             * arm/disarm the soft-fail contract ( see the member note )
             */
            void
            set_soft_fail( const bool aFlag )
            {
                mSoftFail = aFlag ;
            }

            /**
             * true if the last factorization/solve failed softly
             */
            bool
            failed() const
            {
                return mFailed ;
            }

            /**
             * clear a recorded soft failure before the next solve
             */
            void
            clear_failure()
            {
                mFailed = false ;
            }

//------------------------------------------------------------------------------

            virtual void
            solve( SpMatrix & aMatrix,
                   Vector <real> & aLHS,
                   Vector <real> & aRHS );

//------------------------------------------------------------------------------

            virtual void
            solve( SpMatrix & aMatrix,
                   Matrix <real> & aLHS,
                   Matrix <real> & aRHS );

//------------------------------------------------------------------------------

            /**
             * tells if the initialitzation routine has already been called
             */
            bool
            is_initialized() const ;

//------------------------------------------------------------------------------

            /**
             * returns the communication rank of this proc
             */
            proc_t
            rank() const ;

//------------------------------------------------------------------------------

            /**
             * returns the communication size of this proc
             */
            proc_t
            comm_size() const ;

//------------------------------------------------------------------------------

            bool
            uses_mpi() const ;

//------------------------------------------------------------------------------

        virtual void
        initialize(
                SpMatrix & aMatrix,
                // UNSYMMETRIC, and it must stay identical to the MUMPS override's
                // default: default arguments bind STATICALLY, so a base/derived
                // mismatch would silently change the mode with the static type
                // of the pointer the call goes through
                const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                const int_t aNumRhsColumns = 1 );

//------------------------------------------------------------------------------

        void
        free() override;

//------------------------------------------------------------------------------

        /**
         * returns the determinant, if supported by the solver
         * and computation was requested
         */
        virtual real
        get_determinant() const ;

//------------------------------------------------------------------------------

        /**
        * returns the conditioning numbner, if supported by the solver
        * and computation was requested
        */
        virtual real
        get_cond1() const ;

//------------------------------------------------------------------------------

        /**
         * returns the conditioning numbner, if supported by the solver
         * and computation was requested
        */
       virtual real
       get_cond2() const ;

//------------------------------------------------------------------------------

        /**
         * estimated FORWARD error of the last solve, || dx || / || x ||, if
         * supported by the solver and computation was requested.
         *
         * Distinct from the condition numbers above: those bound the worst
         * case the matrix admits, this reports how many digits the solve
         * actually lost. That makes it the cheaper thing to steer on -- a
         * timestepper can read it directly -- while the condition number is
         * the more descriptive figure to show a user
         */
        virtual real
        get_forward_error() const ;

//------------------------------------------------------------------------------

        /**
         * componentwise BACKWARD error of the last solve, if supported by the
         * solver and computation was requested: the smallest relative
         * perturbation of A and b for which the computed x is exact.
         * Near machine precision means the solver did its job and any lost
         * accuracy is the matrix's doing, not the factorization's
         */
        virtual real
        get_backward_error() const ;

//------------------------------------------------------------------------------

        /**
         * the SECOND half of that backward error on its own, if supported.
         * Zero means the term it weights contributes nothing to the forward
         * error, so the matching second condition number is inert for this
         * solve. Deliberately stated in terms of the ERROR BOUND and not of
         * any one solver's internals -- what makes a row of the split empty
         * is the implementation's business, and MUMPS documents its own rule
         * on the override
         */
        virtual real
        get_omega2() const ;

//------------------------------------------------------------------------------

        /**
         * true if this wrapper can solve against an EXISTING factorization
         * rather than rebuilding it. Default false: a wrapper that has not
         * implemented the frozen path must never be asked to take it
         */
        virtual bool
        supports_factorization_reuse() const ;

//------------------------------------------------------------------------------

        /**
         * arm the frozen state: subsequent solves reuse the factorization
         * that the LAST SUCCESSFUL solve built, instead of rebuilding it.
         *
         * This is a scope, not a hint, and the caller owns the whole of it:
         *
         *   solve( A, x, b )          <- builds the factorization
         *   freeze_factorization( A ) <- arm, only after that succeeded
         *   solve( A, y, c ) ...      <- reuse, values of A must NOT change
         *   unfreeze_factorization()  <- always, including on every failure
         *
         * Why a scope rather than a flag: if the values of A change while
         * frozen, the solver returns a valid-LOOKING solution of the wrong
         * matrix, with no error anywhere. Pointer identity cannot detect that
         * -- the wrapper already treats the same pointer as "same structure,
         * different values" and deliberately re-factorizes for exactly this
         * reason. So the contract is that NO assembly may happen inside the
         * scope, the caller guarantees it, and the wrapper checks everything
         * about the matrix it cheaply can ( structure and identity, not an
         * O( nnz ) values scan on every solve ).
         *
         * Checks are ALWAYS-ACTIVE errors rather than assertions: a silent
         * stale-factor solve is precisely the failure that must not survive
         * into a release build.
         *
         * Rank behaviour: only the rank that owns the matrix can record its
         * identity -- on every other rank the matrix argument is a local
         * submatrix or empty, and the distributed solvers never read it. So a
         * worker records only that it is frozen. Uniformity comes from the
         * caller arming and dropping the scope at the same point of the same
         * code path on every rank, never from a rank-local decision.
         *
         * The default implementation throws: a wrapper that does not support
         * reuse must be asked via supports_factorization_reuse() first
         */
        virtual void
        freeze_factorization( const SpMatrix & aMatrix );

//------------------------------------------------------------------------------

        /**
         * drop the frozen state and return to ordinary re-factorizing solves.
         * Safe to call when not frozen, so a scope guard can call it
         * unconditionally on the way out
         */
        virtual void
        unfreeze_factorization();

//------------------------------------------------------------------------------

        /**
         * true while the frozen scope is armed
         */
        virtual bool
        factorization_is_frozen() const ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            virtual void
            initialize();

//------------------------------------------------------------------------------

            /**
             * print a warning ( the turtle ) when the MPI ranks sharing a node ask
             * for more OpenMP threads than the node has PHYSICAL cores. Each rank
             * spawns omp_get_max_threads() BLAS/solver threads, and the runtimes
             * cannot see each other: with OMP_NUM_THREADS unset each one defaults
             * to its affinity mask, which counts every hyperthread.
             *
             * The cost is time and memory, not a hang: measured on this code,
             * 4 ranks x 4 threads on 10 physical cores ran the factorization ~15 %
             * slower than 4 x 2 and used ~13 GiB more ( doc/parallel_execution.md ).
             * BELFEM's element loop carries no omp pragmas, so the extra threads
             * cannot help the assembly at all -- only STRUMPACK and MKL consume them.
             *
             * Self-checking: a no-op unless oversubscribed, a no-op when the
             * physical core count cannot be established ( an unknown budget must
             * never become a recommendation ), and a no-op entirely when built
             * without OpenMP.
             *
             * The budget comes from THIS rank's affinity mask, which no rank can
             * compare against its siblings' without a node-local communicator.
             * The two configurations separate arithmetically instead: disjoint
             * per-rank slices must fit inside the machine, a shared mask cannot,
             * so `cores_in_my_mask * ranks_on_node <= cores_on_machine` classifies
             * them. Ambiguity is resolved toward NOT warning -- a false alarm on a
             * correctly bound run teaches the reader to ignore the box.
             *
             * LIMITATIONS: a shared cgroup smaller than the machine can look like
             * a set of slices and be under-warned; `cpu_set_t` is fixed at
             * CPU_SETSIZE, so a mask wider than that yields no budget at all. In
             * that last case the warning still fires when OMP_NUM_THREADS was
             * never set and several ranks share the node -- it just names no
             * number, because a budget it could not measure must never become a
             * recommendation.
             */
            void
            hatch_turtle();

//------------------------------------------------------------------------------

            void
            mat2vec( const Matrix< real > & aM,
                           Vector< real > & aV );

//------------------------------------------------------------------------------

            void
            vec2mat( const Vector< real > & aV,
                           Matrix< real > & aM );


//------------------------------------------------------------------------------

            SpMatrix *
            matrix() ;

            Vector< real > &
            x() ;

            Vector< real > &
            y() ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline bool
        Wrapper::is_initialized() const
        {
            return mIsInitialized ;
        }

//------------------------------------------------------------------------------

        inline proc_t
        Wrapper::rank() const
        {
            return mCommRank ;
        }

//------------------------------------------------------------------------------

        inline proc_t
        Wrapper::comm_size() const
        {
            return mCommSize ;
        }

//------------------------------------------------------------------------------

        inline bool Wrapper::uses_mpi() const
        {
            return mUsesMPI ;
        }

//------------------------------------------------------------------------------

        inline SpMatrix *
        Wrapper::matrix()
        {
            return mMatrix ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Wrapper::x()
        {
            return *mX ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Wrapper::y()
        {
            return *mY ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_SOLVERWRAPPER_HPP
