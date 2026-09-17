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

#ifndef BELFEM_CL_SOLVER_HPP
#define BELFEM_CL_SOLVER_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SolverParameters.hpp"
#include "cl_SpMatrix.hpp"

#include "en_SolverEnums.hpp"
#include "cl_SolverWrapper.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Unified interface to the sparse direct solvers.
     *
     * @ingroup grp_sparse
     * @see @ref sparse_sparse_usage_guide
     */
    class Solver
    {
        // the type of the chosen solver
        const SolverType mType ;

        SolverParameters mParams ;

        // symmetry mode
        SymmetryMode mSymmetryMode = SymmetryMode::Unsymmetric ;

        solver::Wrapper * mWrapper ;

//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        /** Collective: the constructor synchronizes the solver settings ( not
         *  the solver type ) from rank 0 to every rank, so all ranks must
         *  construct the Solver together. Both constructors share this
         *  contract. */
        Solver( const SolverType aSolverType = gDefaultSolver ) ;

        Solver( const SolverParameters aParams ) ;

//------------------------------------------------------------------------------

        ~Solver() ;

        // NON-COPYABLE, NON-MOVABLE. This class owns the raw pointer
        // mWrapper and deletes it in its destructor, so the implicit copy
        // would be a shallow pointer copy and the second destructor a
        // double free. Copy assignment was already implicitly deleted by
        // the const mType member, and the user-declared destructor
        // suppressed the implicit moves -- which is why a move request
        // used to bind to the implicit COPY constructor, the same hazard
        // under another name. Nothing copies a Solver, which is exactly
        // why this was latent rather than loud.
        //
        // Same shape and same remedy as EigenValues
        // ( cl_FEM_DofMgr_EigenValues.hpp ). NOT the SpMatrix pattern:
        // that one deletes the constructors but keeps assignment on
        // purpose, because an empty SpMatrix is a usable target
        Solver( const Solver & ) = delete ;
        Solver( Solver && ) = delete ;
        Solver & operator=( const Solver & ) = delete ;
        Solver & operator=( Solver && ) = delete ;

//------------------------------------------------------------------------------

         SolverType
         type() const ;

//------------------------------------------------------------------------------

        /**
         * read-only view of the parameters this solver was built with
         * ( e.g. for setup-time consistency checks against the nonlinear
         * tolerances )
         */
         const SolverParameters &
         parameters() const ;

//------------------------------------------------------------------------------

        /**
         * symmetry mode, needed for some solvers
         * @param aMode
         */
        void
        set_symmetry_mode( const SymmetryMode & aMode ) ;

//------------------------------------------------------------------------------

        /**
         * Solves A x = b. It is collective for the MPI backends ( MUMPS,
         * STRUMPACK, PETSc ) and local for the others. The first call
         * initializes the wrapper. It does not return a failure status: when
         * the wrapper's soft-fail contract is armed on a backend that has one
         * ( MUMPS, STRUMPACK, PETSc ), a failed factorization or solve is
         * reported by wrapper()->failed(); every other failure aborts
         * ( BELFEM_ERROR ).
         */
        void
        solve(  SpMatrix       & aMatrix,
                Vector< real > & aLHS,
                Vector< real > & aRHS );

//------------------------------------------------------------------------------

        /**
         * Same contract as the vector overload, with one right-hand side per
         * column of aRHS; the backend resizes aLHS to that shape. The first
         * solve initializes the wrapper with aRHS's column count. PETSc and
         * STRUMPACK have no matrix overload and abort here.
         */
        void
        solve(  SpMatrix       & aMatrix,
                Matrix< real > & aLHS,
                Matrix< real > & aRHS );

//------------------------------------------------------------------------------

        /**
         * this does only do something if PETSc is used
         *
         * @param aPreconditioner
         * @param aKrylovMethod
         * @param aEpsilon
         */
        void
        set_petsc(
                const Preconditioner aPreconditioner,
                const KrylovMethod   aKrylovMethod,
                const real           aEpsilon = 1e-8 );

//------------------------------------------------------------------------------

        /**
         * Silently ignored unless the solver is MUMPS. The three setters below
         * share this contract.
         */
        void
        set_mumps_reordering(
                const MumpsSerialReodrdering   aSerial,
                const MumpsParallelReodrdering aParallel );

        void
        set_mumps_blr(
                const MumpsBlockLowRanking aBlr,
                const real            aEpsilon );

        void
        set_mumps_error_analysis(
                const MumpsErrorAnalysis aSetting );

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

        /**
         * Releases the factorization ( the destructor does it too ); the next
         * solve() re-initializes from scratch, including the RHS column count.
         */
         void
         free() ;

//------------------------------------------------------------------------------

        /**
         * Borrowed: the Solver owns and deletes the wrapper; do not keep it past
         * the Solver's lifetime.
         */
        solver::Wrapper *
        wrapper() ;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        create_wrapper();

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------

    inline solver::Wrapper *
    Solver::wrapper()
    {
        return mWrapper ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_SOLVER_HPP
