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

        /**
         * return the type of the solver
         */
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
         * solves the system
         * @param aMatrix
         * @param aLHS     Left Hand Side
         * @param aRHS     Right Hand Side
         */
        void
        solve(  SpMatrix       & aMatrix,
                Vector< real > & aLHS,
                Vector< real > & aRHS );

//------------------------------------------------------------------------------

        /**
         * solves the system
         * @param aMatrix
         * @param aLHS     Left Hand Side
         * @param aRHS     Right Hand Side
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
         * mumps only
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

        /**
         * this does only do something if PARDISO is used
         *
         */
/*        void
        set_pardiso(
                const PardisoMode aMode ) ; */

//------------------------------------------------------------------------------

        /**
         * tidy up solver manually ( also done by destructor )
         */
         void
         free() ;

//------------------------------------------------------------------------------

        /**
         * expose wrapper object
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
