//
// Created by Christian Messe on 10.07.20.
//

#ifndef BELFEM_CL_SOLVER_HPP
#define BELFEM_CL_SOLVER_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SpMatrix.hpp"

#include "en_SolverEnums.hpp"
#include "cl_SolverWrapper.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    class Solver
    {
        // the type of the chosen solver
        const SolverType mType ;

        // symmetry mode
        SymmetryMode mSymmetryMode = SymmetryMode::Unsymmetric ;

#ifdef BELFEM_PETSC
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // PETSc
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // the preconditioner ( PETSc only )
        Preconditioner mPreconditioner = Preconditioner::UNDEFINED ;

        // the Krylov subspace method ( PETSc only )
        KrylovMethod mKrylovMethod = KrylovMethod::UNDEFINED ;

        // epsilon environment for iterative solver
        real mEpsilon = 1e-8 ;
#endif
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        solver::Wrapper * mWrapper ;

//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        Solver( const SolverType aSolverType = gDefaultSolver,
                const proc_t aMasterRank = 0 ) ;

//------------------------------------------------------------------------------

        ~Solver() ;

//------------------------------------------------------------------------------

        /**
         * return the type of the solver
         */
         SolverType
         type() const ;

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
                const SerialReodrdering   aSerial,
                const ParallelReodrdering aParallel );

        void
        set_mumps_blr(
                const BlockLowRanking aBlr,
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
