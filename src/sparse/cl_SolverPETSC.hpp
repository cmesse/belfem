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

#ifndef BELFEM_CL_SOLVERPETSC_HPP
#define BELFEM_CL_SOLVERPETSC_HPP

#include "petsctools.hpp"
#include "cl_SolverWrapper.hpp"
#include "st_SolverPetscData.hpp"
#include "cl_SolverParameters.hpp"

namespace belfem
{

    namespace solver
    {
        class PETSC : public Wrapper
        {
            const SolverParameters * mParams ;

            // selected preconditioner
            Preconditioner mPreconditioner =mParams->preconditioner() ;

            // selected krylov method
            KrylovMethod mKrylovMethod = mParams->krylov_method() ;

            // relative tolerance for iterative solver
            real mEpsilon = mParams->relative_tolerance() ;

#ifdef BELFEM_PETSC
            // number of colums for RHS and LHS
            PetscInt mNumCols;

            PetscData mData ;

            sparse::PETScAIJ * mDistMatrix = nullptr ;

#endif
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            PETSC( const SolverParameters * aParams );

//------------------------------------------------------------------------------

            ~PETSC() override;

//------------------------------------------------------------------------------

            void
            set(
                    const Preconditioner aPreconditioner,
                    const KrylovMethod aKrylovMethod,
                    const real aEpsilon = 1e-8 );

//------------------------------------------------------------------------------

            void
            solve(
                    SpMatrix      & aMatrix,
                    Vector <real> & aLHS,
                    Vector <real> & aRHS ) override ;

//------------------------------------------------------------------------------

            void
            free() override ;

//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            initialize(
                    SpMatrix & aMatrix,
                    const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                    const int_t aNumRhsColumns = 1 ) override ;


//------------------------------------------------------------------------------
        private :
//------------------------------------------------------------------------------

            /*
             * populates the index vector
             */
            void
            create_indices( const PetscInt & aLength, const PetscInt aOffset = 0 );

//------------------------------------------------------------------------------

            /**
             * create the preconditioner context and the solver object
             */
            void
            create_pc_and_ksp();

//------------------------------------------------------------------------------

            /**
             * link the member matrix to a sparse matrix
             * values are only linked but not comped
             * @param aMatrix
             */
            void
            link_matrix( SpMatrix & aMatrix );

//------------------------------------------------------------------------------

            PetscErrorCode
            set_preconditioner( const Preconditioner aPreconditioner );

//------------------------------------------------------------------------------

            PetscErrorCode
            set_krylovmethod( const KrylovMethod aKrylovMethod );

//------------------------------------------------------------------------------

            PetscErrorCode
            set_initial_guess_flag( const bool aSwitch );

//------------------------------------------------------------------------------

            PetscErrorCode
            set_matrix_ordering( const ReorderingMethod aReorderingMethod );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_SOLVERPETSC_HPP
