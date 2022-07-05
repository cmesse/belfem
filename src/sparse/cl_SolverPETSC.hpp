//
// Created by Christian Messe on 11.07.20.
//

#ifndef BELFEM_CL_SOLVERPETSC_HPP
#define BELFEM_CL_SOLVERPETSC_HPP

#include "petsctools.hpp"
#include "cl_SolverWrapper.hpp"
#include "st_SolverPetscData.hpp"
#include "cl_SolverPetscDistributor.hpp"

namespace belfem
{

    namespace solver
    {
        class PETSC : public Wrapper
        {
            // selected preconditioner
            Preconditioner mPreconditioner = Preconditioner::ASM;

            // selected krylov method
            KrylovMethod mKrylovMethod = KrylovMethod::GMRES;

            // convergence criterion
            real mEpsilon = 1e-8;

#ifdef BELFEM_PETSC
            // number of colums for RHS and LHS
            PetscInt mNumCols;

            PetscData mData ;

            // distributor. Needed only in parallel mode
            PetscDistributor * mDistributor = nullptr ;
#endif
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            PETSC();

//------------------------------------------------------------------------------

            ~PETSC();

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
                    Vector <real> & aRHS );

//------------------------------------------------------------------------------

            void
            free();

//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            initialize(
                    SpMatrix & aMatrix,
                    const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                    const int aNumRhsColumns = 1 );


//------------------------------------------------------------------------------
        private :
//------------------------------------------------------------------------------

            /*
             * populates the index vector
             */
            void
            create_indices( const PetscInt & aLength );

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
        };

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_SOLVERPETSC_HPP
