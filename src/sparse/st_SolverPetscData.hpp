//
// Created by Christian Messe on 13.07.20.
//

#ifndef BELFEM_ST_SOLVERPETSCDATA_HPP
#define BELFEM_ST_SOLVERPETSCDATA_HPP
#include "petsctools.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        /**
         * this struct contains the variables needed by PETSc
         *
         * variables are initialized and destroyed by the PETSc wrapper
         */
        struct PetscData
        {
            // communicator that is used
            MPI_Comm mComm;

            // index arrays for solver
            Vector <PetscInt> mMyPointers ;
            Vector <PetscInt> mMyColumns ;
            Vector <PetscInt> mMyRows ;

            // indices for vectors
            Vector <PetscInt> mVectorIndices ;

            // global number of colums for RHS and LHS
            PetscInt mNumRows ;
            PetscInt mNumCols ;

            // proc local number of colums for RHS and LHS
            PetscInt mMyNumRows ;
            PetscInt mMyNumCols;

            // row offset for indices
            PetscInt mMyRowOffset ;

            // number of nonzeros on this proc
            PetscInt mMyNumNnz ;

            // vector for left hand side
            Vec mLHS ;

            // vector for right hand side
            Vec mRHS ;

            // matrix wrapper
            Mat mMat;

            // solver object
            KSP mKSP;

            // preconditioner context
            PC mPC;
        };

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_ST_SOLVERPETSCDATA_HPP
