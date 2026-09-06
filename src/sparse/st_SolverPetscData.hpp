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
