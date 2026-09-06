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

#ifndef BELFEM_PETSCTOOLS_HPP
#define BELFEM_PETSCTOOLS_HPP


#ifdef BELFEM_PETSC
#include <petscksp.h>
#else

typedef int                 PetscInt;
typedef double              PetscReal ;
typedef int                 PetscErrorCode ;
// vector class
class Vec { public: Vec() = default ; ~Vec() = default ; };

// matrix class
class Mat { public: Mat() = default ; ~Mat() = default ; };

// Preconditioner context
class PC { public: PC() = default ; ~PC() = default ; };

// solver class
class KSP { public: KSP() = default ; ~KSP() = default ; };

#define PETSC_DECIDE 0
#endif


#ifndef BELFEM_MPI
typedef int                 MPI_Comm ;
#endif

#include "typedefs.hpp"

#include "cl_Vector.hpp"
#include "cl_SolverDistMatrix.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    namespace sparse
    {
#ifdef BELFEM_PETSC
        typedef DistMatrixAIJ< PetscInt > PETScAIJ ;
#else
        typedef void PETScAIJ ;
#endif
    }
//------------------------------------------------------------------------------

    string
    petsctools_error_message( const PetscErrorCode aStatus );

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_allocate_vector(
            MPI_Comm         aComm,
            Vec            & aVec,
            const PetscInt   aGlobalLength,
            const PetscInt   aLocalLength=PETSC_DECIDE
    );

//------------------------------------------------------------------------------

    /**
     * convert a BELFEM vector to a PETSC vector
     * aVac must have been initialized.
     * @param aVector
     * @param aIndices
     * @param aVec
     */
    PetscErrorCode
    petsctools_set_vector(
            const Vector <PetscReal> & aVector,
            const Vector <PetscInt> & aIndices,
            Vec & aVec );

//------------------------------------------------------------------------------

    /**
     * convert a PETSC vector to a BELFEM vector
     *
     * aVector must have been initialized
     * @param aVec
     * @param aIndices
     * @param aVector
     *
     */
    PetscErrorCode
    petsctools_get_vector(
            Vec & aVec,
            const Vector <PetscInt> & aIndices,
            Vector< PetscReal > & aVector );

//------------------------------------------------------------------------------

    /**
     * Extract the local portion of a distributed PETSC vector to a BELFEM vector
     *
     * aVector must have been initialized with the correct local size
     * @param aVec
     * @param aVector
     *
     */
    PetscErrorCode
    petsctools_get_local_vector(
            Vec & aVec,
            Vector< PetscReal > & aVector );

//------------------------------------------------------------------------------

    /**
     * Set the local portion of a distributed PETSC vector from a BELFEM vector
     *
     * aVector must contain the local values to set
     * @param aVector
     * @param aVec
     *
     */
    PetscErrorCode
    petsctools_set_local_vector(
            const Vector< PetscReal > & aVector,
            Vec & aVec );

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_create_matrix(
        sparse::PETScAIJ * aDistMatrix,
        Mat & aMat );

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_update_matrix(
        sparse::PETScAIJ * aDistMatrix,
        Mat & aMat );

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_link_matrix(
        COMM_TYPE   aComm,
        SpMatrix & aMatrix,
        Mat      & aMat );

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_notify_matrix_update( Mat & aMat );

//------------------------------------------------------------------------------

}
#endif //BELFEM_PETSCTOOLS_HPP
