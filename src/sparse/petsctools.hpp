//
// Created by Christian Messe on 12.07.20.
//

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
#include "en_SolverEnums.hpp"

namespace belfem
{

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
}
#endif //BELFEM_PETSCTOOLS_HPP
