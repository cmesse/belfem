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

#include "petsctools.hpp"

#include "assert.hpp"

namespace belfem
{
    string
    petsctools_error_message( const PetscErrorCode aStatus )
    {
#ifdef BELFEM_PETSC
        switch( aStatus )
        {
            case( PETSC_ERR_MIN_VALUE ) :
                return "should always be one less then the smallest value";
            case( PETSC_ERR_MEM ) :
                return "unable to allocate requested memory";
            case( PETSC_ERR_SUP ) :
                return "no support for requested operation";
            case( PETSC_ERR_SUP_SYS ) :
                return "no support for requested operation on this computer system";
            case( PETSC_ERR_ORDER ) :
                return "operation done in wrong order";
            case( PETSC_ERR_SIG ) :
                return "signal received";
            case( PETSC_ERR_FP ) :
                return "floating point exception";
            case( PETSC_ERR_COR ) :
                return "corrupted PETSc object";
            case( PETSC_ERR_LIB ) :
                return "error in library called by PETSc";
            case( PETSC_ERR_PLIB ) :
                return "PETSc library generated inconsistent data";
            case( PETSC_ERR_MEMC ) :
                return "memory corruption";
            case( PETSC_ERR_CONV_FAILED ) :
                return "iterative method (KSP or SNES) failed";
            case( PETSC_ERR_USER ) :
                return "user has not provided needed function";
            case( PETSC_ERR_SYS ) :
                return "error in system call";
            case( PETSC_ERR_POINTER ) :
                return "pointer does not point to valid address";
            case( PETSC_ERR_MPI_LIB_INCOMP ) :
                return "MPI library at runtime is not compatible with MPI user compiled with";
            case( PETSC_ERR_ARG_SIZ ) :
                return "nonconforming object sizes used in operation";
            case( PETSC_ERR_ARG_IDN ) :
                return "two arguments not allowed to be the same";
            case( PETSC_ERR_ARG_WRONG ) :
                return "wrong argument (but object probably ok)";
            case( PETSC_ERR_ARG_CORRUPT ) :
                return "null or corrupted PETSc object as argument";
            case( PETSC_ERR_ARG_OUTOFRANGE ) :
                return "input argument, out of range";
            case( PETSC_ERR_ARG_BADPTR ) :
                return "invalid pointer argument";
            case( PETSC_ERR_ARG_NOTSAMETYPE ) :
                return "two args must be same object type";
            case( PETSC_ERR_ARG_NOTSAMECOMM ) :
                return "two args must be same communicators";
            case( PETSC_ERR_ARG_WRONGSTATE ) :
                return "object in argument is in wrong state, e.g. unassembled mat";
            case( PETSC_ERR_ARG_TYPENOTSET ) :
                return "the type of the object has not yet been set";
            case( PETSC_ERR_ARG_INCOMP ) :
                return "two arguments are incompatible";
            case( PETSC_ERR_ARG_NULL ) :
                return "argument is null that should not be";
            case( PETSC_ERR_ARG_UNKNOWN_TYPE ) :
                return "type name doesn't match any registered type";
            case( PETSC_ERR_FILE_OPEN ) :
                return "unable to open file";
            case( PETSC_ERR_FILE_READ ) :
                return "unable to read from file";
            case( PETSC_ERR_FILE_WRITE ) :
                return "unable to write to file";
            case( PETSC_ERR_FILE_UNEXPECTED ) :
                return "unexpected data in file";
            case( PETSC_ERR_MAT_LU_ZRPVT ) :
                return "detected a zero pivot during LU factorization";
            case( PETSC_ERR_MAT_CH_ZRPVT ) :
                return "detected a zero pivot during Cholesky factorization";
            case( PETSC_ERR_INT_OVERFLOW ) :
                return " int overflow " ;
            case( PETSC_ERR_FLOP_COUNT ) :
                return " flop count error " ;
            case( PETSC_ERR_NOT_CONVERGED ) :
                return "solver did not converge";
            case( PETSC_ERR_MISSING_FACTOR ) :
                return "MatGetFactor() failed";
            case( PETSC_ERR_OPT_OVERWRITE ) :
                return "attempted to over write options which should not be changed";
            case( PETSC_ERR_WRONG_MPI_SIZE ) :
                return "example/application run with number of MPI ranks it does not support";
            case( PETSC_ERR_USER_INPUT ) :
                return "missing or incorrect user input";
            case( PETSC_ERR_MAX_VALUE ) :
                return "this is always the one more than the largest error code";
            default :
                return "unknown error" ;
        }
#else
        return "We are not linked against PETSc";
#endif
    }

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_allocate_vector(
            MPI_Comm aComm,
            Vec      & aVec,
            const PetscInt aGlobalLength,
            const PetscInt aLocalLength
    )
    {
#ifdef BELFEM_PETSC
        // create the vector
        PetscErrorCode aStatus = VecCreate( aComm, & aVec );

        BELFEM_ASSERT( aStatus==0,
                      "PETSc has thrown error %i during VecCreate: %s",
                      ( int ) aStatus,
                      petsctools_error_message( aStatus ).c_str() );

        // set the site of the vector
        aStatus = VecSetSizes(
                aVec,
                aLocalLength,
                aGlobalLength );

        BELFEM_ASSERT( aStatus==0,
                      "PETSc has thrown error %i during VecSetSizes: %s",
                      ( int ) aStatus,
                      petsctools_error_message( aStatus ).c_str() );

        // set other parameters from passed options
        aStatus = VecSetFromOptions( aVec );

        BELFEM_ASSERT( aStatus==0,
                      "PETSc has thrown error %i during VecSetFromOptions: %s",
                      ( int ) aStatus,
                      petsctools_error_message( aStatus ).c_str() );

        return aStatus ;
#else
        return 0 ;
#endif
    }

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_set_vector(
            const Vector <PetscReal> & aVector,
            const Vector <PetscInt>  & aIndices,
            Vec & aVec )
    {
#ifdef BELFEM_PETSC
        // number of values to be set
        PetscInt tLength = aIndices.length() ;

        // populate values
        PetscErrorCode aStatus = VecSetValues( aVec,
                tLength,
                aIndices.data(),
                aVector.data(),
                INSERT_VALUES );

        BELFEM_ASSERT( aStatus==0,
              "PETSc has thrown error %i during VecSetValues(): %s",
              ( int ) aStatus,
              petsctools_error_message( aStatus ).c_str() );

        aStatus = VecAssemblyBegin( aVec );

        BELFEM_ASSERT( aStatus==0,
              "PETSc has thrown error %i during VecAssemblyBegin(): %s",
              ( int ) aStatus,
              petsctools_error_message( aStatus ).c_str() );

        aStatus = VecAssemblyEnd( aVec );

        BELFEM_ASSERT( aStatus==0,
              "PETSc has thrown error %i during VecAssemblyEnd(): %s",
              ( int ) aStatus,
              petsctools_error_message( aStatus ).c_str() );

        return aStatus ;
#else
        return 0 ;
#endif
    }

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_get_vector( Vec & aVec, const Vector< PetscInt > & aIndices, Vector< PetscReal > & aVector )
    {
#ifdef BELFEM_PETSC
        // number of values to be read
        PetscInt tLength = aIndices.length() ;

        // get vector data
        PetscErrorCode aStatus = VecGetValues( aVec, tLength, aIndices.data(), aVector.data() ) ;

        BELFEM_ASSERT( aStatus==0,
             "PETSc has thrown error %i during VecGetValues(): %s",
             ( int ) aStatus,
             petsctools_error_message( aStatus ).c_str() );

        return aStatus ;
#else
        return 0 ;
#endif
    }

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_get_local_vector( Vec & aVec, Vector< PetscReal > & aVector )
    {
#ifdef BELFEM_PETSC
        // get the local portion of the vector
        const PetscReal * tArray;
        PetscErrorCode aStatus = VecGetArrayRead( aVec, &tArray );

        BELFEM_ASSERT( aStatus==0,
             "PETSc has thrown error %i during VecGetArrayRead(): %s",
             ( int ) aStatus,
             petsctools_error_message( aStatus ).c_str() );

        // copy local values
        PetscInt tLocalSize = aVector.length();
        for( PetscInt i = 0; i < tLocalSize; ++i )
        {
            aVector( i ) = tArray[ i ];
        }

        // restore the array
        aStatus = VecRestoreArrayRead( aVec, &tArray );

        BELFEM_ASSERT( aStatus==0,
             "PETSc has thrown error %i during VecRestoreArrayRead(): %s",
             ( int ) aStatus,
             petsctools_error_message( aStatus ).c_str() );

        return aStatus ;
#else
        return 0 ;
#endif
    }

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_set_local_vector( const Vector< PetscReal > & aVector, Vec & aVec )
    {
#ifdef BELFEM_PETSC
        // get writable access to the local portion of the vector
        PetscReal * tArray;
        PetscErrorCode aStatus = VecGetArray( aVec, &tArray );

        BELFEM_ASSERT( aStatus==0,
             "PETSc has thrown error %i during VecGetArray(): %s",
             ( int ) aStatus,
             petsctools_error_message( aStatus ).c_str() );

        // copy local values
        PetscInt tLocalSize = aVector.length();
        for( PetscInt i = 0; i < tLocalSize; ++i )
        {
            tArray[ i ] = aVector( i );
        }

        // restore the array
        aStatus = VecRestoreArray( aVec, &tArray );

        BELFEM_ASSERT( aStatus==0,
             "PETSc has thrown error %i during VecRestoreArray(): %s",
             ( int ) aStatus,
             petsctools_error_message( aStatus ).c_str() );

        return aStatus ;
#else
        return 0 ;
#endif
    }

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_create_matrix(
        sparse::PETScAIJ * aDistMatrix,
        Mat & aMat )
    {
#ifdef BELFEM_PETSC

        PetscErrorCode aStatus = MatCreateMPIAIJWithSplitArrays(
            gComm.world(),
            aDistMatrix->n_rows(),
            aDistMatrix->n_rows(),
            aDistMatrix->size(),
            aDistMatrix->size(),
            aDistMatrix->diagonal_pointers(),
            aDistMatrix->diagonal_indices(),
            aDistMatrix->diagonal_values(),
            aDistMatrix->offdiagonal_pointers(),
            aDistMatrix->offdiagonal_indices(),
            aDistMatrix->offdiagonal_values(),
            & aMat );

        BELFEM_ERROR( aStatus==0,
                 "PETSc has thrown error %i during MatCreateMPIAIJWithSplitArrays: %s",
                 ( int ) aStatus,
                 petsctools_error_message( aStatus ).c_str() );

        // NOTE: We do NOT set ISLocalToGlobalMapping because we're using
        // global column indices in the off-diagonal block of MatCreateMPIAIJWithSplitArrays.
        // The garray is only used internally by BELFEM to track which global columns
        // appear in each rank's off-diagonal block.
        // If we were to set the mapping, PETSc would expect local indices instead.

        aStatus =  MatAssemblyBegin( aMat, MAT_FINAL_ASSEMBLY );
        BELFEM_ERROR( aStatus==0,
                    "PETSc has thrown error %i during MatAssemblyBegin: %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

        aStatus =  MatAssemblyEnd( aMat, MAT_FINAL_ASSEMBLY );
        BELFEM_ERROR( aStatus==0,
                    "PETSc has thrown error %i during MatAssemblyEnd: %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

        return aStatus ;
#else
        return 0 ;
#endif
    }

//------------------------------------------------------------------------------
    PetscErrorCode
    petsctools_update_matrix(
        sparse::PETScAIJ * aDistMatrix,
        Mat & aMat )
    {
#ifdef BELFEM_PETSC

        PetscErrorCode aStatus;
        // For matrices created with MatCreateMPIAIJWithSplitArrays, PETSc
        // directly uses the user-provided arrays. The values are already
        // updated in-place by distribute_values(), so we just need to
        // notify PETSc that the values have changed.
        //
        // Since PETSc uses our arrays directly, we only need to signal
        // that values changed by calling MatAssembly.


        // Notify PETSc that matrix values have been modified
        // This triggers preconditioner rebuild on next KSPSolve
        aStatus = MatAssemblyBegin( aMat, MAT_FINAL_ASSEMBLY );

        BELFEM_ERROR( aStatus==0,
                    "PETSc has thrown error %i during MatAssemblyBegin: %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

        aStatus = MatAssemblyEnd( aMat, MAT_FINAL_ASSEMBLY );

        BELFEM_ERROR( aStatus==0,
                    "PETSc has thrown error %i during MatAssemblyEnd: %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

        return aStatus;

#else
    return 0 ;
#endif
    }

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_link_matrix(
        COMM_TYPE   aComm,
        SpMatrix & aMatrix,
        Mat      & aMat )
    {
#ifdef BELFEM_PETSC
        // Ensure zero-based indexing
        aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

        // Create matrix that shares memory with the SpMatrix arrays
        PetscErrorCode aStatus = MatCreateSeqAIJWithArrays(
                aComm,
                aMatrix.n_rows(),
                aMatrix.n_cols(),
                aMatrix.pointers(),
                aMatrix.indices(),
                aMatrix.data(),
                & aMat );

        BELFEM_ERROR( aStatus == 0,
                "PETSc has thrown error %i during MatCreateSeqAIJWithArrays: %s",
                ( int ) aStatus,
                petsctools_error_message( aStatus ).c_str() );

        // Assembly is REQUIRED even when using MatCreateSeqAIJWithArrays
        // This finalizes the matrix structure for use by PETSc
        aStatus = MatAssemblyBegin( aMat, MAT_FINAL_ASSEMBLY );
        BELFEM_ERROR( aStatus == 0,
                "PETSc has thrown error %i during MatAssemblyBegin: %s",
                ( int ) aStatus,
                petsctools_error_message( aStatus ).c_str() );

        aStatus = MatAssemblyEnd( aMat, MAT_FINAL_ASSEMBLY );
        BELFEM_ERROR( aStatus == 0,
                "PETSc has thrown error %i during MatAssemblyEnd: %s",
                ( int ) aStatus,
                petsctools_error_message( aStatus ).c_str() );

        return aStatus;
#else
        return 0;
#endif
    }

//------------------------------------------------------------------------------

    PetscErrorCode
    petsctools_notify_matrix_update( Mat & aMat )
    {
#ifdef BELFEM_PETSC
        // Calling MatAssembly notifies PETSc that matrix values changed,
        // triggering preconditioner rebuild on next KSPSolve
        PetscErrorCode aStatus = MatAssemblyBegin( aMat, MAT_FINAL_ASSEMBLY );
        BELFEM_ERROR( aStatus == 0,
                "PETSc has thrown error %i during MatAssemblyBegin: %s",
                ( int ) aStatus,
                petsctools_error_message( aStatus ).c_str() );

        aStatus = MatAssemblyEnd( aMat, MAT_FINAL_ASSEMBLY );
        BELFEM_ERROR( aStatus == 0,
                "PETSc has thrown error %i during MatAssemblyEnd: %s",
                ( int ) aStatus,
                petsctools_error_message( aStatus ).c_str() );

        return aStatus;
#else
        return 0;
#endif
    }

//------------------------------------------------------------------------------
}