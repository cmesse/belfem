//
// Created by Christian Messe on 13.07.20.
//

#ifndef BELFEM_CL_SOLVERPETSCDISTRIBUTOR_HPP
#define BELFEM_CL_SOLVERPETSCDISTRIBUTOR_HPP



#include "typedefs.hpp"
#include "petsctools.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"
#include "st_SolverPetscData.hpp"

namespace belfem
{
    namespace solver
    {

        /**
         * this class is responsible for distributing a Sparse CSR matrix
         * on several procs so that it can be solved in parallel by PETSc
         */
        class PetscDistributor
        {
            // rank of this proc
            const proc_t   mMyRank ;

            // number of running procs
            const PetscInt  mCommSize ;


            // the data object as owned by the PETSc wrapper
            PetscData & mData ;

#ifdef BELFEM_PETSC
            // the sparse matrix we want to solve
            SpMatrix  & mMatrix ;
#endif
            // data container for vector
            Vector< PetscReal > mRHS ; //?
            Vector< PetscReal > mLHS ; //?

            // communication list, used my master only
            // computed by create_communication_list()
            Vector< proc_t > mCommunicationList ;

#ifdef BELFEM_PETSC

            // how many nonzeros each proc has, master only,
            // computed by compute_rows_and_nnz_per_proc()
            Vector< PetscInt > mNumNnzPerProc ;

            // number of rows per proc, master only
            // computed by compute_rows_and_nnz_per_proc()
            Vector< PetscInt > mNumRowsPerProc ;

            // data container for matrix
            PetscReal * mMatrixData ; //?

            Vector< PetscReal > mSwap ;
#endif

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            PetscDistributor( PetscData & aData, SpMatrix & aMatrix ) ;

//------------------------------------------------------------------------------

            ~PetscDistributor() ;

//------------------------------------------------------------------------------

            // synchronize data for all matrices and vectors
            void
            synchronize(
                    Vector< PetscReal > & aLHS,
                    Vector< PetscReal > & aRHS,
                    bool aHaveLHS ); //?

//------------------------------------------------------------------------------

            void
            recover_lhs( Vector< PetscReal > & aLHS );

            const Vector< proc_t > &
            comm_list() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            // to be called from master only during distribute_indices
            void
            create_communication_list();

//----------------------------------------------------------------------------

            void
            split_ownership();

//------------------------------------------------------------------------------

            // to be called from master only during distribute_indices
            void
            compute_memory();

//------------------------------------------------------------------------------

            /**
             * initialize the matrix on the data container
             */
            PetscErrorCode
            initialize_matrix() ;

            void
            initialize_vectors() ;

//------------------------------------------------------------------------------

            /**
             * distribute the data of the matrix from the master proc
             * and copy the values into the Mat container
             */
            PetscErrorCode
            update_matrix() ;

//------------------------------------------------------------------------------

            void
            distribute_vector( Vector< PetscReal > & aGlobal, Vector< PetscReal > & aLocal );

//------------------------------------------------------------------------------

            void
            collect_vector( Vector< PetscReal > & aGlobal, Vector< PetscReal > & aLocal );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const Vector< proc_t > &
        PetscDistributor::comm_list() const
        {
            return mCommunicationList;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_SOLVERPETSCDISTRIBUTOR_HPP
