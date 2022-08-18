//
// Created by christian on 8/17/22.
//

#ifndef BELFEM_CL_SOLVERSTRUMPACKDISTRIBUTOR_HPP
#define BELFEM_CL_SOLVERSTRUMPACKDISTRIBUTOR_HPP

#include "typedefs.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace solver
    {

        /**
         * this class is responsible for distributing a Sparse CSR matrix
         * on several procs so that it can be solved in parallel by PETSc
         */
        class StrumpackDistributor
        {
            // rank of this proc
            const proc_t mMyRank;

            // number of running procs
            const proc_t mCommSize;

            // communication table for master proc
            Vector< proc_t > mCommTable ;

            // for master
            int mNumRowsOther = 0 ;
            Vector< int > mNumNnzPerProc ;

            // local values of the matrix
            int mMyNNZ = 0 ;
            int mMyNumRows = 0 ;

            Vector< int >  mDist ;
            Vector< int >  mMyIndices ;
            Vector< int >  mMyPointers ;
            Vector< real > mMyValues ;
            Vector< real > mMyRhs ;
            Vector< real > mMyLhs ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            StrumpackDistributor( SpMatrix & aMatrix );

//------------------------------------------------------------------------------

            StrumpackDistributor();

//------------------------------------------------------------------------------

            ~StrumpackDistributor() = default ;

//----------------------------------------------------------------------------

            void
            distribute();

//------------------------------------------------------------------------------

            int
            n_rows() const ;

//------------------------------------------------------------------------------

            int
            nnz() const ;

//------------------------------------------------------------------------------

            const int *
            dist() const ;

//------------------------------------------------------------------------------

            const int *
            pointers() const ;

//------------------------------------------------------------------------------

            const int *
            indices() const ;

//------------------------------------------------------------------------------

            const real *
            values() const ;

//------------------------------------------------------------------------------

            void
            send_values( const SpMatrix & aMatrix );

//------------------------------------------------------------------------------

            void
            receive_values();

//------------------------------------------------------------------------------

            real *
            rhs() ;

//------------------------------------------------------------------------------

            real *
            lhs() ;

//------------------------------------------------------------------------------

            void
            distribute_rhs( Vector< real> & aRHS );

//------------------------------------------------------------------------------

            void
            collect_lhs( Vector< real> & aLHS );

//------------------------------------------------------------------------------

            void
            send_sparsity_pattern( SpMatrix & aMatrix );

//------------------------------------------------------------------------------

            void
            receive_sparsity_pattern();

//------------------------------------------------------------------------------
        };


//------------------------------------------------------------------------------

        inline int
        StrumpackDistributor::n_rows() const
        {
            return mMyNumRows ;
        }

//------------------------------------------------------------------------------

        inline int
        StrumpackDistributor::nnz() const
        {
            return mMyNNZ ;
        }

//------------------------------------------------------------------------------

        inline const int *
        StrumpackDistributor::dist() const
        {
            return mDist.data() ;
        }

//------------------------------------------------------------------------------

        inline const int *
        StrumpackDistributor::pointers() const
        {
            return mMyPointers.data() ;
        }

//------------------------------------------------------------------------------

        inline const int *
        StrumpackDistributor::indices() const
        {
            return mMyIndices.data() ;
        }

//------------------------------------------------------------------------------

        inline const real *
        StrumpackDistributor::values() const
        {
            return mMyValues.data() ;
        }

//------------------------------------------------------------------------------

        inline real *
        StrumpackDistributor::rhs()
        {
            return mMyRhs.data() ;
        }

//------------------------------------------------------------------------------

        inline real *
        StrumpackDistributor::lhs()
        {
            return mMyLhs.data() ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_SOLVERSTRUMPACKDISTRIBUTOR_HPP
