//
// Created by Christian Messe on 10.07.20.
//

#ifndef BELFEM_CL_SOLVERMUMPS_HPP
#define BELFEM_CL_SOLVERMUMPS_HPP

#include "cl_SolverWrapper.hpp"

namespace belfem
{
    namespace solver
    {
        class MUMPS : public Wrapper
        {
#ifdef BELFEM_MUMPS
            // Rank of HOST
            const proc_t mMasterRank ;
#endif
            // vector containing user parameters
            Vector< int >  mIParameters ;
            Vector< real > mRParameters ;

            // vector containing debug information
            Vector< int > mInfo ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            MUMPS( const proc_t aMasterRank = 0 );

//------------------------------------------------------------------------------

            ~MUMPS() ;

//------------------------------------------------------------------------------

            void
            solve(
                    SpMatrix & aMatrix,
                    Vector< real > & aLHS,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            solve(
                    SpMatrix & aMatrix,
                    Matrix< real > & aLHS,
                    Matrix< real > & aRHS );

//------------------------------------------------------------------------------

            real
            get_determinant() const ;

//------------------------------------------------------------------------------

            void
            free();

//------------------------------------------------------------------------------

            void
            set_reordering(
                    const SerialReodrdering   aSerial,
                    const ParallelReodrdering aParallel
            );

//------------------------------------------------------------------------------

            void
            set_block_low_ranking(
                    const BlockLowRanking aBLK,
                    const real            aEpsilon
            );


//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            initialize(
                    SpMatrix & aMatrix,
                    const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                    const int aNumRhsColumns=1 );


//------------------------------------------------------------------------------

            string
            error_message(
                    const int   * aInfo,
                    const int   & aN,
                    const int   & aNNZ ) ;
//------------------------------------------------------------------------------
        private :
//------------------------------------------------------------------------------

            void
            init_defaults();

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_SOLVERMUMPS_HPP
