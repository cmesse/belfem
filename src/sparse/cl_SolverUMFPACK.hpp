//
// Created by Christian Messe on 10.07.20.
//

#ifndef BELFEM_CL_SOLVERDATAUMFPACK_HPP
#define BELFEM_CL_SOLVERDATAUMFPACK_HPP

#include "cl_SolverWrapper.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        class UMFPACK : public Wrapper
        {
            int mTransposedFlag ;

            // symbolic factorization
            void * mSymbolic = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            UMFPACK() ;

//------------------------------------------------------------------------------

            ~UMFPACK() ;

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

            void
            free();

//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            initialize( SpMatrix & aMatrix,
                        const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                        const int aNumRhsColumns = 1 );

//------------------------------------------------------------------------------

            string
            error_message( const int aStatus ) const;

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_SOLVERDATAUMFPACK_HPP
