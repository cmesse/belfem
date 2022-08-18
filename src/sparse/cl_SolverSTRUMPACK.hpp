//
// Created by christian on 8/16/22.
//

#ifndef BELFEM_CL_SOLVERSTRUMPACK_HPP
#define BELFEM_CL_SOLVERSTRUMPACK_HPP

#include "cl_StringList.hpp"
#include "cl_SolverWrapper.hpp"
#ifdef BELFEM_STRUMPACK
#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#elif BELFEM_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#endif
#include <StrumpackSparseSolver.hpp>
#include <StrumpackSparseSolverMPIDist.hpp>
#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_CLANG
#pragma clang diagnostic pop
#endif
#endif

#include "cl_SolverStrumpackDistributor.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        class STRUMPACK : public Wrapper
        {
            int           mArgC ;
            StringList  * mArgV = nullptr ;

            StrumpackDistributor * mDistributor = nullptr ;


#ifdef BELFEM_STRUMPACK
            strumpack::StrumpackSparseSolver<real,int> * mSolver = nullptr ;
            strumpack::StrumpackSparseSolverMPIDist<real, int> * mDistSolver = nullptr ;
#endif
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            STRUMPACK() ;

//------------------------------------------------------------------------------

            ~STRUMPACK() ;

//------------------------------------------------------------------------------

            void
            solve(
                    SpMatrix & aMatrix,
                    Vector< real > & aLHS,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            /*void
            solve(
                    SpMatrix & aMatrix,
                    Matrix< real > & aLHS,
                    Matrix< real > & aRHS );

//------------------------------------------------------------------------------

            void
            free(); */

//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            initialize( SpMatrix & aMatrix,
                        const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                        const int aNumRhsColumns = 1 );

//------------------------------------------------------------------------------

            /*string
            error_message( const int aStatus ) const;  */

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_SOLVERSTRUMPACK_HPP
