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

#ifndef BELFEM_CL_SOLVERSTRUMPACK_HPP
#define BELFEM_CL_SOLVERSTRUMPACK_HPP

#include "cl_StringList.hpp"
#include "cl_SolverWrapper.hpp"
#include "cl_SolverParameters.hpp"

#include "strumpacktools.hpp"


namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        class STRUMPACK : public Wrapper
        {
            const SolverParameters * mParams ;
#ifdef BELFEM_STRUMPACK

            int           mArgC = 0 ;
            StringList  * mArgV = nullptr ;

            sparse::DistMatrix * mDistMatrix = nullptr ;

            strumpack::StrumpackSparseSolver<real,int> * mSolver = nullptr ;
            strumpack::StrumpackSparseSolverMPIDist<real, int> * mDistSolver = nullptr ;
#endif
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /** aParams is borrowed and must outlive the wrapper. */
            STRUMPACK( const SolverParameters * aParams ) ;

//------------------------------------------------------------------------------

            ~STRUMPACK() override ;

//------------------------------------------------------------------------------

            /** Collective. A success/failure verdict is reduced across ranks; a
             *  failure is recorded through flag_failure() when soft-fail is armed
             *  and aborts otherwise. The inherited matrix overload aborts. */
            void
            solve(
                    SpMatrix & aMatrix,
                    Vector< real > & aLHS,
                    Vector< real > & aRHS ) override ;

//------------------------------------------------------------------------------

            void
            free() override;

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
                        const int_t aNumRhsColumns = 1 ) override ;

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_SOLVERSTRUMPACK_HPP
