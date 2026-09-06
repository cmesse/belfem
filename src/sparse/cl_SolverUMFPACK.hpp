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

            ~UMFPACK() override ;

//------------------------------------------------------------------------------

            void
            solve(
                    SpMatrix & aMatrix,
                    Vector< real > & aLHS,
                    Vector< real > & aRHS ) override ;

//------------------------------------------------------------------------------

            void
            solve(
                    SpMatrix & aMatrix,
                    Matrix< real > & aLHS,
                    Matrix< real > & aRHS ) override ;


//------------------------------------------------------------------------------

            void
            free() override ;

//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            initialize( SpMatrix & aMatrix,
                        const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                        const int_t aNumRhsColumns = 1 ) override ;

//------------------------------------------------------------------------------

            string
            error_message( const int aStatus ) const;

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_SOLVERDATAUMFPACK_HPP
