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

#ifndef BELFEM_CL_SOLVERPARDISO_HPP
#define BELFEM_CL_SOLVERPARDISO_HPP

#include "cl_SolverWrapper.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        class PARDISO : public Wrapper
        {
            /**
            * List of Parameters
            * 0: Matrix Type      : 0 - CSR
            *                     : 1 - CSC
            *
            * 1: Indexing Base    :  0 - C++
             *                       1 - Fortran ( SpMatrixIndexingBase )
            *
            * 2: Symmetry Mode    : 11 - Unsymmetric,
            *                        2 - PositiveDefiniteSymmetric,
            *                       -2 - GeneralSymmetric
            * 3: Info Level       : 0 - Silent
            * 4: Precon Exponent
            * 5: Max Number of
            *    Refinement steps : 0 - auto
            * 6: Reordering       : fill-in reducing ordering,
            *                       passed to iparm(2)
            *
            * 7: Compute Determinant : 0 - off
            *                          1 - on
            */
            Vector< int_t > mParameters ;

            /**
             * 0 : phase in which error has occured ( if it was an error )
             * 1 : Number of performed iterative refinement steps ( #7 )
             * 2 : number of nonzeros in the factor LU            ( #18 )
             * 3 : Output: Mflops for LU factorization            ( #19 )
             * 4 : CGS diagnostic                                 ( #20 )
             * 5 : Number of positive eigenvalues                 ( #22 )
             * 6 : Number of negative eigenvalues                 ( #23 )
             * 7 : compute-determinant flag as seen by PARDISO      ( #33 )
             *     the value itself is dparm( 33 ), see get_determinant()
             */
            Vector< int_t > mInfo ;

#ifdef BELFEM_PARDISO
            // default indexing base for solver
            const SpMatrixIndexingBase mIndexingBase = SpMatrixIndexingBase::Fortran ;
#endif

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            PARDISO() ;

//------------------------------------------------------------------------------

            ~PARDISO() override ;

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

            real
            get_determinant() const override;

//------------------------------------------------------------------------------

            void
            free() override ;

//------------------------------------------------------------------------------
        protected :
//------------------------------------------------------------------------------

            void
            initialize( SpMatrix & aMatrix,
                        const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                        const int_t aNumRhsColumns = 1) override ;

//------------------------------------------------------------------------------

            string
            error_message( const int_t aStatus ) const;

//------------------------------------------------------------------------------

            void
            check_status( const int_t aStatus ) ;

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_SOLVERPARDISO_HPP
