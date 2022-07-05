//
// Created by Christian Messe on 11.07.20.
//

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
            * 1: Indexing Base    :  0 - Fortran
             *                       1 - C++
            *
            * 2: Symmetry Mode    : 11 - Unsymmetric,
            *                        2 - PositiveDefiniteSymmetric,
            *                       -2 - GeneralSymmetric
            * 3: Info Level       : 0 - Silent
            * 4: Precon Exponent
            * 5: Max Number of
            *    Refinement steps : 0 - auto
            * 6: Solver Flag      : 0 - direct
            *                       1 - iterative
            *
            * 7: Compute Determinant : 0 - off
            *                          1 - on
            */
            Vector< int > mParameters ;

            /**
             * 0 : phase in which error has occured ( if it was an error )
             * 1 : Number of performed iterative refinement steps ( #7 )
             * 2 : number of nonzeros in the factor LU            ( #18 )
             * 3 : Output: Mflops for LU factorization            ( #19 )
             * 4 : CGS diagnostic                                 ( #20 )
             * 5 : Number of positive eigenvalues                 ( #22 )
             * 6 : Number of negative eigenvalues                 ( #23 )
             * 7 : ! natural log of determinant ( if computed )   ( #33 )
             */
            Vector< int > mInfo ;

#ifdef BELFEM_PARDISO
            // default indexing base for solver
            const SpMatrixIndexingBase mIndexingBase = SpMatrixIndexingBase::Fortran ;
#endif

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            PARDISO() ;

//------------------------------------------------------------------------------

            ~PARDISO() ;

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
        protected :
//------------------------------------------------------------------------------

            void
            initialize( SpMatrix & aMatrix,
                        const SymmetryMode aSymmetryMode = SymmetryMode::Unsymmetric,
                        const int aNumRhsColumns = 1);

//------------------------------------------------------------------------------

            string
            error_message( const int aStatus ) const;

//------------------------------------------------------------------------------

            void
            check_status( const int aStatus ) ;

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_SOLVERPARDISO_HPP
