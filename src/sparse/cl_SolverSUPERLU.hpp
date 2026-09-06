/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_SOLVERSUPERLU_HPP
#define BELFEM_CL_SOLVERSUPERLU_HPP

#ifdef BELFEM_SUPERLU
// <slu_ddefs.h> declares the global extern "C" BLAS prototypes ( dtrsv_,
// dgemm_, ... ) with signatures that clash with the matrix backend, which
// declares the same symbols differently ( e.g. Blaze uses blas_int_t and
// Fortran char-length arguments ). The wrapper never calls BLAS itself and
// libsuperlu.a links the real symbols, so SuperLU's unused prototypes are
// renamed for the duration of this include only. The #undef restores the
// names immediately so the backend's own declarations are untouched.
#define dcopy_ belfem_slu_dcopy_
#define daxpy_ belfem_slu_daxpy_
#define dgemm_ belfem_slu_dgemm_
#define dgemv_ belfem_slu_dgemv_
#define dtrsm_ belfem_slu_dtrsm_
#define dtrsv_ belfem_slu_dtrsv_
#include <slu_ddefs.h>
#undef dcopy_
#undef daxpy_
#undef dgemm_
#undef dgemv_
#undef dtrsm_
#undef dtrsv_
#else
    typedef void SuperMatrix ;
    typedef void superlu_options_t ;
    typedef void GlobalLU_t ;
#endif


#include "cl_SolverWrapper.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        class SUPERLU : public Wrapper
        {

            superlu_options_t * mOptions = nullptr ;

            SuperMatrix * mA = nullptr;

            // column permuted matrix A * Pc ( SLU_NCP ).
            // sp_preorder() aliases the value array of mA, so in-place value
            // updates of the SpMatrix propagate automatically. The SpMatrix
            // must NOT be reassigned or resized between symbolic() and free(),
            // as that would reallocate the value array and dangle this alias.
            SuperMatrix * mAC = nullptr ;

            // factors, owned by SuperLU
            SuperMatrix * mL = nullptr ;    // SLU_SC
            SuperMatrix * mU = nullptr ;    // SLU_NC

            // per-call LU workspace for dgstrf. Must be a valid struct for
            // every call; under plain SamePattern dgstrf reuses perm_c/etree,
            // not this workspace ( that is SamePattern_SameRowPerm ).
            GlobalLU_t * mGlu = nullptr;

            // column permutation ( fill reduction ), length n_cols
            int * mPermC = nullptr ;

            // row permutation ( pivoting ), length n_rows
            int * mPermR = nullptr ;

            // column elimination tree, length n_cols
            int * mEtree = nullptr ;

            bool mHaveSymbolic = false ;
            bool mHaveNumeric = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            SUPERLU() ;

//------------------------------------------------------------------------------

            ~SUPERLU() override ;

//------------------------------------------------------------------------------
            void
            solve( SpMatrix & aMatrix,
                   Vector <real> & aLHS,
                   Vector <real> & aRHS ) override ;

//------------------------------------------------------------------------------

            void
            solve( SpMatrix & aMatrix,
                   Matrix <real> & aLHS,
                   Matrix <real> & aRHS ) override ;

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

            void
            symbolic();

            void
            numeric();

//------------------------------------------------------------------------------

            // tear down the SuperLU factorization only ( no base-class reset );
            // shared by the public free() and by symbolic() re-entry
            void
            free_factorization();

//------------------------------------------------------------------------------

            void
            solve(  Vector< real > & aLHS,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            solve(  Matrix< real > & aLHS,
                    Matrix< real > & aRHS ) ;

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_SOLVERSUPERLU_HPP
