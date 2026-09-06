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

#include "typedefs.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_SolverSUPERLU.hpp"

#ifdef BELFEM_SUPERLU
// SpMatrix index arrays ( belfem::int_t ) are handed to SuperLU without a
// copy, so SuperLU's ::int_t must have the same width. If this fails, rebuild
// SuperLU with a matching XSDK_INDEX_SIZE.
static_assert( sizeof( belfem::int_t ) == sizeof( ::int_t ),
    "BELFEM int_t and SuperLU int_t differ in size" );
#endif

namespace belfem
{
    namespace solver
    {
        SUPERLU::SUPERLU() :
            Wrapper( "SUPERLU  ", false )
        {

        }

        SUPERLU::~SUPERLU()
        {
            this->free();

#ifdef BELFEM_SUPERLU
            // free() releases the factorization; the options struct lives for
            // the whole object lifetime, so it is released here.
            delete mOptions ; mOptions = nullptr ;
#endif
        }

        void
        SUPERLU::solve( SpMatrix  & aMatrix, Vector< real > & aLHS, Vector< real > & aRHS )
        {
            if ( & aMatrix != mMatrix || ! mHaveSymbolic )
            {
                mMatrix = & aMatrix ;
                this->symbolic();
            }

            this->numeric();
            this->solve( aLHS, aRHS );
        }

        void
        SUPERLU::solve( SpMatrix  & aMatrix, Matrix< real > & aLHS, Matrix< real > & aRHS )
        {
            if ( & aMatrix != mMatrix || ! mHaveSymbolic )
            {
                mMatrix = & aMatrix ;
                this->symbolic();
            }

            this->numeric();
            this->solve( aLHS, aRHS );
        }

        void
        SUPERLU::initialize( SpMatrix & aMatrix,
                       const SymmetryMode aSymmetryMode,
                       const int_t aNumRhsColumns  )
        {
#ifdef BELFEM_SUPERLU

            // call initialize function from parent
            Wrapper::initialize();

            mMatrix = & aMatrix ;

            if ( mOptions != nullptr ) delete mOptions ;
            mOptions = new superlu_options_t ;
            set_default_options( mOptions ) ;

            // fill reducing ordering. COLAMD is the SuperLU default and a
            // good general purpose choice for unsymmetric FE matrices.
            // For structurally symmetric Jacobians, MMD_AT_PLUS_A is
            // sometimes faster; expose this as a setting if needed.
            mOptions->ColPerm = COLAMD ;

            switch (  aSymmetryMode )
            {
                case SymmetryMode::GeneralSymmetric:
                case SymmetryMode::PositiveDefiniteSymmetric :
                    mOptions->SymmetricMode = YES ;
                    break ;
                default:
                    mOptions->SymmetricMode = NO ;
                    break ;
            }

            mHaveSymbolic = false ;
            mHaveNumeric  = false ;

#endif
        }

        void SUPERLU::symbolic()
        {
#ifdef BELFEM_SUPERLU
            if ( mHaveSymbolic )
            {
                this->free_factorization();
            }


            BELFEM_ERROR( mMatrix->type() == SpMatrixType::CSC || mMatrix->type() == SpMatrixType::CSR,
                "SuperLU requires a matrix in CSC or CSR format" );

            BELFEM_ERROR( mMatrix->n_rows() > 0 && mMatrix->number_of_nonzeros() > 0
                              && mMatrix->data() != nullptr
                              && mMatrix->indices() != nullptr
                              && mMatrix->pointers() != nullptr,
                "SuperLU: matrix is empty or not fully constructed" );

            mOptions->Trans = mMatrix->type() == SpMatrixType::CSC ? NOTRANS : TRANS ;

            BELFEM_ERROR( mMatrix->n_rows() == mMatrix->n_cols(),
            "SuperLU driver expects a square matrix ( is %lu x %lu )",
                ( long unsigned int ) mMatrix->n_rows(),
                ( long unsigned int ) mMatrix->n_cols() );


            BELFEM_ERROR( mMatrix->n_rows() <= ( index_t ) BELFEM_INT_MAX,
                "SuperLU: matrix dimension %lu exceeds int range",
                ( long unsigned int ) mMatrix->n_rows() );

            BELFEM_ERROR( mMatrix->number_of_nonzeros()
                              <= ( index_t ) std::numeric_limits< ::int_t >::max(),
                "SuperLU: nnz %lu exceeds SuperLU int_t range",
                ( long unsigned int ) mMatrix->number_of_nonzeros() );

            int   tN   = static_cast< int > ( mMatrix->n_rows() );
            int_t tNNZ = static_cast< int_t >( mMatrix->number_of_nonzeros() ) ;

            mA = new SuperMatrix() ;

            mMatrix->set_indexing_base( SpMatrixIndexingBase::Cpp );

            // wrap the SpMatrix arrays, no copy is made
            dCreate_CompCol_Matrix(
                mA,
                tN,
                tN,
                tNNZ,
                mMatrix->data(),
                mMatrix->indices(),
                mMatrix->pointers(),
                SLU_NC,
                SLU_D,
                SLU_GE );

            mPermC = static_cast< int * >( malloc( tN * sizeof( int ) ) ) ;
            mPermR = static_cast< int * >( malloc( tN * sizeof( int ) ) ) ;
            mEtree = static_cast< int * >( malloc( tN * sizeof( int ) ) ) ;

            BELFEM_ERROR( mPermC != nullptr && mPermR != nullptr && mEtree != nullptr,
                "SuperLU: failed to allocate permutation/etree arrays" );

            // compute the column permutation ...
            get_perm_c( mOptions->ColPerm, mA, mPermC );

            // ... and the column elimination tree. Creates mAC = A * Pc,
            // ( AT * Pc in CSR mode ). sp_preorder() writes through mAC and
            // allocates mAC->Store, so the struct must exist first.
            mAC = new SuperMatrix() ;
            sp_preorder( mOptions, mA, mPermC, mEtree, mAC );

            // fresh pattern: the first numeric() must factor from scratch.
            // numeric() switches to SamePattern for subsequent value updates.
            mOptions->Fact = DOFACT ;

            mHaveSymbolic = true ;

#endif
        }

        void
        SUPERLU::numeric()
        {
#ifdef BELFEM_SUPERLU
            BELFEM_ERROR( mHaveSymbolic,
                      "SuperLU::numeric() called before symbolic()" );

            if ( mHaveNumeric )
            {
                Destroy_SuperNode_Matrix( mL );     // SLU_SC
                Destroy_CompCol_Matrix( mU );       // SLU_NC

                mOptions->Fact = SamePattern;

                mHaveNumeric  = false ;
            }

            if ( mL == nullptr ) mL = new SuperMatrix ;
            if ( mU == nullptr ) mU = new SuperMatrix ;

            // persistent supernode workspace, must be a valid struct for dgstrf
            if ( mGlu == nullptr ) mGlu = new GlobalLU_t ;

            SuperLUStat_t tStat ;
            StatInit( & tStat );

            int tPanelSize = sp_ienv( 1 );
            int tRelax     = sp_ienv( 2 );

            int_t tInfo = 0 ;
            dgstrf( mOptions,
               mAC,
               tRelax,
               tPanelSize,
               mEtree,
               nullptr,            // work array: let SuperLU allocate
               0,                  // lwork
               mPermC,
               mPermR,
               mL,
               mU,
               mGlu,
               & tStat,
               & tInfo );

            StatFree( & tStat );

            // dgstrf info: < 0 illegal argument, 0 success,
            // 0 < info <= ncol singular pivot ( 1-based ), > ncol out of memory
            BELFEM_ERROR( tInfo >= 0,
                "SuperLU dgstrf: argument %ld had an illegal value",
                ( long int ) -tInfo );

            // info >= 0: dgstrf has built the factors ( valid when singular,
            // partial on out-of-memory ). Mark them owned before the remaining
            // checks so that a throw in debug ( BELFEM_ERROR ) still lets
            // free_factorization() release mL/mU/mGlu.
            mHaveNumeric = true ;

            BELFEM_ERROR( tInfo <= ( int_t ) mMatrix->n_cols(),
                "SuperLU dgstrf ran out of memory ( %ld bytes allocated when failure occurred )",
                ( long int )( tInfo - ( int_t ) mMatrix->n_cols() ) );

            BELFEM_ERROR( tInfo == 0,
                "SuperLU dgstrf: matrix is singular, U( %ld, %ld ) is exactly zero ( 1-based )",
                ( long int ) tInfo,
                ( long int ) tInfo );
#endif
        }

        void
        SUPERLU::free()
        {
            this->free_factorization();

            // Only the public free() resets the base-class initialized state.
            // symbolic() re-entry calls free_factorization() directly, so
            // reusing the solver across different matrices does not trigger a
            // driver-level re-initialize ( which would orphan the factors ).
            Wrapper::free();
        }

        void
        SUPERLU::free_factorization()
        {
#ifdef BELFEM_SUPERLU

            if ( mHaveNumeric )
            {
                Destroy_SuperNode_Matrix( mL ); delete mL ; mL = nullptr ;
                Destroy_CompCol_Matrix( mU );   delete mU ; mU = nullptr ;
                delete mGlu ; mGlu = nullptr ;
                mHaveNumeric  = false ;
            }

            if ( mHaveSymbolic )
            {
                // IMPORTANT: Store only. The index and value arrays
                // belong to the SpMatrix; Destroy_CompCol_Matrix()
                // would free them and cause a double free in the
                // SpMatrix destructor.
                Destroy_SuperMatrix_Store( mA ); delete mA ; mA = nullptr ;


                Destroy_CompCol_Permuted( mAC ); delete mAC ; mAC = nullptr ;

                ::free( mPermC ); mPermC = nullptr ;
                ::free( mPermR ); mPermR = nullptr ;
                ::free( mEtree ); mEtree = nullptr ;

                mHaveSymbolic = false ;
            }

#endif
        }

        void
        SUPERLU::solve( Vector< real > & aLHS, Vector< real > & aRHS )
        {
#ifdef BELFEM_SUPERLU
            int tN = mMatrix->n_rows() ;
            BELFEM_ERROR( aRHS.length() == ( index_t ) tN,
                       "RHS length %lu does not match matrix size %d",
                       ( long unsigned int ) aRHS.length(), tN );

            // dgstrs overwrites the right hand side with the solution,
            // so we solve in place on the LHS vector
            if( & aLHS != & aRHS )
            {
                aLHS = aRHS ;
            }

            // dense wrapper around the solution vector, no copy
            SuperMatrix tB ;

            dCreate_Dense_Matrix(
                    & tB,
                    tN,
                    1,                  // one right hand side
                    aLHS.data(),
                    tN,                 // leading dimension
                    SLU_DN,
                    SLU_D,
                    SLU_GE );

            SuperLUStat_t tStat ;
            StatInit( & tStat );

            int tInfo = 0 ;

            dgstrs( mOptions->Trans,
                    mL,
                    mU,
                    mPermC,
                    mPermR,
                    & tB,
                    & tStat,
                    & tInfo );

            StatFree( & tStat );

            // only frees the wrapper struct, not the vector data
            Destroy_SuperMatrix_Store( & tB );

            BELFEM_ERROR( tInfo == 0,
                "SuperLU dgstrs returned error = %d", tInfo );
#endif

        }

        void
        SUPERLU::solve( Matrix< real > & aLHS, Matrix< real > & aRHS )
        {
#ifdef BELFEM_SUPERLU
            int tN = mMatrix->n_rows() ;
            int tM = aRHS.n_cols() ;

            BELFEM_ERROR( aRHS.n_rows() == ( index_t ) tN,
                       "RHS row count %lu does not match matrix size %d",
                       ( long unsigned int ) aRHS.n_rows(), tN );

            // make sure the solution matrix has the right shape
            if ( aLHS.n_rows() != ( index_t ) tN || aLHS.n_cols() != ( index_t ) tM )
            {
                aLHS.set_size( tN, tM );
            }

            // pack the right hand side into a contiguous, column-major buffer.
            // Matrix storage may be padded ( Blaze ), so dgstrs cannot operate
            // on aRHS.data() with leading dimension tN directly.
            real * tData = static_cast< real * >( malloc( tN * tM * sizeof( real ) ) );
            BELFEM_ERROR( tData != nullptr, "SuperLU: failed to allocate RHS buffer" );

            index_t tCount = 0 ;
            for ( int j=0; j<tM; ++j )
            {
                for ( int i=0; i<tN; ++i )
                {
                    tData[ tCount++ ] = aRHS( i, j ) ;
                }
            }

            // dense wrapper around the packed buffer; dgstrs solves in place
            SuperMatrix tB ;

            dCreate_Dense_Matrix(
                    & tB,
                    tN,
                    tM,
                    tData,
                    tN,                 // leading dimension
                    SLU_DN,
                    SLU_D,
                    SLU_GE );

            SuperLUStat_t tStat ;
            StatInit( & tStat );

            int tInfo = 0 ;

            dgstrs(  mOptions->Trans,
                    mL,
                    mU,
                    mPermC,
                    mPermR,
                    & tB,
                    & tStat,
                    & tInfo );

            StatFree( & tStat );

            // only frees the wrapper struct, not tData
            Destroy_SuperMatrix_Store( & tB );

            // scatter the solution into the LHS matrix on success
            if ( tInfo == 0 )
            {
                tCount = 0 ;
                for ( int j=0; j<tM; ++j )
                {
                    for ( int i=0; i<tN; ++i )
                    {
                        aLHS( i, j ) = tData[ tCount++ ] ;
                    }
                }
            }

            ::free( tData ) ;

            BELFEM_ERROR( tInfo == 0,
                "SuperLU dgstrs returned error = %d", tInfo );
#endif

        }
    }
}
