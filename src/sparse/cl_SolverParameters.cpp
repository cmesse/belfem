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

#include <cmath>

#include <limits>
#include "cl_SolverParameters.hpp"

#include "commtools.hpp"

namespace belfem
{
    SolverParameters::SolverParameters( const SolverType aType ) :
        mCommRank( comm_rank() ),
        mSolverType( aType )
    {
        // parallel PETSc requires the distributed AIJ format; the CSR
        // default only works on a single proc ( always-active BELFEM_ERROR
        // in PETSC::initialize ). Callers can still override via
        // set_distributed_matrix_type().
        if ( mSolverType == SolverType::PETSc && comm_size() > 1 )
        {
            mDistributedMatrixType = DistributedMatrixType::AIJ ;
        }
    }

    SolverParameters::SolverParameters( const input::Section * aInput ) :
        mCommRank( comm_rank() ),
        mSolverType( this->get_solver_type_from_input( aInput ) )
    {
        // parallel PETSc requires the distributed AIJ format ( see above );
        // an explicit "matrix format" key below still wins
        if ( mSolverType == SolverType::PETSc && comm_size() > 1 )
        {
            mDistributedMatrixType = DistributedMatrixType::AIJ ;
        }

        if ( aInput->key_exists( "matrix format" ) )
        {
            mDistributedMatrixType = belfem::distributed_matrix_type( aInput->get_string( "matrix format" ) );
        }
        if ( aInput->key_exists( "krylov method" ) )
        {
            mKrylovMethod = belfem::krylov_method( aInput->get_string( "krylov method" ) );
        }
        if ( aInput->key_exists( "preconditioner" ) )
        {
            mPreconditioner = belfem::preconditioner( aInput->get_string( "preconditioner" ) );
        }
        if ( aInput->key_exists( "reordering scheme" ) )
        {
            mReorderingMethod = belfem::reordering_method( aInput->get_string( "reordering scheme" ) );
        }
        if ( aInput->key_exists( "compression scheme" ) )
        {
            mCompressionMethod = belfem::compression_method( aInput->get_string( "compression scheme" ) );
        }
        if ( aInput->key_exists( "relative tolerance" ) )
        {
            mRelativeTolerance = aInput->get_real( "relative tolerance" );
            mHaveRelativeTolerance = true ;
        }
        if ( aInput->key_exists( "absolute tolerance" ) )
        {
            // one validation path: the setter checks > 0 and finite
            this->set_absolute_tolerance(
                aInput->get_real( "absolute tolerance" ) );
        }
        if ( aInput->key_exists( "initial guess" ) )
        {
            mUseInitialGuess = aInput->get_bool( "initial guess" );
        }
        if ( aInput->key_exists( "matching" ) )
        {
            mUseMatrixMatching = aInput->get_bool( "matching" );
        }
        if ( aInput->key_exists( "metis nodendp" ) )
        {
            mUseMetisNodeNDP = aInput->get_bool( "metis nodendp" );
        }
        if ( aInput->key_exists( "compression cutoff" ) )
        {
            // one validation path: the setter checks > 0 and finite
            this->set_compression_cutoff(
                aInput->get_real( "compression cutoff" ) );
        }
        if ( aInput->key_exists( "memory budget" ) )
        {
            // an integer number of MB, deliberately without a unit token:
            // unit_to_si() knows lengths and volumes, not bytes, and
            // get_value() would refuse "MB". Read as a real and checked
            // for integrality here: get_int() rounds, and a deck saying
            // 1.6 would silently become 2. Range is checked in the setter
            const real tBudget = aInput->get_real( "memory budget" );
            BELFEM_ERROR( std::isfinite( tBudget ) && tBudget > 0.0
                          && tBudget == std::floor( tBudget ),
                "memory budget in a linear solver section must be a positive whole number of MB, but is %g",
                ( double ) tBudget );
            this->set_memory_budget( ( uint ) tBudget );
        }
        if ( aInput->key_exists( "max iterations" ) )
        {
            // get_int() rounds a real: validate BEFORE the uint store,
            // where a negative value would wrap to a huge budget
            int tMaxIter = aInput->get_int( "max iterations" );
            BELFEM_ERROR( tMaxIter > 0,
                "max iterations in a linear solver section must be positive, but is %i",
                tMaxIter );
            mMaxNumIterations = ( uint ) tMaxIter ;
        }
    }

    SolverType
    SolverParameters::get_solver_type_from_input( const input::Section * aInput ) const
    {
        if ( aInput->key_exists( "library" ) )
        {
            return belfem::solver_type( aInput->get_string( "library" ) );
        }
        else
        {
            return gDefaultSolver;
        }
    }

    SolverParameters::SolverParameters( const SolverParameters & aOther ) :
      mCommRank( comm_rank() ),
      mSolverType(aOther.mSolverType),
      mDistributedMatrixType( aOther.mDistributedMatrixType ),
      mReorderingMethod(aOther.mReorderingMethod),
      mCompressionMethod(aOther.mCompressionMethod),
      mCompressionCutoff( aOther.mCompressionCutoff ),
      mHaveCompressionCutoff( aOther.mHaveCompressionCutoff ),
      mMemoryBudget( aOther.mMemoryBudget ),
      mHaveMemoryBudget( aOther.mHaveMemoryBudget ),
      mUseMatrixMatching( aOther.mUseMatrixMatching ),
      // Solver takes its parameters BY VALUE, so a member missing here is
      // silently reset to its default before the wrapper sees it — that
      // is exactly how a deck's "metis nodendp : false" was being lost
      // ( found by review 2026-08-15 )
      mUseMetisNodeNDP( aOther.mUseMetisNodeNDP ),
      mPreconditioner(aOther.mPreconditioner),
      mKrylovMethod(aOther.mKrylovMethod),
      mRelativeTolerance(aOther.mRelativeTolerance ),
      mHaveRelativeTolerance( aOther.mHaveRelativeTolerance ),
      mAbsoluteTolerance( aOther.mAbsoluteTolerance ),
      mHaveAbsoluteTolerance( aOther.mHaveAbsoluteTolerance ),
      mUseInitialGuess( aOther.mUseInitialGuess ),
      mMaxNumIterations( aOther.mMaxNumIterations )
    {

    }

    void
    SolverParameters::set_distributed_matrix_type( const DistributedMatrixType aType )
    {
        mDistributedMatrixType = aType ;
    }

    void
    SolverParameters::set_preconditioner( const Preconditioner aPreconditioner )
    {
        mPreconditioner = aPreconditioner ;
    }

    void
    SolverParameters::set_krylov_method( const KrylovMethod aKrylovMethod )
    {
        mKrylovMethod = aKrylovMethod ;
    }

    void
    SolverParameters::set_reordering_method( const ReorderingMethod aReorderingMethod )
    {
        mReorderingMethod = aReorderingMethod ;
    }

    void
    SolverParameters::set_compression_method( const CompressionMethod aCompressionMethod )
    {
        mCompressionMethod = aCompressionMethod ;
    }

    void
    SolverParameters::set_relative_tolerance( const real aEpsilon )
    {
        mRelativeTolerance = aEpsilon ;
        mHaveRelativeTolerance = true ;
    }

    void
    SolverParameters::set_absolute_tolerance( const real aEpsilon )
    {
        BELFEM_ERROR( std::isfinite( aEpsilon ) && aEpsilon > 0.0,
            "absolute tolerance must be positive and finite, but is %g",
            ( double ) aEpsilon );
        mAbsoluteTolerance = aEpsilon ;
        mHaveAbsoluteTolerance = true ;
    }

    void
    SolverParameters::set_use_initial_guess( const bool aUse )
    {
        mUseInitialGuess = aUse ;
    }

    void
    SolverParameters::set_matrix_matching( const bool aUse )
    {
        mUseMatrixMatching = aUse ;
    }

    SolverType
    SolverParameters::type() const
    {
        return mSolverType ;
    }


    DistributedMatrixType
    SolverParameters::distributed_matrix_type() const
    {
        return mDistributedMatrixType ;
    }

    Preconditioner
    SolverParameters::preconditioner() const
    {
        return mPreconditioner ;
    }

    KrylovMethod
    SolverParameters::krylov_method() const
    {
        return mKrylovMethod ;
    }

    ReorderingMethod
    SolverParameters::reordering_method() const
    {
        return mReorderingMethod ;
    }

    CompressionMethod
    SolverParameters::compression_method() const
    {
        return mCompressionMethod ;
    }

    real
    SolverParameters::relative_tolerance() const
    {
        return mRelativeTolerance ;
    }

    bool
    SolverParameters::have_relative_tolerance() const
    {
        return mHaveRelativeTolerance ;
    }

    real
    SolverParameters::absolute_tolerance() const
    {
        return mAbsoluteTolerance ;
    }

    bool
    SolverParameters::have_absolute_tolerance() const
    {
        return mHaveAbsoluteTolerance ;
    }

    uint
    SolverParameters::max_iterations() const
    {
        return mMaxNumIterations ;
    }

    real
    SolverParameters::compression_cutoff() const
    {
        return mCompressionCutoff ;
    }

    void
    SolverParameters::set_compression_cutoff( const real aCutoff )
    {
        // > 0: a zero cutoff is MUMPS's lossless-BLR niche, deliberately
        // unreachable from the deck; finite: +inf passes
        // a bare > 0 check and would become rel_tol = inf / CNTL(7) = inf
        BELFEM_ERROR( std::isfinite( aCutoff ) && aCutoff > 0.0,
            "compression cutoff must be positive and finite, but is %g",
            aCutoff );
        mCompressionCutoff = aCutoff ;
        mHaveCompressionCutoff = true ;
    }

    bool
    SolverParameters::have_compression_cutoff() const
    {
        return mHaveCompressionCutoff ;
    }

    void
    SolverParameters::set_memory_budget( const uint aMegaBytes )
    {
        // 0 is the "not stated" value and must not be reachable as a
        // statement: MUMPS reads ICNTL(23) = 0 as "size from the estimate",
        // which is the very behaviour a stated budget replaces
        BELFEM_ERROR( aMegaBytes > 0,
            "memory budget must be a positive number of MB, but is %u",
            ( unsigned int ) aMegaBytes );

        // the value crosses into the Fortran shim as int_t; above its
        // maximum the narrowing turns negative and the shim's guarded
        // write silently leaves ICNTL(23) untouched ( audit finding ).
        // 2^31 MB is 2 PB per process, so nothing real is refused
        BELFEM_ERROR( aMegaBytes <= ( uint ) std::numeric_limits< int_t >::max(),
            "memory budget of %u MB exceeds what the solver interface can carry ( %li MB )",
            ( unsigned int ) aMegaBytes,
            ( long ) std::numeric_limits< int_t >::max() );
        mMemoryBudget = aMegaBytes ;
        mHaveMemoryBudget = true ;
    }

    bool
    SolverParameters::have_memory_budget() const
    {
        return mHaveMemoryBudget ;
    }

    uint
    SolverParameters::memory_budget() const
    {
        return mMemoryBudget ;
    }

    bool
    SolverParameters::use_initial_guess() const
    {
        return mUseInitialGuess ;
    }

    bool
    SolverParameters::use_matrix_matching() const
    {
        return mUseMatrixMatching ;
    }

    bool
    SolverParameters::use_metis_nodendp() const
    {
        return mUseMetisNodeNDP ;
    }

    void
    SolverParameters::synchronize()
    {
        // the have-flags travel with their tolerances so every rank holds
        // the same record of what the deck stated. Since 2026-08-18 the
        // relative have-flag gates nothing ( STRUMPACK always applies the
        // shared value ), but it keeps its slot: the payload is positional
        // and shrinking it from the middle is how pack/unpack desyncs are
        // born ( width history: 7 -> 8 -> 9 -> 12 -> 13 )
        Vector< uint > tIData( 13 );
        Vector< real > tRData( 3 );

        if ( mCommRank == 0 )
        {
            uint tCount = 0 ;
            tIData( tCount++ ) = static_cast< uint >( mDistributedMatrixType );
            tIData( tCount++ ) = static_cast< uint >( mPreconditioner );
            tIData( tCount++ ) = static_cast< uint >( mKrylovMethod );
            tIData( tCount++ ) = static_cast< uint >( mReorderingMethod );
            tIData( tCount++ ) = static_cast< uint >( mCompressionMethod );
            tIData( tCount++ ) = static_cast< uint >( mUseInitialGuess );
            tIData( tCount++ ) = static_cast< uint >( mUseMatrixMatching );
            tIData( tCount++ ) = static_cast< uint >( mHaveRelativeTolerance );
            tIData( tCount++ ) = static_cast< uint >( mUseMetisNodeNDP );
            tIData( tCount++ ) = mMaxNumIterations ;
            tIData( tCount++ ) = static_cast< uint >( mHaveCompressionCutoff );
            tIData( tCount++ ) = static_cast< uint >( mHaveAbsoluteTolerance );
            tIData( tCount++ ) = mMemoryBudget ;

            tRData( 0 ) = mRelativeTolerance ;
            tRData( 1 ) = mCompressionCutoff ;
            tRData( 2 ) = mAbsoluteTolerance ;

            share( tIData );
            share( tRData );
        }
        else
        {
            receive( tIData );
            receive( tRData );
            uint tCount = 0 ;

            mDistributedMatrixType = static_cast< DistributedMatrixType >( tIData( tCount++ ) );
            mPreconditioner = static_cast< Preconditioner >( tIData( tCount++ ) );
            mKrylovMethod = static_cast< KrylovMethod >( tIData( tCount++ ) );
            mReorderingMethod = static_cast< ReorderingMethod >( tIData( tCount++ ) );
            mCompressionMethod = static_cast< CompressionMethod >( tIData( tCount++ ) );
            mUseInitialGuess =  tIData( tCount++ ) != 0 ;
            mUseMatrixMatching = tIData( tCount++ ) != 0 ;
            mHaveRelativeTolerance = tIData( tCount++ ) != 0 ;
            mUseMetisNodeNDP = tIData( tCount++ ) != 0 ;
            mMaxNumIterations = tIData( tCount++ ) ;
            mHaveCompressionCutoff = tIData( tCount++ ) != 0 ;
            mHaveAbsoluteTolerance = tIData( tCount++ ) != 0 ;
            mMemoryBudget = tIData( tCount++ ) ;
            // the have-flag is implied: the setter refuses 0
            mHaveMemoryBudget = mMemoryBudget > 0 ;
            mRelativeTolerance = tRData( 0 );
            mCompressionCutoff = tRData( 1 );
            mAbsoluteTolerance = tRData( 2 );
        }
    }

}
