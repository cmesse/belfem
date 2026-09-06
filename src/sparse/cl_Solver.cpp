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

#include "cl_Solver.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_Timer.hpp"

#include "cl_SolverUMFPACK.hpp"
#include "cl_SolverSUPERLU.hpp"
#include "cl_SolverMUMPS.hpp"
#include "cl_SolverPARDISO.hpp"
#include "cl_SolverPETSC.hpp"
#include "cl_SolverSTRUMPACK.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Solver::Solver( const SolverType aSolverType ) :
        mType( aSolverType ),
        mParams( SolverParameters( aSolverType ) )
    {
        this->create_wrapper();
    }

    Solver::Solver( const SolverParameters aParams ) :
        mType( aParams.type() ),
        mParams( aParams )
    {
        this->create_wrapper();
    }

    void Solver::create_wrapper()
    {
        // make sure that the user settings are consistent over all procs
        mParams.synchronize();

        // the BELFEM-side parallel orderings ( PETSc path, DistMatrix ) need
        // the wrapper libraries linked; MUMPS and STRUMPACK bring their own
        // ParMETIS / PT-Scotch and map these values themselves, so the check
        // is scoped to the solver that uses ours. A missing library is a
        // configuration error, not a fallback. Checked here, after
        // synchronize(), so every rank holds the same value and every rank
        // takes the same reaction
        if ( mParams.type() == SolverType::PETSc )
        {
#ifndef BELFEM_PARMETIS
            BELFEM_ERROR( mParams.reordering_method() != ReorderingMethod::PARMETIS,
                "reordering scheme 'parmetis' requested for PETSc, but we are not linked against ParMETIS" );
#endif
#ifndef BELFEM_PTSCOTCH
            BELFEM_ERROR( mParams.reordering_method() != ReorderingMethod::PTSCOTCH,
                "reordering scheme 'ptscotch' requested for PETSc, but we are not linked against PT-Scotch" );
#endif
        }

        switch ( mParams.type() )
        {
            case ( SolverType::UMFPACK ) :
            {
#ifdef BELFEM_SUITESPARSE
                mWrapper = new solver::UMFPACK();

#else
                BELFEM_ERROR( false,
                    "You are trying to create an UMFPACK solver.\nHowever, we are not linked against UMFPACK." );
#endif
                break;
            }
                case ( SolverType::SUPERLU ) :
            {
#ifdef BELFEM_SUPERLU
                mWrapper = new solver::SUPERLU();

#else
                BELFEM_ERROR( false,
                    "You are trying to create a SUPERLU solver.\nHowever, we are not linked against SUPERLU." );
#endif
                break;
            }
            case ( SolverType::MUMPS ) :
            {
#ifdef BELFEM_MUMPS
                mWrapper = new solver::MUMPS( &mParams );
#else
                BELFEM_ERROR( false,
                             "You are trying to create a MUMPS solver.\nHowever, we are not linked against MUMPS." );
#endif
                break;
            }
            case ( SolverType::STRUMPACK ) :
            {
#ifdef BELFEM_STRUMPACK
                mWrapper = new solver::STRUMPACK( &mParams );
#else
                BELFEM_ERROR( false,
                             "You are trying to create a STRUMPACK solver.\nHowever, we are not linked against STRUMPACK." );
#endif
                break;
            }
            case( SolverType::PARDISO ) :
            {
#ifdef BELFEM_PARDISO
                mWrapper = new solver::PARDISO();
#else
                BELFEM_ERROR( false,
                             "You are trying to create a PARDISO solver.\nHowever, we are not linked against PARDISO." );
#endif
                break;
            }
            case( SolverType::PETSc ) :
            {
#ifdef BELFEM_PETSC
                mWrapper = new solver::PETSC( &mParams );
#else
                BELFEM_ERROR( false,
                             "You are trying to create a PETSc solver.\nHowever, we are not linked against PETSc." );
#endif
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "unknown solver type" );
            }
        }
    }

//------------------------------------------------------------------------------

    Solver::~Solver()
    {
        delete mWrapper ;
    }

//------------------------------------------------------------------------------

    SolverType
    Solver::type() const
    {
        return mType ;
    }

//------------------------------------------------------------------------------

    const SolverParameters &
    Solver::parameters() const
    {
        return mParams ;
    }

//------------------------------------------------------------------------------

    void
    Solver::set_symmetry_mode( const SymmetryMode & aMode )
    {
        mSymmetryMode = aMode ;
    }

//------------------------------------------------------------------------------

    void
    Solver::solve(
            SpMatrix & aMatrix,
            Vector< real > & aLHS,
            Vector< real > & aRHS )
    {
        // make sure that the wrapper has been initialized
        if ( !mWrapper->is_initialized() )
        {
            mWrapper->initialize( aMatrix, mSymmetryMode, 1 );
        }


        if ( gLog.info_level() >= 4 && mWrapper->rank() == 0 )
        {
            Timer tTimer ;

            // solve the system
            mWrapper->solve( aMatrix, aLHS, aRHS );

            // we use umfpack or superlu for small problems, this message is becoming annoying
            if ( this->type() != SolverType::UMFPACK &&  this->type() != SolverType::SUPERLU  )
                message( InfoLevel::Verbose ,
                        "    ... time for solving system using %s  : %u ms\n",
                         mWrapper->label().c_str(),
                         ( unsigned int ) tTimer.stop() );
        }
        else
        {
            // solve the system
            mWrapper->solve( aMatrix, aLHS, aRHS );
        }
    }

//------------------------------------------------------------------------------

    void
    Solver::solve(  SpMatrix       & aMatrix,
            Matrix< real > & aLHS,
            Matrix< real > & aRHS )
    {
        // make sure that the wrapper has been initialized
        if( ! mWrapper->is_initialized() )
        {
            mWrapper->initialize(
                    aMatrix,
                    mSymmetryMode,
                    aRHS.n_cols() );
        }

        // solve the system
        mWrapper->solve( aMatrix, aLHS, aRHS );

    }

//------------------------------------------------------------------------------

    void
    Solver::free()
    {
        mWrapper->free() ;
    }

//------------------------------------------------------------------------------

    void
    Solver::set_petsc(
            const Preconditioner aPreconditioner,
            const KrylovMethod   aKrylovMethod,
            const real           aEpsilon )
    {
        if( mType == SolverType::PETSc )
        {
            // get wrapper as petsc
            solver::PETSC * tPETSC =
                    reinterpret_cast< solver::PETSC * >( mWrapper ) ;

            // write information
            tPETSC->set(
                    aPreconditioner,
                    aKrylovMethod,
                    aEpsilon ) ;
        }
    }

    void
    Solver::set_mumps_reordering(
            const MumpsSerialReodrdering   aSerial,
            const MumpsParallelReodrdering aParallel )
    {
        if( mType == SolverType::MUMPS )
        {
            // get weapper
            solver::MUMPS * tMUMPS =
                    reinterpret_cast< solver::MUMPS * >( mWrapper );

            // write information
            tMUMPS->set_reordering( aSerial, aParallel );
        }
    }

//------------------------------------------------------------------------------

    void
    Solver::set_mumps_blr(
            const MumpsBlockLowRanking aBlr,
            const real            aEpsilon )
    {
        // CAVEAT: initialize() rewrites the CNTL(7) slot
        // from SolverParameters, so a call BEFORE initialize is clobbered
        // — aBlr survives, aEpsilon becomes 0.0 unless the parameters
        // also say blr. The only reliable uses are: configure the
        // PARAMETERS ( method + cutoff ) before constructing the Solver,
        // or call this AFTER initialize ( and again after any
        // re-initialize ). Mutating the parameters post-construction does
        // NOT work either — ICNTL(35) is written once at construction.
        // Zero callers today; this note exists so the API is not mistaken
        // for a working escape hatch.
        if( mType == SolverType::MUMPS )
        {
            // get weapper
            solver::MUMPS * tMUMPS =
                    reinterpret_cast< solver::MUMPS * > ( mWrapper );

            // write information
            tMUMPS->set_block_low_ranking( aBlr, aEpsilon );
        }
    }

//------------------------------------------------------------------------------

    void
    Solver::set_mumps_error_analysis(
            const MumpsErrorAnalysis aSetting )
    {
        if( mType == SolverType::MUMPS )
        {
            // get weapper
            solver::MUMPS * tMUMPS =
                    reinterpret_cast< solver::MUMPS * > ( mWrapper );

            // write information
            tMUMPS->set_error_analysis( aSetting );
        }
    }

//------------------------------------------------------------------------------

    /**
     * this does only do something if PARDISO is used
     *
     */
/*    void
    Solver::set_pardiso( const PardisoMode aMode )
    {
        if( mType == SolverType::PARDISO )
        {
            // get wrapper as petsc
            solver::PARDISO * tPARDISO =
                    reinterpret_cast< solver::PARDISO * >( mWrapper ) ;

            // write information
            tPARDISO->set_mode( aMode ) ;
        }
    } */

//------------------------------------------------------------------------------

}