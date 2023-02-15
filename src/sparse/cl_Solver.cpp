//
// Created by Christian Messe on 10.07.20.
//

#include "cl_Solver.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_Timer.hpp"

#include "cl_SolverUMFPACK.hpp"
#include "cl_SolverMUMPS.hpp"
#include "cl_SolverPARDISO.hpp"
#include "cl_SolverPETSC.hpp"
#include "cl_SolverSTRUMPACK.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Solver::Solver( const SolverType aSolverType,
        const proc_t aMasterRank ) :
        mType( aSolverType )
    {
        switch ( aSolverType )
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
            case ( SolverType::MUMPS ) :
            {
#ifdef BELFEM_MUMPS
                mWrapper = new solver::MUMPS();
#else
                BELFEM_ERROR( false,
                             "You are trying to create a MUMPS solver.\nHowever, we are not linked against MUMPS." );
#endif
                break;
            }
            case ( SolverType::STRUMPACK ) :
            {
#ifdef BELFEM_STRUMPACK
                mWrapper = new solver::STRUMPACK();
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
            case( SolverType::PETSC ) :
            {
#ifdef BELFEM_PETSC
                mWrapper = new solver::PETSC();
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

            message( 4,
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
        if( mType == SolverType::PETSC )
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
            const SerialReodrdering   aSerial,
            const ParallelReodrdering aParallel )
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
            const BlockLowRanking aBlr,
            const real            aEpsilon )
    {
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