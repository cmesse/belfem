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

#include "cl_Logger.hpp"
#include <cstdlib>
#include <string>
#if defined( __GLIBC__ ) && ! defined( __APPLE__ )
#include <malloc.h>          // malloc_trim, see the call in finalize()
#endif

#include <ctime>
#include <cstring>
#include <fstream>
#include <iomanip>

#include "cl_FEM_Controller.hpp"
#include "fn_FEM_ghost_switch.hpp"
#include "cl_Map.hpp"
#include "constants.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_IWG_Timestep.hpp"
#include "fn_check_unit.hpp"
#include "fn_max.hpp"
#include "fn_sum.hpp"

#include "cl_Profiler.hpp"

extern belfem::Logger gLog;
namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        namespace
        {
            /**
             * wall-clock stamp for the per-step summary box, so a log read
             * days later dates itself ( the file mtime only dates the LAST
             * write ). Fixed-width fields throughout — %d and %I are
             * zero-padded and %Z is 3-4 characters here — but the caller
             * still pads to the box width, so a long zone name cannot skew
             * the frame
             */
            std::string
            wallclock_stamp()
            {
                std::time_t tSystemTime = std::time( NULL );
                std::tm * tCalendarTime = std::localtime( &tSystemTime );

                char tBuffer[ 64 ];
                if ( tCalendarTime == nullptr
                    || std::strftime( tBuffer, sizeof( tBuffer ),
                                      "%a %b %d %I:%M:%S %p %Z %Y",
                                      tCalendarTime ) == 0 )
                {
                    return "" ;
                }
                return std::string( tBuffer );
            }
        }

//------------------------------------------------------------------------------

        Controller::Controller( Kernel * aKernel,  Kernel * aKernel2 ) :
            mCommRank( comm_rank() ),
            mCommSize( comm_size() ),
            mKernel( aKernel ),
            mEquation( reinterpret_cast< IWG_Timestep * >( aKernel->dofmgr()->iwg() ) ),
            mMesh( aKernel->mesh() ),
            mKernel2( aKernel2 ),
            mEquation2( aKernel2 == nullptr ? nullptr : reinterpret_cast< IWG_Timestep * >( aKernel2->dofmgr()->iwg() ) ),
            mMesh2( aKernel2 == nullptr ? nullptr : aKernel2->mesh() ),
            mResidualHistory( mStallWindow )
        {
            // internal restart value for the relaxation, deliberately NOT an
            // input parameter: the divergence-counter reset ( see the AIMD
            // blocks ) jumps back to this value, and exposing it invites
            // configurations that detonate mid-divergence
            mOmega0 = std::min( mOmegaMax, mOmegaPicard );

            mKernel->set_controller( this );
            if ( mKernel2 != nullptr )
            {
                mKernel2->set_controller( this );

                // the two-argument constructor is a SECOND attach path: it
                // reaches a thermal kernel without ever calling
                // set_thermal_kernel, so the eigen setup has to be repeated
                // here or it never runs on this path. Idempotent, so a caller
                // that also uses the setter pays nothing. Mirrors the
                // set_soft_fail dual-path below
                this->setup_thermal_eigen();
            }

            mTimer = new Timer();
            mRunningTimeStep = 0 ;

            // NOTE: the MUMPS error analysis is NOT armed here. It is armed
            // per timestep, around the first solve only ( see
            // arm_conditioning_magnetic ); arming it once for the whole run
            // paid for the analysis on every solve of every iteration and
            // read a single value per step. The constructor could not do it
            // anyway — set_params, which parses the key, runs later

            // arm the soft-fail contract: under the controller, a failed
            // factorization ( singular matrix at a strained iterate, e.g.
            // MUMPS INFOG(1) = -10 at quench, 2026-07-27 ) is a failed trial
            // that cuts the timestep, not a run killer
            if ( mKernel->dofmgr()->solver() != nullptr )
            {
                mKernel->dofmgr()->solver()->wrapper()->set_soft_fail( true );
            }
            if ( mKernel2 != nullptr && mKernel2->dofmgr()->solver() != nullptr )
            {
                mKernel2->dofmgr()->solver()->wrapper()->set_soft_fail( true );
            }

            // need consistent checksums for saving and loading fields
            if ( mMesh2 != nullptr && mCommRank == 0 )
            {
                mMesh2->force_checksum( mMesh->checksum() );
            }
        }

        Controller::~Controller()
        {
            delete mTimer;
        }

        const real &
        Controller::time() const
        {
            return mTime ;
        }

        const real &
        Controller::time_thermal() const
        {
            return mTime2 ;
        }

        real
        Controller::epsilon() const
        {
            return mEpsilon ;
        }

        uint
        Controller::iteration() const
        {
            return mIteration ;
        }

        void
        Controller::initialize_timestep()
        {
            // reset timer
            mTimer->reset();
            mOmegaPicard = mOmega0 ;
            mOmegaNewton = mOmega0 ;
            mOmegaPicard2 = mOmega0 ;
            mOmegaNewton2 = mOmega0 ;
            mJustPicard = false ;
            mJustPicard2 = false ;
            mForceNewton = false ;
            mNewtonEscalated = false ;
            mForceNewton2 = false ;
            mNewtonEscalated2 = false ;
            // Always start a timestep with Picard (Messe et al. 2023 §4).
            // The `mAlgorithm` / `mAlgorithmThermal` fields are preserved as
            // the terminal algorithm the promotion logic (when present) will
            // switch to once the residual is small enough.
            mEquation->set_algorithm( SolverAlgorithm::Picard ) ;
            if ( mEquation2 != nullptr )
            {
                // the thermal equation must not enter a fresh attempt with
                // the algorithm of a failed one ( leaked across retries )
                mEquation2->set_algorithm( SolverAlgorithm::Picard ) ;
            }
            // retry hygiene: these also leaked across attempts
            mThermalFrozen = false ;
            mNumIterationsDiv = 0 ;
            mDivergenceStrikes = 0 ;
            mThermalFlipCount = 0 ;
            mReset = false ;

            // reset iteration counter
            mIteration0 = mIteration ;
            mIteration02 = mIteration2 ;
            mIteration = 0 ;
            mIteration2 = 0 ;

            // start the stagnation window fresh for this timestep
            mResidualHistory.clear() ;

            // fresh watchdog trackers for this attempt
            mBestEpsilon  = BELFEM_REAL_MAX ;
            mBestEpsilon2 = BELFEM_REAL_MAX ;
            mWatchdogOmegaPrev  = BELFEM_QUIET_NAN ;
            mWatchdogOmegaPrev2 = BELFEM_QUIET_NAN ;
            mBestEpsilonIteration  = 0 ;
            mBestEpsilonIteration2 = 0 ;

            // a new attempt invalidates the Anderson mixing history
            this->anderson_clear_magnetic() ;
            this->anderson_clear_thermal() ;

            // instrument the first solve of this timestep only; a step that
            // never captures ( frozen thermal, failed first solve ) must
            // print n/a, not the PREVIOUS step's kappa
            mConditionNumbers1.fill( BELFEM_QUIET_NAN );
            mConditionNumbers2.fill( BELFEM_QUIET_NAN );

            this->arm_conditioning_magnetic() ;
            this->arm_conditioning_thermal() ;

            ++mRunningTimeStep ;

            // save point for the thermal state, taken BEFORE any field is
            // touched this step: a rejected timestep restores it across all
            // shifts, including the thermal sub-steps ( see reset_timestep )
            if ( mEquation2 != nullptr )
            {
                mEquation2->make_savepoint() ;
            }

            // update the fields
            mEquation->shift_fields() ;

            // remember old time
            mTime0 = mTime ;

            // update the time
            mTime += mDeltaTime ;

            mKernel->mesh()->time_stamp() = mTime ;
            mEquation->delta_time() = mDeltaTime ;
            mKernel->mesh()->time_step() = mMeshTimeStep ;

            //Solve the circuit problem
            // The circuit MNA solve runs on rank 0 only. Its success/failure must be
            // broadcast so that every rank takes the reset/return branch together —
            // otherwise rank 0 resets while the others march into the collective solve and
            // the run deadlocks. mCircuit is allocated on every rank, so the guard below is
            // rank-consistent and the broadcast is a matched collective.
            int tCircuitOK = 1 ;
            if ( mCircuit != nullptr && mCommRank == 0)
            {
                mCircuit->set_timestep(mDeltaTime) ;
                mCircuit->shift() ;
                mCircuit->compute_MNA_matrix() ;
                tCircuitOK = this->solve_circuit() ? 1 : 0 ;
                if ( tCircuitOK )
                {
                    index_t tCount = 0 ;
                    for ( PhysicalBoundaryCondition * tBC : mCircuitCurrentBCs )
                    {
                        tBC->fix( mCircuit->current( tCount++ ) );
                    }
                    for ( PhysicalBoundaryCondition * tBC : mCircuitVoltageBCs )
                    {
                        tBC->fix( mCircuit->voltage( tCount++ ) );
                    }
                }
            }
            if ( mCircuit != nullptr )
            {
                broadcast( tCircuitOK ) ;
            }
            if ( ! tCircuitOK )
            {
                this->reset_timestep() ;
                return ;
            }

            //update the boundary conditions
            mKernel->compute_boundary_conditions( mTime ) ;

            if( mKernel2 != nullptr )
            {
                mMesh2->set_time_step( mMeshTimeStep );
                mMesh2->time_stamp() = mTime ;

                // shift BEFORE assigning the new delta_time, so that mH( 0 )
                // records the size of the completed step ( same ordering as
                // the magnetic equation above )
                mEquation2->shift_fields() ;
                mEquation2->delta_time() = mDeltaTime ;
                mKernel2->compute_boundary_conditions( mTime ) ;
            }

            // reset epsilon
            mEpsilon = BELFEM_REAL_MAX ;
            mEpsilonAbs = BELFEM_REAL_MAX ;
            mEpsilon2 = BELFEM_REAL_MAX ;

            // print header
            if( mCommRank == 0 && gLog.info_level() > 0 )
            {
                this->print_header() ;
            }
            mFirstFlip = true ;
            mFirstFlip2 = true ;
            mNewtonTrustStreak = 0 ;      // MIT-3a: evidence does not cross timesteps

            // re-arm the thermal stagnation exit
            mThermalFlatCount = 0 ;
            mThermalStalled   = false ;
            mEpsilonAbs2      = BELFEM_REAL_MAX ;

            // an accepted step's magnetic-hit-target latch must not leak
            // into the next attempt
            mMagneticHitTarget = false ;
        }

        void
        Controller::initialize_magnetic()
        {
            // reset timer
            mTimer->reset();
            mOmegaPicard = mOmega0 ;
            mOmegaNewton = mOmega0 ;
            mJustPicard = false ;
            mForceNewton = false ;
            mNewtonEscalated = false ;
            // Always start a timestep with Picard (Messe et al. 2023 §4).
            mEquation->set_algorithm( SolverAlgorithm::Picard ) ;
            mReset = false ;

            // reset iteration counter
            mIteration0 = mIteration ;
            mIteration = 0 ;

            // start the stagnation window fresh for this timestep
            mResidualHistory.clear() ;
            mDivergenceStrikes = 0 ;

            // fresh watchdog tracker for this attempt
            mBestEpsilon = BELFEM_REAL_MAX ;
            mWatchdogOmegaPrev = BELFEM_QUIET_NAN ;
            mBestEpsilonIteration = 0 ;

            // a new attempt invalidates the Anderson mixing history
            this->anderson_clear_magnetic() ;

            // instrument the first solve of this timestep only ( stale-kappa
            // hygiene, cf. initialize_timestep )
            mConditionNumbers1.fill( BELFEM_QUIET_NAN );
            this->arm_conditioning_magnetic() ;

            ++mRunningTimeStep;

            // save point for the thermal state: the thermal problem sub-steps
            // several times within this magnetic step, and a rejected step
            // must restore all of them ( see reset_timestep )
            if ( mEquation2 != nullptr )
            {
                mEquation2->make_savepoint() ;
            }

            // update the fields
            mEquation->shift_fields() ;

            // remember old time
            mTime0 = mTime ;

            // update the time
            mTime += mDeltaTime ;

            mKernel->mesh()->time_stamp() = mTime ;
            mEquation->delta_time() = mDeltaTime ;
            mKernel->mesh()->time_step() = mMeshTimeStep ;

            //Solve the circuit problem
            // The circuit MNA solve runs on rank 0 only. Its success/failure must be
            // broadcast so that every rank takes the reset/return branch together —
            // otherwise rank 0 resets while the others march into the collective solve and
            // the run deadlocks. mCircuit is allocated on every rank, so the guard below is
            // rank-consistent and the broadcast is a matched collective.
            int tCircuitOK = 1 ;
            if ( mCircuit != nullptr && mCommRank == 0)
            {
                mCircuit->set_timestep(mDeltaTime) ;
                mCircuit->shift() ;
                mCircuit->compute_MNA_matrix() ;
                tCircuitOK = this->solve_circuit() ? 1 : 0 ;
                if ( tCircuitOK )
                {
                    index_t tCount = 0 ;
                    for ( PhysicalBoundaryCondition * tBC : mCircuitCurrentBCs )
                    {
                        tBC->fix( mCircuit->current( tCount++ ) );
                    }
                    for ( PhysicalBoundaryCondition * tBC : mCircuitVoltageBCs )
                    {
                        tBC->fix( mCircuit->voltage( tCount++ ) );
                    }
                }
            }
            if ( mCircuit != nullptr )
            {
                broadcast( tCircuitOK ) ;
            }
            if ( ! tCircuitOK )
            {
                this->reset_timestep() ;
                return ;
            }

            //update the boundary conditions
            mKernel->compute_boundary_conditions( mTime ) ;

            // reset epsilon
            mEpsilon = BELFEM_REAL_MAX ;
            mEpsilonAbs = BELFEM_REAL_MAX ;

            // print header
            if( mCommRank == 0 && gLog.info_level() > 0 )
            {
                this->print_header() ;
            }
            mFirstFlip = true ;
            mNewtonTrustStreak = 0 ;      // MIT-3a

        }

        void
        Controller::initialize_thermal()
        {

            mDeltaTime2 = mDeltaTime/mCouplingFactor ;

            // reset timer
            mTimer->reset();
            mOmegaPicard2 = mOmega0 ;
            mOmegaNewton2 = mOmega0 ;
            mJustPicard2 = false ;
            mForceNewton2 = false ;
            mNewtonEscalated2 = false ;
            // Always start a timestep with Picard (Messe et al. 2023 §4).
            mEquation2->set_algorithm( SolverAlgorithm::Picard ) ;
            mResetThermal = false ;

            // reset iteration counter
            mIteration02 = mIteration2 ;
            mIteration2 = 0 ;

            // every thermal sub-step is a fresh fixed-point problem ( the
            // magnetic state moved ): the mixing history must not carry over
            this->anderson_clear_thermal() ;

            // fresh watchdog tracker for this sub-step
            mBestEpsilon2 = BELFEM_REAL_MAX ;
            mWatchdogOmegaPrev2 = BELFEM_QUIET_NAN ;
            mBestEpsilonIteration2 = 0 ;

            // instrument the first solve of this sub-step only ( stale-kappa
            // hygiene, cf. initialize_timestep )
            mConditionNumbers2.fill( BELFEM_QUIET_NAN ) ;

            this->arm_conditioning_thermal() ;

            // update the fields. Shifting at EVERY thermal sub-step is
            // required: the BDF history must hold the previous sub-step at
            // delta_time2 spacing, otherwise the time derivative is wrong by
            // up to the coupling factor. A rejected sub-step is undone by
            // reset_fields ( one shift, one un-shift, see reset_thermal );
            // a rejected MAGNETIC step spans several sub-step shifts and is
            // undone by restoring the savepoint taken in initialize_timestep
            // / initialize_magnetic ( see reset_timestep ).
            mEquation2->shift_fields() ;

            // remember old time
            mTime02 = mTime2 ;

            // update the time
            mTime2 += mDeltaTime2 ;

            mKernel2->mesh()->time_stamp() = mTime2 ;
            mEquation2->delta_time() = mDeltaTime2 ;

            //update the boundary conditions
            mKernel2->compute_boundary_conditions( mTime2 ) ;

            // reset epsilon
            mEpsilon2 = BELFEM_REAL_MAX ;
            mEpsilonAbs2 = BELFEM_REAL_MAX ;

            mFirstFlip2 = true ;

        }

        void
        Controller::impose_voltage_bcs()
        {
            // Voltage/current boundary conditions are patched onto the RHS on rank 0
            // only, after assembly and before the solve. The current rows always come
            // first and are already fixed ( so they are skipped ); the voltage rows are
            // free and get value*dt added.
            if ( mCommRank == 0 )
            {
                // get the abstract nodes
                Cell< Dof * > & tAbstractDofs = mKernel->dofmgr()->abstract_dofs();

                // walk the CURRENT conditions first: the factory creates them
                // ahead of the voltage ones, and they are already fixed, so this
                // loop only advances tCount past them to reach the voltage rows
                index_t tCount = 0;
                for (PhysicalBoundaryCondition * tBC : mKernel->boundary_conditions())
                {
                    if(tBC->type()==BoundaryConditionType::Current || tBC->type() == BoundaryConditionType::CircuitCurrent)
                    {
                        BELFEM_ASSERT( tAbstractDofs( tCount )->is_fixed(), "Dof %lu of node %lu is expected to be fixed but is free",
                            ( long unsigned int ) tAbstractDofs( tCount )->id(), ( long unsigned int ) tAbstractDofs( tCount )->node()->id() ) ;

                        ++tCount ; //Current already fixed, we skip it
                    }
                }

                // get the RHS of the system
                Vector< real > & tRHS = mKernel->dofmgr()->rhs_vector() ;

                for (PhysicalBoundaryCondition * tBC : mKernel->boundary_conditions())
                {
                    if(tBC->type()==BoundaryConditionType::Voltage || tBC->type() == BoundaryConditionType::CircuitVoltage)
                    {
                        BELFEM_ASSERT( ! tAbstractDofs( tCount )->is_fixed(), "Dof %lu of node %lu is expected to be free but is fixed",
                            ( long unsigned int ) tAbstractDofs( tCount )->id(), ( long unsigned int ) tAbstractDofs( tCount )->node()->id() ) ;

                        tRHS( tAbstractDofs( tCount++ )->index()) += tBC->value()*mDeltaTime ;
                    }
                }
            }
        }

        void
        Controller::dump_system_if_requested(
                DofManager * aDofMgr,
                const char * aTag,
                int        & aCount )
        {
            if ( aDofMgr == nullptr ) return ;

            // env-gated assembly dump. Off unless BELFEM_DUMP_SYSTEM is set,
            // so production paths are untouched. Writes the ASSEMBLED system
            // before the solve, which is what the kernel-collapse equivalence
            // check compares.
            if ( std::getenv( "BELFEM_DUMP_SYSTEM" ) == nullptr )
            {
                return ;
            }

            // optional step target, so a LIVE arm and a RESTORE arm can be
            // diffed at the SAME timestep ( BELFEM_DUMP_SYSTEM_STEP=<n> ).
            // <n> counts like the step headers and the warm-restart banner's
            // "next step" line: the value the first assembly after a restart
            // carries. Unset = dump the first four assemblies.
            const char * tWanted = std::getenv( "BELFEM_DUMP_SYSTEM_STEP" );
            bool tStepOk = true ;
            if ( tWanted != nullptr )
            {
                // a mistyped probe must fail loudly, not silently dump
                // nothing ( atoi would fold empty, non-numeric and negative
                // text into a wrong-but-valid step number )
                BELFEM_ERROR( *tWanted != '\0' &&
                    std::strspn( tWanted, "0123456789" ) == std::strlen( tWanted ),
                    "invalid BELFEM_DUMP_SYSTEM_STEP '%s' - expect a plain nonnegative integer",
                    tWanted );

                tStepOk = mRunningTimeStep
                    == ( uint ) std::strtoul( tWanted, nullptr, 10 );
            }

            // one budget PER FIELD. It used to be a single function-local
            // static shared by every path, so whichever field assembled first
            // spent all four dumps and the other never appeared in the output
            if ( aCount < 4 && tStepOk )
            {
                aDofMgr->save_system(
                    "sysdump_" + string( aTag ) + "_"
                        + std::to_string( aCount++ ) + ".hdf5" );
            }
        }

        void
        Controller::anderson_commit_magnetic()
        {
            // only a Picard solve stages a pair; committing during a Newton
            // iterate would advance the counter over an empty history and
            // wrongly suppress the Newton omega growth ( the algorithm state
            // is rank-consistent, so this guard is too )
            if ( mAndersonDepth == 0
                 || mEquation->algorithm() != SolverAlgorithm::Picard ) return ;
            mKernel->dofmgr()->solver_data()->anderson_commit() ;

            // rank-consistent window-fill counter, see the member comment.
            // Upper bound only: a failed mixing solve does not stage ( O3 ),
            // so the counter may over-report — the growth hold then errs on
            // the conservative side, identically on every rank.
            mAndersonCommits = std::min( mAndersonCommits + 1, mAndersonDepth );
        }

        void
        Controller::anderson_discard_magnetic()
        {
            if ( mAndersonDepth == 0 ) return ;
            mKernel->dofmgr()->solver_data()->anderson_discard() ;
        }

        void
        Controller::anderson_clear_magnetic()
        {
            if ( mAndersonDepth == 0 ) return ;
            mKernel->dofmgr()->solver_data()->anderson_clear() ;
            mAndersonCommits = 0 ;
        }

        void
        Controller::anderson_commit_thermal()
        {
            // see anderson_commit_magnetic for the algorithm guard rationale
            if ( mAndersonDepth2 == 0 || mKernel2 == nullptr
                 || mEquation2->algorithm() != SolverAlgorithm::Picard ) return ;
            mKernel2->dofmgr()->solver_data()->anderson_commit() ;
            mAndersonCommits2 = std::min( mAndersonCommits2 + 1, mAndersonDepth2 );
        }

        void
        Controller::anderson_clear_thermal()
        {
            if ( mAndersonDepth2 == 0 || mKernel2 == nullptr ) return ;
            mKernel2->dofmgr()->solver_data()->anderson_clear() ;
            mAndersonCommits2 = 0 ;
        }

//------------------------------------------------------------------------------

        //! format the residual into a fixed 8-character column: plain fixed
        //! point below 10, scientific above. The old min( eps, 9.0 ) clamp
        //! masked blown-up residuals as a stall at 9.000000 ( ts9 trace:
        //! the true residual was 11-23, above the divergence guard )
        static string
        residual_string( const real aEpsilon )
        {
            return aEpsilon < 10.0 ?
                sprint( "%8.6f", aEpsilon ) : sprint( "%8.2e", aEpsilon );
        }

        //! condition numbers into a fixed SIX-character cell, in compact
        //! scientific notation: "1.43e8" rather than C's "1.43e+08". The
        //! reader compares kappa against 1/eps_machine, and that arithmetic
        //! is instant in this form. The mantissa loses a digit once the
        //! exponent needs two, so the footer column never shifts; kappa >= 1
        //! by definition, so a negative exponent only appears if a solver
        //! reports nonsense. NaN prints as "n/a" ( field absent, or its
        //! solver cannot supply kappa cheaply )
        static string
        conditioning_string( const real aKappa )
        {
            if ( ! std::isfinite( aKappa ) || aKappa <= 0.0 )
            {
                return "   n/a" ;
            }

            int  tExponent = ( int ) std::floor( std::log10( aKappa ) );
            real tMantissa = aKappa / std::pow( 10.0, ( real ) tExponent );

            // the printf precision depends on the exponent's digit count, and
            // the rounding-carry threshold depends on the precision ( 9.96
            // under %.1f prints "10.0" and breaks the 6-char cell — Codex+
            // Grok blind agreement ); a carry can also grow the exponent into
            // the next digit class, hence the second pass
            uint tPrecision = 0 ;
            for ( uint tPass = 0; tPass < 2; ++tPass )
            {
                tPrecision = ( tExponent >= 0 && tExponent < 10 )   ? 2 :
                             ( tExponent > -10 && tExponent < 100 ) ? 1 : 0 ;
                const real tCarry = tPrecision == 2 ? 9.995 :
                                    tPrecision == 1 ? 9.95  : 9.5 ;
                if ( tMantissa >= tCarry )
                {
                    tMantissa = 1.0 ;
                    ++tExponent ;
                    continue ;
                }
                break ;
            }

            string tValue =
                tPrecision == 2 ? sprint( "%.2fe%d", tMantissa, tExponent ) :
                tPrecision == 1 ? sprint( "%.1fe%d", tMantissa, tExponent ) :
                                  sprint( "%.0fe%d", tMantissa, tExponent );

            return sprint( "%6s", tValue.c_str() );
        }

        bool
        Controller::try_escalate_to_newton()
        {
            // Picard-breakdown rescue. The lagged-conductivity Picard map loses contraction
            // through the HTS E-J transition: the magnetic residual floors far above
            // mEpsilonSwitch, so the normal Picard->Newton promotion ( iterate_coupled /
            // iterate_magnetic, Messe et al. 2023 Eq. 13 ) never fires and the controller
            // would otherwise cut Δt indefinitely without the residual responding. Before
            // resetting, give the configured Newton tangent one chance at the current step
            // size -- it is often contractive where Picard is not. Escalate at most once per
            // timestep attempt ( mNewtonEscalated ); if Newton also breaks down the caller
            // resets as before.
            if ( mAlgorithm == SolverAlgorithm::NewtonRaphson
                 && ! mNewtonEscalated
                 && mEquation->algorithm() == SolverAlgorithm::Picard
                 && ! std::isnan( mEpsilon )
                 // never hand a blown-up iterate to Newton: escalation exists
                 // for Picard stalling at a LOW floor ( the -55 dB freeze
                 // class ), not for divergence. Above the 1E1 guard the
                 // timestep is too large and must be cut -- and since the
                 // escalation sets mForceNewton, which disables the +10 dB
                 // cut, escalating there would leave Newton flailing with no
                 // exit ( ts9 trace: 75 iterations at +11 dB )
                 && mEpsilon < 1E1 )
            {
                mNewtonEscalated = true ;
                mForceNewton     = true ;
                mJustPicard      = false ;
                mFirstFlip       = false ;
                // MIT-3a: escalation bypasses the normal promotion reset,
                // so the contraction evidence must be cleared here too
                mNewtonTrustStreak = 0 ;

                // Start Newton from the configured relaxation, not the collapsed omega the
                // Picard line search leaves behind ( which would freeze the Newton step ).
                mOmegaNewton = mOmega0 ;
                mEquation->set_algorithm( SolverAlgorithm::NewtonRaphson ) ;

                // the abandoned Picard residuals must not pollute the stall window
                mResidualHistory.clear() ;

                // the algorithm changed: the Anderson history is void
                this->anderson_clear_magnetic() ;

                // restart the watchdog clock: the escalated Newton phase gets
                // a full window to prove itself against the CURRENT residual,
                // not against a stale best from before the breakdown
                mBestEpsilon = mEpsilon ;
                mBestEpsilonIteration = mIteration ;

                if( mCommRank == 0 && gLog.info_level() > 0 )
                {
                    std::cout << sprint( "   │%-71s│",
                        " Picard stalled: escalating to Newton before cutting the timestep" )
                              << std::endl ;
                }
                return true ;
            }
            return false ;
        }

        bool
        Controller::try_escalate_thermal_to_newton()
        {
            // Thermal twin of try_escalate_to_newton, for the Picard freeze: the
            // thermal Picard residual can floor bit-flat ABOVE mEpsilonSwitch2
            // ( observed at -39.9 dB for four iterates on the tapestack3d
            // reproducer, and earlier at -39.5 dB for 100+ ), and since the
            // coupled Picard->Newton promotion fires only below the switch,
            // Newton is gated behind progress that only Newton can make. Give
            // the configured tangent one chance per timestep attempt. Note the
            // reproducer also showed the switch itself is roundoff-sensitive:
            // 8x2 ranks landed at 1.00e-4 and promoted, 4x4 at 1.02e-4 and
            // froze -- the decomposition selected the algorithm.
            if ( mAlgorithmThermal == SolverAlgorithm::NewtonRaphson
                 && ! mNewtonEscalated2
                 && mEquation2->algorithm() == SolverAlgorithm::Picard
                 && ! std::isnan( mEpsilon2 )
                 // same guard as the magnetic twin: escalation exists for a
                 // LOW flat floor, not for divergence -- above the 1E1 guard
                 // the timestep must be cut instead, and mForceNewton2 would
                 // disable that exit
                 && mEpsilon2 < 1E1 )
            {
                mNewtonEscalated2 = true ;
                mForceNewton2     = true ;
                mJustPicard2      = false ;
                mFirstFlip2       = false ;

                // start Newton from the configured relaxation, not the omega
                // the flat Picard phase left behind, cf. the magnetic twin
                mOmegaNewton2 = mOmega0 ;
                mEquation2->set_algorithm( SolverAlgorithm::NewtonRaphson ) ;

                // the algorithm changed: the thermal mixing history is void
                this->anderson_clear_thermal() ;

                // restart the thermal watchdog on the CURRENT residual, cf.
                // the magnetic twin
                mBestEpsilon2 = mEpsilon2 ;
                mBestEpsilonIteration2 = mIteration2 ;

                if( mCommRank == 0 && gLog.info_level() > 0 )
                {
                    std::cout << sprint( "   │%-71s│",
                        " thermal Picard is flat: escalating to Newton" )
                              << std::endl ;
                }
                return true ;
            }
            return false ;
        }

        bool
        Controller::watchdog_magnetic( const real aOmega )
        {
            // the spare clause below compares against the omega of the PREVIOUS
            // watchdog call; every return path must refresh it AFTER the test,
            // or one spared iterate blinds the detector for the rest of the
            // attempt ( three-AI round 2026-08-21, Grok R3 )
            const real tOmegaPrev = mWatchdogOmegaPrev ;
            mWatchdogOmegaPrev = aOmega ;

            if ( mEpsilon < mBestEpsilon )
            {
                mBestEpsilon = mEpsilon ;
                mBestEpsilonIteration = mIteration ;
                return false ;
            }
            if ( mWatchdogWindow == 0
                 || mEpsilon <= 10.0 * mRelativeEpsilonTarget
                 || mIteration <= mMinNumIterations
                 || mIteration - mBestEpsilonIteration < mWatchdogWindow )
            {
                return false ;
            }

            // SPARE a stalled-best iterate while the line search is regrowing
            // its relaxation: a no-new-best window that merely spans an
            // overshoot/backtrack/recovery cycle is not a stall. Measured
            // ( coarse tapestack3d, 2026-08-21 ): the AIMD recovery from a
            // backtracked omega ~0.1 takes ~17-22 iterates at beta = 1.1,
            // window W = 8 cut it mid-climb and CASCADED ( halving improved
            // the best residual and was cut again, 2.10 -> 1.05 -> 0.53 ms ).
            // A genuine stall has omega pinned or shrinking ( the hd floor
            // grinds, and 383@2.10 with omega frozen at 0.053 ), so it still
            // fires. The motivation parallels Chamberlain-Powell-Lemarechal's
            // watchdog technique -- do not punish a nonmonotone step that is
            // still contracting -- though theirs relaxes a line search and
            // this spares a timestep cut. NaN prev ( cold start ) must not
            // spare.
            if ( ! std::isnan( tOmegaPrev ) && aOmega > tOmegaPrev )
            {
                return false ;
            }

            if ( mCommRank == 0 && gLog.info_level() > 0 )
            {
                std::cout << sprint( "   │%-71s│",
                    " Watchdog: no residual improvement, cutting the timestep" )
                          << std::endl ;
            }
            return true ;
        }

        bool
        Controller::watchdog_thermal( const real aOmega )
        {
            // mirror of watchdog_magnetic, cf. the comments there
            const real tOmegaPrev = mWatchdogOmegaPrev2 ;
            mWatchdogOmegaPrev2 = aOmega ;

            if ( mEpsilon2 < mBestEpsilon2 )
            {
                mBestEpsilon2 = mEpsilon2 ;
                mBestEpsilonIteration2 = mIteration2 ;
                return false ;
            }
            if ( mWatchdogWindow2 == 0
                 || mEpsilon2 <= 10.0 * mRelativeEpsilonTarget2
                 || mIteration2 <= mMinNumIterations2
                 || mIteration2 - mBestEpsilonIteration2 < mWatchdogWindow2 )
            {
                return false ;
            }
            if ( ! std::isnan( tOmegaPrev ) && aOmega > tOmegaPrev )
            {
                return false ;
            }
            if ( mCommRank == 0 && gLog.info_level() > 0 )
            {
                std::cout << sprint( "   │%-71s│",
                    " Watchdog: no thermal improvement, cutting the timestep" )
                          << std::endl ;
            }
            return true ;
        }

//------------------------------------------------------------------------------

        bool
        Controller::magnetic_stagnation_forces_reset()
        {
            // Once the window is full and we are still above tolerance, a flat residual
            // ( mean absolute deviation below mStallBand dB over the last mStallWindow
            // steps ) means the iteration is stuck. If Newton has stalled, fall back to
            // Picard for the rest of the timestep ( Messe 2023 hybrid ); if Picard -- the
            // last resort -- has stalled, the timestep is too large and the caller cuts it.
            if ( mResidualHistory.full() && mEpsilon > mRelativeEpsilonTarget )
            {
                real tMean = 0.0 ;
                for ( const real & tValue : mResidualHistory ) tMean += tValue ;
                tMean /= mResidualHistory.size() ;

                real tDeviation = 0.0 ;
                for ( const real & tValue : mResidualHistory ) tDeviation += std::abs( tValue - tMean ) ;
                tDeviation /= mResidualHistory.size() ;

                if ( tDeviation < mStallBand )
                {
                    if ( mEquation->algorithm() == SolverAlgorithm::NewtonRaphson )
                    {
                        // Newton stalled -> latch to Picard for the rest of the timestep,
                        // dropping the Newton-era samples. After a PROMOTION-Newton
                        // stall ( near the residual floor ), Picard resumes from ITS
                        // OWN last relaxation: the collapsed omega of a stalled Newton
                        // reflects the tangent's failure, not Picard's contraction
                        // history ( ts39: importing 0.03 into a Picard that had earned
                        // omega 1.0 cost ~15 crawl iterations ). The resume enters at
                        // HALF throttle like every other algorithm switch — a full-
                        // omega resume from the Newton-stalled iterate can bounce and
                        // burn the backtracking budget into a timestep cut ( ts66 );
                        // the growth rule earns the other half back. After an
                        // ESCALATED-Newton stall the stale Picard omega predates a
                        // genuine Picard breakdown, so the conservative min-merge
                        // stays ( Codex DQ2, uncoupled path has no line search )
                        if ( mForceNewton )
                        {
                            mOmegaPicard = std::max( std::min( mOmegaNewton, mOmegaPicard ), mOmegaMin );
                        }
                        else
                        {
                            mOmegaPicard = std::max( 0.5 * mOmegaPicard, mOmegaMin );
                        }
                        mJustPicard = true ;
                        mForceNewton = false ;
                        mNewtonTrustStreak = 0 ;   // MIT-3a: demotion voids the evidence
                        mEquation->set_algorithm( SolverAlgorithm::Picard ) ;
                        mResidualHistory.clear() ;
                        this->anderson_clear_magnetic() ;
                    }
                    else
                    {
                        // Picard stalled. Try the Newton tangent once before declaring the
                        // timestep too large ( see try_escalate_to_newton ); only if Newton
                        // is unavailable or already spent does the caller cut the timestep.
                        return ! this->try_escalate_to_newton() ;
                    }
                }
            }
            return false ;
        }

        void
        Controller::solve_coupled()
        {
            // certified exit: one call runs the complete nonlinear
            // loop of a timestep. Each trip measures the committed state
            // FIRST ( head ), only then decides whether to solve ( body ) —
            // so the state a timestep commits is always the state whose
            // residual was measured, per field. The driver retries the
            // timestep when mReset comes back set, exactly as before.
            mTripExit = false ;
            mMagBodyRanLastTrip = true ;

            while ( true )
            {
                this->iterate_coupled() ;

                if ( mReset || mTripExit )
                {
                    return ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Controller::iterate_coupled()
        {
            // backup fields
            const Cell< string > & tFieldLabels = mKernel->dofmgr(  )->iwg()->dof_fields() ;

            Cell< Dof * > & tDofs = mKernel->dofmgr(  )->dofs() ;

            uint tNumFields = tFieldLabels.size();
            index_t tNumDofs = tDofs.size() ;

            if ( mBackupFields.size() == 0 )
            {
                mBackupFields.set_size(  tNumFields, {} );
                mBackupDofValues.set_size( tNumDofs );
            }

            // Picard <-> Newton handoff (Messe et al. 2023, Eq. 13): once the
            // residual falls below mEpsilonSwitch, hand off from Picard to the configured
            // terminal Newton algorithm (and back if it rises again), carrying the
            // relaxation across the switch. The mJustPicard latch (set by the stagnation
            // guard below when Newton stalls) pins the solver to Picard for the rest of
            // the timestep.
            const SolverAlgorithm tAlgorithm0 = mEquation->algorithm() ;

            // skip gating: on a trip whose predecessor ran NO magnetic
            // body ( certified-magnetic, thermal still iterating ) the
            // residual did not come from a magnetic step — a promotion here
            // would wipe Anderson history and re-anchor the watchdog on
            // coupling drift ( round-3 C5 ). The handoff idles until a
            // magnetic body runs again.
            if ( mMagBodyRanLastTrip ) {
            // mForceNewton ( set by try_escalate_to_newton on a Picard breakdown ) pins the
            // solver to Newton regardless of residual magnitude; the normal promotion still
            // fires on its own once the residual drops below mEpsilonSwitch. While
            // mJustPicard is latched ( a stalled or line-search-exhausted Newton fell back )
            // the promotion must not even be ATTEMPTED: the flip below would re-anchor the
            // watchdog on every latched iterate -- refreshing the best-residual clock and
            // disabling the very guard that has to end a fruitless Picard grind
            if ( mAlgorithm == SolverAlgorithm::NewtonRaphson && ! mJustPicard
                 && ( mForceNewton || ( mEpsilon < mEpsilonSwitch && mIteration > 1 ) ) )
            {
                if ( mEquation->algorithm() == SolverAlgorithm::Picard )
                {
                    // enter Newton damped, but only on the FIRST promotion of
                    // this attempt: the first full step after a promotion
                    // routinely overshoots from a nearly-converged iterate ( the
                    // reproducible it-3 kick, ts34 trace ), while re-damping on
                    // every re-promotion would ratchet omega down x0.5 per
                    // promote/demote cycle ( Codex+Grok audit ). omega recovers
                    // through the growth rule once the tangent proves
                    // contractive. The escalation path is untouched --
                    // try_escalate_to_newton sets the algorithm directly and
                    // never passes through this flip
                    real tDamp = mFirstFlip ? 0.5 : 1.0 ;
                    mOmegaNewton = std::clamp( tDamp * mOmegaPicard, mOmegaMin, mOmegaMax );
                    mFirstFlip = false ;
                    mNewtonTrustStreak = 0 ;   // MIT-3a: a promotion restarts the evidence

                    // re-anchor the progress watchdog on the promotion, like
                    // the escalation path already does: Newton must prove
                    // itself against ITS OWN entry residual, not against the
                    // Picard best it starts above ( the it-3 entry kick would
                    // otherwise pre-charge the watchdog clock — sidecoatings
                    // cut at exactly 2 + window )
                    mBestEpsilon = mEpsilon ;
                    mBestEpsilonIteration = mIteration ;
                }
                mEquation->set_algorithm( SolverAlgorithm::NewtonRaphson );
            }
            else
            {
                if ( mEquation->algorithm() == SolverAlgorithm::NewtonRaphson && ! mFirstFlip )
                {
                    mOmegaPicard = std::clamp( mOmegaNewton, mOmegaMin, mOmegaMax );
                    mNewtonTrustStreak = 0 ;   // MIT-3a: ordinary demotion voids it
                }
                mEquation->set_algorithm( SolverAlgorithm::Picard );
            }
            if ( mJustPicard )
            {
                mEquation->set_algorithm( SolverAlgorithm::Picard ) ;
            }
            // the stagnation window must hold samples from a single algorithm only
            if ( mEquation->algorithm() != tAlgorithm0 )
            {
                mResidualHistory.clear() ;
                this->anderson_clear_magnetic() ;
            }
            } // end handoff gate ( mMagBodyRanLastTrip )

            // capture the algorithm that runs THIS iteration's solves: a same-call
            // escalation ( try_escalate_to_newton below ) can flip mEquation to
            // Newton after the solve, which must not re-route the omega store or
            // the trust test in the adaptation tail ( Codex+Grok audit )
            const bool tRanNewton = mEquation->algorithm() == SolverAlgorithm::NewtonRaphson ;
            real & tOmega = tRanNewton ? mOmegaNewton : mOmegaPicard ;
            tOmega = std::clamp( tOmega, mOmegaMin, mOmegaMax ) ;

            // AIMD / line-search reference pair: the PREVIOUS trip's
            // residual, captured before the head can overwrite the member on
            // a skip trip. For Picard the body result equals the head value
            // ( pre-update semantics ), so pairing against the head would
            // compare a number with itself and freeze the omega adaptation
            const real tEpsilonPrev    = mEpsilon ;
            const real tEpsilonAbsPrev = mEpsilonAbs ;

            // line-search trial counter: declared at function scope — the
            // Newton trust streak in the adaptation tail consumes it
            uint tBacktracks = 0 ;

            // thermal relaxation store; selected per algorithm in the gated
            // thermal block below ( set there, consumed by the adaptation
            // at the tail and the print )
            real * tOmega2Ptr = nullptr ;
            real tOmega2 = mEquation2 != nullptr &&
                mEquation2->algorithm() == SolverAlgorithm::NewtonRaphson ?
                    mOmegaNewton2 : mOmegaPicard2 ;

            // ---------- MAGNETIC HEAD ----------
            // assemble at the COMMITTED state and measure its residual under
            // the fresh operator, WITHOUT solving. This is the number the
            // exit certifies: the state a timestep commits is the state whose
            // residual was measured. The assembly is not wasted — the first
            // body trial consumes it.
            mEquation->set_omega( tOmega );

            this->reset_dotQ();
            mKernel->dofmgr()->compute_jacobian_and_rhs();
            this->collect_dotQ();

            this->dump_system_if_requested( mKernel->dofmgr(), "magnetic",
                                            mDumpCountMagnetic );

            this->impose_voltage_bcs();

            mKernel->dofmgr()->compute_residual();

            // rank-uniform after the broadcast inside residual()
            const real tEpsilonHead    = mKernel->dofmgr()->residual( mIteration );
            const real tEpsilonHeadAbs = mKernel->dofmgr()->absolute_residual();

            // a NaN head means the committed state itself is broken — cut
            // the timestep before anything consumes it ( round-3 N1: the
            // decision uses the BROADCAST scalar only, never mRhsNorm )
            if ( std::isnan( tEpsilonHead ) )
            {
                this->reset_timestep() ;
                return ;
            }

            const bool tConvergedM = tEpsilonHead    <= mRelativeEpsilonTarget
                                  || tEpsilonHeadAbs <= mAbsoluteEpsilonTarget ;

            // per-attempt latch for the thermal budget cut: the
            // head measurement is an honest feed — sticky, so the post-body
            // latch below stays as belt and braces
            if ( tConvergedM )
            {
                mMagneticHitTarget = true ;
            }

            // the body runs while the state is uncertified or the deck's
            // minimum body count is unmet. min iterations counts BODY solves:
            // a no-body trip can only certify, never iterate
            const bool tSolveM = ( ! tConvergedM )
                              || mIteration < mMinNumIterations ;

            mMagBodyRanLastTrip = tSolveM ;

            if ( ! tSolveM )
            {
                // the certificate becomes the outer residual: the thermal
                // gate, the stall test and the budget latch all judge the
                // honest committed-state value
                mEpsilon    = tEpsilonHead ;
                mEpsilonAbs = tEpsilonHeadAbs ;
            }

            // ---------- MAGNETIC BODY ----------
            const real tEpsilonAbs0 = tEpsilonAbsPrev ;

            if ( tSolveM )
            {

            // backup fields ( restore target of the Newton line search )
            for ( uint f=0; f<tNumFields; ++f )
            {
                mBackupFields( f ) = mMesh->field_data( tFieldLabels( f ) ) ;
            }
            for ( index_t k=0; k<tNumDofs; ++k )
            {
                mBackupDofValues( k ) = tDofs( k )->value() ;
            }

            // the AIMD / exhausted-restore reference stays the PREVIOUS
            // trip's residual ( today's pairing ); Newton's line search
            // replaces its reference with the fresh pre-update value below
            mEpsilon0 = tEpsilonPrev ;

            // the first trial consumes the HEAD assembly and residual —
            // retrials reassemble at the restored state as before
            bool tFirstTrial = true ;

            bool tRun = true ;

            // baseline for the backtracking line search. Seeded from the last
            // accepted iterate's epsilon ( BELFEM_REAL_MAX on the first
            // iteration, reset in initialize_*, so trial 1 always passes ) and
            // REPLACED for Newton trials by the pre-update residual of the
            // CURRENT assembly the moment it is known ( lagged-A metric fix,
            // 2026-08-10 ): mEpsilon0 was measured under the PREVIOUS
            // assembly, and whenever the two assemblies disagree by more than
            // the acceptance band, the stale reference rejected all eight
            // trials flat in omega -- the omega -> 0 limit of any trial is
            // exactly the pre-update residual, so against the honest
            // reference a small enough step is always acceptable and the
            // search degrades gracefully instead of deadlocking ( greg5 out2;
            // matfix survived the same tangent only because its unconditional
            // moved-baseline rescue accepted these flat rejects )
            real tLogEpsilonRef = std::log10( mEpsilon0 ) ;
            real tLogEpsilon ;

            // residual and omega of the previously rejected trial and the count
            // of consecutive flat reject pairs, for the moved-baseline detector
            real tLogEpsilonPrev = BELFEM_QUIET_NAN ;
            real tOmegaTrialPrev = BELFEM_QUIET_NAN ;
            uint tFlatRejects = 0 ;

            while ( tRun )
            {
                // set relaxation parameters
                mEquation->set_omega( tOmega );

                if ( ! tFirstTrial )
                {
                    this->reset_dotQ();
                    mKernel->dofmgr()->compute_jacobian_and_rhs();
                    this->collect_dotQ();

                    this->dump_system_if_requested( mKernel->dofmgr(), "magnetic",
                                                    mDumpCountMagnetic );

                    this->impose_voltage_bcs();

                    mKernel->dofmgr()->compute_residual();
                }
                tFirstTrial = false ;

                mKernel->dofmgr()->solve_from_residual();

                // a failed factorization ( singular matrix, soft-fail
                // contract ) cannot be cured by relaxation: the matrix is
                // assembled at the ACCEPTED iterate, not the trial, and the
                // fields were not touched. Cut the timestep -- but a
                // persistently singular system must not cut forever
                if ( mKernel->dofmgr()->solve_failed() )
                {
                    BELFEM_ERROR( ++mSolverFailCount < 8,
                        "linear solver failed on %u consecutive timestep attempts - the cause is in the solver messages above ( e.g. a singular system, non-convergence, or repeated workspace exhaustion )",
                        ( unsigned int ) mSolverFailCount );
                    this->reset_timestep() ;
                    return ;
                }

                // get the residual
                mEpsilon = mKernel->dofmgr()->residual( mIteration );
                mEpsilonAbs = mKernel->dofmgr()->absolute_residual();
                tLogEpsilon = std::log10( mEpsilon ) ;

                // COND1 / COND2 of the magnetic system -- the Arioli-Demmel-Duff
                // pair, NOT kappa: sampled at the FIRST SOLVE of
                // the timestep and disarmed right after ( no-op unless
                // "mumps error analysis" is on and the solver is MUMPS --
                // NOT the eigen flag; the two diagnostics split 2026-08-30 ).
                // Iteration 0 is structurally Picard ( promotion requires
                // mIteration > 1 ) and Picard accepts its first trial, so
                // this is exactly one instrumented solve — the tBacktracks
                // guard hardens the contract against future promotion-rule
                // changes rather than covering a reachable path today
                if ( mIteration == 0 && tBacktracks == 0 )
                {
                    this->capture_conditioning_magnetic() ;
                }

                // honest Newton reference ( see the declaration comment ):
                // the pre-update residual is identical across the trials of
                // one iteration -- every retry reassembles at the restored
                // entry state -- so overwriting it per trial is idempotent
                if ( tRanNewton )
                {
                    const real tEpsPre = mKernel->dofmgr()->pre_update_residual() ;
                    if ( std::isfinite( tEpsPre ) && tEpsPre > 0.0 )
                    {
                        tLogEpsilonRef = std::log10( tEpsPre ) ;
                    }
                }

                // Picard runs WITHOUT an intra-iterate line search ( Messe
                // et al. 2023 §4, Eq. 14: damp on regression and continue ).
                // Within one frozen assembly a Picard trial cannot be judged:
                // the pre-update epsilon measures the ENTRY state -- the very
                // state the backup holds -- and the lagged post-update
                // residual of a relaxed Picard step is identically
                // (1-omega)*r, pure omega arithmetic. A reject loop therefore
                // re-solves the same system blind, restores the state it just
                // condemned, and exits by copying the stale residual ( the
                // frozen bit-identical print + omega collapse that froze the
                // greg3 run at the flux front ). Regression handling for
                // Picard is the alpha branch of the adaptation tail plus the
                // divergence / stagnation / watchdog guards.
                if ( ! tRanNewton )
                {
                    this->anderson_commit_magnetic() ;
                    tRun = false ;
                }
                // Newton keeps the line search: its post-update residual is
                // non-degenerate ( J != A ). Accept the trial unless it is
                // significantly worse than the last accepted iterate: within
                // 0.3 decade ( ~2x ) always, and below ~+8 dB absolute a
                // regression of up to one decade is tolerated ( the
                // early/mid-settling wander ). What must NOT pass is a
                // multi-decade kick from a deep reference: a full Newton step
                // from a -50 dB iterate that jumps to -30 dB used to sail
                // through the plain absolute clause and could spiral the whole
                // timestep ( ts34 trace, dl20260727 )
                else if ( tLogEpsilon < tLogEpsilonRef + 0.3
                     || ( tLogEpsilon < 0.8 && tLogEpsilon < tLogEpsilonRef + 1.0 ) )
                {
                    // accepted: the staged Anderson pair joins the history
                    this->anderson_commit_magnetic() ;
                    tRun = false ;
                }
                else if ( mKernel2 != nullptr
                          && ! mThermalFrozen
                          && tFlatRejects >= 1
                          && tOmega < 0.6 * tOmegaTrialPrev
                          && std::abs( tLogEpsilon - tLogEpsilonPrev ) < 0.05 )
                {
                    // the rescue below only makes sense when a thermal
                    // baseline EXISTS and moved since the reference was
                    // recorded ( a thermal solve ran last iteration ): on a
                    // magnetic-only run an omega-independent flat regression
                    // is a defect signature, not baseline drift, and this
                    // branch was observed laundering a 44 dB kick there
                    // ( Garber R9 trace, dl20260806 )
                    //
                    // the regression did NOT shrink over two consecutive omega
                    // halvings. A genuine step overshoot vanishes as omega -> 0;
                    // an omega-independent regression means the residual
                    // BASELINE moved under the magnetic system -- near quench
                    // the staggered thermal update shifts rho(T) enough that
                    // the stale reference is a decade off at the SAME state
                    // ( ts1727 cascade: all 8 trials rejected down to
                    // omega = 0.004, then the budget cut the timestep ).
                    // Accept the trial and iterate against the honest current
                    // residual; if the drift outruns the solver, the divergence
                    // and stagnation guards still cut -- through honest
                    // mechanisms instead of a reject storm. TWO flat pairs are
                    // required because a single one could also be a saturated
                    // divergence plateau ( Codex EQ1 ); within one line search
                    // the thermal state is fixed, so genuine drift stays flat
                    // across every halving and only pays one extra solve here
                    //
                    // the STATE is accepted, but the staged Anderson pair is
                    // not: its fixed-point residual mixes two thermal
                    // baselines and would poison the window ( O1 )
                    this->anderson_discard_magnetic() ;
                    this->anderson_clear_magnetic() ;

                    // the rescue must also restart the watchdog clock: the
                    // best-residual tracker references the OLD baseline and
                    // would otherwise cut the step the detector just rescued
                    mBestEpsilon = mEpsilon ;
                    mBestEpsilonIteration = mIteration ;

                    tRun = false ;
                }
                else
                {
                    // reject: restore the last accepted iterate ( fields and dof values )
                    for ( uint f=0; f<tNumFields; ++f )
                    {
                        mMesh->field_data( tFieldLabels( f ) )  = mBackupFields( f ) ;
                    }
                    for ( index_t k=0; k<tNumDofs; ++k )
                    {
                        tDofs( k )->value() = mBackupDofValues( k ) ;
                    }

                    // the rejected trial must not enter the mixing history,
                    // and the restored state invalidates what is in it
                    this->anderson_discard_magnetic() ;
                    this->anderson_clear_magnetic() ;

                    ++tBacktracks ;

                    // flat-pair bookkeeping for the moved-baseline detector:
                    // count consecutive rejects whose residual matches the
                    // previous trial within 0.05 decade AND whose omega
                    // genuinely decreased — at the mOmegaMin floor the halving
                    // clamps to a no-op and identical trials would count as
                    // flat without any physics content ( Grok R1c )
                    tFlatRejects = ( tBacktracks > 1
                        && tOmega < 0.6 * tOmegaTrialPrev
                        && std::abs( tLogEpsilon - tLogEpsilonPrev ) < 0.05 ) ?
                            tFlatRejects + 1 : 0 ;
                    tLogEpsilonPrev = tLogEpsilon ;
                    tOmegaTrialPrev = tOmega ;

                    // A large residual early in the timestep is expected -- the iteration is
                    // still settling from a big initial residual ( e.g. the first Picard
                    // step ), and a transient over-relaxed spike is recoverable by damping.
                    // So only cut the timestep once past the grace period mMinNumIterations.
                    // There, a residual above ~9 dB ( tLogEpsilon > 0.9, i.e. the dB column
                    // = 10*log10(eps) exceeds 9 ) cannot be rescued by under-relaxation at
                    // this timestep size, and an exhausted backtracking budget says the same:
                    // the timestep is too large, so cut it ( the last good iterate is already
                    // restored; reset_timestep() re-runs from the previous converged step at
                    // half the timestep ).
                    // A PROMOTED Newton tangent that rejects every trial from
                    // omega = 0.5 down to 0.5/2^7 is a tangent problem, not a
                    // timestep problem: the step direction is unusable at any
                    // scale ( the greg5 trace: net-current deck with a point
                    // bearing, the tangent's near-null phi-gauge mode -- and
                    // Delta t cuts made it WORSE, since promotion re-fires at
                    // every size ). Restore the last accepted iterate and
                    // LATCH PICARD for the rest of this timestep: Picard was
                    // contracting when it promoted, let it finish the step.
                    // The escalated tangent ( mForceNewton ) is exempt --
                    // there Picard has already broken down, so falling back
                    // would ping-pong and the timestep cut below is the only
                    // remaining lever.
                    if ( tBacktracks >= 8 && ! mForceNewton )
                    {
                        mEpsilon = mEpsilon0 ;
                        mEpsilonAbs = tEpsilonAbs0 ;

                        mJustPicard = true ;
                        mNewtonTrustStreak = 0 ;   // MIT-3a
                        mEquation->set_algorithm( SolverAlgorithm::Picard ) ;

                        // resume at half throttle like every other algorithm
                        // switch; the growth rule earns the rest back
                        mOmegaPicard = std::max( 0.5 * mOmegaPicard, mOmegaMin );

                        // the abandoned Newton samples must not pollute the
                        // stall window, and the watchdog judges the Picard
                        // continuation from here, not from the pre-promotion
                        // best
                        mResidualHistory.clear() ;
                        mBestEpsilon = mEpsilon ;
                        mBestEpsilonIteration = mIteration ;

                        if ( mCommRank == 0 && gLog.info_level() > 0 )
                        {
                            std::cout << sprint( "   │%-71s│",
                                " Newton line search exhausted: finishing this timestep on Picard" )
                                      << std::endl ;
                        }

                        tRun = false ;
                    }
                    // While the solver is a freshly escalated Newton tangent ( mForceNewton,
                    // set by try_escalate_to_newton ), do NOT cut the timestep on the first
                    // over-relaxed overshoot: a full Newton step from a stalled Picard iterate
                    // routinely overshoots, and the whole point of the escalation is to damp it
                    // through the line search. Let the backtracking budget run its course and
                    // only reset once it is exhausted ( tBacktracks >= 8 ); a non-escalated
                    // ( promoted ) Newton trial still resets on the +9 dB threshold as before.
                    // Picard never enters this branch -- its divergence exit is the strike
                    // rule further down.
                    else if ( mIteration > mMinNumIterations
                         && ( ( tLogEpsilon > 0.9 && ! mForceNewton ) || tBacktracks >= 8 ) )
                    {
                        this->reset_timestep() ;
                        return ;
                    }

                    // Within the grace period ( or before the +9 dB threshold ) keep damping.
                    // If the ESCALATED tangent's budget is spent here, stop retrying and accept
                    // the restored last-good iterate ( mEpsilon set to match so the post-loop
                    // bookkeeping stays consistent ) rather than looping forever.
                    else if ( tBacktracks >= 8 )
                    {
                        mEpsilon = mEpsilon0 ;
                        mEpsilonAbs = tEpsilonAbs0 ;
                        tRun = false ;
                    }
                    else
                    {
                        tOmega = std::clamp( 0.5*tOmega, mOmegaMin, mOmegaMax ) ;
                    }
                }
            }

            } // end magnetic body ( tSolveM )

            // freeze the thermal update while the magnetic residual is above
            // the gate: solving the heat equation from a badly unconverged
            // magnetic state feeds a garbage Joule source, and the poisoned
            // temperature poisons rho(T) right back -- the ts17 mutual
            // divergence spiral. With the update frozen, the magnetic solver
            // recovers against a FIXED temperature ( a contracting map ) and
            // the thermal solve resumes once the state is trustworthy.
            //
            // a magnetic iterate the post-solve checks are about to discard
            // ( +10 dB divergence or exhausted budget -- evaluated there
            // with the incremented counter; NaN is caught by the gate
            // comparison itself ) must not feed the thermal solve either:
            // since Picard iterates always accept, this is where a doomed
            // iterate would otherwise reach thermal and charge the shared
            // solver-failure counter for a state that is already condemned
            // gated on tSolveM like tMagneticReset ( round-3 R6 ): a
            // certified skip trip renders no budget verdict — at
            // mIteration == max an ungated prediction would freeze the
            // thermal side of a healthy step and trip the deadlock error
            const bool tMagneticDoomed = tSolveM
                && ( mIteration + 1 > mMaxNumIterations
                     || ( ! mForceNewton
                          && ( mEpsilon > 1E1 ? mDivergenceStrikes + 1 : 0 ) >= 3
                          && mIteration + 1 > mMinNumIterations ) ) ;

            // trip flags: did a thermal body run, and was the
            // thermal state co-certified at the committed magnetic state
            bool tSolvedThermal    = false ;
            bool tThermalCertified = false ;

            // a latched stall makes further thermal work pointless — without
            // this clause the stalled floor would keep taking bodies and the
            // no-body exit below could never fire ( round-3 C1 )
            const bool tUpdateThermal = mKernel2 != nullptr
                && mEpsilon < mThermalUpdateGate
                && ! tMagneticDoomed
                && ! mThermalStalled ;

            if( tUpdateThermal )
            {
                const SolverAlgorithm tAlgorithmT0 = mEquation2->algorithm() ;

                // Picard <-> Newton handoff for the thermal solver ( ported
                // from iterate_thermal ): Newton only once the thermal
                // residual is below its tolerance switch, always Picard on
                // the first iterates of a timestep, and always a Picard
                // restart after a freeze ( the state moved while thermal was
                // frozen; a cold Newton start from there diverges, cf. ts12 )
                // mForceNewton2 ( set by try_escalate_thermal_to_newton on a
                // flat Picard floor ) pins the solver to Newton
                // regardless of residual magnitude, cf. mForceNewton in
                // iterate_magnetic. Without it in this condition the handoff
                // would demote an escalated Newton right back to Picard on
                // the next iterate, since the residual is still above the
                // switch -- which is the whole reason it was escalated.
                if ( mAlgorithmThermal == SolverAlgorithm::NewtonRaphson
                     && ! mJustPicard2
                     && ( mForceNewton2
                          || ( mEpsilon2 < mEpsilonSwitch2 && mIteration2 > 1 ) )
                     && ! mThermalFrozen )
                {
                    if ( mEquation2->algorithm() == SolverAlgorithm::Picard )
                    {
                        // damped Newton entry on the first promotion only,
                        // cf. the magnetic flip in iterate_coupled
                        real tDamp2 = mFirstFlip2 ? 0.5 : 1.0 ;
                        mOmegaNewton2 = std::clamp( tDamp2 * mOmegaPicard2, mOmegaMin2, mOmegaMax2 );
                        mFirstFlip2 = false ;

                        // re-anchor the thermal watchdog on the promotion,
                        // cf. the magnetic flip
                        mBestEpsilon2 = mEpsilon2 ;
                        mBestEpsilonIteration2 = mIteration2 ;
                    }
                    mEquation2->set_algorithm( SolverAlgorithm::NewtonRaphson );
                }
                else
                {
                    // demotion hysteresis: the promotion edge
                    // compares a floating residual against a sharp threshold,
                    // and the reproducer showed partition-order roundoff alone
                    // deciding the comparison. A Newton that is already
                    // running therefore keeps running unless the residual
                    // retreats a full decade above the switch; the fresh-
                    // timestep reset of mEpsilon2 to BELFEM_REAL_MAX sits far
                    // above any band, so the first iterates of a step stay on
                    // Picard exactly as before. The freeze restart rule keeps
                    // priority: after a freeze the state has moved and a warm
                    // Newton continuation from it diverges ( ts12 ).
                    const bool tKeepNewton2 =
                           mEquation2->algorithm() == SolverAlgorithm::NewtonRaphson
                        && mAlgorithmThermal == SolverAlgorithm::NewtonRaphson
                        && ! mJustPicard2
                        && ! mThermalFrozen
                        && mEpsilon2 < 10.0 * mEpsilonSwitch2 ;

                    if ( ! tKeepNewton2 )
                    {
                        if ( mEquation2->algorithm() == SolverAlgorithm::NewtonRaphson && ! mFirstFlip2 )
                        {
                            mOmegaPicard2 = std::clamp( mOmegaNewton2, mOmegaMin2, mOmegaMax2 );
                        }
                        mEquation2->set_algorithm( SolverAlgorithm::Picard );
                    }
                }
                if ( mJustPicard2 )
                {
                    mEquation2->set_algorithm( SolverAlgorithm::Picard ) ;
                }
                mThermalFrozen = false ;

                if ( mEquation2->algorithm() != tAlgorithmT0 )
                {
                    // the algorithm changed: the thermal mixing history is void
                    this->anderson_clear_thermal() ;

                    // chatter guard, behind the demotion band above: a
                    // residual that still crosses the band every few iterates
                    // flips the algorithm. The third demotion within one
                    // attempt latches Picard for the rest of it.
                    if ( tAlgorithmT0 == SolverAlgorithm::NewtonRaphson
                         && ++mThermalFlipCount >= 3 )
                    {
                        mJustPicard2 = true ;
                    }
                }

                // bind the relaxation store of the ACTIVE thermal algorithm
                tOmega2Ptr = mEquation2->algorithm() == SolverAlgorithm::Picard ?
                    & mOmegaPicard2 : & mOmegaNewton2 ;
                *tOmega2Ptr = std::clamp( *tOmega2Ptr, mOmegaMin2, mOmegaMax2 ) ;
                mEquation2->set_omega( *tOmega2Ptr );
                tOmega2 = *tOmega2Ptr ;

                mEpsilon20 = mEpsilon2 ;

                // ---------- THERMAL HEAD ----------
                // Gauss-Seidel position preserved: assembled at the CURRENT
                // state, AFTER any magnetic body this trip — never reused
                // from before it ( round-3 Codex 1 )
                mKernel2->dofmgr()->compute_jacobian_and_rhs();

                this->dump_system_if_requested( mKernel2->dofmgr(), "thermal",
                                                mDumpCountThermal );

                mKernel2->dofmgr()->compute_residual();

                const real tEpsilonTHead    = mKernel2->dofmgr()->residual( mIteration2 );
                const real tEpsilonTHeadAbs = mKernel2->dofmgr()->absolute_residual();

                if ( std::isnan( tEpsilonTHead ) )
                {
                    this->reset_timestep() ;
                    return ;
                }

                const bool tConvergedT = tEpsilonTHead    <= mRelativeEpsilonTarget2
                                      || tEpsilonTHeadAbs <= mAbsoluteEpsilonTarget2 ;

                const bool tSolveT = ( ! tConvergedT )
                                  || mIteration2 < mMinNumIterations2 ;

                if ( ! tSolveT )
                {
                    // co-certificate: thermal measured at the very state the
                    // magnetic head certified this trip
                    mEpsilon2    = tEpsilonTHead ;
                    mEpsilonAbs2 = tEpsilonTHeadAbs ;
                    tThermalCertified = true ;
                }
                else
                {

                tSolvedThermal = true ;

                mKernel2->dofmgr()->solve_from_residual();

                // soft thermal solver failure: cut the timestep
                // ( cf. the magnetic line search above )
                if ( mKernel2->dofmgr()->solve_failed() )
                {
                    BELFEM_ERROR( ++mSolverFailCount < 8,
                        "thermal solver failed on %u consecutive timestep attempts - the cause is in the solver messages above ( e.g. a singular system, non-convergence, or repeated workspace exhaustion )",
                        ( unsigned int ) mSolverFailCount );
                    this->reset_timestep() ;
                    return ;
                }

                mEpsilon2 = mKernel2->dofmgr()->residual( mIteration2 );
                mEpsilonAbs2 = mKernel2->dofmgr()->absolute_residual();

                // COND1 / COND2 of the thermal system, cf. the magnetic capture above
                if ( mIteration2 == 0 )
                {
                    this->capture_conditioning_thermal() ;
                }

                // thermal stagnation exit: with the magnetic system converged
                // and the thermal residual bit-flat, further iterations change
                // nothing -- the loop would otherwise grind to the iteration
                if ( ( mEpsilon <= mRelativeEpsilonTarget
                       || mEpsilonAbs <= mAbsoluteEpsilonTarget )
                     && mEpsilon2 > mRelativeEpsilonTarget2
                     && std::abs( std::log10( mEpsilon2 )
                                - std::log10( mEpsilon20 ) ) < 0.001 )
                {
                    ++mThermalFlatCount ;

                    // two flat iterates with Newton configured but not
                    // running are already conclusive -- Picard has found a
                    // fixed point that is not a root, and no further Picard
                    // iterate will move it. Escalate instead of grinding to
                    // the stall warning. The >= 5 warning below remains the
                    // backstop for when Newton is unavailable, already
                    // running, or was already tried this attempt.
                    if ( mThermalFlatCount >= 2 )
                    {
                        this->try_escalate_thermal_to_newton() ;
                    }

                    if ( mThermalFlatCount >= 5 && ! mThermalStalled )
                    {
                        mThermalStalled = true ;
                        this->print_thermal_stall_warning() ;
                    }
                }
                else
                {
                    mThermalFlatCount = 0 ;
                }

                // no accept/reject loop on the thermal solve: every iterate
                // is accepted, so the staged pair commits right away
                this->anderson_commit_thermal() ;

                // thermal iterates in lockstep with the magnetic kernel here,
                // so its counter must advance too -- otherwise the tThermalReset
                // guard below compares a frozen mIteration2 == 0 against
                // mMinNumIterations2 and never fires ( reset per timestep in
                // initialize_timestep / reset_timestep )
                ++mIteration2 ;

                // progress watchdog on the thermal residual: catches the slow
                // monotone creep at the relaxation floor ( ts15 ) that never
                // trips the eps2 > 1E1 divergence rule
                if ( this->watchdog_thermal( tOmega2 ) )
                {
                    this->reset_timestep() ;
                    return ;
                }

                } // end thermal body ( tSolveT )
            }
            else if ( mKernel2 != nullptr )
            {
                // remember the freeze so the next thermal solve restarts
                // on Picard; the state keeps moving while thermal is frozen,
                // so the thermal mixing history dies with the freeze
                if ( ! mThermalFrozen )
                {
                    this->anderson_clear_thermal() ;
                }
                mThermalFrozen = true ;

                // flat evidence must not carry across a freeze: the state
                // moves while thermal is frozen ( Codex RQ2 )
                mThermalFlatCount = 0 ;
            }

            // a no-body trip whose head(s) met tolerance commits exactly the
            // measured state: magnetic certified above, thermal co-certified
            // this trip, absent, or explicitly accepted at a latched stall
            // ( the stall acceptance is UNCERTIFIED and says so )
            const bool tThermalDone =
                   mKernel2 == nullptr
                || tThermalCertified
                || mThermalStalled ;

            if ( ( ! tSolveM ) && tThermalDone )
            {
                if ( mCommRank == 0 && gLog.info_level() > 0 )
                {
                    // same 20/25/24 column geometry as print_header
                    std::cout << "   ├────────────────────┬─────────────────────────┬────────────────────────┤" << std::endl ;
                    if ( mKernel2 != nullptr && mThermalStalled && ! tThermalCertified )
                    {
                        std::cout << sprint( "   │ accepted at stall  │  magnetic : %10.3e  │  thermal : %10.3e  │",
                            mEpsilon, mEpsilon2 ) << std::endl ;
                    }
                    else if ( mKernel2 != nullptr )
                    {
                        std::cout << sprint( "   │ timestep succeeded │  magnetic : %10.3e  │  thermal : %10.3e  │",
                            mEpsilon, mEpsilon2 ) << std::endl ;
                    }
                    else
                    {
                        std::cout << sprint( "   │ timestep succeeded │  magnetic : %10.3e  │                        │",
                            mEpsilon ) << std::endl ;
                    }

                    // the section is NOT closed here: print_footer() adds the
                    // physics row to it, which it can only do once finalize()
                    // has run the postprocessor
                    mBoxSectionOpen = true ;
                }
                mTripExit = true ;
                return ;
            }

            // deadlock guard ( round-3 R-A ): magnetic certified but thermal
            // neither runnable nor stalled — nothing can advance, and a
            // silent spin here would be the exact silent-failure class
            BELFEM_ERROR( tSolveM || mKernel2 == nullptr
                          || tUpdateThermal || mThermalStalled,
                "coupled loop deadlock: magnetic certified ( eps = %g ) but the thermal update is gated ( thermal update gate = %g ) - review the deck's gate",
                ( double ) mEpsilon,
                ( double ) mThermalUpdateGate );

            if ( tSolveM )
            {
                // increment iteration counter ( a BODY count )
                ++mIteration ;
            }

            // print the line ( body trips and thermal-only trips both moved
            // state; pure certification trips exited above with their own
            // certificate line )
            if( mCommRank == 0 && gLog.info_level() > 0 )
            {
                this->print_line( tOmega, tOmega2 );
            }

            if ( tSolveM )
            {
                // remember this residual (in dB) for the stagnation guard
                mResidualHistory.push( 10. * std::log10( mEpsilon ) ) ;
            }

            // While the escalated Newton tangent is active ( mForceNewton ), do not cut the
            // timestep just because the residual is still above 1E1: recovering from a
            // stalled Picard iterate takes several damped Newton steps to bring the residual
            // down past +10 dB. Termination then rests on the stagnation guard ( a flat
            // residual hands Newton back to Picard, then a second Picard stall resets ) and
            // the mMaxNumIterations ceiling -- both still active below -- not the +10 dB rule.
            //
            // The rule is PERSISTENT, not instantaneous: three consecutive
            // iterates above the bar. A flux-front Picard overshoot from a
            // near-converged state recovers within two alpha halvings when
            // given the chance; an instantaneous cut selected attempts by
            // overshoot height instead of by trend and pinned the timestep
            // at the front ( greg3, t = 6.28 s )
            if ( tSolveM )
            {
                mDivergenceStrikes = mEpsilon > 1E1 ? mDivergenceStrikes + 1 : 0 ;
            }
            const bool tMagneticDiverged = ! mForceNewton
                && mDivergenceStrikes >= 3 && mIteration > mMinNumIterations ;
            // a skip trip runs no magnetic body: its counters and residual
            // describe the last body, so magnetic reset verdicts are only
            // rendered on body trips ( round-3 R6 — the magnetic budget must
            // not burn on idle certification trips )
            const bool tMagneticReset = tSolveM
                && ( mIteration > mMaxNumIterations
                     || tMagneticDiverged
                     || std::isnan( mEpsilon ) ) ;
            // thermal iteration budget: the deck's thermal max iterations,
            // enforced only while the thermal field is the reason the loop
            // is still alive. The magnetic evidence is a per-attempt LATCH:
            // once the magnetic field has reached its target ( rel or abs
            // escape, matching the certified exit ) within this attempt, a later
            // excursion back above target is the continuing thermal updates
            // perturbing a converged field ( the step-267 shape ) — it must
            // not disarm the budget. The certified exit still refuses a degraded
            // magnetic iterate, so the cut can never lose an acceptable
            // state. The thermal solves in lockstep with the magnetic
            // kernel, so mIteration2 tracks the coupled iterate count while
            // both fields converge; an unconditional ceiling here would cap
            // the whole coupled loop at the thermal budget. A latched
            // flat-stall is accepted at the exit, a flat streak in
            // progress defers the cut until it resolves, and tUpdateThermal
            // guards against judging a freeze-stale residual
            if ( mEpsilon <= mRelativeEpsilonTarget
                 || mEpsilonAbs <= mAbsoluteEpsilonTarget )
            {
                mMagneticHitTarget = true ;
            }
            const bool tThermalBudgetSpent = tUpdateThermal
                && mIteration2 > mMaxNumIterations2
                && mMagneticHitTarget
                && mEpsilon2    > mRelativeEpsilonTarget2
                && mEpsilonAbs2 > mAbsoluteEpsilonTarget2
                && ! mThermalStalled
                && mThermalFlatCount == 0 ;
            const bool tThermalReset = mKernel2 != nullptr
                && ( tThermalBudgetSpent
                     || ( mEpsilon2 > 1E1 && mIteration2 > mMinNumIterations2 )
                     || std::isnan( mEpsilon2 ) ) ;
            if( tMagneticReset || tThermalReset )
            {
                // Magnetic Picard breakdown: try the Newton tangent once before cutting Δt.
                // A NaN, an exhausted iteration budget, or a thermal-driven reset cannot be
                // rescued this way, so those fall straight through to reset_timestep().
                if ( tMagneticReset && ! tThermalReset
                     && mIteration <= mMaxNumIterations
                     && ! std::isnan( mEpsilon )
                     && this->try_escalate_to_newton() )
                {
                    // escalated to Newton; keep iterating this timestep at the current Δt
                }
                else
                {
                    //tOmega *= mAlpha ; // reduce the relaxation
                    this->reset_timestep() ;
                    return;
                }
            }

            // stagnation guard: cut the timestep if Picard has stalled ( Newton falls back
            // to Picard internally; see magnetic_stagnation_forces_reset ).
            // Both guards judge magnetic BODY progress — idle certification
            // trips must neither feed nor trigger them ( round-3 R6 )
            if ( tSolveM && this->magnetic_stagnation_forces_reset() )
            {
                this->reset_timestep() ;
                return;
            }

            // progress watchdog: no new residual minimum within the window
            // while still far from tolerance ( catches omega-sawtooth limit
            // cycles the flat-band guard cannot see, cf. ts15 )
            if ( tSolveM && this->watchdog_magnetic( tOmega ) )
            {
                this->reset_timestep() ;
                return;
            }

            // ( the mEpsilonFirst capture that lived here was a write-only
            // relic of the pre-iteration-count timestep controller )
            if ( tSolveM && mIteration > 1 )
            {
                // MIT-3a: the streak is judged for EVERY accepted iterate,
                // not only improving ones. Nesting the reset inside the
                // improvement branch would let a qualifier, an accepted
                // wobble, and another qualifier count as "consecutive"
                // ( Codex, 2026-08-18 ) -- exactly the single-lucky-step
                // case the streak exists to exclude
                if ( tRanNewton
                     && tBacktracks == 0
                     && mEpsilon < 0.5 * mEpsilon0 )
                {
                    ++mNewtonTrustStreak ;
                }
                else
                {
                    mNewtonTrustStreak = 0 ;
                }

                if( mEpsilon < mEpsilon0 )
                {
                    mNumIterationsDiv = 0;
                    // O2 amendment ( ts16 log ): the growth branch stays
                    // ACTIVE while Anderson is on. Holding omega made it a
                    // one-way ratchet to the floor -- alpha-decay and the
                    // backtracking kept shrinking it and nothing grew it
                    // back, starving the mixed step of its beta*r content.
                    // Overshoot is guarded by the flush-on-reject path.
                    real tGrowth = mBeta  + mGamma * 2./ constant::pi * std::atan( ( mEpsilon0 - mEpsilon ) / mEpsilon0 ) ;

                    // a Newton step accepted on the first trial that at least
                    // halves the residual has earned trust: recover omega
                    // geometrically instead of the <= 1.3x adaptation crawl
                    // ( ts34 spent 11 endgame iterations below omega 0.42 on a
                    // healthily contracting Newton ). tRanNewton, not the live
                    // algorithm: a same-call escalation must not classify the
                    // accepted Picard step as Newton trust
                    // MIT-3a ( todo/timestep_collapse_mitigation_design.md
                    // sec 4a ): the qualifying condition builds a STREAK, and
                    // only the second consecutive qualifying step earns the
                    // geometric recovery. One first-trial contraction is not
                    // evidence near a residual floor: on tapestack3d the
                    // doubling fired from omega 0.483 to 0.966 and the very
                    // next iterate lost 1.7 dB, then the punishment cascade
                    // burned the iteration budget and the step was rejected
                    // 8 dB from target. Deliberately NOT gated on the
                    // residual's distance to target: a healthy deep crawl
                    // ( the ts34 case this rule was written for ) sits at the
                    // same distance as the pathological one, so a scale test
                    // cannot separate them -- contraction evidence can
                    if ( mNewtonTrustStreak >= 2 )
                    {
                        tGrowth = 2.0 ;
                    }

                    tOmega *= tGrowth ;
                }
                else if ( mEpsilon < ( 1.0 + mOmegaNoiseBand ) * mEpsilon0 )
                {
                    // neutral band: a wobble at the residual noise floor is
                    // not divergence — hold omega and the divergence counter.
                    // Without this, improvement at a floored residual is a
                    // coin flip and the strict-decrease rule collapses omega
                    // to the floor ( see mOmegaNoiseBand )
                }
                else
                {
                    mNumIterationsDiv+=1 ;
                    if (mNumIterationsDiv >= mMaxNumIterationsDiv)
                    {
                        //Reset the relaxation if we are diverging for more than a certain number of iterations
                        tOmega = mOmega0 ;
                    }
                    else
                    {
                        tOmega *= mAlpha ;
                    }
                }
                tOmega = std::max( std::min( tOmega, mOmegaMax ), mOmegaMin );

            }

            //Same idea for thermal if exists; a frozen thermal update
            //( magnetic residual above the gate ) must not adapt omega2
            //against a stale residual pair. The update goes through the
            //pointer so it lands in the store of the algorithm that
            //actually ran ( mOmegaPicard2 or mOmegaNewton2 ). Moved out of
            //the magnetic adaptation gate: a thermal body on a
            //magnetic-skip trip still adapts omega2, and a thermal HEAD that
            //certified without solving must not ( round-3 C4/R6 )
            if ( tSolvedThermal && tOmega2Ptr != nullptr )
            {
                real & tOmegaT = *tOmega2Ptr ;
                if( mEpsilon2 < mEpsilon20 )
                {
                    // O2 amendment ( ts16 log ): growth stays active,
                    // see the magnetic block for the rationale
                    tOmegaT *= mBeta  + mGamma * 2./ constant::pi * std::atan( ( mEpsilon20 - mEpsilon2 ) / mEpsilon20 ) ;
                }
                else if ( mEpsilon2 < ( 1.0 + mOmegaNoiseBand ) * mEpsilon20 )
                {
                    // neutral band, cf. the magnetic block: a converged
                    // thermal residual wobbling on roundoff must not
                    // grind omega2 to the floor
                }
                else
                {
                    tOmegaT *= mAlpha ;
                }
                tOmegaT = std::max( std::min( tOmegaT, mOmegaMax2 ), mOmegaMin2 );
            }
        }

        void
        Controller::solve_magnetic()
        {
            // certified exit, segregated magnetic loop: one call per
            // timestep. Each trip head-measures the committed state; the
            // trip exits certified, resets, or takes exactly one body
            mTripExit = false ;

            while ( true )
            {
                this->iterate_magnetic() ;

                if ( mReset || mTripExit )
                {
                    return ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Controller::solve_thermal()
        {
            // certified exit, segregated thermal loop ( cf.
            // solve_magnetic ). reset_thermal() sets mResetThermal only —
            // the smaller sub-step exists only after the driver's next
            // initialize_thermal() call, which recomputes mDeltaTime2 and
            // clears the flag. Returning here hands control back to the
            // driver's outer time2 loop, exactly as the old
            // while( run_thermal() ) predicate did ( code-audit P0 )
            mTripExit = false ;

            while ( true )
            {
                this->iterate_thermal() ;

                if ( mReset || mResetThermal || mTripExit )
                {
                    return ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Controller::iterate_magnetic()
        {
            // remember the old epsilon value
            mEpsilon0 = mEpsilon ;

            // Picard <-> Newton handoff (Messe et al. 2023, Eq. 13) — see
            // iterate_coupled for the rationale.
            const SolverAlgorithm tAlgorithm0 = mEquation->algorithm() ;
            // mForceNewton ( set by try_escalate_to_newton on a Picard breakdown ) pins the
            // solver to Newton regardless of residual magnitude; the normal promotion still
            // fires on its own once the residual drops below mEpsilonSwitch. The mJustPicard
            // exclusion keeps a latched fallback from re-anchoring the watchdog every
            // iterate, cf. iterate_coupled
            if ( mAlgorithm == SolverAlgorithm::NewtonRaphson && ! mJustPicard
                 && ( mForceNewton || ( mEpsilon < mEpsilonSwitch && mIteration > 1 ) ) )
            {
                if ( mEquation->algorithm() == SolverAlgorithm::Picard )
                {
                    // damped Newton entry on the first promotion only,
                    // cf. iterate_coupled
                    real tDamp = mFirstFlip ? 0.5 : 1.0 ;
                    mOmegaNewton = std::clamp( tDamp * mOmegaPicard, mOmegaMin, mOmegaMax );
                    mFirstFlip = false ;

                    // re-anchor the progress watchdog on the promotion,
                    // cf. iterate_coupled
                    mBestEpsilon = mEpsilon ;
                    mBestEpsilonIteration = mIteration ;
                }
                mEquation->set_algorithm( SolverAlgorithm::NewtonRaphson );
            }
            else
            {
                if ( mEquation->algorithm() == SolverAlgorithm::NewtonRaphson && ! mFirstFlip )
                {
                    mOmegaPicard = std::clamp( mOmegaNewton, mOmegaMin, mOmegaMax );
                    mNewtonTrustStreak = 0 ;   // MIT-3a: ordinary demotion voids it
                }
                mEquation->set_algorithm( SolverAlgorithm::Picard );
            }
            if ( mJustPicard )
            {
                mEquation->set_algorithm( SolverAlgorithm::Picard ) ;
            }
            if ( mEquation->algorithm() != tAlgorithm0 )
            {
                mResidualHistory.clear() ;
                this->anderson_clear_magnetic() ;
            }

            real & tOmega = mEquation->algorithm() == SolverAlgorithm::Picard ? mOmegaPicard : mOmegaNewton ;

            // set relaxation parameters ( clamp BEFORE handing to the
            // equation, like the coupled and thermal paths )
            tOmega = std::clamp( tOmega, mOmegaMin, mOmegaMax ) ;
            mEquation->set_omega( tOmega );



            this->reset_dotQ();
            mKernel->dofmgr()->compute_jacobian_and_rhs();
            this->collect_dotQ();

            this->dump_system_if_requested( mKernel->dofmgr(), "magnetic",
                                            mDumpCountMagnetic );

            this->impose_voltage_bcs();

            mKernel->dofmgr()->compute_residual();

            const real tEpsilonHead    = mKernel->dofmgr()->residual( mIteration );
            const real tEpsilonHeadAbs = mKernel->dofmgr()->absolute_residual();

            if ( std::isnan( tEpsilonHead ) )
            {
                this->reset_timestep() ;
                return ;
            }

            const bool tConvergedM = tEpsilonHead    <= mRelativeEpsilonTarget
                                  || tEpsilonHeadAbs <= mAbsoluteEpsilonTarget ;

            if ( tConvergedM && mIteration >= mMinNumIterations )
            {
                // certified exit: the committed state is the measured state
                mEpsilon    = tEpsilonHead ;
                mEpsilonAbs = tEpsilonHeadAbs ;

                if ( mCommRank == 0 && gLog.info_level() > 0 )
                {
                    // same 20/25/24 column geometry as print_header
                    std::cout << "   ├────────────────────┬─────────────────────────┬────────────────────────┤" << std::endl ;
                    // a magnetic-only deck really is done here; with a
                    // thermal kernel the sub-steps still follow
                    std::cout << sprint( "   │ %-18s │  magnetic : %10.3e  │                        │",
                        mKernel2 == nullptr ? "timestep succeeded"
                                            : "magnetic converged",
                        mEpsilon ) << std::endl ;

                    // magnetic-only: the footer follows immediately, so the
                    // section stays open for its physics row. With a thermal
                    // kernel the sub-step boxes come in between and this one
                    // has to be closed now
                    if ( mKernel2 == nullptr )
                    {
                        mBoxSectionOpen = true ;
                    }
                    else
                    {
                        std::cout << "   ├────────────────────┴─────────────────────────┴────────────────────────┤" << std::endl ;
                        mBoxSectionOpen = false ;
                    }
                }
                mTripExit = true ;
                return ;
            }

            // ---------- BODY ( the head assembly is consumed here ) ----------
            mKernel->dofmgr()->solve_from_residual();

            // soft solver failure: cut the timestep ( cf. iterate_coupled )
            if ( mKernel->dofmgr()->solve_failed() )
            {
                BELFEM_ERROR( ++mSolverFailCount < 8,
                    "linear solver failed on %u consecutive timestep attempts - the cause is in the solver messages above ( e.g. a singular system, non-convergence, or repeated workspace exhaustion )",
                    ( unsigned int ) mSolverFailCount );
                this->reset_timestep() ;
                return ;
            }

            // get the residual
            mEpsilon = mKernel->dofmgr()->residual( mIteration );
            mEpsilonAbs = mKernel->dofmgr()->absolute_residual();

            // COND1 / COND2 of the magnetic system, cf. iterate_coupled
            if ( mIteration == 0 )
            {
                this->capture_conditioning_magnetic() ;
            }

            // no accept/reject loop here: every iterate is accepted, so the
            // staged Anderson pair commits right away
            this->anderson_commit_magnetic() ;

            // increment iteration counter
            ++mIteration ;

            // print the line
            if( mCommRank == 0 && gLog.info_level() > 0 )
            {
                this->print_line_magnetic( tOmega );
            }

            // remember this residual (in dB) for the stagnation guard
            mResidualHistory.push( 10. * std::log10( mEpsilon ) ) ;

            // persistent divergence rule, cf. iterate_coupled: three
            // consecutive iterates above the bar before the cut. The
            // mForceNewton exemption mirrors the coupled path -- escalation
            // must disable the +10 dB exit or the escalated tangent never
            // gets its damped recovery ( cf. try_escalate_to_newton )
            mDivergenceStrikes = mEpsilon > 1E1 ? mDivergenceStrikes + 1 : 0 ;
            if( mIteration > mMaxNumIterations || ( ( ! mForceNewton && mDivergenceStrikes >= 3 && mIteration > mMinNumIterations ) || std::isnan( mEpsilon ) ) )
            {
                // Magnetic Picard breakdown: try the Newton tangent once before cutting Δt
                // ( see try_escalate_to_newton ). NaN / exhausted-budget cases fall through.
                if ( mIteration <= mMaxNumIterations
                     && ! std::isnan( mEpsilon )
                     && this->try_escalate_to_newton() )
                {
                    // escalated to Newton; keep iterating this timestep at the current Δt
                }
                else
                {
                    //tOmega *= mAlpha ; // reduce the relaxation
                    this->reset_timestep() ;
                    return;
                }
            }

            // stagnation guard: cut the timestep if Picard has stalled ( Newton falls back
            // to Picard internally; see magnetic_stagnation_forces_reset ).
            if ( this->magnetic_stagnation_forces_reset() )
            {
                this->reset_timestep() ;
                return;
            }

            // progress watchdog: no new residual minimum within the window
            // while still far from tolerance ( catches omega-sawtooth limit
            // cycles the flat-band guard cannot see, cf. ts15 )
            if ( this->watchdog_magnetic( tOmega ) )
            {
                this->reset_timestep() ;
                return;
            }

            // ( dead mEpsilonFirst write removed )
            if ( mIteration > 1 )
            {
                if( mEpsilon < mEpsilon0 )
                {
                    mNumIterationsDiv = 0;
                    // O2 amendment ( ts16 log ): the growth branch stays
                    // ACTIVE while Anderson is on. Holding omega made it a
                    // one-way ratchet to the floor -- alpha-decay and the
                    // backtracking kept shrinking it and nothing grew it
                    // back, starving the mixed step of its beta*r content.
                    // Overshoot is guarded by the flush-on-reject path.
                    tOmega *= mBeta  + mGamma * 2./ constant::pi * std::atan( ( mEpsilon0 - mEpsilon ) / mEpsilon0 ) ;
                }
                else if ( mEpsilon < ( 1.0 + mOmegaNoiseBand ) * mEpsilon0 )
                {
                    // neutral band, cf. iterate_coupled: a noise-floor wobble
                    // holds omega and the divergence counter
                }
                else
                {
                    mNumIterationsDiv+=1 ;
                    if (mNumIterationsDiv >= mMaxNumIterationsDiv)
                    {
                        //Reset the relaxation if we are diverging for more than a certain number of iterations
                        tOmega = mOmega0 ;
                    }
                    else
                    {
                        tOmega *= mAlpha ;
                    }
                }
                tOmega = std::max( std::min( tOmega, mOmegaMax ), mOmegaMin );
            }
        }

        void
        Controller::reset_dotQ()
        {
            mMesh->global_variable( "dotQ" )->value() = 0.0 ;
        }

        void
        Controller::collect_dotQ()
        {
            // collect dotQ from the other procs
            if ( mCommSize > 1 )
            {
                if ( mCommRank == 0 )
                {
                    Vector< real > tDotQ( mCommSize, 0.0 );
                    collect( tDotQ, mMesh->global_variable( "dotQ" )->value() );
                    mMesh->global_variable( "dotQ" )->value() = sum ( tDotQ );
                }
                else
                {
                    send( mMesh->global_variable( "dotQ" )->value() );
                }
                comm_barrier();
            }
        }

        void
        Controller::iterate_thermal()
        {
            // remember the old epsilon value
            mEpsilon20 = mEpsilon2 ;

            // Picard <-> Newton handoff (Messe et al. 2023, Eq. 13) — thermal.
            // (No ShiftRegister stagnation guard on the thermal path yet; see iterate_coupled
            // for the magnetic version.)
            const SolverAlgorithm tAlgorithmT0 = mEquation2->algorithm() ;
            if ( mAlgorithmThermal == SolverAlgorithm::NewtonRaphson && ! mJustPicard2
                 && mEpsilon2 < mEpsilonSwitch2 && mIteration2 > 1 )
            {
                if ( mEquation2->algorithm() == SolverAlgorithm::Picard )
                {
                    // damped Newton entry on the first promotion only,
                    // cf. the magnetic flip in iterate_coupled
                    real tDamp2 = mFirstFlip2 ? 0.5 : 1.0 ;
                    mOmegaNewton2 = std::clamp( tDamp2 * mOmegaPicard2, mOmegaMin2, mOmegaMax2 );
                    mFirstFlip2 = false ;

                    // re-anchor the thermal watchdog on the promotion,
                    // cf. iterate_coupled
                    mBestEpsilon2 = mEpsilon2 ;
                    mBestEpsilonIteration2 = mIteration2 ;
                }
                mEquation2->set_algorithm( SolverAlgorithm::NewtonRaphson );
            }
            else
            {
                if ( mEquation2->algorithm() == SolverAlgorithm::NewtonRaphson && ! mFirstFlip2 )
                {
                    mOmegaPicard2 = std::clamp( mOmegaNewton2, mOmegaMin2, mOmegaMax2 );
                }
                mEquation2->set_algorithm( SolverAlgorithm::Picard );
            }
            if ( mJustPicard2 )
            {
                mEquation2->set_algorithm( SolverAlgorithm::Picard ) ;
            }
            if ( mEquation2->algorithm() != tAlgorithmT0 )
            {
                // the algorithm changed: the thermal mixing history is void
                this->anderson_clear_thermal() ;
            }

            real & tOmega = mEquation2->algorithm() == SolverAlgorithm::Picard ? mOmegaPicard2 : mOmegaNewton2 ;
            tOmega = std::clamp( tOmega, mOmegaMin2, mOmegaMax2 ) ;

            // set relaxation parameters
            mEquation2->set_omega( tOmega );
            mKernel2->dofmgr()->compute_jacobian_and_rhs();

            this->dump_system_if_requested( mKernel2->dofmgr(), "thermal",
                                            mDumpCountThermal );

            mKernel2->dofmgr()->compute_residual();

            const real tEpsilonTHead    = mKernel2->dofmgr()->residual( mIteration2 );
            const real tEpsilonTHeadAbs = mKernel2->dofmgr()->absolute_residual();

            if ( std::isnan( tEpsilonTHead ) )
            {
                if ( mDeltaTime2 < mDeltaTime/20 )
                {
                    mCouplingFactor = 20 ;
                    this->reset_timestep() ;
                }
                else
                {
                    this->reset_thermal() ;
                }
                return ;
            }

            const bool tConvergedT = tEpsilonTHead    <= mRelativeEpsilonTarget2
                                  || tEpsilonTHeadAbs <= mAbsoluteEpsilonTarget2 ;

            if ( tConvergedT && mIteration2 >= mMinNumIterations2 )
            {
                // certified exit: the committed state is the measured state
                mEpsilon2    = tEpsilonTHead ;
                mEpsilonAbs2 = tEpsilonTHeadAbs ;

                if ( mCommRank == 0 && gLog.info_level() > 0 )
                {
                    // same 20/25/24 column geometry as print_header
                    std::cout << "   ├────────────────────┬─────────────────────────┬────────────────────────┤" << std::endl ;
                    std::cout << sprint( "   │ substep succeeded  │                         │  thermal : %10.3e  │",
                        mEpsilon2 ) << std::endl ;
                    std::cout << "   ├────────────────────┴─────────────────────────┴────────────────────────┤" << std::endl ;
                    mBoxSectionOpen = false ;
                }
                mTripExit = true ;
                return ;
            }

            // ---------- BODY ( the head assembly is consumed here ) ----------
            mKernel2->dofmgr()->solve_from_residual();

            // soft thermal solver failure: mirror the divergence reset logic
            if ( mKernel2->dofmgr()->solve_failed() )
            {
                BELFEM_ERROR( ++mSolverFailCount < 8,
                    "thermal solver failed on %u consecutive timestep attempts - the cause is in the solver messages above ( e.g. a singular system, non-convergence, or repeated workspace exhaustion )",
                    ( unsigned int ) mSolverFailCount );
                if ( mDeltaTime2 < mDeltaTime/20 )
                {
                    mCouplingFactor = 20 ;
                    this->reset_timestep() ;
                }
                else
                {
                    this->reset_thermal() ;
                }
                return ;
            }

            // get the residual
            mEpsilon2 = mKernel2->dofmgr()->residual( mIteration2 );
            mEpsilonAbs2 = mKernel2->dofmgr()->absolute_residual();

            // COND1 / COND2 of the thermal system, cf. iterate_coupled
            if ( mIteration2 == 0 )
            {
                this->capture_conditioning_thermal() ;
            }

            // no accept/reject loop here: the staged pair commits right away
            this->anderson_commit_thermal() ;

            // increment iteration counter
            ++mIteration2 ;

            // print the line
            if( mCommRank == 0 && gLog.info_level() > 0 )
            {
                this->print_line_thermal( tOmega );
            }

            if( mIteration2 > mMaxNumIterations2 || ((mEpsilon2 > 1E1 && mIteration2 > mMinNumIterations2) || std::isnan( mEpsilon2 ) ) )
            {
                //tOmega *= mAlpha ; // reduce the relaxation
                if (mDeltaTime2 < mDeltaTime/20)
                {
                    mCouplingFactor = 20 ;
                    this->reset_timestep() ;
                }
                else
                {
                    this->reset_thermal() ;
                }
                return;
            }
            // progress watchdog: same consequence as the divergence rule above
            if ( this->watchdog_thermal( tOmega ) )
            {
                if (mDeltaTime2 < mDeltaTime/20)
                {
                    mCouplingFactor = 20 ;
                    this->reset_timestep() ;
                }
                else
                {
                    this->reset_thermal() ;
                }
                return;
            }

            if ( mIteration2 > 2  && std::abs( std::log10(mEpsilon20)-std::log10(mEpsilon2) ) < BELFEM_EPSILON ) //If convergence stagnates, latch to Picard for the rest of the timestep
            {
                mJustPicard2 = true ;
                if ( mEquation2->algorithm() != SolverAlgorithm::Picard )
                {
                    // demoting outside the handoff: the history is void
                    this->anderson_clear_thermal() ;
                }
                mEquation2->set_algorithm( SolverAlgorithm::Picard ) ;
            }

            // ( dead mEpsilonFirst2 write removed )
            if ( mIteration2 > 1 )
            {
                if( mEpsilon2 < mEpsilon20 )
                {
                    mNumIterationsDiv = 0;
                    // O2 amendment ( ts16 log ): growth stays active,
                    // see iterate_coupled for the rationale
                    tOmega *= mBeta  + mGamma * 2./ constant::pi * std::atan( ( mEpsilon20 - mEpsilon2 ) / mEpsilon20 ) ;
                }
                else if ( mEpsilon2 < ( 1.0 + mOmegaNoiseBand ) * mEpsilon20 )
                {
                    // neutral band, cf. iterate_coupled: a noise-floor wobble
                    // holds omega and the divergence counter
                }
                else
                {
                    mNumIterationsDiv+=1 ;
                    if (mNumIterationsDiv >= mMaxNumIterationsDiv)
                    {
                        //Reset the relaxation if we are diverging for more than a certain number of iterations
                        tOmega = mOmega0 ;
                    }
                    else
                    {
                        tOmega *= mAlpha ;
                    }
                }
                tOmega = std::max( std::min( tOmega, mOmegaMax2 ), mOmegaMin2 );
            }
        }

        bool
        Controller::solve_circuit()
        {
            BELFEM_ASSERT( mCircuit != nullptr, "Circuit solver is not initialized" ) ;

            uint tIt_Circuit = 0 ;
            real tOmega_Circuit = 1.0 ;
            real tEpsilon_Circuit = BELFEM_REAL_MAX ;
            real tEpsilon0_Circuit = BELFEM_REAL_MAX ;

            std::cout << std::endl ;
            std::cout <<     "   ┌──────────────────────────────────────────────────────────────────┐"<< std::endl ;
            string tFormat = "   │                          Circuit solver                          │" ;
            string tMessage = sprint( tFormat.c_str(), mRunningTimeStep, mTime * 1000, mDeltaTime * 1000 );
            std::cout << tMessage << std::endl ;
            std::cout <<     "   ├──────────────────────────────────────────────────────────────────┤"<< std::endl ;

            while ( tIt_Circuit < mMaxNumIterations && tEpsilon_Circuit > mRelativeEpsilonTarget  )
            {
                tEpsilon0_Circuit = tEpsilon_Circuit ;

                tIt_Circuit++;

                mCircuit->set_omega(tOmega_Circuit) ;
                mCircuit->compute_jacobian_and_rhs();
                mCircuit->solve() ;
                tEpsilon_Circuit = mCircuit->residual() ;

                real dB = 10. * std::log10( tEpsilon_Circuit );
                tFormat =  "   │ Newton " ;
                tFormat += "Step %2u, residual= %s (%7.2f dB ), relax = %6.4f │" ;

                tMessage = sprint( tFormat.c_str(), tIt_Circuit,
                    residual_string( tEpsilon_Circuit ).c_str(), dB, tOmega_Circuit );

                std::cout << tMessage << std::endl ;

                if(tEpsilon_Circuit < tEpsilon0_Circuit)
                {
                    tOmega_Circuit *= mBeta + mGamma * 2./ constant::pi * std::atan( ( tEpsilon0_Circuit - tEpsilon_Circuit ) / tEpsilon0_Circuit ) ;
                }
                else
                {
                    tOmega_Circuit*=mAlpha ;
                }
                tOmega_Circuit = std::min(tOmega_Circuit, 1.0) ;

            }

            std::cout << "   └──────────────────────────────────────────────────────────────────┘"<< std::endl ;

            if (tIt_Circuit == mMaxNumIterations )
            {
                std::cout << "Time step not convergent, reducing the time step" << std::endl ;
                return false ;
            }


            return true ;
        }

        void Controller::reset_timestep()
        {
            --mRunningTimeStep ;
            mReset = true ;


            // MIT-3a: a rejected step's Newton contractions are not evidence
            // for the retry -- the retry solves a DIFFERENT problem ( smaller
            // Delta t ), and carrying the streak would let the first
            // qualifying step of the retry double omega immediately
            mNewtonTrustStreak = 0 ;

            // Omega reset policy: when the magnetic iteration was diverging ( residual
            // above 0 dB, i.e. mEpsilon > 1 ), carry the under-relaxation over instead of
            // restoring the full step. Restoring omega = 1 here just re-diverges once the
            // timestep is cut. The recovery branch in iterate_magnetic() ( omega grows as
            // the residual drops ) restores omega on its own once the smaller timestep
            // starts converging, so this cannot get stuck at a tiny omega. Only restore
            // the initial omega when the residual was already below 0 dB ( a clean
            // accuracy / timestep-size reset, not a divergence ).
            if ( mEpsilon > 1.0 )
            {
                mOmegaPicard = std::max( std::min( mOmegaPicard, mOmegaNewton ), mOmegaMin );
                mOmegaNewton = mOmegaPicard ;
            }
            else
            {
                mOmegaPicard = mOmega0 ;
                mOmegaNewton = mOmega0 ;
            }
            mEquation->reset_fields() ;
            if ( mEquation2 != nullptr )
            {
                mOmegaPicard2 = 1.0 ;
                mOmegaNewton2 = 1.0 ;

                // the thermal equation may have shifted several sub-steps
                // within this magnetic step; reset_fields would only reverse
                // the last one, so restore the step-start savepoint instead
                mEquation2->restore_savepoint() ;
            }

            // the FIELDS are now restored, but the DOFS still hold the failed
            // attempt's last iterate, and the Picard residual reads the dofs
            // ( r = A x - b through mFieldValues ), not the fields. A
            // non-finite iterate survives x -= omega * delta unchanged, so
            // without this seed a poisoned attempt repeats on every retry —
            // the warm-restart path seeds for exactly this reason
            // ( load_memdump: "the restored state would be invisible to the
            // Picard residual" )
            mKernel->dofmgr()->seed_dof_values() ;
            if ( mKernel2 != nullptr )
            {
                mKernel2->dofmgr()->seed_dof_values() ;
            }
            mIteration = 0 ;
            mIteration2 = 0 ;

            // re-arm the first-promotion latches so the retry damps its first
            // Newton entry again ( currently redundant with the drivers'
            // initializers, but reset must be self-contained — Codex DQ1 )
            mFirstFlip  = true ;
            mFirstFlip2 = true ;

            // re-arm the thermal stagnation exit
            mThermalFlatCount = 0 ;
            mThermalStalled   = false ;
            mEpsilonAbs       = BELFEM_REAL_MAX ;
            mEpsilonAbs2      = BELFEM_REAL_MAX ;

            // the rejected attempt's magnetic-hit-target evidence does not
            // carry to the retry ( it solves a different problem )
            mMagneticHitTarget = false ;

            // the abandoned residuals must not leak into the next attempt
            mResidualHistory.clear() ;
            mDivergenceStrikes = 0 ;

            // neither must the Anderson mixing history and the divergence
            // bookkeeping of the failed attempt
            this->anderson_clear_magnetic() ;
            this->anderson_clear_thermal() ;
            mNumIterationsDiv = 0 ;
            mThermalFlipCount = 0 ;
            mBestEpsilon  = BELFEM_REAL_MAX ;
            mBestEpsilon2 = BELFEM_REAL_MAX ;
            mWatchdogOmegaPrev  = BELFEM_QUIET_NAN ;
            mWatchdogOmegaPrev2 = BELFEM_QUIET_NAN ;
            mBestEpsilonIteration  = 0 ;
            mBestEpsilonIteration2 = 0 ;

            if( mCommRank == 0 && gLog.info_level() > 0 )
            {
                std::cout << "   └───────────────────────────────────────────────────────────────────────┘"<< std::endl ;
            }

            // the minimum timestep is a hard floor in BOTH control modes:
            // never run below it. A cut demanded while Delta t already sits
            // AT the floor cannot change the timestep — instead of aborting,
            // WIDEN the iteration budgets ( doubling per floor retry, capped
            // at mFloorEscalationCap x the deck values ) and retry at the
            // minimum: a slowly converging attempt gets the iterations it
            // needs, and the budgets return to the deck values with the next
            // accepted step ( finalize ). A residual that cannot reach the
            // tolerance at ANY budget keeps retrying at the floor — the
            // deck-side accept for that case is the "absolute tolerance"
            // escape ( sidecoatings: a Delta-t-independent residual floor
            // above the relative tolerance )
            if ( 0.5 * mDeltaTime < mDeltaTimeMin )
            {
                if ( mDeltaTime <= mDeltaTimeMin )
                {
                    // escalation factor 2^retries, capped; identical member
                    // math on every rank ( only the print is rank-guarded )
                    ++mFloorRetries ;

                    // backstop: escalating forever is not a strategy. A
                    // residual that ignores both the timestep and a 4x
                    // iteration budget will not converge by repeating the
                    // attempt, and an unattended job would spend its whole
                    // allocation here. Everything up to the last accepted
                    // step is already saved
                    BELFEM_ERROR( mMaxFloorRetries == 0
                                  || mFloorRetries <= mMaxFloorRetries,
                        "the floor retry limit ( %u ) is exhausted: the timestep stayed\n"
                        "       pinned at the minimum ( %.3e s ) with the iteration budgets\n"
                        "       raised up to %ux the deck values, and the nonlinear residual\n"
                        "       still does not respond to the timestep size.\n"
                        "       Results up to the last accepted step are saved.\n"
                        "       Diagnose with 'compute conditioning : true' in the solver\n"
                        "       section ( a residual floor near kappa * machine epsilon is a\n"
                        "       conditioning limit, not a timestep one ), accept a floored\n"
                        "       residual with 'absolute tolerance' in the nonlinear section,\n"
                        "       or set 'floor retries : 0' to keep retrying indefinitely.",
                        ( unsigned int ) mMaxFloorRetries,
                        ( double ) mDeltaTimeMin,
                        ( unsigned int ) mFloorEscalationCap );
                    uint tFactor = 1 ;
                    for ( uint k = 0; k < mFloorRetries
                                      && tFactor < mFloorEscalationCap; ++k )
                    {
                        tFactor *= 2 ;
                    }
                    mMaxNumIterations  = mMaxNumIterationsDeck  * tFactor ;
                    mMaxNumIterations2 = mMaxNumIterations2Deck * tFactor ;
                    mWatchdogWindow    = mWatchdogWindowDeck    * tFactor ;
                    mWatchdogWindow2   = mWatchdogWindow2Deck   * tFactor ;

                    if ( mCommRank == 0 && gLog.info_level() > 0 )
                    {
                        std::cout << "   ┌───────────────────────────────────────────────────────────────────────┐" << std::endl ;
                        std::cout << sprint( "   │%-71s│",
                            " WARNING: timestep is pinned at the configured minimum" ) << std::endl ;
                        string tLine = sprint(
                            " raising the iteration budget to %u ( %u x deck value )",
                            ( unsigned int ) mMaxNumIterations,
                            ( unsigned int ) tFactor );
                        std::cout << sprint( "   │%-71s│", tLine.c_str() ) << std::endl ;

                        tLine = mMaxFloorRetries == 0 ?
                            sprint( " floor retry %u ( no limit set )",
                                ( unsigned int ) mFloorRetries ) :
                            sprint( " floor retry %u of %u before the run stops",
                                ( unsigned int ) mFloorRetries,
                                ( unsigned int ) mMaxFloorRetries );
                        std::cout << sprint( "   │%-71s│", tLine.c_str() ) << std::endl ;
                        std::cout << sprint( "   │%-71s│",
                            " hint: an 'absolute tolerance' accepts a floored residual" ) << std::endl ;
                        std::cout << "   └───────────────────────────────────────────────────────────────────────┘" << std::endl ;
                    }
                }
                mDeltaTime = mDeltaTimeMin ;
            }
            else
            {
                mDeltaTime *= 0.5 ;
            }

            // A rejection is one of the most diagnostically important events
            // in a run, and it used to be SILENT — the message here was
            // commented out, so the only trace of a rejected step was the
            // same step number reappearing with a smaller Δt in the next
            // header ( 2026-08-13 quench night: 19 rejections, every one
            // reconstructed after the fact by diffing step headers ). Printed
            // AFTER the Δt update so the value shown is the one the retry
            // actually uses — the cut is not always a plain halving, the
            // floor branch above clamps instead.
            if( mCommRank == 0 && gLog.info_level() > 0 )
            {
                // +1: reset_timestep decremented mRunningTimeStep at entry,
                // so the ATTEMPTED step — the one whose header the reader
                // just saw — is one above the counter here ( observed on the
                // first live print: the message named 683 for an attempt
                // headed 684 )
                // same 20/25/24 column geometry as print_header, so a
                // rejection reads as the twin of the header it follows
                std::cout << "   ┌────────────────────┬─────────────────────────┬────────────────────────┐" << std::endl ;
                std::cout << sprint(
                    "   │ timestep rejected  │  retrying step %6u   │   Δ t = %10.4f ms  │",
                    ( unsigned int ) mRunningTimeStep + 1u,
                    mDeltaTime * 1000.0 ) << std::endl ;
                std::cout << "   └────────────────────┴─────────────────────────┴────────────────────────┘" << std::endl ;
            }

            if ( ! mUseLegacyTimestepControl )
            {
                // the failed attempt's cost history is meaningless for the
                // next accepted step, and stale samples would kick the PID
                // through the P/D terms; hold growth for the following
                // accepted steps ( see mPostFailureHold )
                mCtrlErr0 = 1.0 ;
                mCtrlErr1 = 1.0 ;
                mCtrlErr2 = 1.0 ;
                mPostFailureHold = mPostFailureHoldSteps ;
            }

            mDeltaTimeTemporary = BELFEM_QUIET_NAN ;
            mTime = mTime0 ;
            mTime2 = mTime0 ;
            mTime02 = mTime0 ;

            if( mKernel2 != nullptr )
            {
                mMesh2->time_stamp() = mTime ;
            }

            // circuit state is rank-0-owned ( only rank 0 shifts,
            // solves and allocates mX in compute_MNA_matrix ), so only
            // rank 0 may roll it back -- on the other ranks shift_back()
            // walked into the never-allocated solution vector
            if (mCircuit != nullptr && mCommRank == 0)
            {
                mCircuit->shift_back() ;
            }
        }

        void Controller::reset_thermal()
        {
            mResetThermal = true ;
            mOmegaPicard2 = mOmega0 ;
            mOmegaNewton2 = mOmega0 ;
            mEquation2->reset_fields() ;

            // same hole as reset_timestep: the fields are restored, but the
            // thermal dofs still hold the failed sub-step's last iterate
            if ( mKernel2 != nullptr )
            {
                mKernel2->dofmgr()->seed_dof_values() ;
            }

            // the restored state invalidates the thermal mixing history
            this->anderson_clear_thermal() ;

            mBestEpsilon2 = BELFEM_REAL_MAX ;
            mWatchdogOmegaPrev2 = BELFEM_QUIET_NAN ;
            mBestEpsilonIteration2 = 0 ;

            mIteration2 = 0 ;

            // re-arm the thermal first-promotion latch ( cf. reset_timestep )
            mFirstFlip2 = true ;

            mCouplingFactor *= 2.0 ;

            mTime2 = mTime02 ;

            if( mKernel2 != nullptr )
            {
                mMesh2->time_stamp() = mTime2 ;
            }

        }

        void Controller::adjust_timestep()
        {
            // Ensure adjustments only happen after a couple of iterations
            if ( mIteration0 > 0 )
            {

                // Size the next step to keep the solver near its target iteration
                // count: cheap steps grow it, expensive ones shrink it ( Messe et al.
                // 2023, Section 4 ). Non-converged steps are halved separately
                // in reset_timestep().
                // both branches divide by these ( guarded by the mIteration0
                // gate in practice )
                BELFEM_ASSERT( mIteration > 0 && mIterationTarget > 0,
                    "adjust_timestep called with a zero iteration count or target" );

                real tPhi ;
                if ( mUseLegacyTimestepControl )
                {
                    // legacy: memoryless I-control with gain 0.5; the sqrt
                    // damps the response to iteration noise
                    tPhi = std::sqrt( static_cast< real >( mIterationTarget ) / mIteration ) ;
                }
                else
                {
                    // shift the cost-error history and load this step's
                    // sample; the max() floors a degenerate zero sample so it
                    // cannot poison the P/D terms of the following steps
                    mCtrlErr2 = mCtrlErr1 ;
                    mCtrlErr1 = mCtrlErr0 ;
                    mCtrlErr0 = static_cast< real >( std::max( mIteration, 1u ) )
                              / static_cast< real >( mIterationTarget ) ;

                    // multiplicative PID ( Valli, Carey & Coutinho 2002,
                    // CNM 18:131 ): phi = (e1/e0)^kP * (1/e0)^kI
                    //                   * (e1^2/(e0*e2))^kD
                    // The memory smooths the quantized iteration-count signal
                    // that the memoryless sqrt rule overreacts to
                    tPhi = std::pow( mCtrlErr1 / mCtrlErr0, mCtrlKp )
                         * std::pow( 1.0 / mCtrlErr0, mCtrlKi )
                         * std::pow( ( mCtrlErr1 * mCtrlErr1 )
                                     / ( mCtrlErr0 * mCtrlErr2 ), mCtrlKd );

                    // a recent cut vetoes growth until the hold expires, so
                    // the controller cannot re-climb the cliff it just fell
                    // off ( note the first accepted step after a cut skips
                    // this method entirely via the mIteration0 gate above )
                    if ( mPostFailureHold > 0 )
                    {
                        tPhi = std::min( tPhi, 1.0 ) ;
                        --mPostFailureHold ;
                    }
                }

                // growth clamp: 1.5x is safe against the BDF2 zero-stability
                // ratio bound ( 1+sqrt(2), Grigorieff 1983; Hairer & Wanner ).
                // The admissible ratios shrink with the BDF order, so at
                // higher ACTIVE orders the step must grow more carefully —
                // shrinking ( 0.5x ) is always benign
                uint tOrderActive = mEquation->order_active() ;
                if ( mEquation2 != nullptr )
                {
                    tOrderActive = std::max( tOrderActive, mEquation2->order_active() );
                }
                const real tPhiMax = tOrderActive >= 4 ? 1.2 :
                                     tOrderActive == 3 ? 1.4 : 1.5 ;

                tPhi = std::clamp( tPhi, 0.5, tPhiMax ) ;

                // Adapt the number of thermal time steps between magnetic ones
                if ( mIteration2 < mMaxNumIterations2/4 )
                {
                    mCouplingFactor = ceil( mCouplingFactor/2) ;
                }

                // Compute next timestep ensuring it doesn't exceed min/max limit
                mDeltaTime = std::clamp( mDeltaTime*tPhi, mDeltaTimeMin,mDeltaTimeMax);
                mCouplingFactor = std::max(mCouplingFactor, 1.0) ;

                //Remember previous mDeltaTime if it was changed due to save point
                if (!std::isnan( mDeltaTimeTemporary ))
                {
                    //Reset the temporary time stepping value
                    mDeltaTime = mDeltaTimeTemporary ;
                    mDeltaTimeTemporary = BELFEM_QUIET_NAN ;
                }

                // Check if we need to adapt the time step to make a save file at the next time step
                if ( mSaveEvery != 0.0 )
                {
                    // Check if we go pass a time to be saved
                    real tValue1 = mTime/mSaveEvery ;
                    real tValue2 = (mTime+mDeltaTime)/mSaveEvery ;
                    real tValue1Floor = std::floor(tValue1) ;
                    real tValue2Floor = std::floor(tValue2) ;

                    //Catch X.99999 floating point cases
                    if (tValue1-tValue1Floor > 1.0 - 1e-10)
                    {
                        tValue1Floor += 1 ;
                    }

                    if (tValue2-tValue2Floor > 1.0 - 1e-10)
                    {
                        tValue2Floor += 1 ;
                    }

                    //Adapt mDeltaTime if we pass a save point
                    if (tValue1Floor != tValue2Floor )
                    {
                        mDeltaTimeTemporary = mDeltaTime ; //Save the time steps for next iteration
                        mDeltaTime = std::clamp(tValue2Floor*mSaveEvery - mTime, mDeltaTimeMin,mDeltaTimeMax) ;
                        mLastSave = mSave ;
                        mSave = true ; //Save next time step
                    }
                    else
                    {
                        mLastSave = mSave ;
                        mSave = false ; //Don't save next time step
                    }
                }
                else
                {
                    mLastSave = mSave ;
                    mSave = true ; //Save next time step
                }

            }
        }


        void
        Controller::compute_circuit_current()
        {
            uint tCountCircuit = 0 ;
            for (PhysicalBoundaryCondition * tBC : mKernel->boundary_conditions())
            {
                if (tBC->type() == BoundaryConditionType::Current ||
                    tBC->type() == BoundaryConditionType::CircuitCurrent ||
                    tBC->type() == BoundaryConditionType::Voltage)
                {
                    tCountCircuit++;
                }
            }

            // get the abstract nodes
            Cell< Dof * > & tAbstractDofs = mKernel->dofmgr()->abstract_dofs();
            uint tCount = 0 ;
            for (PhysicalBoundaryCondition * tBC : mKernel->boundary_conditions())
            {
                if(tBC->type() == BoundaryConditionType::CircuitVoltage)
                {
                    // fix current in the circuit
                    real tV = mCircuit->voltage( tCount );

                    real tI = mKernel->dofmgr()->dof(tAbstractDofs(tCountCircuit++)->id())->value() ;
                    mCircuit->set_current_and_voltage(tCount++, tI, tV ) ;
                }
            }
        }

        void
        Controller::finalize( const bool aPostProcess )
        {
            mDeltaTime0 = mDeltaTime ;
            mIterationTime = mTimer->stop() ;


            // conditioning, per timestep ( Christian, 2026-08-29 ): the MUMPS
            // COND pair, if asked for, was already sampled on the FIRST
            // iterate of this step by the arm/capture pair ; the eigen
            // estimate, if asked for, is taken HERE, at the END of the step,
            // against the converged matrix -- which is still intact after the
            // last solve, and is the operator a user intuitively expects the
            // number to describe
            if ( mComputeConditioning || mComputeConditioning2 )
            {
                this->compute_conditioning() ;
            }

            // the J/Jc field is written by the Maxwell postprocessor and
            // nowhere else, so the footer must know whether THIS step
            // refreshed it or is about to report the last frame it saw
            mPostProcessed = mKernel->dofmgr()->postprocessors().size() > 0
                          && aPostProcess ;

            if ( mPostProcessed )
            {
                Timer tTimer ;
                mKernel->dofmgr()->postprocess() ;

                mPostprocesingTime = tTimer.stop();
            }

            // Update the FEM component in the circuit
            if ( mCircuit != nullptr && mCommRank == 0)
            {
                mCircuit->save_timestep() ;

                this->compute_circuit_current() ;

            }

            // Get the full LHS for voltage computation
            if ( mCommRank == 0 )
            {
                //Get previous LHS
                mLHS0 = mLHS ;

                // Get the full LHS of current time step
                mKernel->dofmgr()->full_lhs( mLHS );

                // Initialize mLHS0 if this is the first call
                if ( mLHS0.length() == 0 )
                {
                    mLHS0 = mLHS ;
                }
            }

            // a completed timestep clears the consecutive-solver-failure
            // count -- unconditionally, not only when the timestep adapts
            // ( with "adapt timestep : off" the count would otherwise
            // accumulate sporadic failures across the whole run, Grok SQ5 )
            if ( ! mReset )
            {
                mSolverFailCount = 0 ;

                // an accepted step ends the floor regime: the iteration
                // budgets return to their deck values ( the escalation that
                // raised them lives in reset_timestep )
                if ( mFloorRetries > 0 )
                {
                    mFloorRetries      = 0 ;
                    mMaxNumIterations  = mMaxNumIterationsDeck ;
                    mMaxNumIterations2 = mMaxNumIterations2Deck ;
                    mWatchdogWindow    = mWatchdogWindowDeck ;
                    mWatchdogWindow2   = mWatchdogWindow2Deck ;
                }
            }

            if ( ! mReset && mAdaptTimestep )
            {
                this->adjust_timestep() ;
            }

            // every linear solve churns hundreds of MiB through the
            // allocator ( the CSR copy, the numeric refactor, and the Krylov
            // workspace ), and glibc's per-thread arenas keep a few percent
            // of it instead of returning it to the OS. Over a long transient
            // that ratchet reached 61 GiB resident and the OOM killer took
            // the run twice. Nothing here is leaked -- three independent
            // sweeps found no unfreed object -- so the cure is to hand the
            // free heap back once per accepted step, where the cost ( a heap
            // walk, milliseconds ) is invisible against a step measured in
            // minutes. Not called on a rejected step: the retry runs
            // immediately and would only re-fault the same pages back in
#if defined( __GLIBC__ ) && ! defined( __APPLE__ )
            if ( ! mReset )
            {
                malloc_trim( 0 );
            }
#endif

            if( mCommRank == 0 && gLog.info_level() > 0 )
            {
                this->print_footer() ;
            }

            // the per-timestep exchange boundary -- everything this
            // step sent must have been consumed by now. This makes finalize
            // collective when assertions are active: every rank must enter it
            comm_drain_check( "Controller::finalize" );
        }

        void
        Controller::compute_conditioning()
        {
            // MUMPS delivers its COND1 / COND2 pair at the first iterate
            // through its own error analysis ( see capture_conditioning_*,
            // gated by "mumps error analysis" ). The eigenvalue estimate is a
            // SEPARATE quantity behind a SEPARATE key: this function serves
            // "compute conditioning" only, and every field with THAT flag set
            // gets the eigen treatment once per timestep.
            //
            // That fallback returns kappa = |lambda_max| / |lambda_min|, and
            // the small end can be OUT OF REACH: ARPACK accepts a Ritz value
            // once its error bound drops below tol * |lambda|, and that bound
            // cannot fall below the backward error of the matvec, about
            // eps_mach * |lambda_max|. So it needs
            //
            //     tol  >  eps_mach * kappa
            //
            // and a mixed h-phi system reaches kappa ~ 1e17, which would ask
            // for a relative tolerance above 1. EigenValues then warns and
            // hands back a NaN rather than aborting -- the footer shows n/a
            // and the run continues. Where kappa is smaller the estimate is
            // perfectly good, so the path stays live.
            //
            // MUMPS remains the recommended solver for this diagnostic ; the
            // deck parser says so. PETSc and SuperLU could also supply the
            // Arioli/Demmel/Duff numbers -- see
            // todo/conditioning_diagnostic_backends.md
            // each field is asked separately, against its own EigenValues --
            // the DofManager builds one per instance, so the two kernels
            // already own one each and nothing is shared between them.
            // mEigenAnalysisTime is the SUM: asking for both numbers costs
            // two estimates, and the footer reports the total rather than
            // growing a second timing line
            mEigenAnalysisTime = 0 ;

            // THIS FUNCTION OWNS SLOT 0 ONLY -- the eigen estimate. Slots 1,
            // 2 and 3 ( the MUMPS COND1 / COND2 pair, and the omega2 that says
            // whether COND2's term contributes at all ) belong to
            // capture_conditioning_*, which samples them on the FIRST iterate
            // and then DISARMS the error analysis. Re-reading them here would
            // return BELFEM_SIGNALING_NAN, because get_cond1 / get_cond2 hand
            // back a NaN unless ICNTL(11) is still armed ( cl_SolverMUMPS.cpp
            // ) -- so a re-read silently overwrites the captured numbers with
            // n/a. That defect was live on 2026-08-29 and is what this
            // comment exists to prevent recurring.
            //
            // THE NUMBERS ARE NOT THE SAME QUANTITY and must never be
            // compared: MUMPS COND1 is a componentwise 1-norm condition
            // estimate for the solved system WITH ITS ACTUAL RIGHT-HAND SIDE,
            // on the ORIGINAL matrix -- NOT the equilibrated one. MUMPS scales
            // for the factorization ( ICNTL(8) = 77 ), but its error analysis
            // hands DMUMPS_SOL_LCOND an identity weight vector, under its own
            // comment "Notice that D is always the identity"
            // ( MUMPS 5.9.1 dsol_driver.F:6101-6111 ), so no scaling enters
            // the estimate. Saying "of the scaled matrix" here was wrong and
            // was corrected 2026-08-29 -- while slot 0
            // returns |lambda_max| / |lambda_min|, a 2-norm property of the
            // matrix alone -- and that ratio is kappa_2 only for a NORMAL
            // matrix. It is therefore kappa_2 for the symmetric thermal
            // Jacobian and merely a spectral ratio for the nonsymmetric h-phi
            // one; since 2026-08-30 the footer prints the honest common name,
            // |lambda|max/|lambda|min, for both fields.
            // Measured on the same tapestack3d thermal operator:
            // kappa_2 = 2.96e7 against a COND1 in the 1e4-1e5 range. Both are
            // right ; they answer different questions
            if ( mComputeConditioning )
            {
                mTimer->reset() ;

                mConditionNumbers1( 0 ) = mKernel->dofmgr()->eigen_values()->compute_conditioning();

                mEigenAnalysisTime += mTimer->stop() ;
            }

            if ( mComputeConditioning2
                 && mKernel2 != nullptr
                 && mKernel2->dofmgr()->solver() != nullptr )
            {
                mTimer->reset() ;

                // R3b tripwire: mSymmetric is written ONCE, at attach, by
                // setup_thermal_eigen(). If that write is ever lost -- a new
                // attach path that skips the helper, or an EigenValues object
                // rebuilt after attach -- the thermal field silently drops to
                // the nonsymmetric driver. Nothing in the OUTPUT would show
                // it: both drivers converge to the same ratio on a symmetric
                // matrix, so only the flag itself can be checked
                BELFEM_ASSERT( mKernel2->dofmgr()->eigen_values()->is_symmetric(),
                    "thermal EigenValues lost mSymmetric after attach" );

                mConditionNumbers2( 0 ) = mKernel2->dofmgr()->eigen_values()->compute_conditioning();

                mEigenAnalysisTime += mTimer->stop() ;
            }
        }

//------------------------------------------------------------------------------

        /**
         * whether the MUMPS ADD COND2 row belongs in the footer.
         *
         * MUMPS splits the componentwise backward error by row -- rows with a
         * trustworthy denominator into omega1, the rest into omega2 -- and
         * builds the forward error as omega1*COND1 + omega2*COND2. The row is
         * printed when its term can contribute, i.e. when omega2 is nonzero,
         * and ALSO when omega2 is not a number, which means nothing was
         * captured this timestep and the honest n/a belongs in the footer
         * exactly as before.
         *
         * The case being dropped is omega2 == 0. Two things live there and it
         * is worth being precise, because an earlier version of this comment
         * conflated them:
         *
         *   - usually, no row fell into the second category at all, MUMPS
         *     skipped the COND2 estimator, and what it returns is the 1.0 both
         *     condition numbers were initialized to ( dsol_aux.F:959-964,
         *     :1017 ). Printing that reads as "the second condition number is
         *     one", a measurement nobody made. This is the case Christian
         *     asked to stop printing;
         *   - rarely, the category was NOT empty and COND2 really was
         *     estimated, but every exceptional row contributed zero -- TAU == 0
         *     or an exactly vanishing residual, since IW( i, 1 ) = 2 is
         *     assigned outside the guard that updates omega2
         *     ( dsol_aux.F:891-900 ). A real number is then suppressed.
         *
         * That second case is accepted rather than overlooked. It costs
         * nothing in meaning: at omega2 == 0 the term is zero either way, so
         * the number cannot move the error bound it exists to explain. And it
         * cannot be separated from the first from outside MUMPS -- the
         * discriminator is IW, which is internal, and neither omega2 nor COND2
         * alone recovers it.
         */
        static bool
        cond2_row_wanted( const real aOmega2 )
        {
            // exact comparison on purpose: MUMPS assigns the literal ZERO and
            // never arithmetic that lands near it, so this detects an
            // untouched initializer, not a small number
            return ! ( std::isfinite( aOmega2 ) && aOmega2 == 0.0 ) ;
        }

//------------------------------------------------------------------------------

        void
        Controller::setup_thermal_eigen()
        {
            if ( mKernel2 == nullptr ) return ;

            // The thermal Jacobian is symmetric by construction, the magnetic
            // h-phi one is not ( Christian, 2026-08-28 ), so the thermal field
            // gets the Lanczos driver and the magnetic field keeps the default
            // nonsymmetric one. Declared here rather than measured: detecting
            // it would mean comparing A against A^T across a distribution that
            // may already be transposed column blocks.
            //
            // Written ONCE per attach rather than once per timestep. That is
            // safe because the object outlives the call: mEigenValues is
            // assigned only in the DofManager constructor and nothing rebuilds
            // it -- EigenValues::reset() clears mMatrixFlag alone,
            // DofManager::reset() does not touch it, and link_matrix() rebinds
            // the matrix rather than the object. EigenValues itself already
            // states the model: mSymmetric "is set once at setup, never
            // derived from matrix values, so it needs no broadcast".
            //
            // Unconditional on purpose: it selects an ARPACK driver and is
            // inert unless the eigen diagnostic actually runs, so it does not
            // belong behind either diagnostic flag. NOT rank-guarded -- the
            // value is rank-identical by construction and a guard here would
            // desynchronize the drivers across ranks
            mKernel2->dofmgr()->eigen_values()->set_symmetric( true );
        }

        // ----------------------------------------------------------------
        // Diagnostic-configuration warnings, one function per field.
        //
        // Split by field rather than composed into one message because the
        // two fields are validated at DIFFERENT TIMES: the magnetic kernel
        // exists when set_params runs, the thermal one usually does not. The
        // pre-2026-08-30 code composed "magnetic and thermal" in a single
        // sprint at parse time, which is why its thermal branch -- guarded on
        // a kernel that is null there -- never fired on any coupled deck.
        // ----------------------------------------------------------------

        void
        Controller::check_magnetic_diagnostics()
        {
            if ( mCommRank != 0 || mKernel->dofmgr()->solver() == nullptr ) return ;

            const bool tIsMumps =
                mKernel->dofmgr()->solver()->type() == SolverType::MUMPS ;

            // the eigen estimate without MUMPS: allowed, but expensive and
            // blind at the small end
            if ( mComputeConditioning && ! tIsMumps )
            {
                message( InfoLevel::Minimal,
                    "compute conditioning is on without MUMPS for the magnetic field, which is not\n"
                    "                 recommended. The number comes from an eigenvalue estimate: it costs an\n"
                    "                 extra solve per timestep, and it cannot resolve the small end once\n"
                    "                 |lambda|max/|lambda|min exceeds about tol / eps_mach, where it reports\n"
                    "                 n/a instead.\n" );
            }

            // the error analysis without MUMPS: a silent no-op, hence the key
            // names MUMPS
            if ( mMumpsErrorAnalysis && ! tIsMumps )
            {
                message( InfoLevel::Minimal,
                    "mumps error analysis is on for the magnetic field, whose solver is not MUMPS.\n"
                    "                 The key does nothing there: the Arioli-Demmel-Duff COND1 / COND2 pair\n"
                    "                 comes from MUMPS's own error analysis ( ICNTL(11) ), and no other\n"
                    "                 library in this build supplies it.\n" );
            }
        }

        void
        Controller::check_thermal_diagnostics()
        {
            // one-shot: both attach paths call this, and a caller that uses
            // the two-argument constructor AND set_thermal_kernel would
            // otherwise warn twice. Both call sites are collective, so the
            // flag stays rank-identical
            // BOTH prerequisites, not first-caller-wins: a call before
            // set_params must not consume the one-shot against default flags
            // ( set_thermal_kernel before set_params is a legal public order,
            // and set_params re-calls this when the kernel is already live )
            if ( mThermalDiagnosticsChecked
                 || mKernel2 == nullptr
                 || ! mParamsSet ) return ;
            mThermalDiagnosticsChecked = true ;

            if ( mCommRank != 0 || mKernel2->dofmgr()->solver() == nullptr ) return ;

            const bool tIsMumps =
                mKernel2->dofmgr()->solver()->type() == SolverType::MUMPS ;

            if ( mComputeConditioning2 && ! tIsMumps )
            {
                message( InfoLevel::Minimal,
                    "compute conditioning is on without MUMPS for the thermal field, which is not\n"
                    "                 recommended. The number comes from an eigenvalue estimate: it costs an\n"
                    "                 extra solve per timestep, and it cannot resolve the small end once\n"
                    "                 |lambda|max/|lambda|min exceeds about tol / eps_mach, where it reports\n"
                    "                 n/a instead.\n" );
            }

            if ( mMumpsErrorAnalysis2 && ! tIsMumps )
            {
                message( InfoLevel::Minimal,
                    "mumps error analysis is on for the thermal field, whose solver is not MUMPS.\n"
                    "                 The key does nothing there: the Arioli-Demmel-Duff COND1 / COND2 pair\n"
                    "                 comes from MUMPS's own error analysis ( ICNTL(11) ), and no other\n"
                    "                 library in this build supplies it.\n" );
            }
        }

        void
        Controller::arm_conditioning_magnetic()
        {
            // re-armed at every initialize_timestep, captured ( and disarmed )
            // on the FIRST solve of the step. A step that never solves -- a
            // certified predictor accept -- simply keeps the arm for the next
            // one, and its footer honestly shows n/a from the per-step reset
            if ( ! mMumpsErrorAnalysis
                 || mKernel->dofmgr()->solver() == nullptr
                 || mKernel->dofmgr()->solver()->type() != SolverType::MUMPS ) return ;

            mKernel->dofmgr()->solver()->set_mumps_error_analysis(
                MumpsErrorAnalysis::Full );
        }

        void
        Controller::arm_conditioning_thermal()
        {
            if ( ! mMumpsErrorAnalysis2 || mKernel2 == nullptr ) return ;

            // NOTHING but MUMPS arming lives here. The thermal eigen driver is
            // selected in setup_thermal_eigen(), called once per attach --
            // until 2026-08-30 the set_symmetric() call sat in THIS function,
            // between the flag guard and the MUMPS guard, so re-gating the
            // function onto the error-analysis flag would have silently sent
            // the thermal field through the nonsymmetric driver on every deck
            // that did not ask for MUMPS error analysis
            if ( mKernel2->dofmgr()->solver() == nullptr
                 || mKernel2->dofmgr()->solver()->type() != SolverType::MUMPS ) return ;

            mKernel2->dofmgr()->solver()->set_mumps_error_analysis(
                MumpsErrorAnalysis::Full );
        }

        void
        Controller::capture_conditioning_magnetic()
        {
            if ( ! mMumpsErrorAnalysis
                 || mKernel->dofmgr()->solver() == nullptr
                 || mKernel->dofmgr()->solver()->type() != SolverType::MUMPS ) return ;

            mConditionNumbers1( 1 ) = mKernel->dofmgr()->solver()->wrapper()->get_cond1() ;
            mConditionNumbers1( 2 ) = mKernel->dofmgr()->solver()->wrapper()->get_cond2() ;
            mConditionNumbers1( 3 ) = mKernel->dofmgr()->solver()->wrapper()->get_omega2() ;

            // one sample per TIMESTEP, taken on the FIRST iterate: arming is
            // not free ( ICNTL(11) = 1 adds the residual/omega statistics AND
            // the Hager reverse-communication solves of the condition
            // estimator, ITMAX = 5, twice if the second row category is
            // non-empty -- it does NOT add iterative refinement, which is
            // ICNTL(10) and already set to 20 for a single RHS
            // [ cl_SolverMUMPS.cpp, mumpstools.f90 ]; the earlier claim that
            // ICNTL(11) runs refinement was refuted 2026-08-29 ), we cannot
            // know which iterate
            // will be the last, and the native estimate's drift across a step
            // is negligible next to the estimate's own error ( Christian,
            // 2026-08-29 ). Disarmed here ; initialize_timestep re-arms on the
            // next step
            mKernel->dofmgr()->solver()->set_mumps_error_analysis(
                MumpsErrorAnalysis::None );
        }

        void
        Controller::capture_conditioning_thermal()
        {
            if ( ! mMumpsErrorAnalysis2
                 || mKernel2 == nullptr
                 || mKernel2->dofmgr()->solver() == nullptr
                 || mKernel2->dofmgr()->solver()->type() != SolverType::MUMPS ) return ;

            mConditionNumbers2( 1 ) = mKernel2->dofmgr()->solver()->wrapper()->get_cond1() ;
            mConditionNumbers2( 2 ) = mKernel2->dofmgr()->solver()->wrapper()->get_cond2() ;
            mConditionNumbers2( 3 ) = mKernel2->dofmgr()->solver()->wrapper()->get_omega2() ;

            // see the magnetic twin: one sample per timestep, on the first
            // iterate, re-armed by the next initialize_timestep
            mKernel2->dofmgr()->solver()->set_mumps_error_analysis(
                MumpsErrorAnalysis::None );
        }

        void
        Controller::check_iterative_solver_headroom(
                DofManager * aDofMgr,
                const real   aNonlinTol,
                const char * aFieldName ) const
        {
            if ( aDofMgr == nullptr ) return ;

            Solver * tSolver = aDofMgr->solver() ;
            if ( tSolver == nullptr ) return ;

            // only PETSc is gated here. The original premise for exempting
            // the direct libraries — "a factorization's delivered accuracy
            // does not depend on the stated tolerance" — was FALSIFIED for
            // STRUMPACK: its outer GMRES stops at its tolerances.
            // The exemption stands anyway, but NOT because a loose relative
            // pairing is safe there — the 1e-8 / 1e-11 pairing this comment
            // once called "proven" collapsed the tapestack3d timestep, and
            // the 2026-08-18 A/B showed its Picard residual was the linear
            // exit test, not the physics. The exemption survives because
            // the WORKING pairing is itself lin > nonlin: with the class
            // default 1e-10 always applied ( since 2026-08-18 ) against a
            // 1e-11 nonlinear target, factorization-preconditioned GMRES
            // overshoots its exit test by decades ( observed ~1e-15 ), so
            // the lin <= nonlin predicate this gate applies to PETSc would
            // refuse a configuration that demonstrably works. A deck whose
            // STRUMPACK residual crawls instead of dives states a tighter
            // 'relative tolerance' — the expert override, not this gate
            if ( tSolver->type() != SolverType::PETSc ) return ;

            const real tLinTol = tSolver->parameters().relative_tolerance() ;

            // deck and defaults are parsed identically on every rank, so
            // this either passes or aborts collectively
            BELFEM_ERROR( tLinTol <= aNonlinTol,
                "impossible tolerance pairing for the %s field: the nonlinear loop demands a\n"
                "residual of %g, but its PETSc linear solver stops once it reaches a relative\n"
                "tolerance of %g. An iterative solve delivers no more accuracy than it is asked\n"
                "for, so the nonlinear residual can never fall below the linear tolerance.\n"
                "Tighten 'relative tolerance' in the linear %s section to at least the nonlinear\n"
                "value ( two to three decades below it is recommended ), or relax the nonlinear\n"
                "'tolerance'.",
                aFieldName, aNonlinTol, tLinTol, aFieldName );
        }

        void
        Controller::check_compression_headroom(
                DofManager * aDofMgr,
                const real   aNonlinTol,
                const char * aFieldName ) const
        {
            if ( aDofMgr == nullptr ) return ;
            Solver * tSolver = aDofMgr->solver() ;
            if ( tSolver == nullptr ) return ;

            // only the libraries that actually consume the compression
            // switch — a petsc block carrying a stray blr key compresses
            // nothing and must not earn this warning
            if ( tSolver->type() != SolverType::MUMPS
                 && tSolver->type() != SolverType::STRUMPACK ) return ;

            const SolverParameters & tParams = tSolver->parameters() ;
            if ( tParams.compression_method() != CompressionMethod::BLR )
                return ;

            // with BLR the drop tolerance IS the delivered linear
            // accuracy. Two tiers:
            // zero-or-negative headroom ( cutoff >= nonlinear tolerance )
            // is the same impossibility the PETSc gate refuses — hard
            // error; less than two decades but positive is the legitimate
            // memory-bound trade — warning. Both computed on the
            // EFFECTIVE cutoff: erroring only on a STATED value would
            // make stating your cutoff riskier than omitting it, so the
            // inherited default is held to the same bar and the message
            // names its provenance. Parameters are synchronized and the
            // tolerance is parsed on every rank — the verdict is uniform,
            // the abort collective
            const real tCutoff = tParams.compression_cutoff() ;
            if ( tCutoff <= 1e-2 * aNonlinTol ) return ;

            BELFEM_ERROR( tCutoff < aNonlinTol,
                "impossible compression pairing for the %s field: the blr cutoff %g%s is not\n"
                "below the nonlinear tolerance %g. The drop tolerance IS the delivered linear\n"
                "accuracy - the nonlinear loop can never reach its target through a factorization\n"
                "this lossy ( the 2026-07-06 failure class ). State a tighter 'compression\n"
                "cutoff' ( 2-3 decades below the tolerance ), or remove 'compression scheme'.",
                aFieldName,
                tCutoff,
                tParams.have_compression_cutoff() ? "" : " ( the inherited default )",
                aNonlinTol );

            // deliberately UNGUARDED rank-0 print: this is a correctness
            // notice, and -v 0 ( Silent ) must not swallow it the way it
            // swallows the info-gated banners
            if ( mCommRank == 0 )
            {
                std::cout << sprint(
                    "\n    WARNING: %s field combines blr compression ( cutoff %g%s ) with a\n"
                    "    nonlinear tolerance of %g - less than two decades of headroom. The drop\n"
                    "    tolerance IS the delivered linear accuracy: the nonlinear loop may floor\n"
                    "    on the compression, invisibly ( this stalled Newton on 2026-07-06 ).\n"
                    "    State a tighter 'compression cutoff', or remove 'compression scheme'.\n",
                    aFieldName,
                    tCutoff,
                    tParams.have_compression_cutoff() ? "" : ", the inherited default",
                    aNonlinTol ) << std::endl ;
            }
        }

        void
        Controller::print_header()
        {
            // node: do not change the width of this box.
            // it must look good in 80 columns!
            std::cout << std::endl ;
            std::cout <<     "   ┌────────────────────┬─────────────────────────┬────────────────────────┐"<< std::endl ;
            // label the box with the order that is ACTUALLY running this step.
            // During the startup ramp a BDF-p run executes BDF1, BDF2, ... until
            // p states of history exist, so printing the configured p would
            // misreport the first p-1 steps. order_active() is the same value
            // the step-growth clamp reads below, so the label cannot drift away
            // from the scheme in use. The order can also drop back mid-run
            // whenever the history is invalidated, which the label now shows.
            string tScheme ;
            switch ( mTimeStepping )
            {
                case( EulerMethod::BackwardDifference1 ) :
                case( EulerMethod::BackwardDifference2 ) :
                case( EulerMethod::BackwardDifference3 ) :
                case( EulerMethod::BackwardDifference4 ) :
                case( EulerMethod::BackwardDifference5 ) :
                {
                    // order_PENDING, not order_active: the coefficients are
                    // recomputed lazily on the first element of the step, so
                    // at header time order_active() would still be the
                    // previous step's value and the label would lag by one.
                    uint tOrder = mEquation->order_pending() ;
                    if ( mEquation2 != nullptr )
                    {
                        tOrder = std::max( tOrder, mEquation2->order_pending() );
                    }
                    tScheme = sprint( "BDF%u", ( unsigned int ) tOrder );
                    break ;
                }
                case( EulerMethod::ForwardExplicit )     : tScheme = "Explicit" ; break ;
                case( EulerMethod::CrankNicolson )       : tScheme = "CN"       ; break ;
                case( EulerMethod::Galerkin )            : tScheme = "Galerkin" ; break ;
                default                                  : tScheme = "Timestep" ; break ;
            }
            string tFormat = "   │  %-8s %8u │  t : %14.4f ms  │   Δ t = %10.4f ms  │" ;
            string tMessage = sprint( tFormat.c_str(), tScheme.c_str(), mRunningTimeStep, mTime * 1000, mDeltaTime * 1000 );
            std::cout << tMessage << std::endl ;
            std::cout <<     "   ├────────────────────┴─────────────────────────┴────────────────────────┤"<< std::endl ;
        }



        void
        Controller::print_line( const real aOmega, const real aOmega2 )
        {

            real dB = 10. * std::log10( mEpsilon );
            string tFormat = "   │ Magnetic" ;
            tFormat += mEquation->algorithm() == SolverAlgorithm::Picard ? " Picard " : " Newton " ;
            tFormat += "%5u, residual %s (%7.2f dB ), relax %7.5f │" ;

            string tMessage = sprint( tFormat.c_str(), mIteration,
                residual_string( mEpsilon ).c_str(), dB, aOmega );

            std::cout << tMessage << std::endl ;

            if (mKernel2 != nullptr)
            {
                dB = 10. * std::log10( mEpsilon2 );
                tFormat = "   │ Thermal " ;
                tFormat += mEquation2->algorithm() == SolverAlgorithm::Picard ? " Picard " : " Newton " ;
                tFormat += "%5u, residual %s (%7.2f dB ), relax %7.5f │" ;

                tMessage = sprint( tFormat.c_str(), mIteration,
                    residual_string( mEpsilon2 ).c_str(), dB, aOmega2 );

                std::cout << tMessage << std::endl ;
            }
        }

        void
        Controller::print_line_magnetic( const real aOmega )
        {

            real dB = 10. * std::log10( mEpsilon );
            string tFormat = "   │ Magnetic" ;
            tFormat += mEquation->algorithm() == SolverAlgorithm::Picard ? " Picard " : " Newton " ;
            tFormat += "%5u, residual %s (%7.2f dB ), relax %7.5f │" ;

            string tMessage = sprint( tFormat.c_str(), mIteration,
                residual_string( mEpsilon ).c_str(), dB, aOmega );

            std::cout << tMessage << std::endl ;
        }

        void
        Controller::print_line_thermal( const real aOmega )
        {

            real dB = 10. * std::log10( mEpsilon2 );
            string tFormat = "   │ Thermal " ;
            tFormat += mEquation2->algorithm() == SolverAlgorithm::Picard ? " Picard " : " Newton " ;
            tFormat += "%5u, residual %s (%7.2f dB ), relax %7.5f │" ;

            string tMessage = sprint( tFormat.c_str(), mIteration2,
                residual_string( mEpsilon2 ).c_str(), dB, aOmega );

            std::cout << tMessage << std::endl ;
        }

        void
        Controller::print_thermal_stall_warning()
        {
            if ( mCommRank != 0 || gLog.info_level() == 0 )
            {
                return ;
            }

            // max T on the master rank's partition ( diagnostic only; in
            // parallel the hot spot may live on another rank )
            real tTmax = BELFEM_QUIET_NAN ;
            if ( mMesh2 != nullptr && mMesh2->field_exists( "T" ) )
            {
                const Vector< real > & tT = mMesh2->field_data( "T" );
                tTmax = 0.0 ;
                for ( index_t k=0; k<tT.length(); ++k )
                {
                    tTmax = std::max( tTmax, tT( k ) );
                }
            }

            // smallest material table ceiling among the thermal blocks:
            // T at the ceiling freezes the properties and zeroes the
            // dT-derivatives, which is exactly the bit-flat stall signature
            real tTceil = BELFEM_REAL_MAX ;
            for ( Block * tBlock : mKernel2->dofmgr()->blocks() )
            {
                const Material * tMat = tBlock->material() ;
                if ( tMat != nullptr && tMat->is_constant( MaterialProperty::T_max ) )
                {
                    tTceil = std::min( tTceil,
                        tMat->constant_property( MaterialProperty::T_max ) );
                }
            }

            std::cout << sprint( "   │%-71s│",
                " WARNING: thermal residual stalled; accepting this timestep" )
                << std::endl ;

            string tLine = sprint( " thermal residual %8.2e ( %.2f dB )",
                mEpsilon2, 10. * std::log10( mEpsilon2 ) );
            std::cout << sprint( "   │%-71s│", tLine.c_str() ) << std::endl ;

            if ( ! std::isnan( tTmax ) )
            {
                tLine = sprint( " max T ( master rank ) = %.2f K", tTmax );
                if ( tTceil < BELFEM_REAL_MAX )
                {
                    tLine += sprint( ", material T_max = %.2f K", tTceil );
                }
                std::cout << sprint( "   │%-71s│", tLine.c_str() ) << std::endl ;

                if ( tTceil < BELFEM_REAL_MAX && tTmax >= tTceil )
                {
                    std::cout << sprint( "   │%-71s│",
                        " T is pinned at the material table ceiling ( clamp active )" )
                        << std::endl ;
                }
            }
        }

        real
        Controller::get_Tmax()
        {
            BELFEM_ASSERT( comm_rank() == 0 , "Only root may call Controller::get_Tmax() " );
            return mKernel2 != nullptr ? max( mMesh->field_data( "T" )) : gTbulk ;
        }

        real
        Controller::get_Imax()
        {
            BELFEM_ASSERT( comm_rank() == 0 , "Only root may call Controller::get_Imax() " );
            real Imax = 0 ;
            for ( auto tBC : mKernel->boundary_conditions() )
            {
                Imax = std::max( Imax, std::abs( tBC->value() ));
            }
            return Imax ;
        }

        string
        Controller::excitation_cell()
        {
            BELFEM_ASSERT( comm_rank() == 0 , "Only root may call Controller::excitation_cell() " );

            // Which number belongs in this cell depends on what drives the deck.
            // A transport current reports amperes, as it always has. A deck with
            // NO current condition -- a magnetization problem driven purely by a
            // background field -- has no current to report, and printing the
            // background's value() under an "A" label states a field strength in
            // amperes: for a 1 T ramp that reads as tens of thousands of amperes
            // through a conductor carrying none.
            //
            // The background condition imposes H, so it is converted back to the
            // flux density the deck asked for, B = mu0*H. Decks that are neither
            // current- nor background-driven keep the legacy behaviour verbatim,
            // so no existing footer changes.
            //
            // NOTE a voltage condition's value() is in VOLTS and has always been
            // folded into this maximum under an "A" label. That mislabel is not
            // this change's to fix silently, and it survives UNCHANGED in the
            // fallback branch -- but be precise about where: a deck with a
            // voltage condition AND a background one now takes the tesla branch,
            // so it reports B where it used to report max( |V|, |H| ) as amperes.
            // Voltage is deliberately not counted as "a current is provided",
            // because value() is not an ampere there either
            bool tHaveCurrent = false ;
            bool tHaveBackground = false ;
            real tImax = 0.0 ;
            real tHmax = 0.0 ;
            real tAnyMax = 0.0 ;

            for ( auto tBC : mKernel->boundary_conditions() )
            {
                const real tValue = std::abs( tBC->value() );

                tAnyMax = std::max( tAnyMax, tValue );

                switch ( tBC->type() )
                {
                    case BoundaryConditionType::Current :
                    case BoundaryConditionType::CircuitCurrent :
                    {
                        tHaveCurrent = true ;
                        tImax = std::max( tImax, tValue );
                        break ;
                    }
                    case BoundaryConditionType::Background :
                    {
                        tHaveBackground = true ;
                        tHmax = std::max( tHmax, tValue );
                        break ;
                    }
                    default :
                    {
                        break ;
                    }
                }
            }

            // 20 characters per branch: 8 label + 9 value + 3 unit. The %9 is a
            // MINIMUM width, not a maximum, so an oversized number still widens
            // the cell and walks the box border -- unchanged for the ampere
            // branches, which already do this above 99999.999, and reached in
            // the tesla branch only at 1000 T. Printing B rather than H is what
            // keeps a normal background deck inside the box: the same field as
            // an H value is six digits wider
            if ( tHaveCurrent )
            {
                return sprint( " I_max :%9.3f A ", tImax );
            }
            else if ( tHaveBackground )
            {
                return sprint( " B_max :%9.5f T ", tHmax * constant::mu0 );
            }
            else
            {
                return sprint( " I_max :%9.3f A ", tAnyMax );
            }
        }

//------------------------------------------------------------------------------

        real
        Controller::get_JJCmax()
        {
            BELFEM_ASSERT( comm_rank() == 0 , "Only root may call Controller::get_JJCmax() " );
            if ( mMesh->field_exists( "JJCz" ))
            {

                if ( mMesh->number_of_dimensions() == 2 )
                {
                    mLastJJcMax = max( mMesh->field_data( "JJCz" ) ) ;
                }
                else
                {
                    const Vector< real > & JJCx = mMesh->field_data( "JJCx" );
                    const Vector< real > & JJCy = mMesh->field_data( "JJCy" );
                    const Vector< real > & JJCz = mMesh->field_data( "JJCz" );

                    real val = 0 ;
                    index_t n = mMesh->number_of_nodes() ;
                    for ( index_t k=0; k<n; k++ )
                    {
                        val = std::max( val, JJCx( k ) * JJCx( k ) + JJCy( k ) * JJCy( k ) + JJCz( k ) * JJCz( k )  );
                    }
                    mLastJJcMax = std::sqrt( val ) ;
                }
                return mLastJJcMax ;
            }
            else
            {
                return 0.0 ;
            }
        }

        void
        Controller::print_footer()
        {
            string tFormat ;
            string tMessage ;

            // the physics summary opens the footer: it is printed HERE, after
            // finalize() has run the postprocessor, because J/Jc only exists
            // once that pass has written it. Printed from the iteration loop
            // it would always have been one save behind, and in the segregated
            // case its T_max would have predated the thermal sub-steps
            this->print_physics_stats();

            double timecount = static_cast< double > ( mIterationTime ) * 0.001 ;
            unsigned int minutes = ( uint ) timecount / 60 ;
            unsigned int  seconds = ( uint ) timecount % 60 ;


            const bool tMagneticEigenFailed = mComputeConditioning  && ! std::isfinite( mConditionNumbers1( 0 ) ) ;

            // slot 0 of the THERMAL vector, negated: the banner fires when the
            // thermal EIGEN estimate is missing. Reading the magnetic vector,
            // or reading it unnegated, made this predicate report the opposite
            // of its own name ( fixed 2026-08-29 )
            const bool tThermalEigenFailed = mComputeConditioning2
                 && mKernel2 != nullptr
                 && mKernel2->dofmgr()->solver() != nullptr
                 && ! std::isfinite( mConditionNumbers2( 0 ) );

            if ( tMagneticEigenFailed || tThermalEigenFailed )
            {
                // names the EIGEN estimate, not "the conditioning number":
                // when the solver is MUMPS the ADD rows below still carry a
                // perfectly good COND1 / COND2, so a banner claiming the
                // conditioning of this system could not be computed would be
                // false at the moment it is printed
                if ( mKernel->dofmgr()->eigen_values()->backend_label() == "PARPACK" )
                {
                    std::cout << "   │ PARPACK failed to compute the eigenvalue conditioning estimate.       │" << std::endl ;
                }
                else
                {
                    std::cout << "   │ ARPACK failed to compute the eigenvalue conditioning estimate.        │" << std::endl ;
                }
                std::cout << "   ├───────────────────────────────────────────────────────────────────────┤" << std::endl ;
            }

            if ( mComputeConditioning
              || ( mComputeConditioning2
                   && mKernel2 != nullptr
                   && mKernel2->dofmgr()->solver() != nullptr ) )
            {
                // the SUM over the fields that needed an estimate
                tFormat  = "   │                    Time for eigenvalue analysis    :  %12u ms │" ;
                tMessage = sprint( tFormat.c_str(), mEigenAnalysisTime );
                std::cout << tMessage << std::endl ;
            }

            // mPostProcessed, not the mSave shadow: the row must report a
            // measurement that was actually taken, and mPostprocesingTime is
            // only written where the postprocessor ran
            if ( mPostProcessed )
            {
                tFormat = "   │                    Time for postprocessing         :  %12u ms │" ;
                tMessage = sprint( tFormat.c_str(), mPostprocesingTime );
                std::cout << tMessage << std::endl ;
            }

            // Each row is gated on ITS OWN diagnostic, and the closing tail
            // is unconditional. The outer if/else that used to wrap all of
            // this was dropped on 2026-08-30: it keyed the whole block on the
            // EIGEN flags, so a deck asking only for the MUMPS numbers got
            // none of its rows, and its two branches carried byte-identical
            // copies of the tail.
            //
            // The kappa glyph would be TWO BYTES in utf-8, so a %-Ns field
            // would pad by bytes and skew the frame: the label text is
            // written literally and only the fixed-width value is substituted.
            if ( mComputeConditioning )
            {
                // ONE label for both matrix classes. |lambda_max|/|lambda_min|
                // is what the code actually computes; it equals kappa_2 only
                // for a NORMAL matrix, so the old symmetric/nonsymmetric
                // branch printed "κ₂" for the thermal field and the ratio for
                // the magnetic one. Naming the ratio everywhere never
                // over-claims, needs no glossary, and is the label INC-213
                // pointed at when it recorded "κ₂" as part of a misreading
                tFormat  = "   │                    |λ|max/|λ|min, Magnetic         :  %12s    │" ;
                tMessage = sprint( tFormat.c_str(),
                    conditioning_string( mConditionNumbers1( 0 ) ).c_str() );
                std::cout << tMessage << std::endl ;
            }

            if ( mMumpsErrorAnalysis
                 && mKernel->dofmgr()->solver() != nullptr
                 && mKernel->dofmgr()->solver()->type() == SolverType::MUMPS )
            {
                tFormat  = "   │                    MUMPS ADD COND1, Magnetic       :  %12s    │" ;
                tMessage = sprint( tFormat.c_str(),
                    conditioning_string( mConditionNumbers1( 1 ) ).c_str() );
                std::cout << tMessage << std::endl ;

                // OMITTED, not printed as n/a, when the second row category
                // is empty: there is no such quantity for this solve, which
                // is a different statement from "we asked and got nothing"
                if ( cond2_row_wanted( mConditionNumbers1( 3 ) ) )
                {
                    tFormat  = "   │                    MUMPS ADD COND2, Magnetic       :  %12s    │" ;
                    tMessage = sprint( tFormat.c_str(),
                        conditioning_string( mConditionNumbers1( 2 ) ).c_str() );
                    std::cout << tMessage << std::endl ;
                }
            }

            if ( mComputeConditioning2 && mKernel2 != nullptr )
            {
                // see the magnetic twin: one label, always the ratio
                tFormat  = "   │                    |λ|max/|λ|min, Thermal          :  %12s    │" ;
                tMessage = sprint( tFormat.c_str(),
                    conditioning_string( mConditionNumbers2( 0 ) ).c_str() );
                std::cout << tMessage << std::endl ;
            }

            if ( mMumpsErrorAnalysis2
                 && mKernel2 != nullptr
                 && mKernel2->dofmgr()->solver() != nullptr
                 && mKernel2->dofmgr()->solver()->type() == SolverType::MUMPS )
            {
                tFormat  = "   │                    MUMPS ADD COND1, Thermal        :  %12s    │" ;
                tMessage = sprint( tFormat.c_str(),
                    conditioning_string( mConditionNumbers2( 1 ) ).c_str() );
                std::cout << tMessage << std::endl ;

                // see the magnetic twin
                if ( cond2_row_wanted( mConditionNumbers2( 3 ) ) )
                {
                    tFormat  = "   │                    MUMPS ADD COND2, Thermal        :  %12s    │" ;
                    tMessage = sprint( tFormat.c_str(),
                        conditioning_string( mConditionNumbers2( 2 ) ).c_str() );
                    std::cout << tMessage << std::endl ;
                }
            }

            tFormat  = "   │                    Time for timestep iteration     : %6u min %02u s  │" ;
            tMessage = sprint( tFormat.c_str(), minutes, seconds );
            std::cout << tMessage << std::endl ;

            // the stamp gets its own section: a divider separates the
            // measured quantities above from the wall clock below, and it
            // spans the FULL interior instead of a 48-column cell
            std::cout << "   ├───────────────────────────────────────────────────────────────────────┤" << std::endl ;
            std::cout << sprint( "   │                    %-51s│",
                                 wallclock_stamp().c_str() ) << std::endl ;

            std::cout << "   └───────────────────────────────────────────────────────────────────────┘" << std::endl ;
        }

        bool Controller::is_fullycoupled() const
        {
            return mIsFullyCoupled ;
        }

        bool
        Controller::reset() const
        {
            return mReset ;
        }

        void
        Controller::save( const string & aFilename )
        {
            if ( mCommRank == 0 )
            {
                // the temperature goes into every Exodus frame under a name
                // that says what the number is: T_max is the live maximum of
                // the nodal field in a coupled run, T_bulk the bulk value a
                // magnetic-only run holds all its materials at ( get_Tmax
                // returns the one that applies ). Created on the first save
                // — or re-adopted from a memdump, which restores it by name —
                // and refreshed on every one. MaxwellFactory refuses a
                // boundary condition label that would take either name
                Mesh * tMesh = mKernel->mesh() ;
                const string tName = mKernel2 != nullptr ? "T_max" : "T_bulk" ;

                if ( tMesh->global_variable_exists( tName ) )
                {
                    tMesh->global_variable_data( tName ) = this->get_Tmax() ;
                }
                else
                {
                    tMesh->create_global_variable( tName, this->get_Tmax() );
                }

                tMesh->save( aFilename );
                ++mMeshTimeStep ;
            }
        }

        void
        Controller::create_iv_names()
        {
            BELFEM_ASSERT( mCommRank == 0,
                "Only root may call Controller::create_iv_names()" );

            // Abstract dofs sit in the order impose_voltage_bcs relies on: the
            // current conditions first, then the voltage ones, each in
            // creation order — so the k-th terminal condition collected the
            // same way drives the k-th dof. A labelled deck section names its
            // pair ( "current : coil1 { }" -> I_coil1 / U_coil1 ); a label
            // shared by several conditions — one section with several bracket
            // groups, or two sections with one label — gets a running suffix
            // ( I_coil1_1, I_coil1_2 ). Unlabelled conditions, and generators
            // beyond the last condition ( the incidence matrix may carry more
            // cuts than the deck drives ), are positional: I_01 / U_01 ...,
            // zero-padded to the dof count
            Cell< Dof * > & tAbstractDofs = mKernel->dofmgr()->abstract_dofs() ;
            const uint tNumDofs = tAbstractDofs.size() ;

            Cell< PhysicalBoundaryCondition * > tTerminals ;
            for ( PhysicalBoundaryCondition * tBC : mKernel->boundary_conditions() )
            {
                if ( tBC->type() == BoundaryConditionType::Current
                  || tBC->type() == BoundaryConditionType::CircuitCurrent )
                {
                    tTerminals.push( tBC );
                }
            }
            for ( PhysicalBoundaryCondition * tBC : mKernel->boundary_conditions() )
            {
                if ( tBC->type() == BoundaryConditionType::Voltage
                  || tBC->type() == BoundaryConditionType::CircuitVoltage )
                {
                    tTerminals.push( tBC );
                }
            }

            Map< string, uint > tRepeats ;
            for ( PhysicalBoundaryCondition * tBC : tTerminals )
            {
                if ( tBC->label().size() > 0 )
                {
                    tRepeats[ tBC->label() ] = tRepeats.key_exists( tBC->label() ) ?
                            tRepeats( tBC->label() ) + 1 : 1 ;
                }
            }

            const string tFormatI = "I_" + format_with_leading_zeros( tNumDofs );
            const string tFormatU = "U_" + format_with_leading_zeros( tNumDofs );

            mIVNamesI.set_size( tNumDofs, "" );
            mIVNamesU.set_size( tNumDofs, "" );

            Map< string, uint > tSeen ;
            for ( uint k = 0; k < tNumDofs; ++k )
            {
                string tBase ;
                if ( k < tTerminals.size() )
                {
                    // a current row is fixed, a voltage row is free: the
                    // pairing this naming rests on
                    BELFEM_ASSERT( tAbstractDofs( k )->is_fixed() ==
                        ( tTerminals( k )->type() == BoundaryConditionType::Current
                       || tTerminals( k )->type() == BoundaryConditionType::CircuitCurrent ),
                        "Abstract dof %u does not pair with terminal condition %u ( %s )",
                        k, k, to_string( tTerminals( k )->type() ).c_str() );

                    tBase = tTerminals( k )->label() ;

                    if ( tBase.size() > 0 && tRepeats( tBase ) > 1 )
                    {
                        uint tOcc = tSeen.key_exists( tBase ) ? tSeen( tBase ) + 1 : 1 ;
                        tSeen[ tBase ] = tOcc ;
                        tBase += "_" + std::to_string( tOcc );
                    }
                }

                if ( tBase.size() > 0 )
                {
                    mIVNamesI( k ) = "I_" + tBase ;
                    mIVNamesU( k ) = "U_" + tBase ;
                }
                else
                {
                    mIVNamesI( k ) = sprint( tFormatI.c_str(), k + 1 );
                    mIVNamesU( k ) = sprint( tFormatU.c_str(), k + 1 );
                }
            }
        }

        void
        Controller::save_IV( const string & aFilename )
        {
            // Compute the full matrices
            mKernel->dofmgr()->compute_full_matrices();

            if ( mCommRank == 0 )
            {
                // Compute deltaLHS
                Vector< real > tDeltaLHS( mLHS.length() )  ;
                tDeltaLHS = ( mLHS - mLHS0 ) / mDeltaTime0 ;

                // Get the full matrices
                SpMatrix & tM = *mKernel->dofmgr()->full_mass() ;
                SpMatrix & tK = *mKernel->dofmgr()->full_stiffness() ;

                // Compute RHS
                Vector< real > tRHS( mLHS.length() ) ;
                tRHS = tM * tDeltaLHS + tK * mLHS ;

                // Get abstract DoFs
                Cell< fem::Dof * > & tAbstractDofs = mKernel->dofmgr()->abstract_dofs() ;

                // one name pair per abstract dof, shared by the csv header and
                // the mesh globals so a column and a global carry the same
                // name. The conditions do not change after setup, so once
                if ( mIVNamesI.size() != tAbstractDofs.size() )
                {
                    this->create_iv_names() ;
                }

                // an empty file name switches the csv off: the mesh globals
                // below are still refreshed, only the ascii trace is skipped.
                // The stream must be opened EXACTLY once - constructing it on a
                // name already opens ( and truncates ) the file, and a second
                // open() on an open stream returns a null buffer and latches
                // failbit, so every subsequent write is silently dropped
                std::ofstream tFile ;

                if ( aFilename.size() > 0 )
                {
                    // first call overwrites any existing file, later ones append
                    tFile.open( aFilename, mFirstIVSave ? std::ios::trunc
                                                        : std::ios::app );

                    BELFEM_ERROR( tFile.is_open(),
                        "Failed to open IV file %s", aFilename.c_str() );

                    tFile << std::setprecision( 16 );
                }

                // Write header only on first call
                if ( mFirstIVSave && tFile.is_open() )
                {
                    tFile << "time";
                    for ( uint k = 0; k < tAbstractDofs.size(); ++k )
                    {
                        tFile << " " << mIVNamesI( k ) << " " << mIVNamesU( k ) ;
                    }
                    tFile << std::endl ;
                    mFirstIVSave = false ;
                }

                // Write current time
                if ( tFile.is_open() ) tFile << mTime ;

                // Write current and voltage for each abstract DoF
                uint tCount = 0 ;

                for ( fem::Dof * tDof : tAbstractDofs )
                {
                    index_t tIndex = tDof->is_fixed() ? tDof->index() + mKernel->dofmgr()->solver_data()->my_number_of_free_dofs() : tDof->index() ;

                    real I = tDof->value() ;
                    real U = tRHS( tIndex ) ;

                    const string & Ilabel = mIVNamesI( tCount );
                    const string & Ulabel = mIVNamesU( tCount++ );

                    if ( ! mKernel->mesh()->global_variable_exists( Ilabel ) )
                    {
                        mKernel->mesh()->create_global_variable( Ilabel, I );
                    }
                    else
                    {
                        mKernel->mesh()->global_variable( Ilabel )->value() = I ;
                    }
                    if ( ! mKernel->mesh()->global_variable_exists( Ulabel ) )
                    {
                        mKernel->mesh()->create_global_variable( Ulabel, U );
                    }
                    else
                    {
                        mKernel->mesh()->global_variable( Ulabel )->value() = U ;
                    }

                    if ( tFile.is_open() ) tFile << " " << I << " " << U ;

                }

                if ( tFile.is_open() )
                {
                    tFile << std::endl ;
                    tFile.close();
                }
            }

            comm_barrier() ;
        }

        void
        Controller::set_params( const input::Section * aSection )
        {
            const input::Section * tNonLinear = aSection->section_exists( "nonlinear magnetic" )?
                                                aSection->section( "nonlinear magnetic" ) :
                                                aSection->section( "nonlinear" );

            if ( tNonLinear->key_exists( "tolerance" ) )
            {
                mRelativeEpsilonTarget = tNonLinear->get_real( "tolerance" );
            }
            else if ( tNonLinear->key_exists( "relative tolerance" ) )
            {
                mRelativeEpsilonTarget = tNonLinear->get_real( "relative tolerance" );
            }

            if ( tNonLinear->key_exists( "absolute tolerance" ) )
            {
                mAbsoluteEpsilonTarget = tNonLinear->get_real( "absolute tolerance" );

                // a negative target silently disables the escape ( eps_abs is
                // a norm and always exceeds it ) — reject the misconfiguration
                BELFEM_ERROR( mAbsoluteEpsilonTarget >= 0.0,
                    "absolute tolerance must not be negative" );
            }

            if ( tNonLinear->key_exists( "tolerance switch" ) )
            {
                mEpsilonSwitch = tNonLinear->get_real( "tolerance switch" );
            }

            if ( tNonLinear->key_exists( "algorithm" ) )
            {
                string tAlg = tNonLinear->get_string( "algorithm" );

                if ( tAlg == "Newton" || tAlg == "Newton-Raphson" )
                {
                    mAlgorithm = SolverAlgorithm::NewtonRaphson ;
                    mEquation->set_algorithm( SolverAlgorithm::NewtonRaphson ) ;
                }
                else if ( tAlg == "Picard" )
                {
                    mAlgorithm = SolverAlgorithm::Picard ;
                    mEquation->set_algorithm( SolverAlgorithm::Picard ) ;
                }
                else
                {
                    BELFEM_ERROR(false, "Unknown solver algorithm, must be Newton or Picard") ;
                }
            }
            else
            {
                mAlgorithm = SolverAlgorithm::Picard ;
                mEquation->set_algorithm( SolverAlgorithm::Picard ) ;
            }

            if ( tNonLinear->key_exists( "max iterations" ) )
            {
                // validate BEFORE the signed value wraps into the uint member
                int tMax = tNonLinear->get_int( "max iterations" );
                BELFEM_ERROR( tMax > 0,
                    "nonlinear: max iterations must be positive ( is %i )", tMax );
                mMaxNumIterations = ( uint ) tMax ;
            }

            if ( tNonLinear->key_exists( "min iterations" ) )
            {
                int tMin = tNonLinear->get_int( "min iterations" );
                BELFEM_ERROR( tMin >= 0,
                    "nonlinear: min iterations must not be negative ( is %i )", tMin );
                mMinNumIterations = ( uint ) tMin ;
            }

            // cross-check after both keys: fires also against the other
            // key's default ( e.g. min : 150 with max unset = default 100 )
            BELFEM_ERROR( mMaxNumIterations >= mMinNumIterations,
                "nonlinear: max iterations ( %u ) must not be smaller than min iterations ( %u )",
                ( unsigned int ) mMaxNumIterations,
                ( unsigned int ) mMinNumIterations );

            if ( tNonLinear->key_exists( "target iterations" ) )
            {
                // validate BEFORE the signed value wraps into the uint member
                int tTarget = tNonLinear->get_int( "target iterations" );
                BELFEM_ERROR( tTarget > 0,
                    "target iterations must be positive ( is %i )", tTarget );
                mIterationTarget = ( uint ) tTarget ;
            }

            if ( tNonLinear->key_exists( "max relaxation" ) )
            {
                mOmegaMax = tNonLinear->get_real( "max relaxation" );
            }

            if ( tNonLinear->key_exists( "min relaxation" ) )
            {
                mOmegaMin = tNonLinear->get_real( "min relaxation" );
            }

            if ( tNonLinear->key_exists( "stall window" ) )
            {
                // reject negatives BEFORE the int-to-uint conversion: a
                // negative promoted to uint is huge, so the clamp below
                // cannot catch it
                int tWindow = tNonLinear->get_int( "stall window" );
                BELFEM_ERROR( tWindow >= 0,
                    "nonlinear: stall window must not be negative ( is %i )", tWindow );

                // need at least two samples for a meaningful deviation (and to avoid a
                // zero-capacity register whose full() is trivially true)
                mStallWindow = std::max< uint >( ( uint ) tWindow, 2 );
                mResidualHistory.reserve( mStallWindow );
            }

            if ( tNonLinear->key_exists( "stall tolerance" ) )
            {
                mStallBand = tNonLinear->get_real( "stall tolerance" );
            }

            if ( tNonLinear->key_exists( "anderson depth" ) )
            {
                int tDepth = tNonLinear->get_int( "anderson depth" );
                BELFEM_ERROR( tDepth >= 0 && tDepth <= 8,
                    "anderson depth must be in 0 .. 8 ( is %i )", tDepth );
                mAndersonDepth = ( uint ) tDepth ;
                mKernel->dofmgr()->solver_data()->set_anderson_depth( mAndersonDepth );
            }

            if ( tNonLinear->key_exists( "watchdog window" ) )
            {
                int tWindow = tNonLinear->get_int( "watchdog window" );
                BELFEM_ERROR( tWindow >= 0,
                    "nonlinear: watchdog window must not be negative ( is %i, 0 disables )",
                    tWindow );
                mWatchdogWindow = ( uint ) tWindow ;
            }

            // penalty slots on the Maxwell IWG: 0 = ghost eta,
            // 1 = ghost k_reg [Ohm], 2 = coulomb gauge chi.
            // set_penalty() is collective ( rank-0 set + broadcast ) —
            // call it on ALL ranks, never inside a rank guard.
            // chi and eta are read through get_value( key, "-" ) so a
            // dimensioned value ( e.g. "eta : 4 mOhm ;" ) is rejected
            // instead of being silently SI-scaled by get_real().
            if ( tNonLinear->section_exists( "coulomb gauge penalty" ) )
            {
                const input::Section * tGauge =
                        tNonLinear->section( "coulomb gauge penalty" );

                BELFEM_ERROR( tGauge->key_exists( "chi" ),
                    "coulomb gauge penalty block requires the key 'chi' ( set it or remove the block )" );

                value tChi = tGauge->get_value( "chi", "-" );
                BELFEM_ERROR( tChi.first >= 0.0, "chi must not be negative" );

                mEquation->set_penalty( tChi.first, 2 );
            }
            else
            {
                // absent block = gauging OFF. Opt-in since 2026-09-01
                // ( Christian ): on the tape decks the default 1e-4 was
                // invisible to the conditioning estimate and kappa tracked
                // the timestep instead. It was on at 1e-4 from 2026-08-27
                // to then, and off before that. Written, not skipped: the
                // IWG constructor default is 1e-4 and would otherwise stand
                mEquation->set_penalty( 0.0, 2 );
            }

            // the ghost coupling of the thin-shell layers, read through the
            // same helper the thin-shell factory and the mesh cache tag use
            // ( fn_FEM_ghost_switch.hpp ), so the three cannot disagree:
            // -1 = no block = off ( opt-in since 2026-09-01 ), 0 = off,
            // > 0 = on with that eta. A block without eta is refused there.
            // ALWAYS written, on every rank ( set_penalty broadcasts ): the
            // IWG constructor default is 4 and would otherwise survive into
            // the assembly and the log ( audit finding ). With the ghost off
            // there are no ghost facets to assemble, so the 0 is documentary
            const real tEta = fem::read_ghost_eta( aSection );
            mEquation->set_penalty( tEta > 0.0 ? tEta : 0.0, 0 );

            if ( tEta >= 0.0 )
            {
                const input::Section * tGhost =
                        tNonLinear->section( "nitsche ghost penalty" );

                if ( tGhost->key_exists( "k_reg" ) )
                {
                    if ( tEta > 0.0 )
                    {
                        value tKreg = tGhost->get_value( "k_reg", "Ohm" );

                        // k_reg = 0 is rejected: with mRhoMin = 0 both layer
                        // stiffnesses of an SC-SC ghost facet are exactly zero
                        // and the regularized harmonic mean becomes 0/0
                        BELFEM_ERROR( tKreg.first > 0.0, "k_reg must be positive" );
                        mEquation->set_penalty( tKreg.first, 1 );
                    }
                    else if ( mCommRank == 0 )
                    {
                        // not an error: a deck flipping eta to 0 for a test
                        // should not have to delete its k_reg line
                        message( InfoLevel::Minimal,
                            "nitsche ghost penalty: k_reg is stated but eta is 0 -- the ghost is off "
                            "and k_reg is ignored" );
                    }
                }
            }

            //Read thermal if exists
            if (aSection->section_exists( "nonlinear thermal" ))
            {
                const input::Section * tNonLinearThermal = aSection->section( "nonlinear thermal" );

                if ( tNonLinearThermal->key_exists( "tolerance" ) )
                {
                    mRelativeEpsilonTarget2 = tNonLinearThermal->get_real( "tolerance" );
                }
                else if ( tNonLinearThermal->key_exists( "relative tolerance" ) )
                {
                    mRelativeEpsilonTarget2 = tNonLinearThermal->get_real( "relative tolerance" );
                }

                if ( tNonLinearThermal->key_exists( "absolute tolerance" ) )
                {
                    mAbsoluteEpsilonTarget2 = tNonLinearThermal->get_real( "absolute tolerance" );

                    // cf. the magnetic twin above
                    BELFEM_ERROR( mAbsoluteEpsilonTarget2 >= 0.0,
                        "thermal absolute tolerance must not be negative" );
                }

                if ( tNonLinearThermal->key_exists( "tolerance switch" ) )
                {
                    mEpsilonSwitch2 = tNonLinearThermal->get_real( "tolerance switch" );
                }

                if ( tNonLinearThermal->key_exists( "update gate" ) )
                {
                    mThermalUpdateGate = tNonLinearThermal->get_real( "update gate" );
                }

                if ( tNonLinearThermal->key_exists( "algorithm" ) )
                {
                    string tAlg = tNonLinearThermal->get_string( "algorithm" );

                    if ( tAlg == "Newton" || tAlg == "Newton-Raphson" )
                    {
                        mAlgorithmThermal = SolverAlgorithm::NewtonRaphson ;
                    }
                    else if ( tAlg == "Picard" )
                    {
                        mAlgorithmThermal =  SolverAlgorithm::Picard  ;
                    }
                    else
                    {
                        BELFEM_ERROR(false, "Unknown solver algorithm, must be Newton or Picard") ;
                    }
                }
                else
                {
                    mAlgorithmThermal =  SolverAlgorithm::Picard  ;
                }

                if ( tNonLinearThermal->key_exists( "max iterations" ) )
                {
                    // validate BEFORE the signed value wraps into the uint member
                    int tMax = tNonLinearThermal->get_int( "max iterations" );
                    BELFEM_ERROR( tMax > 0,
                        "nonlinear thermal: max iterations must be positive ( is %i )", tMax );
                    mMaxNumIterations2 = ( uint ) tMax ;
                }

                if ( tNonLinearThermal->key_exists( "min iterations" ) )
                {
                    int tMin = tNonLinearThermal->get_int( "min iterations" );
                    BELFEM_ERROR( tMin >= 0,
                        "nonlinear thermal: min iterations must not be negative ( is %i )", tMin );
                    mMinNumIterations2 = ( uint ) tMin ;
                }

                // cross-check after both keys: fires also against the other
                // key's default ( see the magnetic mirror above )
                BELFEM_ERROR( mMaxNumIterations2 >= mMinNumIterations2,
                    "nonlinear thermal: max iterations ( %u ) must not be smaller than min iterations ( %u )",
                    ( unsigned int ) mMaxNumIterations2,
                    ( unsigned int ) mMinNumIterations2 );

                if ( tNonLinearThermal->key_exists( "max relaxation" ) )
                {
                    mOmegaMax2 = tNonLinearThermal->get_real( "max relaxation" );
                }

                if ( tNonLinearThermal->key_exists( "min relaxation" ) )
                {
                    mOmegaMin2 = tNonLinearThermal->get_real( "min relaxation" );
                }

                if ( tNonLinearThermal->key_exists( "anderson depth" ) )
                {
                    int tDepth = tNonLinearThermal->get_int( "anderson depth" );
                    BELFEM_ERROR( tDepth >= 0 && tDepth <= 8,
                        "anderson depth must be in 0 .. 8 ( is %i )", tDepth );
                    mAndersonDepth2 = ( uint ) tDepth ;

                    // the thermal kernel may be linked after set_params;
                    // set_thermal_kernel repeats this forward in that case
                    if ( mKernel2 != nullptr )
                    {
                        mKernel2->dofmgr()->solver_data()->set_anderson_depth( mAndersonDepth2 );
                    }
                }

                if ( tNonLinearThermal->key_exists( "watchdog window" ) )
                {
                    int tWindow = tNonLinearThermal->get_int( "watchdog window" );
                    BELFEM_ERROR( tWindow >= 0,
                        "nonlinear thermal: watchdog window must not be negative ( is %i, 0 disables )",
                        tWindow );
                    mWatchdogWindow2 = ( uint ) tWindow ;
                }

                if ( tNonLinearThermal->key_exists( "coupling" ) )
                {
                    string tCoupling = tNonLinearThermal->get_string( "coupling" );

                    if ( tCoupling == "fully coupled" )
                    {
                        mIsFullyCoupled = true ;
                    }
                    else if ( tCoupling == "segregated" )
                    {
                        mIsFullyCoupled = false ;
                        if ( tNonLinearThermal->key_exists( "coupling factor" ) )
                        {
                            // divisor of the thermal timestep ( see the
                            // "initial timestep" parse ): 0 divides by zero,
                            // a negative flips the timestep sign
                            int tFactor = tNonLinearThermal->get_int( "coupling factor" );
                            BELFEM_ERROR( tFactor > 0,
                                "nonlinear thermal: coupling factor must be positive ( is %i )",
                                tFactor );
                            mCouplingFactor = tFactor ;
                        }
                    }
                    else
                    {
                        BELFEM_ERROR(false, "Unknown coupling type, must be fully coupled or segregated") ;
                    }
                }
            }

            const input::Section * tTimestep = aSection->section( "timestep" );
            {
                if ( tTimestep->key_exists( "initial timestep" ) )
                {
                    value tTime = tTimestep->get_value( "initial timestep", "s");
                    mDeltaTime = tTime.first ;
                    mDeltaTime2 = tTime.first/mCouplingFactor ;

                    // remember it — a warm restart keeps a dumped
                    // delta_time only up to this value ( see load_memdump )
                    mDeltaTimeInitial = tTime.first ;
                }
                if ( tTimestep->key_exists( "maximum timestep" ) )
                {
                    value tTime = tTimestep->get_value( "maximum timestep", "s");
                    mDeltaTimeMax = tTime.first ;
                }

                if ( tTimestep->key_exists( "minimum timestep" ) )
                {
                    value tTime = tTimestep->get_value( "minimum timestep", "s");
                    mDeltaTimeMin = tTime.first ;
                }

                // setup-tier sanity: the controller divides by these and the
                // cut path clamps against them; a NaN or inverted window must
                // fail here, not thousands of timesteps later
                BELFEM_ERROR( mDeltaTime > 0.0,
                    "initial timestep must be set and positive" );
                BELFEM_ERROR( mDeltaTimeMin > 0.0,
                    "minimum timestep must be positive" );
                BELFEM_ERROR( mDeltaTimeMax >= mDeltaTimeMin,
                    "maximum timestep ( %.3e s ) must not be smaller than the minimum ( %.3e s )",
                    ( double ) mDeltaTimeMax, ( double ) mDeltaTimeMin );
                BELFEM_ERROR( mDeltaTime >= mDeltaTimeMin && mDeltaTime <= mDeltaTimeMax,
                    "initial timestep ( %.3e s ) must lie within the minimum/maximum window",
                    ( double ) mDeltaTime );

                value tSimulationTime = tTimestep->get_value( "simulation time", "s");
                mSimulationTime = tSimulationTime.first ;

                if ( tTimestep->key_exists( "adapt timestep" ) )
                {
                    mAdaptTimestep = tTimestep->get_bool( "adapt timestep" );
                }

                // "scheme" is the preferred key, "method" the legacy alias;
                // without either the default ( bdf1 ) applies
                if ( tTimestep->key_exists( "scheme" ) || tTimestep->key_exists( "method" ) )
                {
                    string tMethod = string_to_lower( tTimestep->key_exists( "scheme" ) ?
                        tTimestep->get_string( "scheme" ) :
                        tTimestep->get_string( "method" ) );

                    if ( tMethod == "bdf1" )
                    {
                        mTimeStepping = EulerMethod::BackwardDifference1 ;
                    }
                    else if ( tMethod == "bdf2" )
                    {
                        mTimeStepping = EulerMethod::BackwardDifference2 ;
                    }
                    else if ( tMethod == "bdf3" )
                    {
                        mTimeStepping = EulerMethod::BackwardDifference3 ;
                    }
                    else if ( tMethod == "bdf4" )
                    {
                        mTimeStepping = EulerMethod::BackwardDifference4 ;
                    }
                    else if ( tMethod == "bdf5" )
                    {
                        mTimeStepping = EulerMethod::BackwardDifference5 ;
                    }
                    else if ( tMethod == "explicit" )
                    {
                        mTimeStepping = EulerMethod::ForwardExplicit ;
                    }
                    else if ( tMethod == "crc" || tMethod == "crank-nicolson" )
                    {
                        mTimeStepping = EulerMethod::CrankNicolson ;
                    }
                    else if ( tMethod == "galerkin" )
                    {
                        mTimeStepping = EulerMethod::Galerkin ;
                    }
                    else
                    {
                        BELFEM_ERROR( false, "unknown time stepping scheme: %s", tMethod.c_str() );
                    }
                }

                // master switch for the Anderson mixing ( opt-in: off by
                // default ). On enables the default depths ( magnetic 3,
                // thermal 1 ) where no explicit per-field depth was given;
                // an explicit depth key always wins, including an explicit
                // 0. Off with an explicit nonzero depth is a contradiction
                if ( tTimestep->key_exists( "anderson stabilization" ) )
                {
                    if ( tTimestep->get_bool( "anderson stabilization" ) )
                    {
                        if ( ! tNonLinear->key_exists( "anderson depth" ) )
                        {
                            mAndersonDepth = 3 ;
                        }
                        if ( ! ( aSection->section_exists( "nonlinear thermal" )
                                 && aSection->section( "nonlinear thermal" )->key_exists( "anderson depth" ) ) )
                        {
                            mAndersonDepth2 = 1 ;
                        }
                    }
                    else
                    {
                        BELFEM_ERROR( ! ( tNonLinear->key_exists( "anderson depth" )
                                          && mAndersonDepth > 0 ),
                            "anderson stabilization : false contradicts anderson depth : %u",
                            ( unsigned int ) mAndersonDepth );

                        if ( aSection->section_exists( "nonlinear thermal" ) )
                        {
                            BELFEM_ERROR( ! ( aSection->section( "nonlinear thermal" )->key_exists( "anderson depth" )
                                              && mAndersonDepth2 > 0 ),
                                "anderson stabilization : false contradicts thermal anderson depth : %u",
                                ( unsigned int ) mAndersonDepth2 );
                        }

                        mAndersonDepth  = 0 ;
                        mAndersonDepth2 = 0 ;
                    }
                }
                if ( tTimestep->key_exists( "save every" ) )
                {
                    string tUnitDef = tTimestep->get_units( "save every" );
                    value tValue = unit_to_si( tUnitDef );
                    BELFEM_ERROR(
                            check_unit( tValue, "s" ),
                            "required unit for the offset of the boundary condition is: s" );
                    mSaveEvery = tTimestep->get_value( "save every", "s").first;
                }
                else
                {
                    mSave = true ;
                }
            }

            // ---------------------------------------------------------------
            // solver banner. Printed here, and not next to the parse of each
            // key, because it spans two sections: the penalties come from
            // "nonlinear", the scheme and the Anderson depths from "timestep",
            // which is read above.
            //
            // Every row is printed, including the off ones. Both penalties are
            // opt-in, and both have changed their default inside a single
            // release ( chi was on at 1e-4 from 2026-08-27 to 2026-09-01 ), so
            // the log has to state what a run actually used -- an absent row
            // cannot be told from a knob nobody looked at.
            //
            // Row labels follow the deck sections that set them
            // ( "nitsche ghost penalty", "coulomb gauge penalty",
            //   "anderson stabilization" ) so a reader can grep the input.
            // The scopes are the assembly reality: the ghost is a facet term
            // between thin-shell layer blocks ( maxwell::h_ghost ), while chi
            // rides the volume routines h_picard / h_newton_mu0 / h_newton_mu,
            // which serve DomainType::Conductor and DomainType::ThinShell
            // alike, plus the side connectors
            if ( mCommRank == 0 )
            {
                const real tEtaP  = mEquation->penalty( 0 );
                const real tKregP = mEquation->penalty( 1 );
                const real tChiP  = mEquation->penalty( 2 );

                string tGhostRow = tEtaP > 0.0 ?
                        sprint( "eta %g, k_reg %g Ohm",
                                ( double ) tEtaP, ( double ) tKregP ) :
                        string( "off" );

                string tGaugeRow = tChiP > 0.0 ?
                        sprint( "chi %g", ( double ) tChiP ) :
                        string( "off" );

                message( InfoLevel::Default,
                    "\n    stabilization penalties :\n"
                    "            * Nitsche ghost ( h-shells )  : %s\n"
                    "            * Coulomb gauge ( h-domains ) : %s",
                    tGhostRow.c_str(),
                    tGaugeRow.c_str() );

                // Anderson is one deck switch driving two depths. The thermal
                // one is only meaningful when the deck has a thermal field,
                // so it is named only then
                string tAndersonRow ;

                if ( mAndersonDepth == 0 && mAndersonDepth2 == 0 )
                {
                    tAndersonRow = "off" ;
                }
                else if ( aSection->section_exists( "nonlinear thermal" ) )
                {
                    tAndersonRow = sprint( "on ( magnetic %u, thermal %u )",
                            ( unsigned int ) mAndersonDepth,
                            ( unsigned int ) mAndersonDepth2 );
                }
                else
                {
                    tAndersonRow = sprint( "on ( depth %u )",
                            ( unsigned int ) mAndersonDepth );
                }

                // ForwardExplicit is the only explicit member of EulerMethod;
                // the BDF family, Crank-Nicolson and Galerkin are all implicit
                message( InfoLevel::Default,
                    "\n    timestepping :\n"
                    "            * %s scheme        : %s\n"
                    "            * Anderson stabilization : %s",
                    mTimeStepping == EulerMethod::ForwardExplicit ? "explicit" : "implicit",
                    to_string( mTimeStepping ).c_str(),
                    tAndersonRow.c_str() );
            }

            // warm-restart policy: opt-out. By default an existing memdump
            // is resumed ( announced by the banner in load_memdump );
            // restart : false ignores it and starts fresh
            if ( tTimestep->key_exists( "restart" ) )
            {
                mAllowRestart = tTimestep->get_bool( "restart" );
            }

            // backstop on the floor escalation; 0 = retry indefinitely
            if ( tTimestep->key_exists( "floor retries" ) )
            {
                int tRetries = tTimestep->get_int( "floor retries" );
                BELFEM_ERROR( tRetries >= 0,
                    "floor retries must not be negative ( is %i )", tRetries );
                mMaxFloorRetries = ( uint ) tRetries ;
            }

            // condition-number diagnostics ( print_footer ). Only the flags
            // are read here. Two INDEPENDENT keys since 2026-08-30:
            // "compute conditioning" runs the eigenvalue estimate
            // ( compute_conditioning() ), "mumps error analysis" arms MUMPS
            // ICNTL(11) per timestep around the first solve
            // ( arm_conditioning_* ). Neither is a fallback for the other.
            //
            // The key belongs to a SOLVE, so it is read per field from
            // "linear magnetic" / "linear thermal" -- the same sections that
            // already choose each field's library. A "compute conditioning"
            // in the solver section itself is the shorthand for both, and the
            // per-field key overrides it, mirroring how "linear" is the
            // fallback for the two specific sections
            if ( aSection->key_exists( "compute conditioning" ) )
            {
                mComputeConditioning  = aSection->get_bool( "compute conditioning" );
                mComputeConditioning2 = mComputeConditioning ;
            }

            if ( aSection->section_exists( "linear magnetic" )
                 && aSection->section( "linear magnetic" )->key_exists( "compute conditioning" ) )
            {
                mComputeConditioning = aSection->section( "linear magnetic" )
                        ->get_bool( "compute conditioning" );
            }

            if ( aSection->section_exists( "linear thermal" )
                 && aSection->section( "linear thermal" )->key_exists( "compute conditioning" ) )
            {
                mComputeConditioning2 = aSection->section( "linear thermal" )
                        ->get_bool( "compute conditioning" );
            }

            // "mumps error analysis" -- the SECOND, independent diagnostic.
            // Same three-site chain as the flag above: the solver section is
            // the shorthand for both fields, each linear block overrides its
            // own. Split out of "compute conditioning" on 2026-08-30 because
            // the two answer different questions at very different prices, and
            // neither could be requested alone.
            if ( aSection->key_exists( "mumps error analysis" ) )
            {
                mMumpsErrorAnalysis  = aSection->get_bool( "mumps error analysis" );
                mMumpsErrorAnalysis2 = mMumpsErrorAnalysis ;
            }

            if ( aSection->section_exists( "linear magnetic" )
                 && aSection->section( "linear magnetic" )->key_exists( "mumps error analysis" ) )
            {
                mMumpsErrorAnalysis = aSection->section( "linear magnetic" )
                        ->get_bool( "mumps error analysis" );
            }

            if ( aSection->section_exists( "linear thermal" )
                 && aSection->section( "linear thermal" )->key_exists( "mumps error analysis" ) )
            {
                mMumpsErrorAnalysis2 = aSection->section( "linear thermal" )
                        ->get_bool( "mumps error analysis" );
            }

            // deliberately NOT cleared when mKernel2 is still null here:
            // set_thermal_kernel() may run after set_params ( see the forward
            // below ), so a null kernel at parse time says nothing about
            // whether a thermal field exists. Every consumer of the flag
            // checks mKernel2 itself
            //
            // THE MAGNETIC FIELD IS CHECKED HERE, THE THERMAL ONE IS NOT.
            // Its kernel does not exist yet on the coupled path
            // ( MaxwellFactory::create_controller builds the controller with
            // the magnetic kernel alone ), so a thermal branch in this block
            // could never fire -- which is exactly what the pre-2026-08-30
            // version did, silently, for its whole life. The thermal half
            // lives in check_thermal_diagnostics(), called once the kernel is
            // actually attached
            mParamsSet = true ;
            this->check_magnetic_diagnostics();

            // the two-argument constructor attaches a thermal kernel BEFORE
            // set_params runs, so on that path this is the second of the two
            // operations and the thermal checks belong here. On the ordinary
            // path mKernel2 is null and set_thermal_kernel does it later.
            // check_thermal_diagnostics() is one-shot, so a caller that
            // reaches both sites still gets one message
            if ( mKernel2 != nullptr )
            {
                this->check_thermal_diagnostics();
            }

            // deck-value snapshot of the iteration budgets: the restore
            // targets of the floor escalation ( reset_timestep / finalize )
            mMaxNumIterationsDeck  = mMaxNumIterations ;
            mMaxNumIterations2Deck = mMaxNumIterations2 ;
            mWatchdogWindowDeck    = mWatchdogWindow ;
            mWatchdogWindow2Deck   = mWatchdogWindow2 ;

            // coupled-mode clamp on the thermal watchdog window: with
            // lockstep counting, the LAST iterate that can call
            // watchdog_thermal is mMaxNumIterations ( the doomed check
            // skips the thermal block on the ceiling iterate ), so a window
            // >= maxIter - min2 leaves the watchdog structurally unable to
            // fire before the magnetic ceiling when the thermal best lands
            // early ( ~min2 + 1 ). Conservative policy calibrated to that
            // early-best case — a late best or the spare clauses can still
            // outlast the ceiling. Clamping the DECK snapshot keeps the
            // floor escalation ( which scales window and budget by the same
            // factor ) and the finalize restore coherent. 0 stays the
            // documented disable; segregated decks keep their window. The
            // ternary is load-bearing: a max<uint> form underflows when
            // maxIterDeck <= min2 + 1
            if ( mIsFullyCoupled && mWatchdogWindow2Deck > 0 )
            {
                const uint tBound =
                    mMaxNumIterationsDeck > mMinNumIterations2 + 1
                    ? mMaxNumIterationsDeck - mMinNumIterations2 - 1
                    : 1 ;
                if ( mWatchdogWindow2Deck > tBound )
                {
                    mWatchdogWindow2Deck = tBound ;
                    mWatchdogWindow2     = tBound ;

                    // a magnetic-only deck also lands here ( fully coupled
                    // is the default and the member is dead without a
                    // thermal kernel ) — only tell the user when the deck
                    // actually has a thermal section
                    if ( mCommRank == 0
                         && aSection->section_exists( "nonlinear thermal" ) )
                    {
                        message( InfoLevel::Default,
                            "\n    thermal watchdog window clamped to %u so it can fire\n"
                            "    before the magnetic iteration ceiling ( %u )",
                            ( unsigned int ) tBound,
                            ( unsigned int ) mMaxNumIterationsDeck );
                    }
                }
            }

            // the depths must reach SolverData even when no anderson key was
            // given ( the opt-in master switch may have filled them );
            // idempotent, and set_thermal_kernel repeats the thermal forward
            // when that kernel links later
            mKernel->dofmgr()->solver_data()->set_anderson_depth( mAndersonDepth );
            if ( mKernel2 != nullptr )
            {
                mKernel2->dofmgr()->solver_data()->set_anderson_depth( mAndersonDepth2 );
            }

            // headroom gate ( see the header ): a nonlinear tolerance the
            // linear solver cannot serve is a config error, caught at setup
            // rather than discovered hours in. Repeated in
            // set_thermal_kernel for a late-linked thermal kernel
            this->check_iterative_solver_headroom(
                mKernel->dofmgr(), mRelativeEpsilonTarget, "magnetic" );
            this->check_compression_headroom(
                mKernel->dofmgr(), mRelativeEpsilonTarget, "magnetic" );
            if ( mKernel2 != nullptr )
            {
                this->check_iterative_solver_headroom(
                    mKernel2->dofmgr(), mRelativeEpsilonTarget2, "thermal" );
                this->check_compression_headroom(
                    mKernel2->dofmgr(), mRelativeEpsilonTarget2, "thermal" );
            }
        }

        void
        Controller::set_thermal_kernel( Kernel * aKernel )
        {
            mKernel2 = aKernel ;
            mMesh2 = aKernel->mesh() ;
            mEquation2 = reinterpret_cast<IWG_Timestep *>(mKernel2->dofmgr()->iwg()) ;
            mEquation2->set_algorithm( mAlgorithmThermal ) ;
            mKernel2->set_controller( this );

            // forward the mixing depth if set_params ran before the thermal
            // kernel was linked ( set_anderson_depth is idempotent )
            if ( mAndersonDepth2 > 0 )
            {
                mKernel2->dofmgr()->solver_data()->set_anderson_depth( mAndersonDepth2 );
            }

            // late-attached thermal kernel: arm the soft-fail contract
            // ( cf. the constructor )
            if ( mKernel2->dofmgr()->solver() != nullptr )
            {
                mKernel2->dofmgr()->solver()->wrapper()->set_soft_fail( true );
            }

            // late-attached thermal kernel: repeat the headroom gate and
            // the compression check that set_params could not run against
            // a null kernel
            this->setup_thermal_eigen();
            this->check_thermal_diagnostics();
            this->check_iterative_solver_headroom(
                mKernel2->dofmgr(), mRelativeEpsilonTarget2, "thermal" );
            this->check_compression_headroom(
                mKernel2->dofmgr(), mRelativeEpsilonTarget2, "thermal" );

            // the thermal kernel never receives a material assignment of its
            // own ( auto_set_materials only works across dof managers of one
            // kernel ): share the maxwell-side material so that phi-type
            // domains ( Buffer, Ferro ) can evaluate density, cp and lambda
            // in T_phi
            for ( Block * tBlock : mKernel2->dofmgr()->blocks() )
            {
                if ( tBlock->material() == nullptr )
                {
                    tBlock->set_material( mKernel->dofmgr()->block( tBlock->id() )->material() );
                }
            }

            // relink the maxwell data helpers on both kernels so that they
            // see the thermal kernel (rebuilds any already-allocated helper)
            for ( Block * tBlock : mKernel->dofmgr()->blocks() )
            {
                if ( tBlock->calculator() != nullptr )
                {
                    tBlock->calculator()->link_maxwell( mKernel, mKernel2 );
                }
            }
            for ( Block * tBlock : mKernel2->dofmgr()->blocks() )
            {
                if ( tBlock->calculator() != nullptr )
                {
                    tBlock->calculator()->link_maxwell( mKernel, mKernel2 );
                }
            }
        }

        real Controller::simulation_time() const
        {
            return mSimulationTime ;
        }

        bool
        Controller::save() const
        {
            // without a save grid, every converged step is stored
            if ( mSaveEvery == 0.0 )
            {
                return mSave ;
            }

            // verify rather than predict: mSave is armed by adjust_timestep
            // when the NEXT step is trimmed to land on the save grid, but the
            // flag is sticky — a timestep cut skips the adjust phase and the
            // retry converges at a halved, off-grid time with the stale flag
            // still set ( observed: a stored frame at 13.08087 ms next to the
            // legitimate 13.0 ms save ). A save happens only when the time
            // actually sits on the grid. Tolerance: FP accumulation plus the
            // minimum-timestep clamp, which can land the trimmed step just
            // past the grid point when the remaining distance is below
            // mDeltaTimeMin
            real tGrid = mTime / mSaveEvery ;
            return std::abs( tGrid - std::round( tGrid ) ) * mSaveEvery
                 < mDeltaTimeMin + 1e-6 * mSaveEvery ;
        }

        void
        Controller::set_circuit( Circuit * aCircuit )
        {
            BELFEM_ASSERT( mCircuit == nullptr, "Circuit already set" ) ;
            BELFEM_ASSERT( mKernel != nullptr, "Kernel not set");

            mCircuit = aCircuit ;
            mCircuitCurrentBCs.clear();
            mCircuitVoltageBCs.clear();

            for ( PhysicalBoundaryCondition * tBC : mKernel->boundary_conditions() )
            {
                if ( tBC->type() == BoundaryConditionType::CircuitCurrent )
                {
                    mCircuitCurrentBCs.push( tBC );
                }
                else if ( tBC->type() == BoundaryConditionType::CircuitVoltage )
                {
                    mCircuitVoltageBCs.push( tBC );
                }
            }
        }

        void
        Controller::synchronize_history_fields(
                IWG_Timestep * aEquation,
                DofManager   * aDofMgr,
                const uint     aStepCount,
                const string & aPath,
                const char   * aEquationName )
        {
            // The resumed step SHIFTS before it assembles, so the fields
            // must hold the PRE-shift levels 0..r-1 with
            // r = min( count, order-1 ). Mid-ramp that is one deeper than
            // the scalar restore's own r ( min( count-1, order-1, 4 ) );
            // once the ramp saturates the two agree. Level order-1 is
            // never read pre-shift — the first shift fills it from
            // order-2 — so it is deliberately outside this set.
            // timestepping_order() via the IWG base pointer: the override
            // is private on IWG_Timestep itself
            uint tOrder = aDofMgr->iwg()->timestepping_order() ;
            uint tRequired = std::min( aStepCount,
                                       tOrder > 0 ? tOrder - 1 : 0 );
            if ( tRequired == 0 ) return ;

            // EXACTLY the levels verified below: handing a deeper,
            // never-filled level to distribute_fields would fault rather
            // than help — FieldData::distribute walks entity indices with
            // no empty guard, and EDGE/FACE history vectors stay empty
            // until something sizes them, so a PARALLEL MID-RAMP restart
            // would reproduce the very crash this function prevents
            // ( Grok — invisible to the serial probe )
            Cell< string > tParents ;
            Cell< string > tHistory ;
            aEquation->history_field_labels( tParents, tHistory, tRequired );
            if ( tHistory.size() == 0 ) return ;

            // rank 0 holds the loaded levels; workers hold the empty
            // shells create_old_dof_fields made during initialize().
            // Same machinery as the per-step field synch — entity
            // multiplicities ( edge/face strides ) included. NOTE: the
            // SAVE side needs no collect twin — DofManager::solve ends
            // every step with distribute( all_fields ) from a master that
            // holds the full solution, so rank 0's fields are
            // solution-fresh and its shift builds the true history
            // ( stale-save hypothesis refuted by both reviewers; a future
            // fully-distributed solver would break this and must revisit )
            aDofMgr->distribute_fields( tHistory );

            for ( const string & tLabel : tParents )
            {
                const index_t tParentLength =
                    aDofMgr->mesh()->field( tLabel )->data().length() ;

                for ( uint s = 0; s < tRequired; ++s )
                {
                    const string tOldLabel = tLabel + std::to_string( s );
                    const index_t tLength =
                        aDofMgr->mesh()->field( tOldLabel )->data().length() ;

                    // a named error here is the difference between this
                    // line and the raw SEGV / bounds throw the first
                    // execution of this path produced: whatever
                    // corrupts a level, THIS names the field and rank
                    BELFEM_ERROR( tLength == tParentLength,
                        "BDF history restore incomplete ( %s ): field '%s' has length %lu on rank %u, expected %lu\n"
                        "( memdump %s promises %u history steps ) - the dump's history is unusable",
                        aEquationName,
                        tOldLabel.c_str(),
                        ( long unsigned int ) tLength,
                        ( unsigned int ) mCommRank,
                        ( long unsigned int ) tParentLength,
                        aPath.c_str(),
                        ( unsigned int ) aStepCount );
                }
            }
        }

        void
        Controller::save_memdump( const string & aPath )
        {
            if ( mCommRank == 0 )
            {

                HDF5 tFile( aPath, FileMode::NEW );

                tFile.create_group( "meta" );
                mMesh->save_meta( tFile.active_group(), mRunningTimeStep );

                // Persist the controller's next-step timestep. The warm-restart
                // work changed the READ side: load_memdump caps a restored value at the
                // deck's initial timestep before the first post-restart solve.
                // The dumped value remains useful for diagnostics and for a
                // future loader that can safely re-enter at the earned step.
                // The BDF history below is still restored verbatim; that is
                // the order-cliff fix.
                // finalize() has already adjusted Delta t for the next step.
                // The helper
                // needs an lvalue: active_group() returns hid_t by value
                hid_t  tGroup  = tFile.active_group() ;
                herr_t tStatus = 0 ;
                hdf5::save_scalar_to_file( tGroup, "delta_time", mDeltaTime, tStatus );

                // BDF integrator state, so a warm start resumes at full
                // order instead of re-anchoring the ramp at BDF1 ( the
                // restart cliff, 2026-08-15 ). bdf_last_dt is the IWG's OWN
                // delta -- the LAST COMPLETED step h_n; the delta_time above
                // is the NEXT step and is a different quantity
                {
                    Vector< real > tH ;
                    uint tStepCount ;
                    real tLastDt ;

                    mEquation->save_history_state( tH, tStepCount, tLastDt );
                    hdf5::save_vector_to_file( tGroup, "bdf_h", tH, tStatus );
                    hdf5::save_scalar_to_file( tGroup, "bdf_step_count",
                        static_cast< index_t >( tStepCount ), tStatus );
                    hdf5::save_scalar_to_file( tGroup, "bdf_last_dt", tLastDt, tStatus );

                    if ( mEquation2 != nullptr )
                    {
                        mEquation2->save_history_state( tH, tStepCount, tLastDt );
                        hdf5::save_vector_to_file( tGroup, "bdf_h2", tH, tStatus );
                        hdf5::save_scalar_to_file( tGroup, "bdf_step_count2",
                            static_cast< index_t >( tStepCount ), tStatus );
                        hdf5::save_scalar_to_file( tGroup, "bdf_last_dt2", tLastDt, tStatus );
                    }
                }

                tFile.close_active_group();

                tFile.create_group( "fields" );
                mMesh->save_fields( tFile.active_group() );
                tFile.close_active_group();

                if ( mMesh->number_of_global_variables() > 0 )
                {
                    tFile.create_group( "globals" );
                    mMesh->save_globals( tFile.active_group() );
                    tFile.close_active_group();
                }

                if ( mMesh2 != nullptr )
                {
                    tFile.create_group( "fields2" );

                    mMesh2->save_fields( tFile.active_group() );
                    tFile.close_active_group();

                    if ( mMesh2->number_of_global_variables() > 0 )
                    {
                        tFile.create_group( "globals2" );
                        mMesh2->save_globals( tFile.active_group() );
                        tFile.close_active_group();
                    }
                }

                if ( mCircuit != nullptr )
                {
                    tFile.create_group( "circuit" );
                    mCircuit->save_state( tFile.active_group() );
                    tFile.close_active_group();
                }
                tFile.close();
            }
            comm_barrier();
        }

        void
        Controller::load_memdump( const string & aPath )
        {
            // BDF integrator state from the dump ( flags + payload live
            // outside the rank-0 block: every rank takes the broadcasts )
            uint tHaveBdf = 0, tHaveBdf2 = 0 ;
            Vector< real > tBdfH, tBdfH2 ;
            uint tBdfStepCount = 0, tBdfStepCount2 = 0 ;
            real tBdfLastDt = 0.0, tBdfLastDt2 = 0.0 ;

            bool tExists = false ;

            if ( mCommRank == 0 )
            {
                tExists = file_exists( aPath );

                // warm restart is opt-OUT: by default an existing dump is
                // resumed ( announced by the banner below ); a deck sets
                // timestep { restart : false ; } to ignore it and start
                // fresh — the dump is then overwritten at the first save
                if ( tExists && ! mAllowRestart )
                {
                    if ( gLog.info_level() > 0 )
                    {
                        std::cout << "   found " << aPath
                            << " — ignoring it ( restart : false ); starting fresh"
                            << std::endl ;
                    }
                    tExists = false ;
                }
                broadcast( tExists );

                if ( tExists )
                {
                    HDF5 tFile( aPath, FileMode::OPEN_RDONLY );
                    tFile.select_group( "meta" );
                    mRunningTimeStep = mMesh->load_meta( tFile.active_group() );

                    // restore the dumped timestep if present ( older dumps
                    // lack the key — they fall back to the deck's initial
                    // value ). Read while the meta group is still selected;
                    // lvalue required, cf. save_memdump
                    hid_t tGroup = tFile.active_group() ;
                    if ( hdf5::dataset_exists( tGroup, "delta_time" ) )
                    {
                        herr_t tStatus = 0 ;
                        hdf5::load_scalar_from_file( tGroup, "delta_time", mDeltaTime, tStatus );

                        // a corrupt dump must fail here, not propagate a NaN
                        // through the clamp into the time loop
                        BELFEM_ERROR( std::isfinite( mDeltaTime ) && mDeltaTime > 0.0,
                            "corrupt delta_time in memdump %s", aPath.c_str() );
                    }

                    // BDF integrator state ( read here while /meta is open;
                    // broadcast + restore happen after the field synch below )
                    // presence = ALL THREE keys; a partial triple is a
                    // corrupt dump and must fail with a NAMED error here,
                    // not through the generic HDF5 helper downstream
                    {
                        const bool tH  = hdf5::dataset_exists( tGroup, "bdf_h" );
                        const bool tC  = hdf5::dataset_exists( tGroup, "bdf_step_count" );
                        const bool tD  = hdf5::dataset_exists( tGroup, "bdf_last_dt" );

                        BELFEM_ERROR( tH == tC && tC == tD,
                            "incomplete BDF history ( magnetic ) in memdump %s",
                            aPath.c_str() );

                        if ( tH )
                        {
                            herr_t tStatus = 0 ;
                            index_t tCount = 0 ;
                            hdf5::load_vector_from_file( tGroup, "bdf_h", tBdfH, tStatus );
                            hdf5::load_scalar_from_file( tGroup, "bdf_step_count", tCount, tStatus );
                            hdf5::load_scalar_from_file( tGroup, "bdf_last_dt", tBdfLastDt, tStatus );
                            tBdfStepCount = static_cast< uint >( tCount );
                            tHaveBdf = 1 ;
                        }
                    }
                    {
                        const bool tH  = hdf5::dataset_exists( tGroup, "bdf_h2" );
                        const bool tC  = hdf5::dataset_exists( tGroup, "bdf_step_count2" );
                        const bool tD  = hdf5::dataset_exists( tGroup, "bdf_last_dt2" );

                        BELFEM_ERROR( tH == tC && tC == tD,
                            "incomplete BDF history ( thermal ) in memdump %s",
                            aPath.c_str() );

                        if ( tH )
                        {
                            herr_t tStatus = 0 ;
                            index_t tCount = 0 ;
                            hdf5::load_vector_from_file( tGroup, "bdf_h2", tBdfH2, tStatus );
                            hdf5::load_scalar_from_file( tGroup, "bdf_step_count2", tCount, tStatus );
                            hdf5::load_scalar_from_file( tGroup, "bdf_last_dt2", tBdfLastDt2, tStatus );
                            tBdfStepCount2 = static_cast< uint >( tCount );
                            tHaveBdf2 = 1 ;
                        }
                    }

                    tFile.close_active_group();
                    tFile.select_group( "fields" );
                    mMesh->load_fields( tFile.active_group() );
                    tFile.close_active_group();

                    if ( hdf5::group_exists( tFile.active_group(), "globals"))
                    {
                        tFile.select_group( "globals" );
                        mMesh->load_globals( tFile.active_group() );
                        tFile.close_active_group();
                    }

                    if ( mMesh2 != nullptr && hdf5::group_exists( tFile.active_group(), "fields2"))
                    {
                        tFile.select_group( "fields2" );
                        mMesh2->load_fields( tFile.active_group() );
                        tFile.close_active_group();
                        if ( hdf5::group_exists( tFile.active_group(), "globals2"))
                        {
                            tFile.select_group( "globals2" );
                            mMesh2->load_globals( tFile.active_group() );
                            tFile.close_active_group();
                        }
                    }

                    if ( mCircuit != nullptr && hdf5::group_exists( tFile.active_group(), "circuit"))
                    {
                        tFile.select_group( "circuit" );
                        mCircuit->load_state( tFile.active_group() );
                        tFile.close_active_group();
                    }
                    tFile.close();
                }
            }
            else
            {
                broadcast( tExists );
            }
            comm_barrier();

            if ( tExists )
            {
                broadcast( mKernel->mesh()->time_step() );
                broadcast( mKernel->mesh()->time_stamp() );

                mMeshTimeStep = mKernel->mesh()->time_step();
                mTime = mKernel->mesh()->time_stamp();

                // slave the thermal clock to the restored magnetic time: a dump
                // is only ever written after the segregated catch-up loop has
                // equalized the two clocks — without this, a
                // segregated warm restart re-integrates thermal from t = 0 and
                // rotates the restored BDF history into garbage
                mTime2 = mTime ;

                broadcast( mRunningTimeStep );
                broadcast( tHaveBdf );
                broadcast( tHaveBdf2 );

                BELFEM_ERROR( tHaveBdf > 0,
                    "memdump %s carries no magnetic BDF history "
                    "( bdf_h / bdf_step_count / bdf_last_dt ). Either the dump "
                    "predates 2026-08-15 or it is truncated; either way it cannot "
                    "be resumed. Start fresh with timestep { restart : false ; }.",
                    aPath.c_str() );

                if ( mEquation2 != nullptr )
                {
                    BELFEM_ERROR( tHaveBdf2 > 0,
                        "memdump %s carries no thermal BDF history "
                        "( bdf_h2 / bdf_step_count2 / bdf_last_dt2 ), but this is a "
                        "coupled run. Either the dump is magnetic-only or it predates "
                        "2026-08-15. Start fresh with timestep { restart : false ; }.",
                        aPath.c_str() );
                }

                // warm-start timestep: rank 0 holds the dumped value ( or the
                // deck initial for older dumps ) — every rank takes this
                // collective, no HDF5 on non-root
                broadcast( mDeltaTime );
                mDeltaTime = std::clamp( mDeltaTime, mDeltaTimeMin, mDeltaTimeMax );

                // cap the RE-ENTRY step at the deck's initial
                // timestep. This is a restart-entry cap, not the normal growth
                // limiter: it prevents the first post-restart solve from
                // initialising at a large dt. Not announced on its own: the
                // WARM RESTART banner below prints the delta t that results
                if ( mDeltaTime > mDeltaTimeInitial )
                {
                    mDeltaTime = mDeltaTimeInitial ;
                }

                // initialize the dof manager if it has not been already.
                // Seed mode: the restored fields are the truth — a full
                // initialize would write the fixed dofs' current values
                // ( factory constants, pre-load BELFEM_EPS currents ) over
                // them. Fixed dof values stay stale until the first
                // compute_boundary_conditions, which runs before assembly
                mKernel->dofmgr()->initialize( true );

                // synchronize the fields ( globals stay rank-0-only by
                // design — the dump restored them on the master mesh and
                // no worker code reads them )
                mKernel->dofmgr()->distribute_fields( mKernel->dofmgr()->iwg()->all_fields() );

                // seed the free dofs from the restored fields: a manager
                // that was initialized before the load ( thermal, in the
                // factory ) still holds the values it was seeded with back
                // then, and the restored state would be invisible to the
                // Picard residual. For the manager initialized
                // just above this is a no-op repeat of its own seeding
                mKernel->dofmgr()->seed_dof_values() ;

                if ( mKernel2 != nullptr )
                {
                    mKernel2->dofmgr()->initialize( true );
                    mKernel2->dofmgr()->distribute_fields( mKernel2->dofmgr()->iwg()->all_fields() );
                    mKernel2->dofmgr()->seed_dof_values() ;
                }

                // a restored state has no valid mixing history ( defensive:
                // the history is never persisted )
                this->anderson_clear_magnetic() ;
                this->anderson_clear_thermal() ;

                // BDF warm start: broadcast the integrator payload and
                // restore it on EVERY rank ( each rank computes its own
                // coefficients ). Present-but-invalid data is a hard error
                // ( delta_time precedent ). The flags were broadcast and
                // REQUIRED above, so the guard here is now structural
                if ( tHaveBdf > 0 )
                {
                    broadcast( tBdfH );
                    broadcast( tBdfStepCount );
                    broadcast( tBdfLastDt );

                    // the integrator state promises history the FIELD side
                    // must actually hold: distribute the numbered levels
                    // ( they are not in all_fields, so the synch above
                    // never moved them — the warm-restart crash ), then verify on
                    // EVERY rank before the integrator is told it has
                    // history. initialize() above has already run
                    // init_qold_table on all ranks, so the fields exist
                    this->synchronize_history_fields( mEquation,
                        mKernel->dofmgr(), tBdfStepCount, aPath,
                        "magnetic" );

                    BELFEM_ERROR( mEquation->restore_history_state(
                            tBdfH, tBdfStepCount, tBdfLastDt ),
                        "invalid BDF history in memdump %s ( magnetic )",
                        aPath.c_str() );
                }
                // flag already broadcast and required above -- re-broadcasting
                // it here would be a second collective against the ranks' one
                if ( tHaveBdf2 > 0 && mEquation2 != nullptr )
                {
                    broadcast( tBdfH2 );
                    broadcast( tBdfStepCount2 );
                    broadcast( tBdfLastDt2 );

                    this->synchronize_history_fields( mEquation2,
                        mKernel2->dofmgr(), tBdfStepCount2, aPath,
                        "thermal" );

                    BELFEM_ERROR( mEquation2->restore_history_state(
                            tBdfH2, tBdfStepCount2, tBdfLastDt2 ),
                        "invalid BDF history in memdump %s ( thermal )",
                        aPath.c_str() );
                }

                // raw cout like the WARM RESTART banner below, and inside
                // the same guard: message( Default ) is silent at Minimal
                // info level, and a cold-start MUST not be silent
                if ( mCommRank == 0 && gLog.info_level() > 0 )
                {
                    if ( tHaveBdf > 0 )
                    {
                        std::cout << sprint(
                            "    warm start resumes the BDF history ( magnetic%s )",
                            ( tHaveBdf2 > 0 && mEquation2 != nullptr ) ? " + thermal" : "" )
                            << std::endl ;

                        // the cold-start arms that used to live here are gone:
                        // both triples are REQUIRED above, so a dump
                        // without them never reaches this banner
                    }
                }

                // a warm restart must be VISIBLE: say where the run resumes
                if ( mCommRank == 0 && gLog.info_level() > 0 )
                {
                    std::cout << "   ┌───────────────────────────────────────────────────────────────────────┐" << std::endl ;
                    string tLine = sprint(
                        " WARM RESTART from %s", aPath.c_str() );
                    std::cout << sprint( "   │%-71s│", tLine.c_str() ) << std::endl ;
                    // print the step the first assembly will carry ( the loop
                    // increments mRunningTimeStep before assembling ), so the
                    // banner matches the step headers and is safe to copy into
                    // BELFEM_DUMP_SYSTEM_STEP
                    tLine = sprint(
                        " resuming at t = %.4f ms, next step %u, delta t = %.4f ms",
                        mTime * 1000.0,
                        ( unsigned int ) ( mRunningTimeStep + 1 ),
                        mDeltaTime * 1000.0 );
                    std::cout << sprint( "   │%-71s│", tLine.c_str() ) << std::endl ;
                    std::cout << "   └───────────────────────────────────────────────────────────────────────┘" << std::endl ;
                }
            }
        }

        void
        Controller::print_physics_stats()
        {
            BELFEM_ASSERT( mCommRank == 0,
                "Only root may call Controller::print_physics_stats() " );

            // same 20/25/24 column geometry as print_header. The divider
            // continues the section the caller left open, or opens a fresh
            // one when the last thing printed was a closed box
            std::cout << ( mBoxSectionOpen ?
                "   ├────────────────────┼─────────────────────────┼────────────────────────┤" :
                "   ├────────────────────┬─────────────────────────┬────────────────────────┤" )
                << std::endl ;

            // the temperature cell says what its number is: T_max is the
            // live maximum of the nodal field in a coupled run, T_bulk the
            // fixed value a magnetic-only run holds its materials at
            // ( get_Tmax returns the one that applies ). The label is
            // right-aligned so the colon sits in the same column either way
            const string tTcell = sprint( "%7s : %10.4f K  ",
                mKernel2 != nullptr ? "T_max" : "T_bulk", this->get_Tmax() );

            if ( ! mMesh->field_exists( "JJCz" ) )
            {
                // no postprocessor has ever written the field: the cell stays
                // EMPTY rather than carrying a zero that would read as a
                // measurement
                std::cout << sprint( "   │%s│                         │%s│",
                    this->excitation_cell().c_str(), tTcell.c_str() ) << std::endl ;
            }
            else if ( mPostProcessed )
            {
                // the postprocessor ran on THIS state, so the field is live
                std::cout << sprint( "   │%s│  J/Jc_max : %10.6f  │%s│",
                    this->excitation_cell().c_str(), this->get_JJCmax(), tTcell.c_str() ) << std::endl ;
            }
            else
            {
                // J/Jc is refreshed only where the postprocessor runs, which
                // is on saved steps. Between two saves the number is the one
                // measured at the last save, and the asterisk says so --
                // recomputing it here would rescan the SAME old field and
                // dress a stale value as a fresh one
                std::cout << sprint( "   │%s│  J/Jc_max : %10.6f ⃰ │%s│",
                    this->excitation_cell().c_str(), mLastJJcMax, tTcell.c_str() ) << std::endl ;
            }

            std::cout << "   ├────────────────────┴─────────────────────────┴────────────────────────┤" << std::endl ;
            mBoxSectionOpen = false ;
        }
    } // namespace fem
}
