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
#ifndef CL_FEM_CONTROLLER_HPP
#define CL_FEM_CONTROLLER_HPP
#include "cl_FEM_DofMgr_SolverData.hpp"
#include "typedefs.hpp"
#include "cl_Timer.hpp"
#include "cl_Input_Section.hpp"
#include "cl_Circuit.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_ShiftRegister.hpp"

namespace belfem
{
    class Mesh ;

    namespace fem
    {
        class Kernel ;
        class IWG_Timestep ;

        /**
         * @brief Nonlinear iteration controller: relaxation, timestep adaptation and convergence policy.
         *
         * @ingroup grp_fem_kernel
         * @see @ref fem_kernel_nonlinear_controller_theory
         */
        class Controller
        {
            const proc_t   mCommRank ;
            const proc_t   mCommSize ;
            Kernel       * mKernel = nullptr ;
            IWG_Timestep * mEquation = nullptr ;
            Mesh         * mMesh = nullptr ;
            Kernel       * mKernel2   = nullptr ;
            IWG_Timestep * mEquation2 = nullptr ;
            Mesh         * mMesh2     = nullptr ;
            Circuit      * mCircuit   = nullptr ;

            real mTime = 0.0 ;
            real mDeltaTime = BELFEM_QUIET_NAN ;
            real mDeltaTimeTemporary = BELFEM_QUIET_NAN ; //Temporary time step to reach the time steps to be saved
            real mDeltaTime0 = BELFEM_QUIET_NAN ;

            //! the deck's "initial timestep", kept as parsed. It is the
            //! re-entry cap for a warm restart: load_memdump keeps a dumped
            //! step only up to this value. MEASURED: the dt of the first
            //! post-restart step decides the restart's outcome and cannot be
            //! repaired afterwards. Working hypothesis, not proof: STRUMPACK
            //! derives solver-entry state ( ordering, matching, equilibration,
            //! tiny-pivot threshold ) from the FIRST Jacobian of the process
            //! and reuses it. Stored rather than snapshotted inside
            //! load_memdump so the invariant does not depend on that being
            //! called before any timestep adaptation
            real mDeltaTimeInitial = BELFEM_QUIET_NAN ;
            uint mRunningTimeStep = 0 ;
            uint mMeshTimeStep = 1 ;

            real  mTime2 =  0.0 ;
            real  mDeltaTime2 = BELFEM_QUIET_NAN ;
            real mDeltaTimeTemporary2 = BELFEM_QUIET_NAN ; //Temporary time step to reach the time steps to be saved

            // - - - - - - begin user settings - - - - - -

            real mAlpha = 0.5 ;
            real mBeta  = 1.1 ;
            real mGamma = 0.4 ;

            uint mMinNumIterations = 2 ;
            uint mMaxNumIterations = 100 ;
            uint mIterationTarget = 20 ;

            //! consecutive iterates with epsilon above the +10 dB divergence
            //! bar; the cut fires at 3, so the Eq.-14 alpha damping gets two
            //! halvings to catch a flux-front Picard overshoot before the
            //! timestep is abandoned ( an instantaneous cut was trend-blind:
            //! the same overshoot recovered in-step whenever it peaked
            //! below the bar ). Reset per attempt and on any iterate at or
            //! below the bar
            uint mDivergenceStrikes = 0 ;

            // stagnation guard: if the mean absolute deviation of the last
            // mStallWindow residuals (in dB) stays below mStallBand while the
            // residual is still above tolerance, the timestep is cut instead of
            // grinding on to mMaxNumIterations (see iterate_coupled).
            // The band must sit ABOVE the noise of a floored residual, or the
            // Newton->Picard demotion this guard exists for never fires: the
            // sidecoatings floor wandered +-0.05..0.2 dB against the old
            // 0.001 dB band and burned the full watchdog window instead.
            // Genuine convergence moves >= O(1) dB per iterate, so 0.2 dB
            // separates the two regimes. Deck override: "stall tolerance"
            uint mStallWindow = 5 ;
            real mStallBand   = 0.2 ; // dB

            // progress watchdog: cut the timestep when the residual has not
            // reached a NEW MINIMUM for this many iterations while still above
            // 10x tolerance. Catches what the flat-band guard cannot see:
            // omega-sawtooth limit cycles ( locally "converging" on any short
            // window ) and slow monotone creep of the thermal residual at the
            // relaxation floor. 0 disables. Input key: "watchdog window".
            uint mWatchdogWindow  = 30 ;
            uint mWatchdogWindow2 = 30 ;

            //! relaxation seen by the previous watchdog_* call, per field: the
            //! spare clause reads "omega grew since last iterate" as a line
            //! search recovering rather than a stall. NaN = no history, never
            //! spare. Reset at ATTEMPT boundaries only ( initialize_*, reset_* )
            //! -- NOT at promotion/escalation re-anchors, which restart the
            //! best-residual clock but do not invalidate last iterate's omega
            real mWatchdogOmegaPrev  = BELFEM_QUIET_NAN ;
            real mWatchdogOmegaPrev2 = BELFEM_QUIET_NAN ;

            uint mMinNumIterations2 = 2 ;
            uint mMaxNumIterations2 = 100 ;
            uint mIterationTarget2 = 20 ;

            uint mMaxNumIterationsDiv = 10 ; //Number of diverging iterations before resetting the relaxation
            uint mNumIterationsDiv = 0; //Number of diverging iterations
            // conditioning diagnostic, per field. It describes a property of
            // ONE solve, and the two fields configure their solvers
            // separately ( "linear magnetic" / "linear thermal", falling back
            // to "linear" ), so the flag lives with them.
            // "compute conditioning" in the solver section still sets both.
            //
            // These gate the EIGENVALUE estimate and nothing else. The MUMPS
            // Arioli/Demmel/Duff numbers have their own flag below -- the two
            // diagnostics are independent, and either may be asked for alone
            bool mComputeConditioning  = false ; // magnetic; always costs ARPACK
            bool mComputeConditioning2 = false ; // thermal

            // "mumps error analysis", per field, read from the same three
            // places as the flag above. Gates MUMPS ICNTL(11) ONLY -- the
            // COND1 / COND2 / omega2 triple that feeds the "MUMPS ADD" footer
            // rows. A no-op on any other library, which is why the key names
            // MUMPS: see the warnings raised for a field that asks for it
            // without one
            bool mMumpsErrorAnalysis  = false ; // magnetic
            bool mMumpsErrorAnalysis2 = false ; // thermal

            // guards check_thermal_diagnostics() against warning twice: both
            // attach paths call it, and a caller can use both
            bool mThermalDiagnosticsChecked = false ;

            // whether set_params has run. The thermal diagnostics must not be
            // judged against the DEFAULT flags, so their check requires both
            // prerequisites -- parsed params AND an attached kernel -- and
            // runs from whichever call arrives second, in any order
            bool mParamsSet = false ;

            //! neutral band of the relaxation adaptation: a residual within
            //! (1 + band) of the previous iterate is a wobble at the noise
            //! floor, not divergence — omega and the divergence counter hold
            //! instead of decaying. Without it, improvement at a floored
            //! residual is a coin flip and the strict-decrease AIMD collapses
            //! omega geometrically to mOmegaMin ( sidecoatings trace:
            //! 0.5 -> 1.3e-4 in 30 iterations at a flat -46 dB floor ).
            //! 0.05 = 0.2 dB, matched to the stall band above
            real mOmegaNoiseBand = 0.05 ;

            // MIT-3a: consecutive first-trial Newton steps that at least
            // halved the residual. The geometric omega recovery
            // ( tGrowth = 2.0 ) requires TWO of them, because a single
            // qualifying step is not evidence of contraction near a noise
            // floor -- on tapestack3d the doubling fired once ( omega
            // 0.483 -> 0.966 ) and regressed on the very next iterate.
            // SCOPE: this addresses the x2 rule ONLY; a later loss in the
            // same trace came from ordinary arctan growth and is NOT fixed
            // here. A healthily
            // contracting Newton earns the streak on its second step and
            // keeps the acceleration; a lucky step never does.
            //
            // Reset on: timestep start ( both drivers ), reject,
            // promotion, escalation, and every demotion ON THE COUPLED
            // PATH -- the evidence belongs to one uninterrupted Newton
            // run, not to the timestep. NOTE this is NOT simply
            // "wherever mFirstFlip is reset": iterate_magnetic's
            // promotion clears mFirstFlip without touching the streak.
            // That twin never increments and never consumes the streak
            // ( it has no x2 rule ), and initialize_magnetic zeroes at
            // every step start, so the omission is inert -- but do not
            // rely on the mFirstFlip locator when adding a path
            uint mNewtonTrustStreak = 0 ;

            real mEpsilonSwitch = 1e-4 ; // residual at which Picard hands off to Newton (Messe et al. 2023, Eq. 13)
            real mRelativeEpsilonTarget = 1e-6 ; // relative error criterion
            //! absolute error criteria: DISABLED by default ( 0.0 ) — the raw
            //! ||A x - b|| is dimensional, so a universal default is unit- and
            //! load-scale dependent ( Codex+Grok RQ3 ). Opt in per deck via
            //! "absolute tolerance" in the nonlinear ( thermal ) sections.
            real mAbsoluteEpsilonTarget = 0.0 ;

            real mEpsilonSwitch2 = 1e-3 ; // thermal residual at which Picard hands off to Newton
            real mRelativeEpsilonTarget2 = 1e-6 ; // relative thermal error criterion
            real mAbsoluteEpsilonTarget2 = 0.0 ;  // absolute thermal error criterion ( see above: opt-in )

            //! coupled runs only: skip the thermal update while the magnetic
            //! residual exceeds this gate. DISABLED by default ( ts17
            //! 2026-07-22 showed it converts a mutual divergence spiral into
            //! a limit cycle pinned at the gate instead of rescuing it — a
            //! lost state needs a timestep cut, not update reordering ).
            //! Kept as an experiment knob: "update gate" in the nonlinear
            //! thermal section.
            real mThermalUpdateGate = BELFEM_REAL_MAX ;

            //! set while the thermal update is frozen by the gate; forces a
            //! Picard restart on the next thermal solve ( the state moved
            //! while thermal was frozen — a cold Newton start from there
            //! diverges, cf. ts12 ). Inert while the gate is disabled.
            bool mThermalFrozen = false ;

            real mSimulationTime = BELFEM_REAL_MAX ;
            bool mAdaptTimestep = true ;
            real mSaveEvery = 0.0 ;
            bool mLastSave = false ;
            bool mSave = true ;
            
            // - - - - - - end user settings - - - - - -

            // todo: add flag if we want to reset omega at each timestep

            // end user settings
            real mOmegaNewton = 1.0 ;
            real mOmegaPicard = 1.0 ;
            real mOmegaNewton2 = 1.0 ;
            real mOmegaPicard2 = 1.0 ;
            bool mJustPicard = false ;
            bool mJustPicard2 = false ;

            //! thermal Picard<->Newton demotions within the current attempt;
            //! the third one latches mJustPicard2 ( chatter guard )
            uint mThermalFlipCount = 0 ;

            //! Anderson mixing depths ( 0 = off ), forwarded to the kernels'
            //! SolverData. Opt-in: off by default; enable via
            //! timestep { anderson stabilization : true ; } ( depths 3 / 1 )
            //! or per field with nonlinear { anderson depth : n ; }
            uint mAndersonDepth  = 0 ;
            uint mAndersonDepth2 = 0 ;

            //! committed-pair counters ( = mixing window fill ). Tracked HERE,
            //! not queried from SolverData: the history lives on the master
            //! rank only, and the omega adaptation below must take the same
            //! branch on every rank. Used to suppress the AIMD growth while
            //! Anderson supplies the step direction ( O2 ).
            uint mAndersonCommits  = 0 ;
            uint mAndersonCommits2 = 0 ;

            // Picard-breakdown rescue ( see try_escalate_to_newton ). When the magnetic
            // Picard iteration diverges at a residual far above mEpsilonSwitch -- so the
            // normal Picard->Newton promotion never fires -- the controller escalates to
            // the configured Newton tangent once before cutting the timestep. mForceNewton
            // pins the solver to Newton for the rest of the attempt; mNewtonEscalated caps
            // the escalation at one per timestep attempt so it cannot ping-pong.
            bool mForceNewton = false ;
            bool mNewtonEscalated = false ;

            // Thermal twin of the pair above ( see try_escalate_thermal_to_newton,
            // ). The thermal Picard residual can floor bit-flat ABOVE
            // mEpsilonSwitch2, and since the coupled Picard->Newton promotion fires
            // only below the switch, Newton is gated behind progress that only
            // Newton can make. Same lifecycle as the magnetic pair: cleared in
            // initialize_timestep / initialize_thermal, at most one escalation
            // per timestep attempt.
            bool mForceNewton2 = false ;
            bool mNewtonEscalated2 = false ;

            SolverAlgorithm mAlgorithm = SolverAlgorithm::Picard ;
            SolverAlgorithm mAlgorithmThermal = SolverAlgorithm::Picard ;
            //! default: BDF1, the validated baseline ( Messe et al. 2023 §4 );
            //! the deck picks another scheme via timestep { scheme : ... ; }
            EulerMethod          mTimeStepping = EulerMethod::BackwardDifference1 ;

            real mOmega0 = BELFEM_QUIET_NAN ;

            real mOmegaMin = 0.001 ;
            real mOmegaMax = 1.0 ;
            real mOmegaMin2 = 0.1 ;
            real mOmegaMax2 = 1.0 ;

            real mTime0 = 0.0 ;
            real mTime02  = 0.0 ;
            Timer * mTimer = nullptr ;

            real mEpsilon  = BELFEM_REAL_MAX ;
            real mEpsilon0 = BELFEM_REAL_MAX ;

            real mEpsilon2 = BELFEM_REAL_MAX ;
            real mEpsilon20 = BELFEM_REAL_MAX ;

            //! absolute thermal residual norm ||A x - b|| of the last solve;
            //! consumed by the absolute-tolerance escape of the certified
            //! exit ( mAbsoluteEpsilonTarget2 was parsed but dead )
            real mEpsilonAbs2 = BELFEM_REAL_MAX ;

            //! absolute magnetic residual norm, the twin of mEpsilonAbs2:
            //! consumed by the certified-exit heads so that a deck-set
            //! "absolute tolerance" in the nonlinear ( magnetic ) section
            //! is live ( it was parsed but dead — Claude+Grok jury
            //! 2026-08-09 ). With the default target 0.0 the escape never
            //! fires and the behaviour is unchanged
            real mEpsilonAbs = BELFEM_REAL_MAX ;

            //! thermal stagnation exit: consecutive coupled iterations with a
            //! bit-flat thermal residual while the magnetic system is
            //! converged; further iterations provably change nothing ( the
            //! T-clamp fixed point, see compute_T_fem ). Trips mThermalStalled,
            //! which ends the timestep through the certified exit with a
            //! warning.
            uint mThermalFlatCount = 0 ;
            bool mThermalStalled = false ;

            //! per-attempt latch: the magnetic field reached its relative or
            //! absolute target at least once in this attempt. Feeds the
            //! thermal budget cut in iterate_coupled — a magnetic excursion
            //! back above target caused by continuing thermal updates must
            //! not disarm the budget ( the step-267 shape ).
            //! Reset per attempt ( initialize_timestep / reset_timestep ).
            bool mMagneticHitTarget = false ;

            //! trip flags: whether the previous trip of the
            //! certified-exit loop ran a magnetic body ( gates the
            //! Picard/Newton handoff — an idle certification trip must not
            //! promote ), and whether the current trip ended in a certified
            //! exit ( consumed by solve_* )
            bool mMagBodyRanLastTrip = true ;

            bool mTripExit = false ;

            //! progress-watchdog trackers: best residual of the current
            //! attempt and the iteration it was seen at
            real mBestEpsilon  = BELFEM_REAL_MAX ;
            real mBestEpsilon2 = BELFEM_REAL_MAX ;
            uint mBestEpsilonIteration  = 0 ;
            uint mBestEpsilonIteration2 = 0 ;

            //! consecutive timestep cuts caused by soft solver failures
            //! ( singular matrix ); a persistently singular system must
            //! eventually abort loudly instead of cutting forever. Cleared
            //! in finalize() on every successful timestep.
            uint mSolverFailCount = 0 ;

            //! timestep-control mode: true restores the legacy behaviour
            //! ( iteration-target sqrt rule, no rejection memory ). Internal
            //! switch by design — no input key until the PID mode is
            //! validated; see src/fem/doc/timestepping_strategy.md
            bool mUseLegacyTimestepControl = false ;

            //! PID gains on the cost error e = N_iterations / N_target,
            //! applied multiplicatively ( log-space PID, Valli, Carey &
            //! Coutinho 2002, CNM 18:131, doi:10.1002/cnm.475 ). The legacy
            //! sqrt rule is the pure-I point ( kP=0, kI=0.5, kD=0 ) of this
            //! controller. Derivative ships disabled: the iteration count is
            //! quantized and D-action amplifies its noise
            real mCtrlKp = 0.15 ;
            real mCtrlKi = 0.30 ;
            real mCtrlKd = 0.0 ;

            //! cost-error history e_n, e_{n-1}, e_{n-2} of ACCEPTED steps
            //! only; 1.0 = on target. Reset to 1.0 on every timestep cut —
            //! a rejected attempt's cost carries no information about the
            //! next accepted one, and stale samples would kick the P/D terms
            real mCtrlErr0 = 1.0 ;
            real mCtrlErr1 = 1.0 ;
            real mCtrlErr2 = 1.0 ;

            //! post-rejection growth hold ( SUNDIALS practice ): number of
            //! adjust_timestep calls after a cut in which the ratio is
            //! capped at 1, so the controller cannot re-climb the cliff it
            //! just fell off. Note the first accepted step after a cut skips
            //! adjust_timestep entirely ( the mIteration0 gate ), so growth
            //! resumes on the fourth accepted step after the cut
            uint mPostFailureHold      = 0 ;
            uint mPostFailureHoldSteps = 2 ;

            //! floor escalation ( Christian 2026-08-10: allow more iterations
            //! rather than aborting ): each cut demanded while Delta t already
            //! sits at mDeltaTimeMin doubles the iteration budgets, up to
            //! mFloorEscalationCap x the deck values, and the attempt retries
            //! AT the floor. The budgets return to the deck values with the
            //! next accepted step ( finalize ).
            uint mFloorRetries       = 0 ;
            uint mFloorEscalationCap = 4 ;

            //! backstop on the escalation: after this many consecutive floor
            //! retries the run stops with a diagnosis instead of spinning
            //! forever on a residual that does not respond to Delta t at any
            //! iteration budget ( unattended batch jobs would otherwise burn
            //! their whole allocation ). Everything up to the last ACCEPTED
            //! step is already on disk. 0 = no backstop, retry indefinitely.
            //! Input key: timestep { floor retries : n ; }
            uint mMaxFloorRetries = 20 ;

            //! deck-configured iteration budgets, the restore targets after
            //! a floor escalation; snapshot taken at the end of set_params
            uint mMaxNumIterationsDeck  = 100 ;
            uint mMaxNumIterations2Deck = 100 ;
            uint mWatchdogWindowDeck    = 30 ;
            uint mWatchdogWindow2Deck   = 30 ;

            //! warm-restart policy ( opt-OUT ): by
            //! default a rerun resumes from an existing memdump -- made
            //! VISIBLE by the warm-restart banner in load_memdump. A deck
            //! sets timestep { restart : false ; } to ignore the dump and
            //! start fresh ( the dump is overwritten at the first save )
            bool mAllowRestart = true ;

            //!
            uint mIteration = 0 ;
            uint mIteration2 = 0 ;
            uint mIterationTime = 0 ;
            uint mIteration0 = 0 ;
            uint mIteration02 = 0 ;
            uint mEigenAnalysisTime = 0 ;
            uint mPostprocesingTime = 0 ;

            //! condition numbers of the magnetic and thermal system matrices,
            //! refreshed EVERY TIMESTEP and printed in the footer. Reset to
            //! NaN at the top of each step so a step that captured nothing
            //! prints n/a rather than the previous step's figure. A field
            //! with "mumps error analysis" samples COND1/COND2 on the FIRST
            //! iterate of the step ( slots 1-3; the native estimate's drift
            //! across the step is negligible next to its own error ) ; a
            //! field with "compute conditioning" runs the eigen path at the
            //! END of the step, against the converged matrix, which is still
            //! intact after the last solve ( slot 0 ). The two keys are
            //! independent

            //! Diagnostic only — nothing in the iteration scheme
            //! consumes them. NaN means "not available": the field has no
            //! kernel, or its solver cannot supply kappa cheaply
            // slot 0 : the eigen estimate ( spectral ratio, or kappa_2 when
            //          the symmetric driver earned the name )
            // slot 1 : MUMPS ADD COND1
            // slot 2 : MUMPS ADD COND2
            // slot 3 : MUMPS omega2 -- NOT reported, it is the discriminator
            //          that says whether slot 2 is a measurement at all
            Vector< real > mConditionNumbers1  = { BELFEM_QUIET_NAN, BELFEM_QUIET_NAN, BELFEM_QUIET_NAN, BELFEM_QUIET_NAN };
            Vector< real > mConditionNumbers2  = { BELFEM_QUIET_NAN, BELFEM_QUIET_NAN, BELFEM_QUIET_NAN, BELFEM_QUIET_NAN };

            real mDeltaTimeMax = BELFEM_REAL_MAX ;
            real mDeltaTimeMin = 1e-10 ;
            bool mReset = false ;
            bool mResetThermal = false ;
            bool mFirstFlip = false ;
            bool mFirstFlip2 = false ;

            bool mIsFullyCoupled = true ;
            real mCouplingFactor = 1 ; //Number of thermal time steps within each magnetic time steps


            Cell < PhysicalBoundaryCondition * > mCircuitCurrentBCs ;
            Cell < PhysicalBoundaryCondition * > mCircuitVoltageBCs ;

            Vector< real > mLHS ;
            Vector< real > mLHS0 ;
            bool mFirstIVSave = true ;

            //! names of the I/U pairs save_IV writes, one per abstract dof,
            //! shared by the csv header and the mesh globals ( create_iv_names )
            Cell< string > mIVNamesI ;
            Cell< string > mIVNamesU ;

            // ring buffer of the last mStallWindow residuals (in dB)
            ShiftRegister< real > mResidualHistory ;

            Cell< Vector< real > > mBackupFields ;
            Vector< real > mBackupDofValues ;

            //! last J/Jc maximum seen, refreshed by get_JJCmax()
            real mLastJJcMax = 0. ;

            //! true when a three-column row has been printed and its section
            //! deliberately left OPEN, so the physics row that print_footer()
            //! adds continues it with a "┼" divider instead of opening a
            //! second section under the first
            bool mBoxSectionOpen = false ;

            //! true when the postprocessor ran on the state of the step that
            //! just finished. The J/Jc field is written there and nowhere else,
            //! so this decides whether print_physics_stats() reports a live
            //! number or the last one it saw. Set in finalize().
            bool mPostProcessed = false ;

        public:

            Controller( Kernel * aKernel,  Kernel * aKernel2 = nullptr ) ;

            ~Controller();

            void
            initialize_timestep();

            void
            initialize_magnetic();

            void
            initialize_thermal();

//------------------------------------------------------------------------------

            //! certified-exit loops: one call per timestep, owning
            //! the complete nonlinear iteration; the committed state is
            //! always head-measured before an exit. Drivers retry on reset()
            void
            solve_coupled();

            void
            solve_magnetic();

            void
            solve_thermal();

            bool
            solve_circuit();

            void
            finalize( const bool aPostProcess = true );

            const real &
            time() const;

            const real &
            time_thermal() const;

            real
            epsilon() const ;

            uint
            iteration() const ;


            bool
            is_fullycoupled() const ;

            bool
            reset() const ;

            void
            save( const std::string & aFilename );

            void
            save_IV( const std::string & aFilename );

            /**
             * names the I/U pair of every abstract dof: after the deck label
             * of the terminal condition that drives it, else positional
             * ( rank 0 only, once; see the implementation for the rules )
             */
            void
            create_iv_names();

            void
            set_params( const input::Section * aSection );

            void
            set_thermal_kernel( Kernel * aKernel );

            real
            simulation_time()  const ;

            bool
            save()  const ;

            void
            set_circuit( Circuit * aCircuit );

            Circuit *
            circuit() ;

            Kernel *
            kernel() ;

            Kernel *
            thermal_kernel() ;

            void
            save_memdump( const string & aPath );

            void
            load_memdump( const string & aPath );

            EulerMethod
            euler_method() const ;


        private:

            // imposes the voltage/current boundary conditions on the RHS ( rank 0 only,
            // after assembly and before the solve ); shared by the magnetic iterate paths
            //! the per-trip bodies and the old outer predicates are
            //! internals of the solve_* loops — the uncertified outer path
            //! must not survive as a callable API ( round-3 R-G )

            void
            iterate_coupled();

            void
            iterate_magnetic();

            void
            iterate_thermal();

            void
            impose_voltage_bcs();

            // env-gated diagnostic dump of the assembled system, shared by both
            // iterate paths; off unless BELFEM_DUMP_SYSTEM is set.
            // aTag names the field in the file name, and aCount is that field's
            // OWN budget -- one shared counter let whichever field assembled
            // first spend all four dumps and starve the other
            void
            dump_system_if_requested( DofManager * aDofMgr,
                                      const char * aTag,
                                      int        & aCount );

            //! per-field dump budgets. Members rather than function-local
            //! statics: a static is shared by every Controller in the process,
            //! which is wrong the moment a second one exists
            int mDumpCountMagnetic = 0 ;
            int mDumpCountThermal  = 0 ;


            // magnetic stagnation guard: if the residual window is flat ( deviation below
            // mStallBand ), Newton falls back to Picard internally; returns true only when
            // Picard itself has stalled, signalling the caller to reset_timestep()
            bool
            magnetic_stagnation_forces_reset();

            // progress watchdogs ( see mWatchdogWindow ): update the best-residual
            // tracker and return true when the attempt has made no new minimum for
            // the configured window while still far from tolerance -- UNLESS the
            // field's relaxation grew since the previous call, which marks a line
            // search recovering from an overshoot rather than a stall ( replayed
            // 2026-08-21: fires on the omega-pinned grinds, spares the AIMD
            // recovery whose cut used to cascade the timestep down ).
            //
            // aOmega must be the relaxation that PRODUCED the current residual --
            // the caller's post-line-search tOmega / tOmega2 -- never re-derived
            // from the live algorithm, which a same-call escalation or a
            // stagnation demotion may already have flipped or mutated.
            bool
            watchdog_magnetic( const real aOmega );

            bool
            watchdog_thermal( const real aOmega );

            // Anderson bookkeeping ( plan §3.2 ): SolverData stages the pair,
            // the controller decides its fate. Helpers keep the commit
            // counters and the per-kernel SolverData in sync; all are cheap
            // no-ops while the respective depth is 0.
            void anderson_commit_magnetic();
            void anderson_discard_magnetic();
            void anderson_clear_magnetic();
            void anderson_commit_thermal();
            void anderson_clear_thermal();

            // Picard-breakdown rescue: if the configured terminal algorithm is Newton and
            // the magnetic Picard iteration has broken down ( residual stuck far above
            // mEpsilonSwitch, so the normal promotion never fires ), switch to the Newton
            // tangent for the rest of this timestep attempt instead of cutting Δt. Returns
            // true if it escalated ( caller keeps iterating ), false otherwise ( caller
            // resets ). Escalates at most once per attempt, guarded by mNewtonEscalated.
            bool
            try_escalate_to_newton();

            // Thermal twin: called from the coupled stagnation
            // detector when the thermal Picard residual sits bit-flat above
            // mEpsilonSwitch2 with Newton configured but not running. Returns
            // true if it escalated. At most once per attempt ( mNewtonEscalated2 ).
            bool
            try_escalate_thermal_to_newton();

            void
            reset_timestep();

            void
            reset_thermal();

            void
            adjust_timestep();

            void
            compute_circuit_current() ;

            void
            compute_conditioning();

            // MUMPS computes its condition estimate during the SOLVE, so the
            // error analysis is armed for exactly ONE solve per timestep --
            // the first iterate -- and disarmed as soon as the value is read.
            // Leaving it armed ( the pre-2026-08-10 behaviour ) paid for the
            // analysis on every solve of every iteration and read one value.
            // The first iterate is also the better sample: it is assembled at
            // the converged previous step, so kappa is measured at a
            // comparable state in every timestep and the series is a trend.
            // All four are no-ops unless that FIELD's error-analysis flag is
            // set ( mMumpsErrorAnalysis / mMumpsErrorAnalysis2 -- NOT the
            // eigen flags, which gate a different diagnostic ) and its own
            // solver is MUMPS. They run on EVERY rank, unguarded:
            // the ICNTL setting must match across the communicator, and the
            // branch is taken identically everywhere
            void
            arm_conditioning_magnetic();

            void
            arm_conditioning_thermal();

            void
            capture_conditioning_magnetic();

            void
            capture_conditioning_thermal();

            //! selects the thermal field's ARPACK driver, once per attach.
            //! Called from BOTH attach paths -- set_thermal_kernel() and the
            //! two-argument constructor -- because either can be the only one
            //! a caller uses. Idempotent, unconditional, and never rank-guarded
            void
            setup_thermal_eigen();

            //! per-field validation of the two diagnostic keys: the eigen
            //! estimate without MUMPS, and the error analysis without MUMPS
            //! ( a silent no-op ). Split by field because the two are
            //! validated at different times -- the magnetic kernel exists at
            //! set_params, the thermal one usually does not
            void
            check_magnetic_diagnostics();

            //! thermal twin, called from BOTH attach paths. One-shot.
            void
            check_thermal_diagnostics();

//-----------------------------------------------------------------------------

            //! setup gate: behind an ITERATIVE linear solver, the nonlinear
            //! loop can never push its residual below the linear relative
            //! tolerance — the Krylov iteration stops working the moment its
            //! target is met, so the linear tolerance IS the delivered
            //! accuracy ( headroom rule, Messe et al. 2023 §2.7 context ).
            //! Hard error when the deck asks for the impossible. The
            //! predicate is deliberately narrow — solver type == PETSc —
            //! NOT because direct libraries deliver machine epsilon
            //! ( refuted by measurement ), but because STRUMPACK's working
            //! default pairing is itself lin > nonlin: its outer GMRES
            //! overshoots the exit test by decades, so this predicate
            //! would refuse a configuration proven to work ( see the
            //! rationale in the implementation ). The lossy
            //! ( blr ) case, where the headroom rule DOES apply, is
            //! handled by check_compression_headroom below — two tiers:
            //! zero-or-negative headroom errors, thin-positive warns
            void
            check_iterative_solver_headroom(
                    DofManager * aDofMgr,
                    const real   aNonlinTol,
                    const char * aFieldName ) const ;

//-----------------------------------------------------------------------------

            //! sibling of the gate above for LOSSY factorizations
            //! ( compression scheme : blr on MUMPS/STRUMPACK ): the drop
            //! tolerance is the delivered linear accuracy. Two tiers:
            //! cutoff >= nonlinear tolerance is a hard ERROR
            //! ( the PETSc-gate impossibility ),
            //! thin-but-positive headroom ( < 2 decades ) a rank-0
            //! WARNING — the legitimate memory-bound trade, eyes open.
            //! Both on the effective cutoff, default included, so the
            //! inherited 1e-8 is held to the same bar as a stated value.
            //! Warning printed unguarded: -v 0 must not swallow it
            void
            check_compression_headroom(
                    DofManager * aDofMgr,
                    const real   aNonlinTol,
                    const char * aFieldName ) const ;

//-----------------------------------------------------------------------------

            //! warm-restart companion of restore_history_state:
            //! distributes the numbered history fields ( the equation's
            //! history_field_labels — deliberately NOT in all_fields, the
            //! per-step synch must not pay for them ) to the workers and
            //! hard-verifies every level the resumed step's shift will
            //! read, on every rank, BEFORE the integrator is told it has
            //! history. COLLECTIVE — call on all ranks
            void
            synchronize_history_fields(
                    IWG_Timestep * aEquation,
                    DofManager   * aDofMgr,
                    const uint     aStepCount,
                    const string & aPath,
                    const char   * aEquationName );

            void
            print_header();

            void
            print_line( const real aOmega, const real aOmega2 );

            void
            print_line_magnetic( const real aOmega );

            void
            print_line_thermal( const real aOmega );

            void
            print_footer();

            //! warning box for the thermal stagnation exit: reports the
            //! accepted thermal residual, the master-rank max T and the
            //! smallest material T_max ( the table-clamp diagnostic )
            void
            print_thermal_stall_warning();

            real
            get_Tmax();

            real
            get_Imax();

            /**
             * the leftmost cell of the timestep footer, 20 characters wide so
             * the box geometry survives. A deck driven by a transport current
             * reports that current in amperes; a deck driven only by a
             * background field reports the field as a flux density in TESLA,
             * whatever unit the deck itself used to state it
             */
            string
            excitation_cell();

            real
            get_JJCmax();

            void
            print_physics_stats();

            void
            reset_dotQ();

            void
            collect_dotQ();

        };

        inline
        Circuit *
        Controller::circuit()
        {
            return mCircuit;
        }

        inline
        Kernel *
        Controller::kernel()
        {
            return mKernel;
        }

        inline
        Kernel *
        Controller::thermal_kernel()
        {
            return mKernel2;
        }

        inline
        EulerMethod
        Controller::euler_method() const
        {
            return mTimeStepping;
        }

    }
}

#endif //CL_FEM_CONTROLLER_HPP
