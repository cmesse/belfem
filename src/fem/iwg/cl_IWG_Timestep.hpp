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

#ifndef BELFEM_CL_IWG_TIMESTEP_HPP
#define BELFEM_CL_IWG_TIMESTEP_HPP


#include "cl_IWG.hpp"
#include "en_SolverEnums.hpp"
#include "cl_TimestepMatrices.hpp"

namespace belfem
{
    namespace fem
    {
        class Element ;

        /**
         * @brief Base class for transient problems.
         *
         * @ingroup grp_fem_iwg
         * @see @ref fem_iwg_iwg_usage_guide
         */
        class IWG_Timestep: public IWG
        {
            // Unit-test seam, no runtime cost and no public API widened.
            //
            // The BDF coefficients ( mAlpha, mBeta ) and the step-size history
            // that determines them ( mH, mStepCount ) are private on purpose.
            // The production writers are shift_fields, reset_fields,
            // restore_savepoint — each needing an initialized DofManager —
            // and, since 2026-08-15, the validated memdump seam
            // restore_history_state ( which needs none ). Without this
            // declaration the
            // variable-step coefficient formulas below are reachable by no
            // test at all: they are pure arithmetic, but their inputs cannot
            // be posed and their outputs cannot be read.
            //
            // The probe is defined in tests/fem/test_BdfTimestepMethod.cpp and
            // exists only there. Nothing in src/ may use it.
            friend class BdfCoefficientProbe ;

            EulerMethod mMethod = EulerMethod::BackwardDifference1 ;



            Cell< Cell< Vector< real > * > > mFieldData ;

            real mAlpha = 1.0 ;
            Vector< real > mBeta ;
            Vector< real > mH ;

            // true if BDF coefficients must be recomputed before the next assembly
            bool mCoeffsDirty = true ;

            // savepoint storage: deep copies of the dof fields and the step-size
            // history, so a rejected coupled timestep can restore the state
            // across any number of shift_fields calls ( see make_savepoint )
            Cell< Cell< Vector< real > > > mFieldSnapshot ;
            Vector< real > mHSnapshot ;
            real mDeltaTimeSnapshot = BELFEM_QUIET_NAN ;
            real mHDroppedSnapshot = BELFEM_QUIET_NAN ;
            uint mStepCountSnapshot = 0 ;
            bool mHaveSavepoint = false ;

            // the step size that falls off the end of mH on a shift: the
            // rotation destroys it, but reset_fields must restore it for a
            // BDF5 retry, which reads mH( 3 )
            real mHDropped = BELFEM_QUIET_NAN ;

            // number of timesteps started ( shifted ), drives the startup order ramp
            uint mStepCount = 0 ;

            // order of the scheme actually running this step ( startup ramp
            // may run below the configured mOrder ); consumed by collect_qhist
            uint mOrderActive = 1 ;

            // remember whether the scheme was configured with a stiffness matrix
            bool mHaveStiffness = true ;

            // timestepping pointer as configured by set_timestepping_method
            void
            ( IWG_Timestep:: * mTimestep )(
                    Matrix< real > & aJ,
                    Vector< real > & aRHS );

            // scheme that actually runs this step; differs from mTimestep only
            // during the startup ramp of a BDF-p run ( step 1: BDF1, step 2: BDF2, ... )
            void
            ( IWG_Timestep:: * mTimestepActive )(
                    Matrix< real > & aJ,
                    Vector< real > & aRHS ) = nullptr ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            //real mDeltaTime = 1.0 ;
            uint mOrder = 1 ;

            TimestepMatrices * mTimeStepMatrices ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Timestep( const IwgType aType,
                          const ModelDimensionality aDimensionality,
                          const IwgMode aMode=IwgMode::Iterative,
                          const SymmetryMode aSymmetryMode=SymmetryMode::Unsymmetric,
                          const DofMode      aDofMode=DofMode::AllBlocksEqual,
                          const SideSetDofLinkMode aSideSetDofLinkMode=SideSetDofLinkMode::FacetOnly );


            ~IWG_Timestep() override ;

            // set the timestepping method
            void
            set_timestepping_method( const EulerMethod aMethod, const bool aHaveStiffness = true ) override;

            // return the timestepping method
            EulerMethod
            method() const override;

            // BDF order actually running this step ( startup ramp )
            uint
            order_active() const;

            // BDF order the UPCOMING assembly will run; ask this, not
            // order_active(), before the step has been assembled
            uint
            order_pending() const;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) override;

//------------------------------------------------------------------------------

            void
            shift_fields() override;

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

            /**
             * compute BDF coefficients (alpha, beta) for variable-step BDF
             * and select the scheme for the upcoming step (startup ramp).
             * Must be called AFTER delta_time is set to the upcoming step
             * and AFTER shift_fields has been called. Called lazily from
             * compute_jacobian_and_rhs when shift_fields/reset_fields have
             * flagged the coefficients as dirty.
             */
            void
            compute_bdf_coefficients();

//------------------------------------------------------------------------------

            /**
             * The beta-weighted dof history of the scheme running this step:
             *
             *     qhist = beta0*q0 - beta1*q1 + beta2*q2 - beta3*q3 + beta4*q4
             *
             * truncated to the order actually running, where q0 is the previous
             * step, q1 the one before it, and so on.
             *
             * WHAT IT IS FOR. A BDF step evaluates the time derivative as
             * ( alpha*q - qhist ) / 1, so the residual contracts the mass matrix
             * against exactly that combination. Anything that differentiates the
             * mass term must therefore differentiate against THIS vector, or the
             * Jacobian silently stops being the tangent of the residual it is
             * paired with. That is why this function exists rather than each
             * producer assembling its own history: it is the single source of
             * truth, used by the bdf2-5 right-hand sides and by every
             * nonlinear-mass producer ( maxwell::phi_ferro_newton and the
             * maxwell / thermal h-kernels ).
             *
             * FOUR THINGS THAT ARE EASY TO GET WRONG:
             *
             * 1. The alternating SIGNS live here, not in mBeta. compute_bdf_
             *    coeffs_* store beta as positive magnitudes; this function
             *    applies the +-+- pattern. Do not re-apply it at a call site.
             *
             * 2. alpha is NOT included. qhist is only the history half; the
             *    caller pairs it with alpha*q for the current state.
             *
             * 3. The order used is mOrderActive, not mOrder. During the startup
             *    ramp a BDF-p run executes at reduced order for its first steps,
             *    so the number of history terms — and the beta values, which also
             *    depend on the variable step sizes — change from step to step.
             *    At order <= 1 the result is exactly q0, with no beta applied.
             *
             * 4. The returned reference aliases the calculator's qswap() scratch
             *    buffer, which belongs to the CURRENTLY LINKED ELEMENT. It stays
             *    valid only until the next collect_qhist() call or the next
             *    element link. Read it, use it, do not store it: caching it
             *    across elements, or holding it across anything that may touch
             *    qswap(), reads another element's history.
             *
             * @return reference to the calculator's qswap() buffer, filled with
             *         the history combination described above
             */
            const Vector< real > &
            collect_qhist();

//------------------------------------------------------------------------------

            void
            reset_fields() override;

//------------------------------------------------------------------------------

            /**
             * save a deep copy of all dof fields, the step-size history and
             * the step counter. Call at the start of a timestep, BEFORE
             * shift_fields. Unlike reset_fields ( which reverses exactly one
             * shift ), restore_savepoint undoes any number of shifts — needed
             * for the thermal equation, which sub-steps several times within
             * one magnetic timestep.
             */
            void
            make_savepoint();

//------------------------------------------------------------------------------

            /**
             * restore the state stored by make_savepoint
             */
            void
            restore_savepoint();

//------------------------------------------------------------------------------

            /**
             * Export the multi-step integrator state for the memdump, so a
             * warm restart can resume at full BDF order instead of
             * re-anchoring the ramp at order 1 ( the restart cliff observed
             * 2026-08-15: a BDF1 re-entry into an active state is a
             * different tangent and detonated the Newton promotion ).
             *
             * The triple is ( mH, mStepCount, mDeltaTime ) and the third
             * member is load-bearing: at save time this object's mDeltaTime
             * still holds the LAST COMPLETED step size h_n — the controller
             * assigns the upcoming step only AFTER shift_fields — and the
             * first post-restore shift pushes exactly this value into
             * mH( 0 ). The controller's own ( already dumped ) delta_time is
             * the NEXT step and must not be confused with it.
             */
            void
            save_history_state(
                Vector< real > & aH,
                uint           & aStepCount,
                real           & aLastDeltaTime ) const ;

//------------------------------------------------------------------------------

            /**
             * the dof labels ( aParents ) and their numbered history
             * labels ( aHistory, label + "0".."order-1" ) — the fields
             * shift_fields() rotates. Single source of the naming
             * convention shared with create_old_dof_fields(); the
             * controller's warm-restart field synch is the consumer
             * ( these labels are deliberately NOT in all_fields )
             */
            void
            history_field_labels(
                Cell< string > & aParents,
                Cell< string > & aHistory,
                const uint       aDepth = 0 ) const ;

//------------------------------------------------------------------------------

            /**
             * Counterpart of save_history_state. Validates and restores;
             * returns false — leaving the cold-start ramp untouched — when
             * the data does not validate:
             *   - aH must match mH in length, all entries finite and >= 0
             *   - aStepCount == 0 is a reject ( it would silently ramp to
             *     BDF1 behind a full history )
             *   - aLastDeltaTime must be finite and > 0; it is HISTORY and
             *     is deliberately NOT clamped to the deck's Δt window
             *   - the slots a resumed step and its possible retry read must
             *     be positive: aH( 0 .. r-1 ) with
             *     r = min( aStepCount - 1, mOrder - 1, 4 ); note mH( 4 ) is
             *     the never-written spare slot and zero there is NORMAL
             * On success mStepCount restores capped at mOrder ( the ramp
             * saturates there ), mDeltaTime takes aLastDeltaTime, and the
             * coefficients are marked dirty for the lazy recompute. mOrder
             * and mMethod stay deck authority. mHDropped stays NAN: it is
             * unreachable before the first shift refills it from mH( 3 ).
             */
            bool
            restore_history_state(
                const Vector< real > & aH,
                const uint             aStepCount,
                const real             aLastDeltaTime ) ;

 //------------------------------------------------------------------------------

            TimestepMatrices *
            matrices() override ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * a default interface to compute the matrices of a timestepping
             * scheme of shape of
             * \f$ M \, \dot q + K \, q = f \f$
             *
             * The results are written into the TimestepMatrices object returned
             * by matrices(): the mass matrix M, the stiffness matrix K, the
             * load vector f, and — for Newton — the contracted derivative blocks
             * dMdX_times_x, dMdX_times_h, dKdX_times_x and dFdX from which
             * assemble_dJdx() builds the Newton correction.
             *
             * @param aElement element the matrices are computed for
             */
            void
            compute_mkf( Element * aElement) override;

//------------------------------------------------------------------------------

            void
            compute_timestep( Matrix< real > & aJ,
                              Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            set_field( DofManagerBase * aField ) override;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            // The startup ramp, in one place. BDF-p needs p history states, so
            // the first steps run the highest order the history supports.
            // Const and side-effect free, so it can be asked before the lazy
            // recompute has fired; compute_bdf_coefficients assigns its result
            // to mOrderActive, and order_pending() reports it to the log.
            inline uint
            ramped_order() const
            {
                switch( mMethod )
                {
                    case( EulerMethod::BackwardDifference2 ) :
                    case( EulerMethod::BackwardDifference3 ) :
                    case( EulerMethod::BackwardDifference4 ) :
                    case( EulerMethod::BackwardDifference5 ) :
                    {
                        return ( mOrder > 1 && mStepCount < mOrder )
                                ? ( mStepCount > 0 ? mStepCount : 1 )
                                : mOrder ;
                    }
                    default :
                    {
                        // every other scheme is single-state
                        return 1 ;
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            create_old_dof_fields();

//------------------------------------------------------------------------------

            void
            static_subfield( Matrix< real > & aM, Vector< real > & aRHS );

//------------------------------------------------------------------------------
            void
            explicit_euler( Matrix< real > & aM, Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            crank_nicolson( Matrix< real > & aM, Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            galerkin( Matrix< real > & aM, Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            mass_matrix_only(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );
//------------------------------------------------------------------------------

            void
            stiffness_matrix_only(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf1(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf2(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf3(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf4(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf5(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf1_nok(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf2_nok(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf3_nok(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf4_nok(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            bdf5_nok(
                Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            compute_bdf_coeffs_2();

//------------------------------------------------------------------------------

            void
            compute_bdf_coeffs_3();

//------------------------------------------------------------------------------

            void
            compute_bdf_coeffs_4();

//------------------------------------------------------------------------------

            void
            compute_bdf_coeffs_5();

//------------------------------------------------------------------------------

            void
            derivative( Matrix< real > & aJ,
                Vector< real > & aRHS  );

//------------------------------------------------------------------------------

            void
            euler_no_stiffness( Matrix< real > & aJ,
                                  Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            derivative_no_stiffness( Matrix< real > & aJ,
                                  Vector< real > & aRHS );


//------------------------------------------------------------------------------

            uint
            timestepping_order() const override ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        // return the timestepping method
        inline EulerMethod
        IWG_Timestep::method() const
        {
            return mMethod ;
        }

//------------------------------------------------------------------------------

        // BDF order actually running this step ( startup ramp may sit below
        // the configured order ); the controller's step-growth clamp reads it
        inline uint
        IWG_Timestep::order_active() const
        {
            return mOrderActive ;
        }

//------------------------------------------------------------------------------

        // BDF order the UPCOMING assembly will run. The coefficients are
        // recomputed lazily, on the first element of a step, so between
        // shift_fields and that first assembly mOrderActive still holds the
        // PREVIOUS step's order. Anything reporting the order before assembly
        // -- the timestep header, above all -- must ask this instead of
        // order_active(), or it lags the ramp by one step.
        inline uint
        IWG_Timestep::order_pending() const
        {
            return mCoeffsDirty ? this->ramped_order() : mOrderActive ;
        }

//------------------------------------------------------------------------------

        inline uint
        IWG_Timestep::timestepping_order() const
        {
            return mOrder ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IWG_TIMESTEP_HPP
