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

#ifndef BELFEM_CL_OPTIMIZER_HPP
#define BELFEM_CL_OPTIMIZER_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Objective.hpp"
#include "en_Opt_Algorithm.hpp"
#include "en_Opt_Status.hpp"

namespace belfem
{
    namespace opt
    {
//------------------------------------------------------------------------------

        /**
         * driver that optimizes an Objective using a chosen Algorithm.
         *
         * This is the counterpart of ode::Integrator: it owns the numerical
         * bookkeeping (bounds, tolerances, evaluation budget) while remaining
         * agnostic of what is being optimized. The concrete solver library
         * (nlopt) is a private implementation detail of cl_Optimizer.cpp.
         *
         * The Objective is borrowed by reference (non-owning) and must outlive
         * every call to optimize(), just as ode::Integrator borrows its ODE.
         */
        class Optimizer
        {
            // objective to be minimized or maximized
            Objective & mObjective;

            // algorithm used for this optimizer
            Algorithm mAlgorithm;

            // lower / upper bounds on the design variables. Set together via
            // set_bounds(); both left empty means an unbounded problem.
            Vector< real > mLowerBounds;
            Vector< real > mUpperBounds;

            // relative tolerance on the objective value
            real mFtolRel = 1.0e-8;

            // relative tolerance on the design vector
            real mXtolRel = 1.0e-8;

            // maximum number of objective evaluations
            uint mMaxEval = 1000;

            // if true, the objective is maximized instead of minimized
            bool mMaximize = false;

            // detailed diagnostic from the solver library for the most recent
            // optimize() call; empty when no detail was reported
            string mErrmsg;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Optimizer( Objective & aObjective, const Algorithm aAlgorithm );

//------------------------------------------------------------------------------

            ~Optimizer() = default;

//------------------------------------------------------------------------------

            /**
             * set the lower and upper bounds of the design variables.
             * Each vector must have length objective.dimension().
             *
             * Setting aLowerBounds( i ) == aUpperBounds( i ) fixes design
             * variable i: nlopt eliminates that dimension internally, so
             * parameters can be pinned without changing the Objective.
             */
            void
            set_bounds(
                    const Vector< real > & aLowerBounds,
                    const Vector< real > & aUpperBounds );

//------------------------------------------------------------------------------

            /**
             * run the optimization, starting from and writing back into aX.
             *
             * @param aX     on entry the initial guess, on exit the optimum
             * @param aValue on exit the objective value at aX
             * @return       the termination status
             *
             * aX and aValue are only meaningful when is_usable( status ) is
             * true; on failure they must not be trusted. When the solver
             * reports a detailed diagnostic, it is available via errmsg().
             */
            Status
            optimize( Vector< real > & aX, real & aValue );

//------------------------------------------------------------------------------

            /**
             * detailed diagnostic from the solver for the most recent
             * optimize() call ( e.g. which bound fails lb <= ub ); empty
             * when no detail was reported. Complements the generic
             * error_message( status ).
             */
            inline const string &
            errmsg() const;

//------------------------------------------------------------------------------

            /**
             * relative tolerance on the objective value ( ftol_rel )
             */
            inline real &
            ftol_rel();

//------------------------------------------------------------------------------

            /**
             * relative tolerance on the design vector ( xtol_rel )
             */
            inline real &
            xtol_rel();

//------------------------------------------------------------------------------

            /**
             * maximum number of objective evaluations; 0 means unlimited
             */
            inline uint &
            max_eval();

//------------------------------------------------------------------------------

            /**
             * flag selecting maximization ( true ) over minimization ( false )
             */
            inline bool &
            maximize();

//------------------------------------------------------------------------------

            /**
             * the algorithm used by this optimizer
             */
            inline const Algorithm &
            algorithm() const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline real &
        Optimizer::ftol_rel()
        {
            return mFtolRel;
        }

//------------------------------------------------------------------------------

        inline real &
        Optimizer::xtol_rel()
        {
            return mXtolRel;
        }

//------------------------------------------------------------------------------

        inline uint &
        Optimizer::max_eval()
        {
            return mMaxEval;
        }

//------------------------------------------------------------------------------

        inline bool &
        Optimizer::maximize()
        {
            return mMaximize;
        }

//------------------------------------------------------------------------------

        inline const Algorithm &
        Optimizer::algorithm() const
        {
            return mAlgorithm;
        }

//------------------------------------------------------------------------------

        inline const string &
        Optimizer::errmsg() const
        {
            return mErrmsg;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_OPTIMIZER_HPP
