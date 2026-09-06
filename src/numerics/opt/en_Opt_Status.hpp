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

#ifndef BELFEM_EN_OPT_STATUS_HPP
#define BELFEM_EN_OPT_STATUS_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace opt
    {
        /**
         * outcome of an optimization run. The success codes mirror the
         * termination reasons reported by the underlying library so the caller
         * can tell why the solver stopped without depending on nlopt headers.
         */
        enum class Status
        {
            // --- success: the optimizer converged on a stopping criterion ---
            SUCCESS          = 0, // generic success
            STOPVAL_REACHED  = 1, // objective dropped below the stop value
            FTOL_REACHED     = 2, // relative change of objective below ftol
            XTOL_REACHED     = 3, // relative change of design vector below xtol
            MAXEVAL_REACHED  = 4, // reached the maximum number of evaluations
            MAXTIME_REACHED  = 5, // reached the maximum wall-clock time

            // --- failure: no usable result was produced -------------------
            FAILURE          = 6, // generic failure
            INVALID_ARGS     = 7, // bounds/tolerances inconsistent
            OUT_OF_MEMORY    = 8, // allocation failed
            ROUNDOFF_LIMITED = 9, // halted by roundoff, result may be inaccurate
            FORCED_STOP      = 10 // objective requested an early stop
        };

//------------------------------------------------------------------------------

        /**
         * true if the status corresponds to a converged / usable result
         */
        inline bool
        is_success( const Status aStatus )
        {
            return aStatus <= Status::MAXTIME_REACHED;
        }

//------------------------------------------------------------------------------

        /**
         * true if the returned design point is usable, even when the solver
         * did not terminate on a clean convergence criterion. This is the
         * MUMPS-style error/warning distinction: ROUNDOFF_LIMITED is the
         * routine "converged as far as floating point allows" exit of the
         * Powell-type algorithms ( BOBYQA/COBYLA ) under tight tolerances,
         * and FORCED_STOP returns the best point found before the stop.
         * Guard fatal aborts with this; use is_success() to detect a clean
         * convergence.
         */
        inline bool
        is_usable( const Status aStatus )
        {
            return is_success( aStatus )
                || aStatus == Status::ROUNDOFF_LIMITED
                || aStatus == Status::FORCED_STOP;
        }

//------------------------------------------------------------------------------

        /**
         * human-readable explanation of an optimizer status, suitable for a
         * BELFEM_ERROR message when is_usable() is false. Defined in
         * cl_Optimizer.cpp. Analogous to MUMPS::error_message().
         */
        string
        error_message( const Status aStatus );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_EN_OPT_STATUS_HPP
