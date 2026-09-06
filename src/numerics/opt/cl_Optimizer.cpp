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

#include "cl_Optimizer.hpp"
#include "assert.hpp"
#include "fn_sprint.hpp"

#ifdef BELFEM_NLOPT
#include <nlopt.h>
#include <type_traits>
#endif

namespace belfem
{
    namespace opt
    {
//------------------------------------------------------------------------------

        Optimizer::Optimizer( Objective & aObjective, const Algorithm aAlgorithm ) :
            mObjective( aObjective ),
            mAlgorithm( aAlgorithm )
        {
        }

//------------------------------------------------------------------------------

        void
        Optimizer::set_bounds(
                const Vector< real > & aLowerBounds,
                const Vector< real > & aUpperBounds )
        {
            BELFEM_ERROR( aLowerBounds.length() == mObjective.dimension(),
                          "lower bound vector has length %u but objective dimension is %u",
                          ( unsigned int ) aLowerBounds.length(),
                          ( unsigned int ) mObjective.dimension() );

            BELFEM_ERROR( aUpperBounds.length() == mObjective.dimension(),
                          "upper bound vector has length %u but objective dimension is %u",
                          ( unsigned int ) aUpperBounds.length(),
                          ( unsigned int ) mObjective.dimension() );

            mLowerBounds = aLowerBounds;
            mUpperBounds = aUpperBounds;
        }

//------------------------------------------------------------------------------

        // Maps the library-agnostic Status onto an explanatory sentence. The
        // wording follows the nlopt result-code documentation (nlopt-2.7.1,
        // doc/docs/NLopt_Reference.md, "Return values") but is phrased in
        // BELFEM terms and carries no nlopt symbols, so it is compiled whether
        // or not NLOPT support is present.
        string
        error_message( const Status aStatus )
        {
            switch ( aStatus )
            {
                case ( Status::SUCCESS ) :
                    return "Optimization converged successfully.";

                case ( Status::STOPVAL_REACHED ) :
                    return "Optimization stopped because the target objective value "
                           "( stopval ) was reached.";

                case ( Status::FTOL_REACHED ) :
                    return "Optimization converged: the relative change of the objective "
                           "value fell below the ftol tolerance.";

                case ( Status::XTOL_REACHED ) :
                    return "Optimization converged: the relative change of the design "
                           "vector fell below the xtol tolerance.";

                case ( Status::MAXEVAL_REACHED ) :
                    return "Optimization stopped because the maximum number of objective "
                           "evaluations ( maxeval ) was reached; the result may not be "
                           "converged.";

                case ( Status::MAXTIME_REACHED ) :
                    return "Optimization stopped because the maximum allowed run time "
                           "( maxtime ) was reached; the result may not be converged.";

                case ( Status::FAILURE ) :
                    return "Generic optimization failure.";

                case ( Status::INVALID_ARGS ) :
                    return "Invalid arguments: a lower bound may exceed its upper bound, "
                           "the initial guess may lie outside the bounds, or the selected "
                           "algorithm may be incompatible with the problem.";

                case ( Status::OUT_OF_MEMORY ) :
                    return "The optimizer ran out of memory.";

                case ( Status::ROUNDOFF_LIMITED ) :
                    return "Optimization halted because roundoff errors limited further "
                           "progress; the returned result is typically still usable.";

                case ( Status::FORCED_STOP ) :
                    return "Optimization was forced to stop from within the objective "
                           "function ( forced termination ).";

                default :
                    return sprint( "Unrecognized optimizer status code %i.",
                                   ( int ) aStatus );
            }
        }

//------------------------------------------------------------------------------
#ifdef BELFEM_NLOPT
//------------------------------------------------------------------------------

        // the bounds and design vectors are handed to nlopt through
        // Vector<real>::data(), so real must be binary-compatible with the
        // double* that the nlopt C API expects.
        static_assert( std::is_same< real, double >::value,
                       "opt::Optimizer passes Vector<real>::data() directly to the "
                       "nlopt C API, which requires real == double" );

//------------------------------------------------------------------------------

        // translate the library-agnostic selector into an nlopt algorithm.
        // kept local so the enum mapping is the only place that knows nlopt.
        static nlopt_algorithm
        to_nlopt_algorithm( const Algorithm aAlgorithm )
        {
            switch ( aAlgorithm )
            {
                case ( Algorithm::BOBYQA )     : return NLOPT_LN_BOBYQA ;
                case ( Algorithm::COBYLA )     : return NLOPT_LN_COBYLA ;
                case ( Algorithm::NELDERMEAD ) : return NLOPT_LN_NELDERMEAD ;
                case ( Algorithm::SBPLX )      : return NLOPT_LN_SBPLX ;
                case ( Algorithm::PRAXIS )     : return NLOPT_LN_PRAXIS ;
                case ( Algorithm::MMA )        : return NLOPT_LD_MMA ;
                case ( Algorithm::SLSQP )      : return NLOPT_LD_SLSQP ;
                case ( Algorithm::LBFGS )      : return NLOPT_LD_LBFGS ;
                default :
                {
                    BELFEM_ERROR( false, "unknown optimization algorithm" );
                    return NLOPT_LN_BOBYQA ;
                }
            }
        }

//------------------------------------------------------------------------------

        // translate an nlopt result into the BELFEM status enum
        static Status
        to_status( const nlopt_result aResult )
        {
            switch ( aResult )
            {
                case ( NLOPT_SUCCESS )          : return Status::SUCCESS ;
                case ( NLOPT_STOPVAL_REACHED )  : return Status::STOPVAL_REACHED ;
                case ( NLOPT_FTOL_REACHED )     : return Status::FTOL_REACHED ;
                case ( NLOPT_XTOL_REACHED )     : return Status::XTOL_REACHED ;
                case ( NLOPT_MAXEVAL_REACHED )  : return Status::MAXEVAL_REACHED ;
                case ( NLOPT_MAXTIME_REACHED )  : return Status::MAXTIME_REACHED ;
                case ( NLOPT_INVALID_ARGS )     : return Status::INVALID_ARGS ;
                case ( NLOPT_OUT_OF_MEMORY )    : return Status::OUT_OF_MEMORY ;
                case ( NLOPT_ROUNDOFF_LIMITED ) : return Status::ROUNDOFF_LIMITED ;
                case ( NLOPT_FORCED_STOP )      : return Status::FORCED_STOP ;
                default :                         return Status::FAILURE ;
            }
        }

//------------------------------------------------------------------------------

        // bundles the objective with preallocated work buffers so the
        // trampoline performs no heap allocation per objective evaluation
        // ( doc/coding_philosophy.md: no hidden allocations on hot paths ).
        struct EvalContext
        {
            Objective &    mObjective;   // borrowed, non-owning
            Vector< real > mX;           // design point buffer, length dimension()
            Vector< real > mGradient;    // gradient buffer, length dimension()
            Vector< real > mEmpty;       // length 0, passed for derivative-free algorithms
        };

//------------------------------------------------------------------------------

        // C-style trampoline matching nlopt_func. It recovers the EvalContext
        // from the void* payload and forwards the call to the virtual method,
        // mirroring the function-pointer indirection in ode::Integrator.
        static double
        objective_trampoline(
                unsigned       aN,
                const double * aX,
                double       * aGradient,
                void         * aData )
        {
            EvalContext * tContext = static_cast< EvalContext * >( aData );

            BELFEM_ASSERT( tContext->mX.length() == aN,
                           "nlopt requested dimension %u but the context is sized for %u",
                           aN, ( unsigned int ) tContext->mX.length() );

            // copy the incoming design point into the preallocated buffer
            for ( unsigned int i = 0; i < aN; ++i )
            {
                tContext->mX( i ) = aX[ i ];
            }

            // derivative-free algorithms ( aGradient == nullptr ) receive an
            // empty vector; gradient-based ones the preallocated buffer
            real tValue = tContext->mObjective.compute_objective(
                    tContext->mX,
                    aGradient == nullptr ? tContext->mEmpty : tContext->mGradient );

            if ( aGradient != nullptr )
            {
                BELFEM_ASSERT( tContext->mGradient.length() == aN,
                               "objective resized the gradient vector to %u, expected %u",
                               ( unsigned int ) tContext->mGradient.length(), aN );

                for ( unsigned int i = 0; i < aN; ++i )
                {
                    aGradient[ i ] = tContext->mGradient( i );
                }
            }

            return tValue ;
        }

//------------------------------------------------------------------------------

        Status
        Optimizer::optimize( Vector< real > & aX, real & aValue )
        {
            const uint tN = mObjective.dimension();

            BELFEM_ERROR( aX.length() == tN,
                          "design vector has length %u but objective dimension is %u",
                          ( unsigned int ) aX.length(), ( unsigned int ) tN );

            // reset the diagnostic of the previous run
            mErrmsg.clear();

            // create the solver instance
            nlopt_opt tOpt = nlopt_create( to_nlopt_algorithm( mAlgorithm ), tN );
            BELFEM_ERROR( tOpt != nullptr, "failed to create nlopt optimizer" );

            // apply bounds if they were provided
            if ( mLowerBounds.length() == tN )
            {
                nlopt_set_lower_bounds( tOpt, mLowerBounds.data() );
            }
            if ( mUpperBounds.length() == tN )
            {
                nlopt_set_upper_bounds( tOpt, mUpperBounds.data() );
            }

            // preallocated evaluation context handed to the trampoline as
            // f_data; nlopt stores this pointer, so tContext must outlive the
            // nlopt_optimize call below ( it does — same scope ).
            EvalContext tContext {
                    mObjective,
                    Vector< real >( tN ),
                    Vector< real >( tN, 0.0 ),
                    Vector< real >() };

            // register the objective through the trampoline
            if ( mMaximize )
            {
                nlopt_set_max_objective( tOpt, & objective_trampoline, & tContext );
            }
            else
            {
                nlopt_set_min_objective( tOpt, & objective_trampoline, & tContext );
            }

            // stopping criteria
            nlopt_set_ftol_rel( tOpt, mFtolRel );
            nlopt_set_xtol_rel( tOpt, mXtolRel );
            nlopt_set_maxeval( tOpt, ( int ) mMaxEval );

            // run the optimization
            double tValue = 0.0;
            nlopt_result tResult = nlopt_optimize( tOpt, aX.data(), & tValue );
            aValue = tValue;

            // salvage nlopt's detailed diagnostic ( e.g. which bound fails
            // lb <= ub ) before the handle is destroyed; NULL when none
            const char * tErrmsg = nlopt_get_errmsg( tOpt );
            if ( tErrmsg != nullptr )
            {
                mErrmsg = tErrmsg;
            }

            nlopt_destroy( tOpt );

            return to_status( tResult );
        }

//------------------------------------------------------------------------------
#else  // BELFEM_NLOPT
//------------------------------------------------------------------------------

        Status
        Optimizer::optimize( Vector< real > & aX, real & aValue )
        {
            BELFEM_ERROR( false,
                          "BELFEM was compiled without NLOPT support. "
                          "Reconfigure with -DUSE_NLOPT=ON to use opt::Optimizer." );

            return Status::FAILURE;
        }

//------------------------------------------------------------------------------
#endif // BELFEM_NLOPT
//------------------------------------------------------------------------------
    }
}
