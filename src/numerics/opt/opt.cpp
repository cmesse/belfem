//
// Minimal example for the opt::Optimizer interface.
// Minimizes the 2D Rosenbrock function with a derivative-free algorithm.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "assert.hpp"

#include "cl_Vector.hpp"
#include "cl_Objective.hpp"
#include "cl_Optimizer.hpp"
#include "en_Opt_Algorithm.hpp"
#include "en_Opt_Status.hpp"

using namespace belfem;
using namespace opt;

Communicator gComm;
Logger       gLog( 5 );

//------------------------------------------------------------------------------

// f(x,y) = (a - x)^2 + b*(y - x^2)^2, minimum 0 at (a, a^2)
class Rosenbrock : public Objective
{
    const real mA = 1.0;
    const real mB = 100.0;

//------------------------------------------------------------------------------
public:
//------------------------------------------------------------------------------

    Rosenbrock() :
        Objective( 2 )
    {
    }

//------------------------------------------------------------------------------

    real
    compute_objective(
            const Vector< real > & aX,
                  Vector< real > & aGradient ) override
    {
        const real tX = aX( 0 );
        const real tY = aX( 1 );

        // derivative-free algorithm: aGradient is empty and left untouched
        return ( mA - tX ) * ( mA - tX )
             + mB * ( tY - tX * tX ) * ( tY - tX * tX );
    }

//------------------------------------------------------------------------------
};

//------------------------------------------------------------------------------

int main( int argc, char * argv[] )
{
    gComm.init( argc, argv );

    // the objective
    Rosenbrock tObjective;

    // the optimizer
    Optimizer tOptimizer( tObjective, Algorithm::BOBYQA );

    // search box
    Vector< real > tLower = { -5.0, -5.0 };
    Vector< real > tUpper = {  5.0,  5.0 };
    tOptimizer.set_bounds( tLower, tUpper );

    tOptimizer.xtol_rel() = 1.0e-10;
    tOptimizer.max_eval() = 10000;

    // initial guess
    Vector< real > tX = { -1.2, 1.0 };
    real tValue = 0.0;

    Status tStatus = tOptimizer.optimize( tX, tValue );

    // intended failure-handling pattern: abort only when the result is not
    // usable ( ROUNDOFF_LIMITED under tight tolerances is a normal BOBYQA
    // exit and still delivers the optimum ), and combine the generic status
    // explanation with nlopt's detailed diagnostic
    BELFEM_ERROR( is_usable( tStatus ), "%s\n%s",
                  error_message( tStatus ).c_str(),
                  tOptimizer.errmsg().c_str() );

    std::cout << "status   : " << ( uint ) tStatus << " ( "
              << error_message( tStatus ) << " )" << std::endl;
    std::cout << "optimum  : ( " << tX( 0 ) << " , " << tX( 1 ) << " )" << std::endl;
    std::cout << "value    : " << tValue << std::endl;

    return gComm.finalize();
}
