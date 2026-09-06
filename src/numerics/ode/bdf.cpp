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

#include "typedefs.hpp"
#include "constants.hpp"

#include "cl_Communicator.hpp"
#include "commtools.hpp"

#include "cl_Logger.hpp"
#include "cl_Vector.hpp"


#include "cl_ShiftRegister.hpp"
#include "cl_BDF.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );


int main( int    argc,
          char * argv[] )
{
    // the stepwith
    ShiftRegister< real > tH( 5 );

    // the values from the last timesteps
    ShiftRegister< real > tY( tH.capacity() + 1 );

    ode::BDF tSolver( tH );

    real tStep = 1.0 ;

    // time value
    real tT = 0.0 ;

    // initial value
    tY.push( 0.0 );

    // constant
    real tC = 0.05 * constant::pi ;

    for ( size_t  k=0; k<10; ++k )
    {
        tH.push( tStep );

        // compute time value
        tT += tH( 0 );

        // compute derivative
        real tF = tC * std::cos( tT * tC );

        real tY_guess = tSolver.eval( tY, tF );

        // tSolver.deval( tY, false );

        real tY_exact = std::sin( tT * tC );

        // example function
        tY.push( tY_guess );

        std::cout << "#ORDER " << tH.size() << std::endl ;

        tSolver.coefficients().print( "Coefficients" );

        std::cout << std::endl ;

        std::cout << "guess : " << tY_guess << " " << "exact: " << tY_exact << std::endl ;

        tStep *= 0.95 ; // <-- uncomment to test with variable stepwidth
    }

}