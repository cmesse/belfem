//
// Created by Christian Messe on 28.12.19.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"

#include "cl_Vector.hpp"
#include "cl_ODE.hpp"
#include "cl_ODE_Integrator.hpp"
#include "en_ODE_Type.hpp"
#include "en_ODE_Status.hpp"

using namespace belfem;
using namespace ode;

Communicator gComm;
Logger       gLog( 5 );

//------------------------------------------------------------------------------

class MyODE : public ODE
{
//------------------------------------------------------------------------------
public:
//------------------------------------------------------------------------------

    MyODE() :
        ODE( 1 )
    {
    }

//------------------------------------------------------------------------------
    void
    compute(
            const real           & aT,
            const Vector< real > & aY,
                  Vector< real > & adYdT )
    {
        adYdT( 0 ) = 3.0 / ( 2.0 * aT * aT ) + aY( 0 ) / ( 2.0 * aT );
    }

//------------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    // create the state
    Vector< real > tY = { 0.0 };

    // the derivative vector
    Vector< real > tdYdT = { 0.0 };

    // create the ODE
    MyODE tMyODE;

    // create the integrator
    Integrator tIntegrator( tMyODE, Type::RK45 );

    real tTime = 1.0;

    tIntegrator.timestep() = 10.0;
    tIntegrator.maxtime() = 0.930 ;

    while ( ( uint ) tIntegrator.step( tTime, tY ) < 2 )
    {
        std::cout << tTime << " " << tY( 0 ) << std::endl;
    }
    return gComm.finalize();
}