//
// Created by Christian Messe on 2019-01-19.
//

#include "fn_Create_Glue_Poly.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    void
    create_glue_poly(
            const real          & aX,
            const real          & aDeltaX,
            const Vector <real> & aF,
            Vector <real>       & aC )
    {
        // get boundary conditions on the left side
        const real & f0  = aF( 0 );
        const real & df0 = aF( 1 );
        const real & f1  = aF( 2 );
        const real & f2  = aF( 3 );
        const real & df2 = aF( 4 );

        // get some shortcuts
        const real & X  = aX;
        const real   X2 = std::pow( aX, 2 );
        const real   X3 = std::pow( aX, 3 );
        const real   X4 = std::pow( aX, 4 );

        const real & dX  = aDeltaX;
        const real   dX2 = std::pow( aDeltaX, 2 );
        const real   dX3 = std::pow( aDeltaX, 3 );
        const real   dX4 = std::pow( aDeltaX, 4 );

        const real b = (f0-2.0*f1+f2);

        // populate data
        aC.set_size( 5 );

        aC( 0 ) = -(df0*dX)+df2*dX-2.0*b;

        aC( 1 ) = dX*((df0+df2)*dX+f0-f2)+4.0*(df0*dX-df2*dX+2.0*b)*X;

        aC( 2 ) =      dX2*(df0*dX-df2*dX+4.0*b)
                      -3.0*dX*((df0+df2)*dX+f0-f2)*X
                      -6.0*(df0*dX-df2*dX+2*b)*X2;

        aC( 3 ) = -(df0*dX*(dX+X)*(dX2+dX*X-4.0*X2))+(dX-X)
                   *(-((dX+X)*(3.0*dX*(f0-f2)+8.0*b*X))
                   +df2*dX*(-dX2+dX*X+4.0*X2));

        aC( 4 ) = -2.0*b*X4+dX3*X*(3.0*f0-3.0*f2+df0*X-df2*X)
                 +dX*X3*(-f0+f2-df0*X+df2*X)
                 +dX4*(4.0*f1+(df0+df2)*X)-dX2*X2*(-4.0*b+(df0+df2)*X);

        aC /= 4.0 * dX4;
    }

//------------------------------------------------------------------------------
}