//
// Created by christian on 10/17/22.
//

#ifndef BELFEM_CL_BHEXTRAPOLATIONODE_HPP
#define BELFEM_CL_BHEXTRAPOLATIONODE_HPP

#include "cl_ODE.hpp"
#include "cl_Vector.hpp"
namespace belfem
{
    class BhExtrapolationOde : public ode::ODE
    {
        // coefficients for polynomial
        Vector< real > mCoeffs ;

        // value where we go into a straight line
        real mBswitch ;

//-----------------------------------------------------------------------------
    public:
//-----------------------------------------------------------------------------

        BhExtrapolationOde( const real aBref,
               const real adHdBref,
               const real ad2HdB2ref,
               const real aBmax );

        ~BhExtrapolationOde() = default ;

//-----------------------------------------------------------------------------

        void
        compute(
                const real           & aT,
                const Vector< real > & aY,
                Vector< real > & adYdT ) ;

    };
}
#endif //BELFEM_CL_BHEXTRAPOLATIONODE_HPP
