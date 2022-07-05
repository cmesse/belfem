//
// Created by christian on 4/21/22.
//

#ifndef BELFEM_NIST_FUNCTIONS_HPP
#define BELFEM_NIST_FUNCTIONS_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"

namespace belfem
{
    namespace nist
    {
//----------------------------------------------------------------------------

        /**
         * this is the specific heat capacity polynomial as used in NASA CEA
         */
        real
        cp_janaf( const Vector< real > & aC, const real aT );

//----------------------------------------------------------------------------

        /**
         * this is the derivative of the specific heat capacity polynomial
         */
        real
        dcpdT_janaf( const Vector< real > & aC, const real aT );

//----------------------------------------------------------------------------

        /**
         * a general polynomial used for specific heat and thermal conductivity
         * @param aCoeffs
         * @param aT
         * @return
         */
        real
        proppoly( const Vector< real > & aC, const real aT  );

//----------------------------------------------------------------------------

        // derivative of nist poly
        real
        dproppoly( const Vector< real > & aC, const real aT );

//----------------------------------------------------------------------------

        real
        rho0( const Vector< real > & aP, const real aRRR, const real aT );

//----------------------------------------------------------------------------

        /**
         * a function for thermal conductivity under zero electric field
         * @param aT
         * @return
         */
         real
         lambda0( const Vector< real > & aP, const real aBeta, const real aT );

//----------------------------------------------------------------------------

        /**
         * This is the derivative of the inverse of lambda0.
         * This function is needed to find the temperature where lambda peaks
         * @param aT
         * @return
         */
        real
        dinvlambda0dT( const Vector< real > & aP, const real aBeta, const real aT );

//----------------------------------------------------------------------------

        void
        extend_resistivity( const Vector< real > & aPoly1,
                            const         real     aTswitch,
                                  Vector< real > & aPoly0 );

//----------------------------------------------------------------------------

         /**
          * Kohler's law needs a cutoff where it is replaced with a
          * cubic polynomial
          *
          * the electric resistivity is corrected so that
          * rho + rho_b0 * (1 + K),
          *
          * where K=10^polyval(A,log10(x))
          *
          * and x = rho_ref/rho_b0 * | b |
          *
          * at the curoff, we replace K with
          * K0 = polyval(B,x)
          *
          * @param aA       coefficients for Kohler correction, let k = 10^
          * @param aKmin    cutoff where polynomial is used
          * @param aXmin0   initial guess of value where cutoff happens
          * @param aB       coefficients for help polynomial
          *
          *
          */
         real
         extend_kohler( const Vector< real >  & aA,
                     const real              aKmin,
                     const real              aXmin0,
                           Vector< real >  & aB );

//----------------------------------------------------------------------------

        inline real
        kohler( const Vector< real >  & aA,
                const Vector< real >  & aB,
                const         real      aXmin,
                const         real      aX )
        {
            return aX < aXmin ?
                   polyval( aB, aX ) :
                   std::pow( 10, polyval( aA, std::log10( aX ) ) );
        }

//----------------------------------------------------------------------------
    }
}

#endif //BELFEM_NIST_FUNCTIONS_HPP
