//
// Created by Christian Messe on 26.08.26.
//

#include <gtest/gtest.h>
#include <cmath>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "constants.hpp"
#include "cl_Gas.hpp"

#include "cl_GM_EoS_Nitrogen.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

/*
 * Span, Lemmon, Jacobsen and Wagner, J. Phys. Chem. Ref. Data 29( 6 ):1361
 * ( 2000 ). Critical point from Table 1, ancillary equations from Sec. 3.
 *
 * These checks exist because the derivative test cannot see an error in the
 * coefficient tables. A single wrong exponent in the residual table -- j of
 * term 20 -- once left the vapor pressure, the critical point and the
 * saturated vapor density all correct while the saturated LIQUID density was
 * 0.5 % too dense. The liquid line is what caught it, so it is tested here.
 */
TEST( GASMODELS, Nitrogen_Vapor )
{
    Gas tRef( "N2" );

    EoS_Nitrogen tGas( tRef );

    // Table 1
    const real tTcrit   = 126.192 ;
    const real tPcrit   = 3.3958e6 ;
    const real tRhocrit = 313.300 ;

//----------------------------------------------------------------------------
// vapor pressure ancillary
//----------------------------------------------------------------------------

    /*
     * Two states the ancillary does not fit by construction. The normal
     * boiling point and the triple point are therefore genuine checks of it,
     * and both are highly sensitive to the sign of the third coefficient.
     */

    // normal boiling point, 101325 Pa at 77.355 K
    EXPECT_NEAR( tGas.p_vap( 77.355 ) / 101325.0, 1.0, 1e-4 );

    // triple point, 12520 Pa at 63.151 K
    EXPECT_NEAR( tGas.p_vap( 63.151 ) / 12520.0, 1.0, 1e-3 );

    // and the inverse
    EXPECT_NEAR( tGas.T_vap( 101325.0 ), 77.355, 1e-2 );

//----------------------------------------------------------------------------
// critical point
//----------------------------------------------------------------------------

    /*
     * The equation of state must return the critical pressure at the critical
     * temperature and volume. This ties rho_crit to the coefficient tables:
     * the paper prints the critical density in molar form, 11.1839 mol/dm^3,
     * and an error in converting it to mass units shows up right here.
     */
    EXPECT_NEAR( tGas.p( tTcrit, 1.0 / tRhocrit ) / tPcrit, 1.0, 1e-4 );

//----------------------------------------------------------------------------
// saturated liquid density
//----------------------------------------------------------------------------

    /*
     * Against the paper's own saturated liquid density ancillary,
     *
     *   rho' / rho_c = exp( sum_i N_i * theta^k_i ) ,  theta = 1 - T / Tc
     *
     * evaluated here rather than taken from a table, so the reference is the
     * paper and not a transcription of it.
     */
    const Vector< real > tN = {  1.48654237, -0.280476066,
                                 0.0894143085, -0.119879866 } ;

    const Vector< real > tK = {  0.3294, 2.0 / 3.0, 8.0 / 3.0, 35.0 / 6.0 } ;

    const Vector< real > tT = { 63.151, 70.0, 77.355, 90.0, 100.0, 110.0 } ;

    for( uint k = 0; k < tT.length(); ++k )
    {
        const real T     = tT( k ) ;
        const real theta = 1.0 - T / tTcrit ;

        real tExponent = 0.0 ;

        for( uint i = 0; i < 4; ++i )
        {
            tExponent += tN( i ) * std::pow( theta, tK( i ) ) ;
        }

        const real tRhoRef = tRhocrit * std::exp( tExponent ) ;

        // just inside the liquid branch of the saturation line
        const real tRhoEoS = 1.0 / tGas.v( T, 1.0001 * tGas.p_vap( T ) ) ;

        // the equation of state is quoted to 0.02 % here; allow 0.1 %
        EXPECT_NEAR( tRhoEoS / tRhoRef, 1.0, 1e-3 );
    }

//----------------------------------------------------------------------------
// ideal gas limit
//----------------------------------------------------------------------------

    const real tR = constant::Rm / tRef.data( 0 )->M() ;

    for( real p : { 1.0e2, 1.0e3, 1.0e4 } )
    {
        const real v = tGas.v( 300.0, p );

        EXPECT_NEAR( p * v / ( tR * 300.0 ), 1.0, 1e-4 );
    }
}
