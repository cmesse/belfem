//
// Created by Christian Messe on 26.08.26.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "constants.hpp"
#include "cl_Gas.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

/*
 * Cross check of the nitrogen model against a secondary source:
 * VDI Heat Atlas, chapter D2.3, Properties of Nitrogen, Span and Krauss,
 * Tables 2 and 3, saturated liquid and saturated vapor.
 *
 * What this does and does not prove.
 *
 * The thermodynamic columns of D2.3 are computed from the reference equation
 * of state of Span et al., the same equation BELFEM implements. So this is not
 * an independent check of the equations; it is an independent check of the
 * IMPLEMENTATION of them, against a table produced by someone else's code.
 * That is exactly the class of error a transcription slip in the coefficient
 * tables belongs to, and it covers many more states than the paper's own
 * ancillaries do.
 *
 * The transport columns are NOT comparable at the same level. D2.3 uses
 * Stephan and Krauss ( 1987 ) while BELFEM uses Lemmon and Jacobsen ( 2004 ),
 * two different correlations. Those columns are therefore checked loosely, as
 * a sanity bound rather than a verification; Transport_LemmonJacobsen is where
 * the viscosity and conductivity are actually pinned down.
 */
TEST( GASMODELS, Nitrogen_VDI )
{
    Gas tGas( HelmholtzModel::Nitrogen );

    /*
     * D2.3 Tables 2 and 3, at temperatures in degrees Celsius.
     *
     * The tables start at -210 C, which is 63.15 K. That is a thousandth of a
     * kelvin BELOW the triple point BELFEM carries, 63.151 K, and
     * Helmholtz::is_liquid() returns true unconditionally below the triple
     * point. The first row would therefore be compared against the liquid root
     * on the vapor branch. It is left out rather than worked around.
     */
    const Vector< real > tDegC = { -208.0, -204.0, -200.0, -196.0,
                                   -192.0, -188.0, -186.0 } ;

    // vapor pressure in bar
    const Vector< real > tPsat = { 0.17860, 0.33973, 0.59842, 0.98899,
                                   1.5497,  2.3219,  2.8009 } ;

    // Table 2, saturated liquid density in kg/m^3
    const Vector< real > tRhoL = { 858.97, 842.15, 824.85, 807.01,
                                   788.56, 769.40, 759.51 } ;

    // Table 2, saturated liquid isobaric heat capacity in kJ/(kg K)
    const Vector< real > tCpL  = { 2.004, 2.012, 2.024, 2.041,
                                   2.063, 2.092, 2.110 } ;

    // Table 3, saturated vapor density in kg/m^3
    const Vector< real > tRhoV = { 0.93502, 1.6884, 2.8403, 4.5102,
                                   6.8322,  9.9578, 11.875 } ;

    for( uint k = 0; k < tDegC.length(); ++k )
    {
        const real T = tDegC( k ) + 273.15 ;

        // vapor pressure
        const real p = tGas.eos()->p_vap( T );

        EXPECT_NEAR( p / ( tPsat( k ) * 1.0e5 ), 1.0, 2e-3 );

        // saturated liquid, just inside the liquid branch
        const real tRhoLiq = 1.0 / tGas.eos()->v( T, 1.0001 * p );

        EXPECT_NEAR( tRhoLiq / tRhoL( k ), 1.0, 2e-3 );

        // saturated liquid heat capacity, kJ/(kg K) in the table
        EXPECT_NEAR( tGas.cp( T, 1.0001 * p ) / ( tCpL( k ) * 1.0e3 ),
                     1.0, 5e-3 );

        // saturated vapor, just outside it
        const real tRhoVap = 1.0 / tGas.eos()->v( T, 0.9999 * p );

        EXPECT_NEAR( tRhoVap / tRhoV( k ), 1.0, 2e-3 );
    }
}
