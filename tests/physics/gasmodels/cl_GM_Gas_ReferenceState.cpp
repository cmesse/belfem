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
 * The Helmholtz fluids and the CEA ideal gas model must sit on the SAME
 * enthalpy and entropy scale.
 *
 * Helmholtz::set_reference_point() forces h and s of the real fluid to the CEA
 * standard state values at 298.15 K and 1 bar, so agreement there is by
 * construction and this test only confirms the wiring. Agreement at any OTHER
 * temperature is not by construction: it holds only as far as the ideal gas
 * part of the Helmholtz energy and the CEA heat capacity polynomial agree, so
 * 273.15 K is the check with something to say.
 *
 * Methane used to be anchored on its normal boiling point instead, which put it
 * alone on a different scale -- its entropy differed from the CEA model by 43 %.
 * This test exists so that cannot come back unnoticed.
 */
TEST( GASMODELS, Gas_ReferenceState )
{
    const Cell< string > tSpecies = { "H2", "O2", "N2", "CH4" } ;

    const real p = gastables::gPref ;

    for( uint k = 0; k < tSpecies.size(); ++k )
    {
        Gas tHelmholtz( tSpecies( k ), GasModel::HELMHOLTZ );
        Gas tIdealGas ( tSpecies( k ), GasModel::IDGAS );

//----------------------------------------------------------------------------
// at the reference state itself
//----------------------------------------------------------------------------

        const real T0 = gastables::gTref ;

        EXPECT_NEAR( tHelmholtz.s( T0, p ) / tIdealGas.s( T0, p ), 1.0, 1e-6 );
        EXPECT_NEAR( tHelmholtz.h( T0, p ) / tIdealGas.h( T0, p ), 1.0, 1e-6 );

//----------------------------------------------------------------------------
// away from it
//----------------------------------------------------------------------------

        /*
         * 273.15 K at 1 bar. All four fluids are well into the gas phase here,
         * so the real gas departure is small and what remains is the
         * difference between the two caloric representations. Measured worst
         * case over the four is 0.006 % on entropy.
         */
        const real T1 = 273.15 ;

        // measured worst case over the four fluids: 0.006 %
        EXPECT_NEAR( tHelmholtz.s( T1, p ) / tIdealGas.s( T1, p ), 1.0, 1e-4 );

        /*
         * Enthalpy is held looser than entropy on purpose. Both inherit the
         * same 0.2 to 0.4 % difference between the two caloric fits, but for
         * oxygen and nitrogen the absolute enthalpy at 273 K is small enough
         * that the same discrepancy is a larger fraction of it. Measured worst
         * case over the four fluids: 0.017 %.
         */
        EXPECT_NEAR( tHelmholtz.h( T1, p ) / tIdealGas.h( T1, p ), 1.0, 5e-4 );

        /*
         * Heat capacity is a looser story and deliberately checked looser: the
         * Helmholtz ideal gas part and the CEA polynomial are independent fits
         * of the same quantity and differ by up to 0.4 % here. That is a real
         * difference between the two models, not an offset that could be
         * removed.
         */
        EXPECT_NEAR( tHelmholtz.cp( T1, p ) / tIdealGas.cp( T1, p ), 1.0, 1e-2 );
    }
}
