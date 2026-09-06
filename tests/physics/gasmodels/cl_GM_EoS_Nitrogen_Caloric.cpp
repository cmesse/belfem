//
// Created by Christian Messe on 26.08.26.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "cl_Gas.hpp"

#include "cl_GM_EoS_Nitrogen.hpp"
#include "fn_GM_Helmholz_DerivTest.hpp"

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

/*
 * Span, Lemmon, Jacobsen and Wagner, A Reference Equation of State for the
 * Thermodynamic Properties of Nitrogen for Temperatures from 63.151 to 1000 K
 * and Pressures to 2200 MPa, J. Phys. Chem. Ref. Data 29( 6 ):1361 ( 2000 ).
 */
TEST( GASMODELS, Nitrogen_Caloric )
{
    Gas tRef( "N2" );

    EoS_Nitrogen tGas( tRef );

//----------------------------------------------------------------------------
// analytic derivatives against finite differences
//----------------------------------------------------------------------------

    /*
     * This catches an error in a derivative function, but NOT an error in the
     * coefficient tables: it only asserts that the derivatives are consistent
     * with the phi they are derived from. Nitrogen_Vapor covers the tables.
     */
    Vector< real > tR2( 7 );

    deriv_test( tGas, tR2 );

    // phi0_t
    EXPECT_NEAR( tR2( 0 ), 1.0, 1e-6 );

    // phi0_tt
    EXPECT_NEAR( tR2( 1 ), 1.0, 1e-6 );

    // phir_t
    EXPECT_NEAR( tR2( 2 ), 1.0, 1e-6 );

    // phir_d
    EXPECT_NEAR( tR2( 3 ), 1.0, 1e-6 );

    // phir_tt
    EXPECT_NEAR( tR2( 4 ), 1.0, 1e-6 );

    // phir_dt
    EXPECT_NEAR( tR2( 5 ), 1.0, 1e-6 );

    // phir_dd
    EXPECT_NEAR( tR2( 6 ), 1.0, 1e-6 );
}
