/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 *
 * Unit tests for the jc/n law angle contract: the caller delivers
 * the UNFOLDED field-to-tape-normal angle theta in [0, pi]; laws that are
 * even in theta must evaluate identically on both lobes.
 */

#include <gtest/gtest.h>

#include "constants.hpp"
#include "cl_JcFunction_ModifiedKim.hpp"

using namespace belfem;
using namespace belfem::material;

//------------------------------------------------------------------------------

// ModifiedKim uses cos^2/sin^2 only, so it must be even under
// theta -> pi - theta: the unfolded angle changes nothing for this law
TEST( JcFunctionModifiedKim, EvenUnderUnfold )
{
    // YBCO-ish parameters: jc0, B0, k, alpha
    JcFunctionModifiedKim tLaw( 3e10, 0.04, 0.3, 0.8 );

    for ( real tB : { 0.01, 0.1, 0.5, 2.0 } )
    {
        for ( uint k = 0; k <= 90; ++k )
        {
            const real tTheta = k * constant::pi / 180.0 ;
            const real tLo = tLaw.eval( tB, tTheta );
            const real tHi = tLaw.eval( tB, constant::pi - tTheta );
            EXPECT_NEAR( tHi / tLo, 1.0, 1e-14 )
                << "theta = " << k << " deg, B = " << tB ;
        }
    }
}

//------------------------------------------------------------------------------

// the historical folded input [0, pi/2] and the unfolded input give the
// same value for an even law: acos-free cross-check of that contract
TEST( JcFunctionModifiedKim, FoldedInputEquivalent )
{
    JcFunctionModifiedKim tLaw( 3e10, 0.04, 0.3, 0.8 );

    const real tB = 0.3 ;
    for ( uint k = 91; k <= 180; ++k )
    {
        const real tTheta  = k * constant::pi / 180.0 ;
        const real tFolded = constant::pi - tTheta ;
        EXPECT_NEAR( tLaw.eval( tB, tTheta ) / tLaw.eval( tB, tFolded ),
                     1.0, 1e-14 );
    }
}

//------------------------------------------------------------------------------
