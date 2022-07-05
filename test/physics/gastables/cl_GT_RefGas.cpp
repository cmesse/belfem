//
// Created by Christian Messe on 30.08.19.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_linspace.hpp"
#include "fn_r2.hpp"

#define protected public
#define private   public
#include "GT_globals.hpp"
#include "cl_GT_RefGasFactory.hpp"
#include "cl_GT_RefGas.hpp"
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;

TEST( GASTABLES, RefGas )
{
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Initialization
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    gastables::RefGasFactory tFactory;

    // create a new reference gas of type Nitrogen
    gastables::RefGas * tGas = tFactory.create_refgas( "N2" );

    for ( uint j=0; j<2; ++j )
    {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// reference data taken from NIST Chemistry WebBook, SRD 69
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // number of samples to be tested
        const uint tN = 19;

        // Temperatures to be tested
        Vector<real> tT = linspace( 200.0, 2000.0, tN );
        const real tEpsilon = 0.001;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Specific Heat
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const Vector<real> tCp = { 29.2305, 29.1714, 29.2727, 29.594, 30.1181,
                                   30.7609, 31.4389, 32.096, 32.7027, 33.2475,
                                   33.7294, 34.1523, 34.5224, 34.8462, 35.1302,
                                   35.38, 35.6006, 35.7963, 35.9707 };


        Vector<real> tValues( tN );


        // Cp test
        for ( uint k = 0; k < tN; ++k )
        {
            tValues( k ) = tGas->Cp( tT( k ));
        }

        EXPECT_NEAR( r2( tValues, tCp ), 1.0, tEpsilon );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Enthalpy
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const Vector<real> tH = { 5798.88, 8717.69, 11638.3, 14579.7, 17563.9,
                                  20607.3, 23717.3, 26894.3, 30134.8, 33432.8,
                                  36782.2, 40176.7, 43610.9, 47079.6, 50578.8,
                                  54104.5, 57653.8, 61223.8, 64812.4 };


        for ( uint k = 0; k < tN; ++k )
        {
            tValues( k ) = tGas->H( tT( k ));
        }

        EXPECT_NEAR( r2( tValues, tH ), 1.0, tEpsilon );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Entropy
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const Vector<real> tS = { 179.935, 191.77, 200.172, 206.734, 212.174,
                                  216.864, 221.015, 224.757, 228.17, 231.313,
                                  234.227, 236.944, 239.489, 241.882, 244.14,
                                  246.278, 248.306, 250.236, 252.077 };

        for ( uint k = 0; k < tN; ++k )
        {
            tValues( k ) = tGas->S( tT( k ));
        }

        EXPECT_NEAR( r2( tValues, tS ), 1.0, tEpsilon );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Dynamic Viscosity
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // mu
        const Vector<real> tMu = { 1.29235E-05, 1.78971E-05, 2.22166E-05, 2.60779E-05, 2.96021E-05,
                                   3.28684E-05, 3.59317E-05, 3.88315E-05, 4.15973E-05, 4.42515E-05,
                                   4.68117E-05, 4.92920E-05, 5.17037E-05, 5.40561E-05, 5.63570E-05,
                                   5.86129E-05, 6.08293E-05, 6.30109E-05, 6.51618E-05 };

        for ( uint k = 0; k < tN; ++k )
        {
            tValues( k ) = tGas->mu( tT( k ));
        }

        EXPECT_NEAR( r2( tValues, tMu ), 1.0, tEpsilon );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Thermal Conductivity
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const Vector<real> tLambda = { 0.0186592, 0.025858, 0.032205, 0.0381425, 0.0439172,
                                       0.0496046, 0.055197, 0.0606663, 0.0659906, 0.0711597,
                                       0.0761733, 0.0810379, 0.0857632, 0.090361, 0.0948433,
                                       0.0992215, 0.103507, 0.107709, 0.111837 };

        for ( uint k = 0; k < tN; ++k )
        {
            tValues( k ) = tGas->lambda( tT( k ));
        }

        EXPECT_NEAR( r2( tValues, tLambda ), 1.0, tEpsilon );

        tGas->set_mode( RefGasMode::SPLINE );

    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // delete the gas object
    delete tGas;
}