/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * The removed 'expk' key is refused by the two boundary-condition factories,
 * on every source shape, and a section without it still builds. The circuit
 * copy of the same refusal is covered in tests/circuit/test_HybridFactory.cpp.
 *
 * Both factories are built straight from a 'boundary conditions' section; no
 * mesh is involved. They hand ownership of their conditions to the kernel, so
 * an accepted deck deletes them here. On a refused deck the factory frees
 * what it built before the throw: every condition, with the source function
 * it already owns, and the function of the condition being refused.
 *
 * Mu0Band: the thin-shell normal-field recovery admits an h-conductor whose
 * constant permeability is mu0 within a relative rounding allowance of 1e-6.
 */

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_InputFile.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"
#include "cl_MaxwellBoundaryConditionFactory.hpp"
#include "cl_ThermalBoundaryConditionFactory.hpp"
#include "cl_MaxwellFactory.hpp"
#include "constants.hpp"

using namespace belfem;

namespace
{
    // a 'boundary conditions' deck in the test temp directory, under an
    // absolute path so the case does not depend on the working directory;
    // removed on scope exit
    class Deck
    {
        const std::string mPath ;

    public:

        Deck( const std::string & aName, const std::string & aBody ) :
            mPath( ::testing::TempDir() + aName )
        {
            std::ofstream tFile( mPath );
            tFile << "boundary conditions\n{\n" << aBody << "\n}\n";
        }

        Deck( const Deck & ) = delete ;
        Deck( Deck && ) = delete ;
        Deck & operator=( const Deck & ) = delete ;
        Deck & operator=( Deck && ) = delete ;

        ~Deck()
        {
            std::remove( mPath.c_str() );
        }

        const std::string &
        path() const
        {
            return mPath ;
        }
    };

    bool
    contains( const std::string & aText, const std::string & aPart )
    {
        return aText.find( aPart ) != std::string::npos ;
    }

    void
    delete_conditions( Cell< fem::PhysicalBoundaryCondition * > & aConditions )
    {
        for ( fem::PhysicalBoundaryCondition * tBC : aConditions )
        {
            delete tBC ;
        }
        aConditions.clear();
    }

    // build the Maxwell factory on a deck and hand back what its constructor
    // threw; an empty string means the deck was accepted. Only
    // std::runtime_error is caught, so any other exception fails the test
    std::string
    maxwell_refusal( const std::string & aName, const std::string & aBody )
    {
        Deck tDeck( aName, aBody );
        InputFile tFile( tDeck.path() );

        try
        {
            fem::MaxwellBoundaryConditionFactory tFactory(
                    tFile.section( "boundary conditions" ), 3 );
            EXPECT_EQ( tFactory.boundary_conditions().size(), 1u );
            delete_conditions( tFactory.boundary_conditions() );
            return "" ;
        }
        catch ( const std::runtime_error & aError )
        {
            return aError.what();
        }
    }

    // the same for the thermal factory
    std::string
    thermal_refusal( const std::string & aName, const std::string & aBody )
    {
        Deck tDeck( aName, aBody );
        InputFile tFile( tDeck.path() );

        try
        {
            fem::ThermalBoundaryConditionFactory tFactory(
                    tFile.section( "boundary conditions" ) );
            EXPECT_EQ( tFactory.boundary_conditions().size(), 1u );
            delete_conditions( tFactory.boundary_conditions() );
            return "" ;
        }
        catch ( const std::runtime_error & aError )
        {
            return aError.what();
        }
    }

    const char * gExpkRefusal = "no longer accepted" ;

    // one section of type aType on sideset 1; aShape is the waveform body.
    // The Maxwell cases use a background field: its amplitude has a unit
    // ( A/m or T ), while a Maxwell dirichlet section checks the amplitude
    // against an empty unit and refuses a dimensioned one
    std::string
    bc_section( const std::string & aType,
                const std::string & aShape,
                const bool          aWithExpk )
    {
        return "    " + aType + "\n    {\n        sideset : 1 ;\n" + aShape
             + ( aWithExpk ? "        expk : 99 ;\n" : "" ) + "    }" ;
    }

    std::string
    background( const std::string & aShape, const bool aWithExpk )
    {
        return bc_section( "background",
                           "        direction : 0, 0, 1 ;\n" + aShape, aWithExpk );
    }

    std::string
    dirichlet( const std::string & aShape, const bool aWithExpk )
    {
        return bc_section( "dirichlet", aShape, aWithExpk );
    }

    const std::string gMaxwellSigmoid =
        "        type : sigmoid ;\n        amplitude : 1 A/m ;\n"
        "        period : 1 s ;\n        offset : 0 s ;\n" ;

    const std::string gMaxwellConstant =
        "        type : constant ;\n        amplitude : 1 A/m ;\n" ;

    const std::string gThermalSigmoid =
        "        type : sigmoid ;\n        amplitude : 80 K ;\n"
        "        period : 1 s ;\n        offset : 0 s ;\n" ;

    const std::string gThermalConstant =
        "        type : constant ;\n        amplitude : 80 K ;\n" ;
}

//------------------------------------------------------------------------------
// Maxwell
//------------------------------------------------------------------------------

// the refusal sits before the waveform dispatch: the key was only ever read
// on the sigmoid branch, so it is refused on every shape, and on a section
// with no 'type' at all
TEST( AdmissionRefusals, ExpkOnAMaxwellBcIsRefused )
{
    std::string tMessage = maxwell_refusal( "expk_mx_sigmoid.conf",
                                            background( gMaxwellSigmoid, true ) );
    EXPECT_TRUE( contains( tMessage, gExpkRefusal ) ) << tMessage ;

    tMessage = maxwell_refusal( "expk_mx_constant.conf",
                                background( gMaxwellConstant, true ) );
    EXPECT_TRUE( contains( tMessage, gExpkRefusal ) ) << tMessage ;

    tMessage = maxwell_refusal( "expk_mx_untyped.conf", background( "", true ) );
    EXPECT_TRUE( contains( tMessage, gExpkRefusal ) ) << tMessage ;
}

//------------------------------------------------------------------------------

TEST( AdmissionRefusals, MaxwellBcWithoutExpkStillBuilds )
{
    std::string tMessage = maxwell_refusal( "noexpk_mx_sigmoid.conf",
                                            background( gMaxwellSigmoid, false ) );
    EXPECT_EQ( tMessage, "" );

    tMessage = maxwell_refusal( "noexpk_mx_constant.conf",
                                background( gMaxwellConstant, false ) );
    EXPECT_EQ( tMessage, "" );
}

//------------------------------------------------------------------------------

// the first section is built with its function before the second is refused,
// so the refusal runs over a condition that already owns its function. The
// freeing itself can only turn red under valgrind
TEST( AdmissionRefusals, ARefusedSecondMaxwellBcFreesTheFirst )
{
    const std::string tMessage = maxwell_refusal( "expk_mx_second.conf",
        background( gMaxwellConstant, false ) + "\n" + background( gMaxwellSigmoid, true ) );
    EXPECT_TRUE( contains( tMessage, gExpkRefusal ) ) << tMessage ;
}

//------------------------------------------------------------------------------
// thermal
//------------------------------------------------------------------------------

TEST( AdmissionRefusals, ExpkOnAThermalBcIsRefused )
{
    std::string tMessage = thermal_refusal( "expk_th_sigmoid.conf",
                                            dirichlet( gThermalSigmoid, true ) );
    EXPECT_TRUE( contains( tMessage, gExpkRefusal ) ) << tMessage ;

    tMessage = thermal_refusal( "expk_th_constant.conf",
                                dirichlet( gThermalConstant, true ) );
    EXPECT_TRUE( contains( tMessage, gExpkRefusal ) ) << tMessage ;

    tMessage = thermal_refusal( "expk_th_untyped.conf", dirichlet( "", true ) );
    EXPECT_TRUE( contains( tMessage, gExpkRefusal ) ) << tMessage ;
}

//------------------------------------------------------------------------------

TEST( AdmissionRefusals, ThermalBcWithoutExpkStillBuilds )
{
    std::string tMessage = thermal_refusal( "noexpk_th_sigmoid.conf",
                                            dirichlet( gThermalSigmoid, false ) );
    EXPECT_EQ( tMessage, "" );

    tMessage = thermal_refusal( "noexpk_th_constant.conf",
                                dirichlet( gThermalConstant, false ) );
    EXPECT_EQ( tMessage, "" );
}

//------------------------------------------------------------------------------

// the thermal twin of ARefusedSecondMaxwellBcFreesTheFirst
TEST( AdmissionRefusals, ARefusedSecondThermalBcFreesTheFirst )
{
    const std::string tMessage = thermal_refusal( "expk_th_second.conf",
        dirichlet( gThermalConstant, false ) + "\n" + dirichlet( gThermalSigmoid, true ) );
    EXPECT_TRUE( contains( tMessage, gExpkRefusal ) ) << tMessage ;
}

//------------------------------------------------------------------------------

TEST( Mu0Band, RoundedMu0IsAdmittedAndARealSusceptibilityIsNot )
{
    using belfem::fem::MaxwellFactory ;
    const belfem::real tMu0 = belfem::constant::mu0 ;

    EXPECT_TRUE( MaxwellFactory::mu_is_mu0( tMu0 ) );
    EXPECT_TRUE( MaxwellFactory::mu_is_mu0( 4.0e-7 * belfem::constant::pi ) ) << "the pre-2019 exact value" ;
    EXPECT_TRUE( MaxwellFactory::mu_is_mu0( 1.256637e-6 ) ) << "seven significant digits" ;

    EXPECT_FALSE( MaxwellFactory::mu_is_mu0( 1.25664e-6 ) ) << "six significant digits" ;
    EXPECT_FALSE( MaxwellFactory::mu_is_mu0( 1.2566e-6 ) ) << "five significant digits" ;
    EXPECT_FALSE( MaxwellFactory::mu_is_mu0( tMu0 * ( 1.0 - 1.0e-5 ) ) ) << "copper's susceptibility" ;
    EXPECT_FALSE( MaxwellFactory::mu_is_mu0( 1000.0 * tMu0 ) );
    EXPECT_FALSE( MaxwellFactory::mu_is_mu0( BELFEM_QUIET_NAN ) );
}
