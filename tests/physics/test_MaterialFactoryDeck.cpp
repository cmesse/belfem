/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * The deck-facing contract of MaterialFactory and material::Alloy.
 *
 * AlloyRRR: an alloy derives its residual resistivity from RRR in the
 * constructor, so a missing or non-physical RRR is refused there, whether the
 * value comes from the factory argument or from the deck; an in-formula
 * rrr<N> token and the RRR key build the same material.
 *
 * BuiltinForms: the three spellings of a builtin section -- the inner key,
 * the section label, and the legacy 'name : builtin ;' header -- each select
 * the material named, with the deck's RRR forwarded into the constructor.
 *
 * AllowList: a key or subsection the resolved shape never reads is refused
 * instead of dropped.
 *
 * ConstantLaws: constant jc and n through the deck leave the power law
 * finite, with ec defaulted by the factory.
 *
 * Every factory here is built with tables off: no rho database is written
 * into the working directory; the temperature splines are still built. On a
 * refused deck whatever was built before the throw leaks: the factory does
 * not own its materials, a throwing constructor runs no destructor, and an
 * alloy refused inside its own constructor leaks its components too. One
 * material per refused deck bounds that to a single case.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_InputFile.hpp"
#include "cl_Material.hpp"
#include "cl_Material_Metal.hpp"
#include "cl_Material_Copper.hpp"
#include "cl_Material_Silver.hpp"
#include "cl_Material_Aluminum.hpp"
#include "cl_Material_Alloy.hpp"
#include "cl_MaterialFactory.hpp"

using belfem::real ;
using belfem::Material ;
using belfem::MaterialFactory ;
using belfem::MaterialProperty ;
using belfem::MaterialType ;
using belfem::ResistivityLaw ;

namespace
{
    // the constant-jc overload of the power law reads rho( gTbulk ): pin it
    // for the duration of a test and hand the previous value back afterwards
    struct TbulkGuard
    {
        const real mSaved ;

        explicit TbulkGuard( const real aT ) : mSaved( belfem::gTbulk )
        {
            belfem::gTbulk = aT ;
        }

        ~TbulkGuard()
        {
            belfem::gTbulk = mSaved ;
        }
    };

    // a 'materials' deck in the test temp directory; the parser needs a
    // file. The path is absolute, so the case does not depend on the working
    // directory. Removed again by the destructor so the test leaves no litter
    struct Deck
    {
        const std::string mPath ;

        Deck( const std::string & aName, const std::string & aBody ) :
            mPath( ::testing::TempDir() + "t51_" + aName + ".conf" )
        {
            std::ofstream tFile( mPath );
            tFile << "materials\n{\n" << aBody << "\n}\n";
        }

        ~Deck()
        {
            std::remove( mPath.c_str() );
        }
    };

    void
    delete_materials( MaterialFactory & aFactory )
    {
        for ( auto & tEntry : aFactory.materials() )
        {
            delete tEntry.second ;
        }
    }

    // build the factory on a deck and hand back what its constructor threw;
    // an empty string means the deck was accepted. Only std::runtime_error is
    // caught, so any other exception fails the test
    std::string
    refusal( const std::string & aName, const std::string & aBody )
    {
        Deck tDeck( aName, aBody );
        belfem::InputFile tFile( tDeck.mPath );

        try
        {
            MaterialFactory tFactory( tFile.section( "materials" ), false );
            delete_materials( tFactory );
            return "" ;
        }
        catch ( const std::runtime_error & aError )
        {
            return aError.what();
        }
    }

    bool
    contains( const std::string & aText, const std::string & aPart )
    {
        return aText.find( aPart ) != std::string::npos ;
    }

    // the residual the regula-falsi in Metal::set_RRR() drives to zero,
    // evaluated on the rho_0 the material stores afterwards
    void
    expect_rrr_identity( Material * aMaterial, const real aRRR )
    {
        belfem::material::Metal * tMetal = dynamic_cast< belfem::material::Metal * >( aMaterial );
        ASSERT_NE( tMetal, nullptr );

        const real tTref = tMetal->constant_property( MaterialProperty::T_ref_rho_i );
        EXPECT_TRUE( std::isfinite( tTref ) );
        EXPECT_GT( tTref, 0.0 );

        const real tRho0 = tMetal->constant_property( MaterialProperty::rho_0 );
        EXPECT_TRUE( std::isfinite( tRho0 ) );
        EXPECT_GT( tRho0, 0.0 );
        EXPECT_NEAR( ( tMetal->rho_i_custom( tTref ) + tRho0 ) / tRho0, aRRR, 1e-9 );
    }

    const char * gAlloyRefusal = "requires a valid RRR > 1" ;
}

//------------------------------------------------------------------------------
// R1: the alloy's RRR gate
//------------------------------------------------------------------------------

TEST( AlloyRRR, AlloyWithoutValidRRRIsRefused )
{
    MaterialFactory tFactory ;

    // the factory argument, missing and at the boundary
    EXPECT_THROW( tFactory.create_material( "pb40sn60", BELFEM_QUIET_NAN, false ),
                  std::runtime_error );
    EXPECT_THROW( tFactory.create_material( "pb40sn60", 1.0, false ),
                  std::runtime_error );

    try
    {
        tFactory.create_material( "pb40sn60", BELFEM_QUIET_NAN, false );
        ADD_FAILURE() << "an alloy without RRR was accepted" ;
    }
    catch ( const std::runtime_error & aError )
    {
        EXPECT_TRUE( contains( aError.what(), gAlloyRefusal ) ) << aError.what();
    }

    // just above the boundary the alloy constructs
    Material * tAlloy = tFactory.create_material( "pb40sn60", 1.5, false );
    ASSERT_NE( tAlloy, nullptr );
    EXPECT_DOUBLE_EQ( tAlloy->constant_property( MaterialProperty::RRR ), 1.5 );
    delete tAlloy ;

    // the deck twins: no RRR key, and RRR at the boundary
    std::string tMessage = refusal( "alloy_no_rrr",
        "    solder\n    {\n        builtin : pb40sn60 ;\n    }" );
    EXPECT_TRUE( contains( tMessage, gAlloyRefusal ) ) << tMessage ;

    tMessage = refusal( "alloy_rrr_one",
        "    solder\n    {\n        builtin : pb40sn60 ;\n        RRR : 1 ;\n    }" );
    EXPECT_TRUE( contains( tMessage, gAlloyRefusal ) ) << tMessage ;
}

//------------------------------------------------------------------------------
// R2: the rrr token and the RRR key
//------------------------------------------------------------------------------

TEST( AlloyRRR, RRRTokenAndArgumentAgree )
{
    // three routes to the same alloy: the RRR key beside the formula, the
    // token inside the formula as a section type ( lowercased by the
    // parser ), and the token in mixed case through the inner key, which
    // keeps the deck's spelling
    Deck tDeck( "alloy_routes",
        "    solder\n    {\n        builtin : pb40sn60 ;\n        RRR : 5 ;\n    }\n"
        "    pb40sn60rrr5 : builtin\n    {\n    }\n"
        "    mixed\n    {\n        builtin : Pb40Sn60RRR5 ;\n    }" );
    belfem::InputFile tFile( tDeck.mPath );
    MaterialFactory tFactory( tFile.section( "materials" ), false );

    ASSERT_EQ( tFactory.materials().size(), 3u );
    ASSERT_TRUE( tFactory.materials().key_exists( "solder" ) );
    ASSERT_TRUE( tFactory.materials().key_exists( "pb40sn60rrr5" ) );
    ASSERT_TRUE( tFactory.materials().key_exists( "mixed" ) );

    Material * tByKey   = tFactory.materials()( "solder" );
    Material * tByToken = tFactory.materials()( "pb40sn60rrr5" );
    Material * tByMixed = tFactory.materials()( "mixed" );

    for ( Material * tMat : { tByKey, tByToken, tByMixed } )
    {
        ASSERT_NE( dynamic_cast< belfem::material::Alloy * >( tMat ), nullptr );
        EXPECT_EQ( tMat->type(), MaterialType::PureMetal );
        EXPECT_DOUBLE_EQ( tMat->constant_property( MaterialProperty::RRR ), 5.0 );
    }

    // the label keeps the name the material was created under; only the
    // values are required to agree
    EXPECT_EQ( tByKey->label(),   "pb40sn60" );
    EXPECT_EQ( tByToken->label(), "pb40sn60rrr5" );
    EXPECT_EQ( tByMixed->label(), "Pb40Sn60RRR5" );

    const real T = 61.3 ;

    for ( Material * tOther : { tByToken, tByMixed } )
    {
        for ( MaterialProperty tProperty : { MaterialProperty::rho_0,
                                             MaterialProperty::M,
                                             MaterialProperty::T_max } )
        {
            const real tReference = tByKey->constant_property( tProperty );
            EXPECT_TRUE( std::isfinite( tReference ) );
            EXPECT_NEAR( tOther->constant_property( tProperty ), tReference,
                         1e-12 * std::abs( tReference ) ) << tOther->label();
        }

        const real tRho = tByKey->rho( T );
        EXPECT_TRUE( std::isfinite( tRho ) );
        EXPECT_NEAR( tOther->rho( T ), tRho, 1e-12 * std::abs( tRho ) ) << tOther->label();
    }

    delete_materials( tFactory );
}

//------------------------------------------------------------------------------
// R3: the three builtin forms
//------------------------------------------------------------------------------

TEST( BuiltinForms, ThreeBuiltinFormsAndAnRRRYieldFourMaterials )
{
    // form 1: the inner key; form 2: the section label; form 3: the legacy
    // header, which the parser stores both as a key of 'materials' and as the
    // section's own label; and an alloy through form 1
    Deck tDeck( "builtin_forms",
        "    cu\n    {\n        builtin : copper ;\n        RRR : 100 ;\n    }\n"
        "    silver : builtin\n    {\n        RRR : 50 ;\n    }\n"
        "    aluminum : builtin ;\n    {\n        RRR : 20 ;\n    }\n"
        "    solder\n    {\n        builtin : pb40sn60 ;\n        RRR : 5 ;\n    }" );
    belfem::InputFile tFile( tDeck.mPath );
    MaterialFactory tFactory( tFile.section( "materials" ), false );

    ASSERT_EQ( tFactory.materials().size(), 4u );

    for ( const char * tKey : { "cu", "silver", "aluminum", "solder" } )
    {
        ASSERT_TRUE( tFactory.materials().key_exists( tKey ) ) << tKey ;
    }

    Material * tCopper   = tFactory.materials()( "cu" );
    Material * tSilver   = tFactory.materials()( "silver" );
    Material * tAluminum = tFactory.materials()( "aluminum" );
    Material * tSolder   = tFactory.materials()( "solder" );

    // each form built the class it names, and no other of the four
    EXPECT_NE( dynamic_cast< belfem::material::Copper * >( tCopper ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Silver * >( tCopper ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Aluminum * >( tCopper ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Alloy * >( tCopper ), nullptr );

    EXPECT_NE( dynamic_cast< belfem::material::Silver * >( tSilver ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Copper * >( tSilver ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Aluminum * >( tSilver ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Alloy * >( tSilver ), nullptr );

    EXPECT_NE( dynamic_cast< belfem::material::Aluminum * >( tAluminum ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Copper * >( tAluminum ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Silver * >( tAluminum ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Alloy * >( tAluminum ), nullptr );

    EXPECT_NE( dynamic_cast< belfem::material::Alloy * >( tSolder ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Copper * >( tSolder ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Silver * >( tSolder ), nullptr );
    EXPECT_EQ( dynamic_cast< belfem::material::Aluminum * >( tSolder ), nullptr );

    // labels are class names; the map key is the deck section type
    EXPECT_EQ( tCopper->label(),   "Copper" );
    EXPECT_EQ( tSilver->label(),   "Silver" );
    EXPECT_EQ( tAluminum->label(), "Aluminum" );

    for ( Material * tMat : { tCopper, tSilver, tAluminum, tSolder } )
    {
        EXPECT_EQ( tMat->type(), MaterialType::PureMetal ) << tMat->label();
    }

    // the deck's RRR reached each constructor
    EXPECT_DOUBLE_EQ( tCopper->constant_property( MaterialProperty::RRR ),   100.0 );
    EXPECT_DOUBLE_EQ( tSilver->constant_property( MaterialProperty::RRR ),    50.0 );
    EXPECT_DOUBLE_EQ( tAluminum->constant_property( MaterialProperty::RRR ),  20.0 );
    EXPECT_DOUBLE_EQ( tSolder->constant_property( MaterialProperty::RRR ),     5.0 );

    {
        SCOPED_TRACE( "copper" );
        expect_rrr_identity( tCopper, 100.0 );
    }
    {
        SCOPED_TRACE( "silver" );
        expect_rrr_identity( tSilver, 50.0 );
    }
    {
        SCOPED_TRACE( "aluminum" );
        expect_rrr_identity( tAluminum, 20.0 );
    }

    // the alloy through the deck is the alloy through the argument
    MaterialFactory tDirect ;
    Material * tReference = tDirect.create_material( "pb40sn60", 5.0, false );
    const real tRho0 = tReference->constant_property( MaterialProperty::rho_0 );
    EXPECT_TRUE( std::isfinite( tRho0 ) );
    EXPECT_GT( tRho0, 0.0 );
    EXPECT_NEAR( tSolder->constant_property( MaterialProperty::rho_0 ), tRho0, 1e-12 * tRho0 );
    delete tReference ;

    delete_materials( tFactory );
}

//------------------------------------------------------------------------------
// R4: the per-shape allow-list
//------------------------------------------------------------------------------

TEST( AllowList, UnusedMaterialKeysAreRefused )
{
    // RRR on a superconductor: read by the factory, dropped by the
    // constructor, so only the allow-list can name it
    std::string tMessage = refusal( "ybco_rrr",
        "    ybco\n    {\n        builtin : ybco ;\n        jc : 3e10 A/m^2 ;\n"
        "        n : 25 ;\n        RRR : 50 ;\n    }" );
    EXPECT_TRUE( contains( tMessage, "never reads the key 'rrr'" ) ) << tMessage ;

    // a superconductor key on a pure metal
    tMessage = refusal( "copper_resistivity_type",
        "    cu\n    {\n        builtin : copper ;\n        resistivity type : riva ;\n    }" );
    EXPECT_TRUE( contains( tMessage, "never reads the key 'resistivity type'" ) ) << tMessage ;

    // a labelled subsection is looked up by nobody
    tMessage = refusal( "ybco_labelled_defect",
        "    ybco\n    {\n        builtin : ybco ;\n        jc : 3e10 A/m^2 ;\n"
        "        n : 25 ;\n        defect : mine\n        {\n        }\n    }" );
    EXPECT_TRUE( contains( tMessage, "carries a labelled subsection" ) ) << tMessage ;
}

//------------------------------------------------------------------------------
// R5: constant jc and n
//------------------------------------------------------------------------------

TEST( ConstantLaws, ConstantJcAndNLeaveThePowerLawFinite )
{
    TbulkGuard tGuard( 77.0 );

    Deck tDeck( "ybco_constants",
        "    ybco\n    {\n        builtin : ybco ;\n        jc : 3e10 A/m^2 ;\n        n : 25 ;\n    }" );
    belfem::InputFile tFile( tDeck.mPath );
    MaterialFactory tFactory( tFile.section( "materials" ), false );

    ASSERT_EQ( tFactory.materials().size(), 1u );
    ASSERT_TRUE( tFactory.materials().key_exists( "ybco" ) );
    Material * tHts = tFactory.materials()( "ybco" );
    EXPECT_EQ( tHts->type(), MaterialType::HTS );

    // the factory defaults ec when the deck omits it
    ASSERT_TRUE( tHts->is_constant( MaterialProperty::jc ) );
    ASSERT_TRUE( tHts->is_constant( MaterialProperty::n ) );
    ASSERT_TRUE( tHts->is_constant( MaterialProperty::ec ) );

    const real ec = tHts->constant_property( MaterialProperty::ec );
    const real jc = tHts->constant_property( MaterialProperty::jc );
    const real n  = tHts->constant_property( MaterialProperty::n );

    EXPECT_DOUBLE_EQ( ec, 1e-4 );
    EXPECT_DOUBLE_EQ( jc, 3e10 );
    EXPECT_DOUBLE_EQ( n, 25.0 );

    // the field accessors fall back to the constant
    EXPECT_DOUBLE_EQ( tHts->jc( 0.7, 0.4 ), jc );
    EXPECT_DOUBLE_EQ( tHts->jc( 0.7, 0.4, 77.0 ), jc );

    // the closed form the constant overload evaluates: the power law in
    // parallel with rho_n at the bulk temperature, floored at zero
    const real tRhoN = tHts->rho( 77.0 );
    ASSERT_TRUE( std::isfinite( tRhoN ) );
    ASSERT_GT( tRhoN, 0.0 );

    for ( real tRatio : { 0.8, 1.0, 1.3 } )
    {
        const real J = tRatio * jc ;
        const real tRhoPL = std::max( ( ec / jc ) * std::pow( J / jc, n - 1.0 ), 0.0 );
        const real tExpected = 1.0 / ( ( 1.0 / tRhoN ) + ( 1.0 / tRhoPL ) );

        const real tRho = tHts->rho_powerlaw( J );
        EXPECT_TRUE( std::isfinite( tRho ) ) << "J / jc = " << tRatio ;
        EXPECT_GT( tRho, 0.0 ) << "J / jc = " << tRatio ;
        EXPECT_NEAR( tRho, tExpected, 1e-12 * tExpected ) << "J / jc = " << tRatio ;
    }

    delete_materials( tFactory );

    // the riva boundary on a constant n: n = 1 is the ohmic resistor and is
    // refused, n = 25 constructs under the riva law
    std::string tMessage = refusal( "ybco_riva_n_one",
        "    ybco\n    {\n        builtin : ybco ;\n        jc : 3e10 A/m^2 ;\n"
        "        n : 1 ;\n        resistivity type : riva ;\n    }" );
    EXPECT_TRUE( contains( tMessage, "greater than 1" ) ) << tMessage ;

    Deck tRivaDeck( "ybco_riva",
        "    ybco\n    {\n        builtin : ybco ;\n        jc : 3e10 A/m^2 ;\n"
        "        n : 25 ;\n        resistivity type : riva ;\n    }" );
    belfem::InputFile tRivaFile( tRivaDeck.mPath );
    MaterialFactory tRivaFactory( tRivaFile.section( "materials" ), false );

    ASSERT_TRUE( tRivaFactory.materials().key_exists( "ybco" ) );
    EXPECT_EQ( tRivaFactory.materials()( "ybco" )->resistivity_law(), ResistivityLaw::Riva );

    delete_materials( tRivaFactory );
}
