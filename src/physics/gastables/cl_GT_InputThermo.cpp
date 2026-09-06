/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <cstdlib>

#include "stringtools.hpp"
#include "filetools.hpp"
#include "cl_Vector.hpp"

#include "cl_GT_HeatPoly.hpp"
#include "cl_GT_HeatPolyCustom.hpp"

#include "cl_GT_InputThermo.hpp"

#include "fn_GT_fix_label.hpp"
#include "fn_GT_fix_capitals.hpp"
#include "fn_GT_parse.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        InputThermo::InputThermo( const string & aPath )
                : Ascii( aPath, FileMode::OPEN_RDONLY )
        {
            // line counter
            uint tLineCount = 0;

            for ( const std::string & tLine : mBuffer )
            {
                if ( tLine.size() > 0 )
                {
                    // get first character of line
                    char tChar = tLine.substr( 0, 1 ).c_str()[ 0 ];

                    // first character must neither be minus, space or exclamation mark
                    if ( tChar != 32 && tChar != 33 && tChar != 45 )
                    {
                        // get name of gas
                        string tLabel = fix_label( tLine.substr( 0, tLine.find_first_of( " " )));

                        // special treatment of "(a)" tag
                        auto tFlag = tLabel.find( "(a)" );

                        if ( tFlag < tLabel.length())
                        {
                            tLabel = tLabel.substr( 0, tFlag );
                        }

                        // skip the file header and the END PRODUCTS / END REACTANTS markers
                        if ( tLabel != "thermo" && tLabel.substr( 0, 3 ) != "END" )
                        {
                            // add word to map
                            mMap[ tLabel ] = tLineCount;
                        }
                    }
                }
                // increment counter
                ++tLineCount;
            }
        }

//------------------------------------------------------------------------------

        bool
        InputThermo::entry_exists( RefGas * aRefgas )
        {
            return mMap.key_exists( aRefgas->label() );
        }

//------------------------------------------------------------------------------

        void
        InputThermo::read_data( RefGas * aRefgas )
        {
            if ( this->entry_exists( aRefgas ) )
            {
                // name of the record, for error messages
                const string & tName = aRefgas->label();

                // get line in buffer
                uint tLineCount = mMap( tName ) + 1;

                BELFEM_ERROR( tLineCount < mBuffer.size(),
                        "thermo.inp truncated at record of %s", tName.c_str() );

                string tLine = mBuffer( tLineCount++ );

                // fields reach to the liquid flag at column 51
                BELFEM_ERROR( tLine.length() >= 52,
                        "Record line of %s in thermo.inp is too short",
                        tName.c_str() );

                // get number of polys
                uint tNumPolys = parse_int( tLine.substr( 0, 2 ),
                        "polynomial count in thermo.inp" );

                // default exponents
                Vector<real> tDefaultExponents = { -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0 };

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // check if components have been generated already, which happens
                // when a species is read from more than one source
                if( ! aRefgas->has_components() )
                {
                    // read the composition: five pairs of A2 symbol + F6.2 amount
                    uint tOff = 10;

                    for ( uint k = 0; k < 5; ++k )
                    {
                        string tComp = fix_capitals( clean_string( tLine.substr( tOff, 2 )));
                        real tValue = parse_real_optional( tLine.substr( tOff + 2, 6 ),
                                0.0, "composition in thermo.inp" );

                        if ( tComp.length() > 0 )
                        {
                            // add component to gas
                            aRefgas->add_component( tComp, tValue );
                        }

                        tOff += 8;
                    }

                    aRefgas->set_component_flag();
                }
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // read the liquid flag
                uint tFlag = parse_int( tLine.substr( 51, 1 ),
                        "liquid flag in thermo.inp" );
                if ( tFlag == 1 )
                {
                    aRefgas->set_liquid_flag();
                }
                else
                {
                    aRefgas->unset_liquid_flag();
                }

                // create words from line
                Cell< string > tWords = string_to_words( tLine );

                // read the molar mass
                aRefgas->set_molar_mass( parse_real( tWords( tWords.size() - 2 ),
                        "molar mass in thermo.inp" ) * 0.001 );

                // read the formation enthalpy
                aRefgas->set_reference_formation_enthalpy( parse_real(
                        tWords( tWords.size() - 1 ), "formation enthalpy in thermo.inp" ) );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                for ( uint k = 0; k < tNumPolys; ++k )
                {
                    BELFEM_ERROR( tLineCount < mBuffer.size(),
                            "thermo.inp truncated in record of %s", tName.c_str() );

                    // get next line
                    tLine = mBuffer( tLineCount++ );

                    // skip commented and empty lines
                    while ( tLine.empty() || tLine.at( 0 ) == '!' )
                    {
                        BELFEM_ERROR( tLineCount < mBuffer.size(),
                                "thermo.inp truncated in record of %s", tName.c_str() );

                        // jump to next line
                        tLine = mBuffer( tLineCount++ );
                    }

                    // the CEA range line is 80 columns wide
                    BELFEM_ERROR( tLine.length() >= 80,
                            "Range line too short in record of %s", tName.c_str() );

                    // read lower temperature
                    real tTmin = parse_real( tLine.substr( 1, 10 ),
                            "Tmin in thermo.inp" );

                    // read upper temperature
                    real tTmax = parse_real( tLine.substr( 12, 10 ),
                            "Tmax in thermo.inp" );

                    // read number of exponents
                    uint tN = parse_int( tLine.substr( 22, 1 ),
                            "exponent count in thermo.inp" );

                    Vector<real> tExponents( tN );


                    // read the exponents
                    uint tOff = 23;
                    for ( uint i = 0; i < tN; ++i )
                    {
                        tExponents( i ) = parse_real( tLine.substr( tOff, 5 ),
                                "exponent in thermo.inp" );
                        tOff += 5;
                    }

                    // only needed for first polynomial
                    if ( k == 0 )
                    {
                        // read the reference enthalpy
                        aRefgas->set_reference_enthalpy( parse_real(
                                tLine.substr( 64, 16 ), "reference enthalpy in thermo.inp" ));
                    }

                    BELFEM_ERROR( tLineCount < mBuffer.size(),
                            "thermo.inp truncated in record of %s", tName.c_str() );

                    // next line
                    tLine = mBuffer( tLineCount++ );

                    BELFEM_ERROR( tLine.length() >= 80,
                            "Coefficient line too short in record of %s", tName.c_str() );

                    Vector<real> tCoefficients( tN );

                    tOff = 0;
                    uint tM = ( tN < 5 ) ? tN : 5;

                    for ( uint i = 0; i < tM; ++i )
                    {
                        tCoefficients( i ) = parse_real(
                                tLine.substr( tOff, 16 ).replace( 12, 1, "e" ),
                                "coefficient in thermo.inp" );
                        tOff += 16;
                    }

                    BELFEM_ERROR( tLineCount < mBuffer.size(),
                            "thermo.inp truncated in record of %s", tName.c_str() );

                    // next line
                    tLine = mBuffer( tLineCount++ );

                    BELFEM_ERROR( tLine.length() >= 80,
                            "Coefficient line too short in record of %s", tName.c_str() );

                    tOff = 0;
                    for ( uint i = 5; i < tN; ++i )
                    {
                        tCoefficients( i ) = parse_real(
                                tLine.substr( tOff, 16 ).replace( 12, 1, "e" ),
                                "coefficient in thermo.inp" );
                        tOff += 16;
                    }

                    // enthalpy constant
                    real tEnthalpyConstant = parse_real(
                            tLine.substr( 48, 16 ).replace( 12, 1, "e" ),
                            "enthalpy constant in thermo.inp" );

                    // entropy constant
                    real tEntropyConstant = parse_real(
                            tLine.substr( 64, 16 ).replace( 12, 1, "e" ),
                            "entropy constant in thermo.inp" );

                    // pointer to polynomial
                    HeatPoly * tHeatPoly = nullptr;

                    if ( tN == 7 && tExponents == tDefaultExponents )
                    {
                        // create a default heat polynomial
                        tHeatPoly = new HeatPoly(
                                tTmin,
                                tTmax,
                                tEnthalpyConstant,
                                tEntropyConstant,
                                tCoefficients );
                    }
                    else
                    {
                        // create a custom polynomial
                        tHeatPoly = new HeatPolyCustom(
                                tTmin,
                                tTmax,
                                tEnthalpyConstant,
                                tEntropyConstant,
                                tCoefficients,
                                tExponents );
                    }
                    // add polynomial to list
                    aRefgas->add_heat_poly( tHeatPoly );
                }
            }
        }

    } /* namespace gastables */
} /* namespace belfem */