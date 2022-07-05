//
// Created by Christian Messe on 2019-08-18.
//




#include <cstdlib>

#include "stringtools.hpp"
#include "filetools.hpp"
#include "cl_Vector.hpp"

#include "cl_GT_HeatPoly.hpp"
#include "cl_GT_HeatPolyCustom.hpp"

#include "cl_GT_InputThermo.hpp"

#include "fn_GT_fix_label.hpp"
#include "fn_GT_fix_capitals.hpp"

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

            for ( std::string tLine : mBuffer )
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

                        // add word to map
                        mMap[ tLabel ] = tLineCount;
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
                // get line in buffer
                uint tLineCount = mMap( aRefgas->label()) + 1;

                string tLine = mBuffer( tLineCount++ );

                // get number of polys
                uint tNumPolys = std::stoi( tLine.substr( 0, 2 ));

                // default exponents
                Vector<real> tDefaultExponents = { -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0 };

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // check if components have been generated already
                // this can be the case if eg cryo data exist in a separate file
                if( ! aRefgas->has_components() )
                {
                    // read the composition
                    uint tOff = 10;

                    for ( uint k = 0; k < 5; ++k )
                    {
                        string tComp = fix_capitals( clean_string( tLine.substr( tOff, 3 )));
                        real tValue = std::stod( clean_string( tLine.substr( tOff + 3, 4 )));

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
                uint tFlag = std::stoi( tLine.substr( 51, 1 ));
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
                aRefgas->set_molar_mass( std::stod( tWords( tWords.size() - 2 ) ) * 0.001 );

                // read the formation enthalpy
                aRefgas->set_reference_formation_enthalpy( std::stod( tWords( tWords.size() - 1 ) ) );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                for ( uint k = 0; k < tNumPolys; ++k )
                {
                    // get next line
                    tLine = mBuffer( tLineCount++ );

                    // skip commented lines
                    while ( tLine.at( 0 ) == '!' )
                    {
                        // jump to next line
                        tLine = mBuffer( tLineCount++ );
                    }

                    // read lower temperature
                    real tTmin = std::stod( tLine.substr( 1, 10 ));

                    // read upper temperature
                    real tTmax = std::stod( tLine.substr( 12, 10 ));

                    // read number of exponents
                    uint tN = std::stoi( tLine.substr( 22, 1 ));

                    Vector<real> tExponents( tN );


                    // read the exponents
                    uint tOff = 23;
                    for ( uint i = 0; i < tN; ++i )
                    {
                        tExponents( i ) = std::stod( tLine.substr( tOff, 5 ));
                        tOff += 5;
                    }

                    // only needed for first polynomial
                    if ( k == 0 )
                    {
                        // read the reference enthalpy
                        aRefgas->set_reference_enthalpy( std::stod( tLine.substr( 64, 16 )));
                    }

                    // next line
                    tLine = mBuffer( tLineCount++ );

                    Vector<real> tCoefficients( tN );

                    tOff = 0;
                    uint tM = ( tN < 5 ) ? tN : 5;

                    for ( uint i = 0; i < tM; ++i )
                    {
                        tCoefficients( i ) = std::stod( tLine.substr( tOff, 16 ).replace( 12, 1, "e" ));
                        tOff += 16;
                    }

                    // next line
                    tLine = mBuffer( tLineCount++ );

                    tOff = 0;
                    for ( uint i = 5; i < tN; ++i )
                    {
                        tCoefficients( i ) = std::stod( tLine.substr( tOff, 16 ).replace( 12, 1, "e" ));
                        tOff += 16;
                    }

                    // enthalpy constant
                    real tEnthalpyConstant = std::stod( tLine.substr( 48, 16 ).replace( 12, 1, "e" ));

                    // entropy constant
                    real tEntropyConstant = std::stod( tLine.substr( 64, 16 ).replace( 12, 1, "e" ));

                    // pointer to polynomial
                    HeatPoly * tHeatPoly = nullptr;

                    if ( tExponents == tDefaultExponents )
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