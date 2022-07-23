//
// Created by Christian Messe on 14.04.20.
//

#include "en_Materials.hpp"

namespace belfem
{
    string
    to_string( const MaterialType & aMaterialType )
    {
        switch ( aMaterialType )
        {
            case( MaterialType::Unity ) :
            {
                return "Unity" ;
            }
            case( MaterialType::Simple ) :
            {
                return "Simple" ;
            }
            case( MaterialType::Altramat80 ) :
            {
                return "Altramat80" ;
            }
            case( MaterialType::Aluminum ) :
            {
                return "Aluminum" ;
            }
            case( MaterialType::Copper ) :
            {
                return "Copper" ;
            }
            case( MaterialType::Inconel718 ) :
            {
                return "Inconel718" ;
            }
            case( MaterialType::TI6AL4V ) :
            {
                return "Ti-6Al-4V" ;
            }
            case( MaterialType::CCSIC ) :
            {
                return "C/C-SiC" ;
            }
            case( MaterialType::Rohacell51 ) :
            {
                return "Rohacell51" ;
            }
            case( MaterialType::CuCrZr ) :
            {
                return "CuCrZr" ;
            }
            case( MaterialType::Zirconia ) :
            {
                return "Zirconia" ;
            }
            case( MaterialType::Inconel750X ) :
            {
                return "Inconel750X" ;
            }
            case( MaterialType::Hastelloy ) :
            {
                return "Hastelloy-C276" ;
            }
            case( MaterialType::Silver ) :
            {
                return "Silver" ;
            }
            default:
            {
                return "undefined" ;
            }
        }
    }

    MaterialType
    string_to_material_type( const string & aString, const bool aNoErrorIfUnknown )
    {
        // make string lower case
        const string tString = string_to_lower( aString );

        if( tString == "unity" )
        {
            return MaterialType::Unity ;
        }
        if( tString == "air" )
        {
            return MaterialType::Air ;
        }
        else if ( tString == "simple" )
        {
            return MaterialType::Simple ;
        }
        else if ( tString == "aluminum" )
        {
            return MaterialType::Aluminum ;
        }
        else if ( tString == "copper" )
        {
            return MaterialType::Copper ;
        }
        else if ( tString == "inconel718" || tString == "inconel-718"  )
        {
            return MaterialType::Inconel718 ;
        }
        else if ( tString == "altramat80" || tString == "altramat"  )
        {
            return MaterialType::Altramat80 ;
        }
        else if ( tString == "ti-6al-4v" || tString == "ti6al4v"  || tString == "ti/6al/4v" )
        {
            return MaterialType::TI6AL4V ;
        }
        else if ( tString == "rohacell51" || tString == "rohacell" )
        {
            return MaterialType::Rohacell51 ;
        }
        else if ( tString == "ccsic" || tString == "c/c-sic"  )
        {
            return MaterialType::CCSIC ;
        }
        else if ( tString == "cucrzr" || tString == "cucr1zr"  )
        {
            return MaterialType::CuCrZr ;
        }
        else if ( tString == "zirconia" || tString == "zro2"  )
        {
            return MaterialType::Zirconia ;
        }
        else if ( tString == "inconel750x" || tString == "inconel-750x"  )
        {
            return MaterialType::Inconel750X ;
        }
        else if ( tString == "hastelloy" || tString == "hastelloyc276" || tString == "hastelloy-c276" || tString == "hastelloy c-276" )
        {
            return MaterialType::Hastelloy ;
        }
        else if ( tString == "silver" || tString == "ag" || tString == "argentum" )
        {
            return MaterialType::Silver ;
        }
        else
        {
            BELFEM_ERROR( aNoErrorIfUnknown, "Unknown material type: %s", aString.c_str() );
            return MaterialType::UNDEFINED ;
        }
    }
}