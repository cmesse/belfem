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
            case( MaterialType::HastelloyC276 ) :
            {
                return "Hastelloy-C276" ;
            }
            case( MaterialType::Silver ) :
            {
                return "Silver" ;
            }
            case( MaterialType::YBCO ) :
            {
                return "YBCO" ;
            }
            case( MaterialType::Pb40Sn60 ) :
            {
                return "Pb40Sn60" ;
            }
            case( MaterialType::SAE301 ) :
            {
                return "SAE-301" ;
            }
            case( MaterialType::SAE316 ) :
            {
                return "SAE-316" ;
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
        if ( tString == "simple" )
        {
            return MaterialType::Simple ;
        }
        if ( tString == "aluminum" || tString == "al" )
        {
            return MaterialType::Aluminum ;
        }
        if ( tString == "copper" || tString == "cu"  )
        {
            return MaterialType::Copper ;
        }
        if ( tString == "inconel718" || tString == "inconel-718" || tString == "2.4668" )
        {
            return MaterialType::Inconel718 ;
        }
        if ( tString == "altramat80" || tString == "altramat"  )
        {
            return MaterialType::Altramat80 ;
        }
        if ( tString == "ti-6al-4v" || tString == "ti6al4v"  || tString == "ti/6al/4v"  || tString == "3.7165" )
        {
            return MaterialType::TI6AL4V ;
        }
        if ( tString == "rohacell51" || tString == "rohacell" )
        {
            return MaterialType::Rohacell51 ;
        }
        if ( tString == "ccsic" || tString == "c/c-sic"  )
        {
            return MaterialType::CCSIC ;
        }
        if ( tString == "cucrzr" || tString == "cucr1zr"  )
        {
            return MaterialType::CuCrZr ;
        }
        if ( tString == "zirconia" || tString == "zro2"  )
        {
            return MaterialType::Zirconia ;
        }
        if ( tString == "inconel750x" || tString == "inconel-750x"  || tString == "2.4669" )
        {
            return MaterialType::Inconel750X ;
        }
        if ( tString == "hastelloy" || tString == "hastelloyc276" || tString == "hastelloy-c276" || tString == "hastelloy c-276" )
        {
            return MaterialType::HastelloyC276 ;
        }
        if ( tString == "silver" || tString == "ag" || tString == "argentum" )
        {
            return MaterialType::Silver ;
        }
        if ( tString == "ybco" || tString == "yb123" )
        {
            return MaterialType::YBCO ;
        }
        if (   tString == "pbsn" || tString == "pb40sn60" || tString == "pb40-sn60" || tString == "pb40/sn60"
                 || tString == "snpb" || tString == "sn60pb40" || tString == "sn60-pb40" || tString == "sn60/pb40" )
        {
            return MaterialType::Pb40Sn60 ;
        }
        if (   tString == "sae301" || tString == "sae-301" || tString == "aisi301" || tString ==  "aisi-301"
            || tString == "301" || tString ==  "1.4310")
        {
            return MaterialType::SAE301 ;
        }
        if (   tString == "sae316" || tString == "sae-316" || tString == "aisi316" || tString == "aisi-316"
                || tString == "314" || tString ==  "1.4401")
        {
            return MaterialType::SAE316 ;
        }
        else
        {
            BELFEM_ERROR( aNoErrorIfUnknown, "Unknown material type: %s", aString.c_str() );
            return MaterialType::UNDEFINED ;
        }
    }
}