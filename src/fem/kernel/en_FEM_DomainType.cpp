//
// Created by christian on 10/4/21.
//


#include "stringtools.hpp"
#include "assert.hpp"
#include "en_FEM_DomainType.hpp"
namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        string
        to_string( const DomainType aDomaintype )
        {
            switch ( aDomaintype )
            {
                case ( DomainType::SuperConductor ) :
                {
                    return "SuperConductor";
                }
                case ( DomainType::Coil ) :
                {
                    return "Coil";
                }
                case ( DomainType::FerroMagnetic ) :
                {
                    return "FerroMagnetic";
                }
                case ( DomainType::Air ) :
                {
                    return "Air";
                }
                case ( DomainType::InterfaceScAir ) :
                case ( DomainType::InterfaceFmAir ) :
                case ( DomainType::InterfaceScFm ) :
                {
                    return "Interface";
                }
                case ( DomainType::Cut ) :
                {
                    return "Cut";
                }
                case ( DomainType::Boundary ) :
                {
                    return "Boundary";
                }
                case ( DomainType::ThinShell ) :
                {
                    return "Shell";
                }
                default:
                {
                    return "Undefined";
                }
            }
        }

//------------------------------------------------------------------------------

        DomainType
        domain_type( const string & aString )
        {
            string tString = string_to_lower( aString );

            if ( tString == "superconductor" )
            {
                return DomainType::SuperConductor;
            }
            else if ( tString == "coil" )
            {
                return DomainType::Coil;
            }
            else if ( tString == "ferro" || tString == "ferromagnetic" )
            {
                return DomainType::FerroMagnetic;
            }
            else if ( tString == "air" || tString == "vacuum" )
            {
                return DomainType::Air;
            }
            else if ( tString == "cut" )
            {
                return DomainType::Cut ;
            }
            else if ( tString == "boundary" )
            {
                return DomainType::Boundary;
            }
            else if ( tString == "thinshell" || tString == "tape" || tString == "shell"  )
            {
                return DomainType::ThinShell;
            }
            else
            {
                BELFEM_ERROR( false, "Unknown Domain Type: %s", aString.c_str());
                return DomainType::UNDEFINED;
            }
        }

//------------------------------------------------------------------------------
    }
}