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
                case ( DomainType::Conductor ) :
                {
                    return "Conductor";
                }
                case ( DomainType::Coil ) :
                {
                    return "Coil";
                }
                case ( DomainType::Ferro ) :
                {
                    return "Ferro";
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
                case ( DomainType::Symmetry ) :
                case ( DomainType::SymmetryAir ) :
                case ( DomainType::SymmetryFerro ) :
                case ( DomainType::SymmetryConductor ) :
                {
                    return "Symmetry";
                }
                case ( DomainType::AntiSymmetry ) :
                case ( DomainType::AntiSymmetryAir ) :
                case ( DomainType::AntiSymmetryFerro ) :
                case ( DomainType::AntiSymmetryConductor ) :
                {
                    return "AntiSymmetry";
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

            if ( tString == "conductor" ||  tString == "superconductor"  )
            {
                return DomainType::Conductor;
            }
            else if ( tString == "coil" )
            {
                return DomainType::Coil;
            }
            else if ( tString == "ferro" || tString == "ferromagnetic" )
            {
                return DomainType::Ferro;
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
            else if ( tString == "symmetry" )
            {
                return DomainType::Symmetry;
            }
            else if ( tString == "antisymmetry" )
            {
                return DomainType::AntiSymmetry;
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