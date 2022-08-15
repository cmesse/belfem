//
// Created by christian on 10/4/21.
//

#ifndef BELFEM_EN_FEM_DOMAINTYPE_HPP
#define BELFEM_EN_FEM_DOMAINTYPE_HPP

namespace belfem
{
    namespace fem
    {
        //
        enum class DomainType
        {
            Default = 0,         // default setting for blocks and sidesets
            Air             =  1,  // Maxwell Specific, do not change order!
            Coil            =  2,  // Maxwell Specific, do not change order!
            FerroMagnetic   =  3,  // Maxwell Specific, do not change order!
            SuperConductor  =  4,  // Maxwell Specific, do not change order!
            InterfaceScAir  =  5,  // Maxwell Specific
            InterfaceFmAir  =  6,  // Maxwell Specific
            InterfaceScFm   =  7,  // Maxwell Specific
            SymmetrySc      =  8,  // Maxwell Specific
            SymmetryFm      =  9,  // Maxwell Specific
            SymmetryAir     = 10,  // Maxwell Specific
            AntiSymmetrySc  = 11,  // Maxwell Specific
            AntiSymmetryFm  = 12,  // Maxwell Specific
            AntiSymmetryAir = 13,  // Maxwell Specific
            Cut             = 14,  // Maxwell Specific
            Boundary        = 15,  // Maxwell Specific
            ThinShell       = 16,  // Maxwell Specific
            Ghost           = 17,  // Maxwell Specific
            UNDEFINED = 18
        };

//------------------------------------------------------------------------------

        string
        to_string( const DomainType aDomaintype );

//------------------------------------------------------------------------------

        DomainType
        domain_type( const string & aString );

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_EN_FEM_DOMAINTYPE_HPP
